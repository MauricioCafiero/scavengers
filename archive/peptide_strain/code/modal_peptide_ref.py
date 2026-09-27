"""Run the apo peptide hydrogen relaxation on a Modal GPU.

`peptide_apo.py fold` produces the apo structures locally, because Boltz is fast on this laptop's GPU.
`peptide_apo.py energy` then protonates each one, relaxes its hydrogens with the heavy atoms fixed, and
takes a single point -- that part is UMA and slow on CPU, so this ships it to a rented GPU and brings
back `peptide_reference.json`.

It is the same script on both sides. `peptide_apo.py energy` picks CUDA when it sees a device and CPU
otherwise, so running it here instead costs time and nothing else; there is no second implementation to
drift. Use the local route when there is no credit to spend -- there are only three apo folds per shell.

Staged through `runs/<ligand>/modal/<tag>/` exactly as `modal_score.py` does, so what crossed the
network is recoverable independently of what was installed.

Usage:
    python code/peptide_apo.py fold runs/octinoxate --shell 3     # first, locally
    modal run code/modal_peptide_ref.py --shell 3
"""
import json
import os
import sys
import time

import modal

PACKAGES = [
    "torch==2.13.0",
    "fairchem-core==2.22.0",
    "ase==3.29.0",
    "rdkit==2026.3.6",
    "numpy==2.4.6",
    "pdbfixer==1.12.0",
]

image = (
    modal.Image.debian_slim(python_version="3.12")
    .pip_install(*PACKAGES)
    .add_local_dir("code", remote_path="/root/code")
)

app = modal.App("peptidebuilder-apo-reference")
cache = modal.Volume.from_name("peptidebuilder-model-cache", create_if_missing=True)


@app.function(image=image, gpu="L4", volumes={"/cache": cache},
              secrets=[modal.Secret.from_name("huggingface-secret")], timeout=5400)
def apo_remote(payload: dict, shell: str, fmax: float, steps: int, tf32: bool = False) -> dict:
    """Relax the apo structures' hydrogens on the GPU and return peptide_reference.json."""
    import subprocess

    os.environ.setdefault("HF_HOME", "/cache/huggingface")
    os.environ.setdefault("TORCH_HOME", "/cache/torch")
    # TF32 off by default, to match the laptop's arithmetic. See modal_score.py for the measurement
    # behind this: it made no difference to the one case tested, but the point of a reference is that
    # it is comparable with everything else, so the cheaper default is the conservative one.
    if not tf32:
        os.environ["NVIDIA_TF32_OVERRIDE"] = "0"
    for var in ("HF_TOKEN", "HUGGINGFACE_TOKEN", "HUGGING_FACE_HUB_TOKEN"):
        tok = os.environ.get("HF_TOKEN") or os.environ.get("HUGGINGFACE_TOKEN")
        if tok:
            os.environ.setdefault(var, tok)

    run = "/tmp/run"
    for rel, blob in payload.items():
        dest = os.path.join(run, rel)
        os.makedirs(os.path.dirname(dest), exist_ok=True)
        with open(dest, "wb") as fh:
            fh.write(blob)

    import torch
    if not torch.cuda.is_available():
        return {"error": "no CUDA device in the container; refusing to run on a rented CPU"}
    device = str(torch.cuda.get_device_name(0))
    print(f"device: {device}", flush=True)

    cmd = [sys.executable, "/root/code/peptide_apo.py", "energy", run,
           "--fmax", str(fmax), "--steps", str(steps)]
    if shell:
        cmd += ["--shell", str(shell)]
    print("running: " + " ".join(cmd), flush=True)
    t0 = time.time()
    proc = subprocess.run(cmd, capture_output=True, text=True, cwd="/root/code")
    elapsed = time.time() - t0
    log = proc.stdout + ("\n--- stderr ---\n" + proc.stderr if proc.stderr else "")
    print(f"peptide_apo.py exited {proc.returncode} after {elapsed/60:.1f} min", flush=True)

    ref_path = os.path.join(run, "peptide_reference.json")
    ref = json.load(open(ref_path)) if os.path.exists(ref_path) else {}
    xyz = {}
    apo = os.path.join(run, "apo")
    if os.path.isdir(apo):
        for fn in sorted(os.listdir(apo)):
            if fn.endswith("_relaxed_h.xyz"):
                xyz[fn] = open(os.path.join(apo, fn)).read()
    cache.commit()
    return {"reference": ref, "xyz": xyz, "log": log, "returncode": proc.returncode,
            "elapsed_s": elapsed, "device": device}


@app.local_entrypoint()
def main(shell: str = "", outdir: str = "runs/octinoxate", fmax: float = 0.10, steps: int = 75,
         tag: str = "", max_minutes: int = 90, dry_run: bool = False, install: bool = True):
    import glob
    import shutil

    apo_dir = os.path.join(outdir, "apo")
    pattern = os.path.join(apo_dir, "boltz_results_apo_s*", "predictions", "*", "*_model_0.cif")
    cifs = sorted(glob.glob(pattern))
    if shell:
        cifs = [c for c in cifs if f"apo_s{shell}_" in os.path.basename(c)]
    if not cifs:
        sys.exit(f"no apo structures under {apo_dir}\n"
                 f"run: python code/peptide_apo.py fold {outdir} --shell {shell or '<N>'}")

    payload = {}
    for c in cifs:
        payload[os.path.relpath(c, outdir)] = open(c, "rb").read()

    tag = tag or f"apo_shell{shell or 'all'}_" + time.strftime("%Y%m%d_%H%M%S")
    stage = os.path.join(outdir, "modal", tag)
    os.makedirs(os.path.join(stage, "sent"), exist_ok=True)
    os.makedirs(os.path.join(stage, "returned"), exist_ok=True)
    for rel, blob in payload.items():
        dest = os.path.join(stage, "sent", rel)
        os.makedirs(os.path.dirname(dest), exist_ok=True)
        open(dest, "wb").write(blob)

    print(f"\n{len(cifs)} apo structures, {sum(len(b) for b in payload.values())/1e6:.2f} MB, "
          f"fmax {fmax}, {steps} steps, hard timeout {max_minutes} min")
    for c in cifs:
        print(f"  {os.path.basename(c).replace('_model_0.cif','')}")
    if dry_run:
        print("\n--dry-run: nothing sent, nothing billed")
        return

    res = apo_remote.with_options(timeout=max_minutes * 60).remote(payload, shell, fmax, steps)
    if "error" in res:
        sys.exit(f"remote refused: {res['error']}")

    open(os.path.join(stage, "modal_run.log"), "w").write(res["log"])
    for fn, text in res["xyz"].items():
        open(os.path.join(stage, "returned", fn), "w").write(text)
    open(os.path.join(stage, "returned", "peptide_reference.json"), "w").write(
        json.dumps(res["reference"], indent=2))
    print(f"\nremote {res['elapsed_s']/60:.1f} min on {res['device']}, "
          f"{len(res['reference'])} references returned")
    if res["returncode"] != 0:
        print(f"WARNING: peptide_apo.py exited {res['returncode']} -- read "
              f"{os.path.join(stage, 'modal_run.log')} before trusting these")

    if not install:
        print("--no-install: left in the staging directory only")
        return

    # Merge rather than overwrite: peptide_reference.json accumulates across shells, and each shell's
    # entries are keyed by apo structure name, so a later shell must not drop an earlier one's.
    out = os.path.join(outdir, "peptide_reference.json")
    store = json.load(open(out)) if os.path.exists(out) else {}
    store.update(res["reference"])
    with open(out, "w") as fh:
        json.dump(store, fh, indent=2)
    for fn in res["xyz"]:
        shutil.copy2(os.path.join(stage, "returned", fn), os.path.join(apo_dir, fn))
    print(f"merged {len(res['reference'])} into {out} ({len(store)} total)")
    for name, d in sorted(res["reference"].items()):
        print(f"  {name:<22} E_apo {d['e_apo_kcal']:>14.3f}  "
              f"({d['n_atoms']} atoms, charge {d['peptide_charge']:+d})")
