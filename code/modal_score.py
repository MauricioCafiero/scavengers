"""Run the UMA scoring step on a Modal GPU instead of this laptop's CPU.

Scoring is the only expensive stage in the pipeline: 40-50 minutes per complex on six CPU cores, so
twelve complexes is most of a night. `binding_energy.py` already picks CUDA up when it is there, so
nothing about the science changes -- this script ships the inputs to a GPU container, runs that same
script unmodified, and brings the results back.

**Nothing local is modified.** The remote side calls `code/binding_energy.py` as a subprocess, so
there is one implementation of the scoring and this file cannot drift from it.

Everything that crosses the network is staged under `runs/<ligand>/modal/<tag>/`, which keeps what was
sent, exactly what came back, the remote log, and a manifest recording the GPU, the timings, the
package versions and a checksum per file. The results are then copied into the places the rest of the
pipeline looks -- `boltz/`, `boltz/structures/`, `uma_logs/` -- so the staging directory is the
provenance record and the canonical locations stay the working copies.

Usage:
    # everything except scoring, locally:
    SKIP_SCORING=1 zsh code/run_shell.sh s3 3 "tryptophan:2,..."

    # then the scoring on a GPU:
    modal run code/modal_score.py --structures-file runs/octinoxate/modal/s3.txt --out binding_shell3.csv
    modal run code/modal_score.py --structures s3_orig_control,s3_orig_f4 --gpu A10G

    # inspect without spending anything:
    modal run code/modal_score.py --structures s3_orig_control --dry-run

Cost. The container is billed per second while it runs, so the structures are scored **sequentially in
one container**: fairchem compiles its model on first use, and one container reuses that compilation
where a fan-out would pay it twelve times over. Rough list prices are L4 ~$0.80/hour and A10G
~$1.10/hour, and a 650-atom complex takes single-digit minutes on either, so twelve structures should
land near half an hour. `--max-minutes` is a hard timeout, not a guess: the container is killed when it
expires, which bounds the spend whatever happens. Check current prices at modal.com/pricing rather
than trusting the figures here.
"""
import json
import os
import sys
import time

import modal

APP_NAME = "peptidebuilder-uma-scoring"

# Pinned to what this repo runs locally, so a number computed on the GPU is comparable with one
# computed on the laptop. fairchem 2.22 refuses the UMA 1.0 checkpoint, hence uma-s-1p2p1 downstream.
PACKAGES = [
    "torch==2.13.0",
    "fairchem-core==2.22.0",
    "ase==3.29.0",
    "rdkit==2026.3.6",
    "numpy==2.4.6",
    "pdbfixer==1.12.0",     # brings openmm, for adding hydrogens to Boltz's heavy-atom output
]

image = (
    modal.Image.debian_slim(python_version="3.12")
    .pip_install(*PACKAGES)
    .add_local_dir("code", remote_path="/root/code")
)

app = modal.App(APP_NAME)

# The UMA checkpoint is around a gigabyte and comes from Hugging Face. Caching it in a Volume means
# only the first run pays for the download; later runs start scoring almost immediately.
cache = modal.Volume.from_name("peptidebuilder-model-cache", create_if_missing=True)


@app.function(
    image=image,
    gpu=os.environ.get("PEPTIDEBUILDER_MODAL_GPU", "L4"),
    volumes={"/cache": cache},
    secrets=[modal.Secret.from_name("huggingface-secret")],
    timeout=int(os.environ.get("PEPTIDEBUILDER_MODAL_TIMEOUT", "5400")),
)
def score_remote(payload: dict, structures: str, terms: str, out_name: str,
                 tf32: bool = False) -> dict:
    """Score the given structures on the GPU and return a gzipped tar of everything produced.

    `payload` maps a path relative to the run directory to its bytes: the ligand geometry and one
    Boltz CIF per structure. That is all `binding_energy.py` reads.
    """
    import io
    import subprocess
    import tarfile

    os.environ.setdefault("HF_HOME", "/cache/huggingface")
    os.environ.setdefault("TORCH_HOME", "/cache/torch")

    # TF32 off by default. binding_energy.py sets the device and nothing else, so torch's CUDA
    # defaults apply, and on Ampere or Ada hardware that can mean TF32 matmuls -- reduced mantissa
    # against the CPU's full float32. Measured on an L4, that costs agreement with a CPU score:
    # +0.9 kcal/mol on the interaction and +4.6 on the ligand strain for the same structure, with the
    # free-ligand relaxation taking 59 steps instead of 25 to reach the same fmax. The criteria are
    # identical; only the arithmetic differs, and a different trajectory finds a different minimum.
    #
    # NVIDIA_TF32_OVERRIDE is read by the CUDA libraries at initialisation, so setting it here
    # reaches the subprocess without any change to the scoring code. Pass --tf32 to allow it, which
    # is faster and fine if nothing is being compared against a CPU-scored result.
    if not tf32:
        os.environ["NVIDIA_TF32_OVERRIDE"] = "0"
    # fairchem reads whichever of these the installed version wants; setting both is harmless.
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
    # Everything returned must be a plain built-in type. `torch.__version__` is a TorchVersion, not a
    # str, so returning it unpickles only where torch is installed -- and torch is deliberately not
    # installed locally, so the whole result came back as a DeserializationError after the GPU had
    # already done the work.
    info = {
        "cuda_available": bool(torch.cuda.is_available()),
        "device": str(torch.cuda.get_device_name(0)) if torch.cuda.is_available() else "cpu",
        "torch": str(torch.__version__),
    }
    print(f"device: {info['device']}  (cuda={info['cuda_available']})", flush=True)
    if not info["cuda_available"]:
        # Better to fail loudly than to quietly bill GPU time for a CPU run that would be slower
        # than doing it locally.
        return {"error": "no CUDA device in the container; refusing to score on a rented CPU",
                "info": info}

    cmd = [sys.executable, "/root/code/binding_energy.py", run,
           "--terms", terms, "--out", out_name, "--structures", structures]
    print("running: " + " ".join(cmd), flush=True)
    t0 = time.time()
    proc = subprocess.run(cmd, capture_output=True, text=True, cwd="/root/code")
    elapsed = time.time() - t0
    log = proc.stdout + ("\n--- stderr ---\n" + proc.stderr if proc.stderr else "")
    print(f"binding_energy.py exited {proc.returncode} after {elapsed/60:.1f} min", flush=True)

    # Collect whatever it produced. Tarred and gzipped because a function's return value is size
    # limited and 36 xyz files plus a csv plus the partials would be wasteful sent raw.
    buf = io.BytesIO()
    produced = []
    with tarfile.open(fileobj=buf, mode="w:gz") as tar:
        for root, _dirs, files in os.walk(os.path.join(run, "boltz")):
            for fn in files:
                full = os.path.join(root, fn)
                rel = os.path.relpath(full, run)
                if rel in payload or rel.endswith("_model_0.cif"):
                    continue          # do not ship the inputs back
                tar.add(full, arcname=rel)
                produced.append(rel)
    cache.commit()
    return {
        "tar": buf.getvalue(),
        "log": log,
        "returncode": proc.returncode,
        "elapsed_s": elapsed,
        "produced": sorted(produced),
        "info": info,
    }


def _collect(outdir, names):
    """The inputs binding_energy.py reads: the ligand geometry and one CIF per structure."""
    payload, missing = {}, []
    lig = os.path.join(outdir, "ligand.xyz")
    if not os.path.exists(lig):
        sys.exit(f"no ligand.xyz under {outdir}")
    payload["ligand.xyz"] = open(lig, "rb").read()
    for n in names:
        rel = os.path.join("boltz", f"boltz_results_{n}", "predictions", n, f"{n}_model_0.cif")
        full = os.path.join(outdir, rel)
        if os.path.exists(full):
            payload[rel] = open(full, "rb").read()
        else:
            missing.append(n)
    return payload, missing


@app.local_entrypoint()
def main(structures: str = "", structures_file: str = "", outdir: str = "runs/octinoxate",
         out: str = "binding_modal.csv", terms: str = "interaction,strain_ligand",
         gpu: str = "L4", tag: str = "", max_minutes: int = 90,
         dry_run: bool = False, install: bool = True, force: bool = False,
         tf32: bool = False):
    """Stage the inputs, score them on a GPU, and copy the results into the usual places."""
    import hashlib
    import shutil
    import tarfile

    if structures_file:
        names = [ln.strip() for ln in open(structures_file) if ln.strip()
                 and not ln.startswith("#")]
    else:
        names = [s.strip() for s in structures.split(",") if s.strip()]
    if not names:
        sys.exit("give --structures a,b,c or --structures-file with one name per line")

    payload, missing = _collect(outdir, names)
    if missing:
        sys.exit(f"no folded structure for: {', '.join(missing)}\n"
                 f"fold them first -- this script only scores")

    # Do not pay to recompute something already on disk. binding_energy.py has no such check itself.
    already = [n for n in names
               if os.path.exists(os.path.join(outdir, "boltz", f"partial_{n}.json"))]
    if already and not force:
        names = [n for n in names if n not in already]
        print(f"already scored, skipping: {', '.join(already)}")
        if not names:
            sys.exit("every requested structure is already scored; pass --force to redo")
        payload, _ = _collect(outdir, names)

    tag = tag or time.strftime("%Y%m%d_%H%M%S")
    stage = os.path.join(outdir, "modal", tag)
    sent = os.path.join(stage, "sent")
    got = os.path.join(stage, "returned")
    os.makedirs(sent, exist_ok=True)
    os.makedirs(got, exist_ok=True)
    for rel, blob in payload.items():
        dest = os.path.join(sent, rel)
        os.makedirs(os.path.dirname(dest), exist_ok=True)
        open(dest, "wb").write(blob)

    mb = sum(len(b) for b in payload.values()) / 1e6
    print(f"\n{len(names)} structures, {mb:.2f} MB of input, gpu={gpu}, "
          f"hard timeout {max_minutes} min")
    print(f"staging under {stage}")
    for n in names:
        print(f"  {n}")
    if dry_run:
        print("\n--dry-run: nothing sent, nothing billed")
        return

    t0 = time.time()
    # .remote() is required to invoke; calling the Function object raises TypeError.
    res = score_remote.with_options(
        gpu=gpu, timeout=max_minutes * 60,
    ).remote(payload, ",".join(names), terms, out, tf32)
    wall = time.time() - t0

    if "error" in res:
        open(os.path.join(stage, "modal_run.log"), "w").write(json.dumps(res, indent=2))
        sys.exit(f"remote refused: {res['error']}")

    open(os.path.join(stage, "modal_run.log"), "w").write(res["log"])
    with tarfile.open(fileobj=__import__("io").BytesIO(res["tar"]), mode="r:gz") as tar:
        tar.extractall(got)

    manifest = {
        "tag": tag, "gpu": gpu, "tf32": tf32, "structures": names, "terms": terms, "out": out,
        "returncode": res["returncode"],
        "remote_scoring_minutes": round(res["elapsed_s"] / 60, 1),
        "wall_minutes_including_startup": round(wall / 60, 1),
        "minutes_per_structure": round(res["elapsed_s"] / 60 / max(len(names), 1), 1),
        "device": res["info"], "packages": PACKAGES,
        "produced": res["produced"],
        "sha256": {},
    }
    for root, _d, files in os.walk(got):
        for fn in files:
            full = os.path.join(root, fn)
            rel = os.path.relpath(full, got)
            manifest["sha256"][rel] = hashlib.sha256(open(full, "rb").read()).hexdigest()[:16]
    open(os.path.join(stage, "manifest.json"), "w").write(json.dumps(manifest, indent=2))

    print(f"\nremote scoring {manifest['remote_scoring_minutes']} min "
          f"({manifest['minutes_per_structure']} min/structure), "
          f"wall {manifest['wall_minutes_including_startup']} min including startup")
    print(f"device: {res['info']['device']}")
    print(f"{len(res['produced'])} files returned into {got}")
    if res["returncode"] != 0:
        print(f"\nWARNING: binding_energy.py exited {res['returncode']} -- "
              f"read {os.path.join(stage, 'modal_run.log')} before trusting the results")

    if not install:
        print("\n--no-install: results left in the staging directory only")
        return

    # Copy into the places the rest of the pipeline looks. The staging copy stays as the record.
    copied, skipped = [], []
    for rel in sorted(manifest["sha256"]):
        src = os.path.join(got, rel)
        dst = os.path.join(outdir, rel)
        if os.path.exists(dst) and not force:
            skipped.append(rel)
            continue
        os.makedirs(os.path.dirname(dst), exist_ok=True)
        shutil.copy2(src, dst)
        copied.append(rel)
    log_dst = os.path.join(outdir, "uma_logs", f"score_{tag}_modal.log")
    os.makedirs(os.path.dirname(log_dst), exist_ok=True)
    shutil.copy2(os.path.join(stage, "modal_run.log"), log_dst)

    print(f"\ninstalled {len(copied)} files into {outdir}")
    if skipped:
        print(f"left alone because they already exist ({len(skipped)}): "
              f"{', '.join(skipped[:4])}{' ...' if len(skipped) > 4 else ''}")
        print("pass --force to overwrite")
    print(f"scoring log -> {log_dst}")
    print(f"provenance  -> {os.path.join(stage, 'manifest.json')}")
