"""Run the MD production leg on a Modal GPU, using MD_openmm's own dynamics code.

Dynamics is the one GPU-bound step in the chain. Everything else stays local: `omd build` needs the
OpenFF/GAFF2 stack to parameterise the ligand, and `omd mmgbsa` needs AmberTools -- both conda-only and
both cheap on CPU. Only `run` wants a real GPU, and this laptop has no CUDA.

Measured, on a 21,225-atom WT-MetaD system in the rotaxanes work: local CPU 27.874 ms/step against
Modal A10G 0.570 ms/step, a 49x speedup. Against this laptop's OpenCL (1.72 ms/step at 8,647 atoms,
so roughly 4.6 ms/step scaled to 23k) an A10G is about 8x faster, which turns a 20 ns run of a
23,121-atom complex from ~14 hours into ~1.6 hours for about $1.74.

**This imports `openmm_md.dynamics.run` rather than reimplementing the protocol.** A Modal trajectory
and a local one are then produced by the same tested code -- same minimisation, same restrained
equilibration, same integrator, same timestep, same reporters. That is not a nicety: these trajectories
are compared as binding energies, and a protocol difference would be indistinguishable from a result.

`dynamics.py` imports only numpy, openmm and `.config`, and `openmm_md/__init__.py` imports nothing, so
that import pulls none of OpenFF, openmmforcefields, ParmEd or AmberTools. The image is correspondingly
light. Going through the `omd` CLI instead would drag in the whole stack, because `cli.py` imports
`build_system` at module level.

Lessons taken from the two prior Modal/GPU ports, each of which cost something to learn:

* **Warm up before timing.** OpenMM's CUDA platform compiles and autotunes lazily, so an un-warmed probe
  reads 2-5x too slow. The rotaxanes Colab notebook records learning this "the expensive way on A10G".
* **Explicit platform, never `auto`.** MD_openmm's `auto` probe omits CUDA (`dynamics.py:55`), so on a
  GPU host it silently selects CPU -- renting a GPU to compute on a processor. Passing `CUDA` explicitly
  takes the `requested != "auto"` branch, which raises if CUDA is missing instead of degrading.
* **Full-system trajectory on the Volume, stripped copy returned.** A solute-only trajectory breaks warm
  restarts, since `setPositions` needs every atom, and it only fails *after* a leg has been paid for.
  The storage it saves is pennies a month against dollars an hour of GPU.
* **The image block must stay byte-identical between scripts** that share it, or the cached build is
  missed and billed builder minutes are spent rebuilding.

Usage:
    modal run code/modal_md.py::push    --structure s2_esm2_control
    modal run code/modal_md.py::probe   --structure s2_esm2_control          # cheap, measures ns/day
    modal run code/modal_md.py::produce --structure s2_esm2_control --steps 10000000
    modal run code/modal_md.py::strip   --structure s2_esm2_control --stride 10
    modal run code/modal_md.py::fetch   --structure s2_esm2_control

`produce` strips the solute itself, so `strip` is only needed when that step failed or when a different
stride is wanted -- it works from the trajectory already on the Volume and costs no GPU.
"""
import json
import os
import sys
import time
from pathlib import Path

import modal

GPU = os.environ.get("PEPTIDEBUILDER_MD_GPU", "A10G")
RATE_USD_PER_S = {"A10G": 1.10, "L4": 0.80, "T4": 0.59, "A100": 2.10}.get(GPU, 1.10) / 3600.0

# MD_openmm's source, imported remotely. Not vendored: it must stay the same code that produced the
# local trajectories.
MD_SRC = Path("~/python_mac/openmm/src").expanduser()

# Modal's own micromamba image, so conda packages land in the DEFAULT environment and the container's
# interpreter *is* the conda interpreter. A debian_slim image with micromamba bolted on creates a
# separate env, and importing its openmm from the system Python fails with `CXXABI_1.3.15 not found`
# because the compiled libOpenMM links conda's libstdc++, not the system one.
#
# `cuda-version=12.9` is the pin that matters and is not optional: conda-forge's unpinned openmm is
# built against a newer toolkit than Modal's host driver accepts, and a Context then dies with
# `CUDA_ERROR_UNSUPPORTED_PTX_VERSION (222)` -- after the GPU has been billed. This combination is the
# one the rotaxanes benchmark actually ran on (0.570 ms/step, A10G).
#
# openmm=8.4 matches the local openmm-md env, so a Modal trajectory and a local one come from the same
# OpenMM as well as the same protocol code.
#
# Keep this block byte-identical if it is ever copied into a sibling script: any difference misses the
# cached build and spends billed builder minutes.
image = (
    modal.Image.micromamba(python_version="3.11")
    .micromamba_install("openmm=8.4", "mdtraj", "numpy", "cuda-version=12.9",
                        channels=["conda-forge"])
    .env({"OMP_NUM_THREADS": "4"})
    .add_local_dir(str(MD_SRC), remote_path="/root/src")
)

app = modal.App("peptidebuilder-md", image=image)
state = modal.Volume.from_name("peptidebuilder-md-state", create_if_missing=True)

# ---------------------------------------------------------------------------------------------
# no GPU on these two: staging and retrieval must never attach one
# ---------------------------------------------------------------------------------------------
@app.function(volumes={"/state": state}, timeout=60 * 20)
def _push(structure: str, payload: dict) -> list:
    d = f"/state/{structure}"
    os.makedirs(d, exist_ok=True)
    for name, blob in payload.items():
        with open(os.path.join(d, name), "wb") as fh:
            fh.write(blob)
    state.commit()
    return sorted(os.listdir(d))


@app.function(volumes={"/state": state}, timeout=60 * 30)
def _fetch(structure: str, names: list) -> dict:
    d = f"/state/{structure}"
    if not os.path.isdir(d):
        return {"error": f"nothing staged for {structure!r}"}
    out = {}
    for n in names or sorted(os.listdir(d)):
        p = os.path.join(d, n)
        if os.path.isfile(p):
            out[n] = open(p, "rb").read()
    return out


@app.function(volumes={"/state": state}, timeout=60 * 10)
def _ls(structure: str) -> list:
    d = f"/state/{structure}"
    if not os.path.isdir(d):
        return []
    return sorted(f"{n}  {os.path.getsize(os.path.join(d,n))/1e6:.1f} MB"
                  for n in os.listdir(d) if os.path.isfile(os.path.join(d, n)))


# ---------------------------------------------------------------------------------------------
# GPU
# ---------------------------------------------------------------------------------------------
@app.function(gpu=GPU, volumes={"/state": state}, timeout=60 * 25)
def _probe(structure: str, steps: int = 4000, warmup: int = 2000) -> dict:
    """Measure the real step rate, discarding a warmup. Never trust an un-warmed number.

    OpenMM's CUDA platform compiles and autotunes lazily, so the first thousands of steps run several
    times slower than steady state. The rotaxanes Colab notebook records reading 2-5x too slow from an
    un-warmed probe, "the expensive way on A10G".
    """
    import openmm as mm
    from openmm import app as oapp, unit, XmlSerializer

    d = f"/state/{structure}"
    try:
        platform = mm.Platform.getPlatformByName("CUDA")
    except Exception as e:
        return {"error": f"no CUDA platform in the container: {e}"}

    pdb = oapp.PDBFile(f"{d}/complex.pdb")
    system = XmlSerializer.deserializeSystem(open(f"{d}/system.xml").read())
    integ = mm.LangevinMiddleIntegrator(300 * unit.kelvin, 1 / unit.picosecond,
                                        2 * unit.femtosecond)
    integ.setConstraintTolerance(1e-6)
    try:
        sim = oapp.Simulation(pdb.topology, system, integ, platform)
    except Exception as e:
        # A Context that dies here is the cuda-version mismatch, not a bad system.
        return {"error": f"CUDA context failed ({e}); check the cuda-version pin in the image"}
    sim.context.setPositions(pdb.positions)
    sim.minimizeEnergy(maxIterations=200)
    print(f"minimised {system.getNumParticles()} particles, warming up {warmup} steps", flush=True)

    sim.step(warmup)
    t0 = time.time()
    sim.step(steps)
    dt = time.time() - t0

    ms = 1000.0 * dt / steps
    return {"gpu": GPU, "openmm": mm.version.version,
            "n_particles": system.getNumParticles(),
            "ms_per_step": round(ms, 4),
            "ns_per_day": round((steps * 2e-6) / dt * 86400, 1),
            "warmup_steps": warmup, "timed_steps": steps,
            "est_20ns_hours": round(10_000_000 * ms / 1000 / 3600, 2),
            "est_20ns_usd": round(10_000_000 * ms / 1000 * RATE_USD_PER_S, 2)}


@app.function(gpu=GPU, volumes={"/state": state}, timeout=60 * 60 * 6)
def _produce(structure: str, steps: int, leg: str = "L1", stride: int = 10) -> dict:
    """The production leg, via MD_openmm's own dynamics.run, explicitly on CUDA."""
    sys.path.insert(0, "/root/src")
    import openmm as mm
    try:
        mm.Platform.getPlatformByName("CUDA")
    except Exception as e:
        return {"error": f"no CUDA; refusing to run on a rented CPU: {e}"}

    from openmm_md.config import Config
    from openmm_md.dynamics import run as run_dynamics

    d = Path(f"/state/{structure}")
    out = d / f"prod_{leg}"
    out.mkdir(parents=True, exist_ok=True)

    cfg = Config()
    cfg.platform = "CUDA"      # explicit: `auto` in this version cannot select CUDA
    t0 = time.time()
    run_dynamics(d / "system.xml", d / "complex.pdb", out, cfg, steps)
    dt = time.time() - t0

    # Strip to the solute for the return trip. The full-system dcd stays on the Volume so a warm
    # restart remains possible -- setPositions needs every atom, and a solute-only frame fails.
    try:
        wrapped = _strip_solute(out / "traj.dcd", d / "complex.pdb", out, stride)
    except Exception as e:
        wrapped = {"error": str(e)[:200]}

    state.commit()
    return {"gpu": GPU, "steps": steps, "leg": leg,
            "minutes": round(dt / 60, 1), "usd": round(dt * RATE_USD_PER_S, 2),
            "ms_per_step": round(1000 * dt / steps, 4), "wrapped": wrapped,
            "files": sorted(p.name for p in out.iterdir() if p.is_file())}


def _strip_solute(dcd: Path, top: Path, out: Path, stride: int = 10) -> dict:
    """Write the wrapped, solute-only trajectory `omd mmgbsa` needs, beside the full one.

    Order matters, and getting it backwards cost a leg's worth of confusion: **image_molecules on the
    FULL solvated system first, then slice to the solute**. mdtraj picks periodic-imaging anchors by
    molecule size, so on a solute-only trajectory (648 atoms here, largest molecule 604) nothing
    qualifies and it raises "Could not find any anchor molecules". With water present the peptide is
    unambiguously the anchor. This is the order `openmm_md.analyze._write_wrapped` uses, which is why
    the local path always worked and the first Modal leg's strip did not.

    Chunked, and strided: 20,000 frames of 23,121 atoms is ~5.5 GB loaded whole, and MM/GBSA converges
    on far fewer frames than that -- the 20 ns local reference moved by 0.04 kcal/mol over its last
    decile.
    """
    import mdtraj as md

    ions = {"NA", "CL", "K", "MG", "CA", "ZN", "SOD", "CLA", "POT"}
    waters = {"HOH", "WAT", "TIP3", "SOL"}
    acc, sel = None, None
    for chunk in md.iterload(str(dcd), top=str(top), chunk=50, stride=stride):
        w = chunk.image_molecules(inplace=False)          # full system: anchors resolvable
        if sel is None:
            sel = [a.index for res in w.topology.residues
                   if res.name not in waters | ions for a in res.atoms]
            if not sel:
                raise ValueError("no solute atoms found to wrap")
        s = w.atom_slice(sel)
        acc = s if acc is None else acc.join(s)
    if acc is None:
        raise ValueError(f"no frames read from {dcd}")
    acc = acc.center_coordinates()
    acc.save(str(out / "traj_wrapped.xtc"))
    acc[0].save(str(out / "traj_wrapped.pdb"))
    return {"frames": int(acc.n_frames), "atoms": int(acc.n_atoms), "stride": stride}


@app.function(volumes={"/state": state}, timeout=60 * 40, memory=8192)
def _strip(structure: str, leg: str = "L1", stride: int = 10) -> dict:
    """Redo the strip from the trajectory already on the Volume. No GPU -- this is mdtraj on CPU."""
    d = Path(f"/state/{structure}")
    out = d / f"prod_{leg}"
    dcd = out / "traj.dcd"
    if not dcd.exists():
        return {"error": f"no trajectory at {dcd}"}
    try:
        r = _strip_solute(dcd, d / "complex.pdb", out, stride)
    except Exception as e:
        return {"error": str(e)[:300]}
    state.commit()
    r["files"] = sorted(p.name for p in out.iterdir() if p.is_file())
    return r


# ---------------------------------------------------------------------------------------------
# local entrypoints
# ---------------------------------------------------------------------------------------------
def _local_dir(structure, outdir="runs/octinoxate"):
    d = Path(outdir) / "md" / structure / "system"
    if not (d / "system.xml").exists():
        sys.exit(f"no built system at {d}\nrun `omd build` locally first")
    return d


@app.local_entrypoint()
def push(structure: str, outdir: str = "runs/octinoxate"):
    """Stage system.xml + complex.pdb onto the Volume. Once per structure."""
    d = _local_dir(structure, outdir)
    payload = {n: (d / n).read_bytes() for n in ("system.xml", "complex.pdb")}
    mb = sum(len(v) for v in payload.values()) / 1e6
    print(f"pushing {mb:.1f} MB for {structure}")
    print("on the volume:", _push.remote(structure, payload))


@app.local_entrypoint()
def probe(structure: str, steps: int = 4000, warmup: int = 2000):
    """Measure ns/day on the GPU before committing to a long run. A few minutes, cents."""
    r = _probe.remote(structure, steps, warmup)
    if "error" in r:
        sys.exit(r["error"])
    print(json.dumps(r, indent=2))
    print(f"\n20 ns would be about {r['est_20ns_hours']} h and ${r['est_20ns_usd']} on {GPU}")


@app.local_entrypoint()
def produce(structure: str, steps: int = 10_000_000, leg: str = "L1",
            outdir: str = "runs/octinoxate", install: bool = True):
    """Run the production leg, then pull the stripped trajectory back for local mmgbsa."""
    t0 = time.time()
    r = _produce.remote(structure, steps, leg)
    if "error" in r:
        sys.exit(r["error"])
    print(json.dumps({k: v for k, v in r.items() if k != "files"}, indent=2))
    print("files on the volume:", r["files"])
    if not install:
        return
    got = _fetch.remote(structure, [f"prod_{leg}/traj_wrapped.xtc",
                                    f"prod_{leg}/traj_wrapped.pdb",
                                    f"prod_{leg}/energy.csv"])
    dest = Path(outdir) / "md" / structure / f"prod_{leg}_modal"
    dest.mkdir(parents=True, exist_ok=True)
    for name, blob in got.items():
        if name == "error":
            print(f"fetch: {blob}")
            continue
        p = dest / Path(name).name
        p.write_bytes(blob)
        print(f"  wrote {p} ({len(blob)/1e6:.1f} MB)")
    print(f"\nwall {(time.time()-t0)/60:.1f} min. Next, locally:")
    m = Path(outdir) / "md" / structure
    print(f"  ~/miniforge3/envs/openmm-md/bin/omd mmgbsa --protein {m}/protein_fixed.pdb \\")
    print(f"      --ligand {m}/ligand_prepped.sdf \\")
    print(f"      --traj {dest}/traj_wrapped.xtc --topology {dest}/traj_wrapped.pdb \\")
    print(f"      --out-dir {dest}/mmgbsa --no-auto-cofactors --run")


@app.local_entrypoint()
def fetch(structure: str, outdir: str = "runs/octinoxate", leg: str = "L1"):
    """Pull the stripped trajectory and energy log for a leg already run."""
    got = _fetch.remote(structure, [f"prod_{leg}/traj_wrapped.xtc",
                                    f"prod_{leg}/traj_wrapped.pdb",
                                    f"prod_{leg}/energy.csv"])
    dest = Path(outdir) / "md" / structure / f"prod_{leg}_modal"
    dest.mkdir(parents=True, exist_ok=True)
    for name, blob in got.items():
        if name == "error":
            sys.exit(blob)
        (dest / Path(name).name).write_bytes(blob)
        print(f"  wrote {dest / Path(name).name} ({len(blob)/1e6:.1f} MB)")


@app.local_entrypoint()
def ls(structure: str):
    """What is on the Volume for this structure. No GPU."""
    for line in _ls.remote(structure):
        print("  " + line)


@app.local_entrypoint()
def strip(structure: str, leg: str = "L1", stride: int = 10, outdir: str = "runs/octinoxate"):
    """Re-do the solute strip on the Volume, then pull it back. Use when produce's strip failed."""
    r = _strip.remote(structure, leg, stride)
    if "error" in r:
        sys.exit(r["error"])
    print(json.dumps(r, indent=2))
    got = _fetch.remote(structure, [f"prod_{leg}/traj_wrapped.xtc", f"prod_{leg}/traj_wrapped.pdb"])
    dest = Path(outdir) / "md" / structure / f"prod_{leg}_modal"
    dest.mkdir(parents=True, exist_ok=True)
    for name, blob in got.items():
        if name == "error":
            sys.exit(blob)
        (dest / Path(name).name).write_bytes(blob)
        print(f"  wrote {dest / Path(name).name} ({len(blob)/1e6:.1f} MB)")
