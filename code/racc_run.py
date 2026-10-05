"""The dynamics step of one docked/co-folded leg, run on the Reading ARC from a SLURM job.

The local tail (`run_dock_pose_md.sh` stages analyze onward) is repo-side, so the cluster runs
exactly one command: `omd run`. This wrapper carries the leg name as an env var so one file serves
every leg, and `STEPS` passes through for the 100k smoke launches (see CLAUDE.md, cluster section).

Staging, submission and the tail are all driven over ssh from the Mac; RACC.md (untracked) holds the
connection recipe and the host-specific facts, CLAUDE.md the workflow. Cluster-side paths assume the
workspace `~/remote_work/openmm` and env `~/.conda/envs/openmm-md` (RACC.md maps them).

Launch (from ~/remote_work/openmm, in an activated env):
    submit_gpu.sh racc_run.py          # full 20 ns
    STEPS=100000 submit_gpu.sh racc_run.py 10   # smoke
"""
import os
import shutil
import subprocess
import sys

LEG = os.environ.get("LEG", "shuffle_control_dock7")
WORKROOT = os.path.expanduser(os.environ.get("WORKROOT", "~/remote_work/openmm"))
# RUNREL is the leg's parent path under WORKROOT; defaults to octinoxate but is overridable for a
# second molecule (e.g. RUNREL=runs/oxybenzone/md). The scratch dcd dir keys on LEG only, so it needs
# no change as long as leg names stay unique across systems.
RUNREL = os.environ.get("RUNREL", "runs/octinoxate/md")
DCD_DIR = "/scratch5/gaussian/io927423/dcd"
STEPS = os.environ.get("STEPS", "10000000")
OMD = os.path.expanduser("~/.conda/envs/openmm-md/bin/omd")
base = os.path.join(WORKROOT, RUNREL, LEG)
prod = os.path.join(base, "prod_20ns")

# Preflight, before the GPU does anything: the 2026-10-04 quota fill killed both in-flight legs 70
# minutes into a launch, which a first-line check turns into an instant SLURM failure with no GPU
# time billed. The particle count comes off the serialized system (seconds on CPU); the dcd carries
# all of it, 3 float coords per particle per strided frame (stride 500).
sysxml = open(f"{base}/system/system.xml").read()
from openmm import XmlSerializer  # noqa: E402  -- the env's openmm, imported after the path is known
n_particles = XmlSerializer.deserialize(sysxml).getNumParticles()
dcd_bytes = int(n_particles * (int(STEPS) / 500 + 1) * 12 * 1.05)


def free_gb(path):
    return shutil.disk_usage(os.path.expanduser(path)).free / 2**30


need_scratch = dcd_bytes / 2**30 * 1.2  # 20% headroom over the projection
if free_gb(DCD_DIR) < need_scratch:
    print(f"ABORTING {LEG}: projected dcd {dcd_bytes / 2**30:.1f} GB; "
          f"scratch {free_gb(DCD_DIR):.1f} GB free (need {need_scratch:.1f})")
    sys.exit(1)
print(f"[preflight] {n_particles} particles, projected dcd {dcd_bytes / 2**30:.1f} GB, OK")

# The trajectory is WRITTEN to scratch, not written home and moved at the end: home's quota killed
# two legs mid-write on 2026-10-04, and an end-of-run move saves nothing when the kill is what
# happens mid-write. omd's `--out-dir` is pointed at a scratch directory via a symlink in the leg
# dir, so the dcd lands on scratch from its first byte while every path reading `<leg>/prod_20ns`
# (racc_run elsewhere, the Mac-side fetch, run_dock_pose_md.sh's tail) stays unchanged. The small
# out-dir residents (energy.csv, checkpoint.chk, final.pdb) ride along to scratch; home stays at
# staging + system size, a few MB.
os.makedirs(DCD_DIR, exist_ok=True)
scratch_prod = os.path.join(DCD_DIR, f"{LEG}_prod")
if os.path.lexists(prod) and not os.path.islink(prod):
    # a real prod dir means a legacy/failed run's leftovers; set them aside rather than delete
    os.replace(prod, prod + "_home_old")
os.makedirs(scratch_prod, exist_ok=True)
if not os.path.islink(prod):
    os.symlink(scratch_prod, prod)

# Same invocation as run_dock_pose_md.sh's dynamics stage, platform OpenCL: the conda-forge OpenMM
# build has no CUDA platform (RACC.md), and the platform string matches the OpenCL legs everywhere
# else in the project.
subprocess.run([OMD, "run",
                "--system", f"{base}/system/system.xml",
                "--topology", f"{base}/system/complex.pdb",
                "--out-dir", prod,
                "--steps", STEPS,
                "--platform", "OpenCL"], check=True)

# Strip to the wrapped solute HERE, on the cluster -- mirroring modal_md.py::produce. The full-system
# dcd (21k+ atoms, ~5 GB) stays on scratch for a possible re-strip/warm restart, and only the ~100 MB
# solute-only traj_wrapped.{xtc,pdb} is fetched to the Mac. That Mac has 8 GB of RAM and mdtraj loads a
# trajectory whole, so fetching the full box and stripping *locally* (as the racc lane did before) OOM'd
# the analyze step at ~4-5 GB. This is the SAME `omd analyze` run_dock_pose_md.sh runs locally
# (its line 93), so the wrapped solute is produced by identical code on both lanes; left unstrided so the
# frame count matches the other racc legs (bg33_1_dock2 etc.) rather than modal's strided 2,000.
subprocess.run([OMD, "analyze",
                "--traj", f"{prod}/traj.dcd",
                "--topology", f"{base}/system/complex.pdb",
                "--out-dir", prod], check=True)

# Fetched from the Mac side as `dcd/<LEG>_prod/traj_wrapped.{xtc,pdb}` (NOT the full dcd, which stays
# here); run_dock_pose_md.sh's whole-run analyze guard then skips, and its window loop slices the
# wrapped solute instead of the full box. Nothing needs moving.
print(f"[racc_run] full trajectory at {scratch_prod}/traj.dcd (kept on scratch)")
print(f"[racc_run] wrapped solute at {scratch_prod}/traj_wrapped.xtc -- fetch this, not the dcd")