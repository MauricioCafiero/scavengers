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
DCD_DIR = "/scratch5/gaussian/io927423/dcd"
STEPS = os.environ.get("STEPS", "10000000")
OMD = os.path.expanduser("~/.conda/envs/openmm-md/bin/omd")
base = os.path.join(WORKROOT, "runs/octinoxate/md", LEG)

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


need_home = 2.0                     # checkpoints, energy.csv and the leg dir also write home
need_dcd = dcd_bytes / 2**30 * 1.2  # the dcd lands on scratch; 20% headroom over projection
if free_gb(WORKROOT) < need_home or free_gb(DCD_DIR) < need_dcd:
    print(f"ABORTING {LEG}: projected dcd {dcd_bytes / 2**30:.1f} GB; free "
          f"home {free_gb(WORKROOT):.1f} (need {need_home:.0f}), "
          f"scratch {free_gb(DCD_DIR):.1f} (need {need_dcd:.1f})")
    sys.exit(1)
print(f"[preflight] {n_particles} particles, projected dcd {dcd_bytes / 2**30:.1f} GB, OK")

# Same invocation as run_dock_pose_md.sh's dynamics stage, platform OpenCL: the conda-forge OpenMM
# build has no CUDA platform (RACC.md), and the platform string matches the OpenCL legs everywhere
# else in the project.
subprocess.run([OMD, "run",
                "--system", f"{base}/system/system.xml",
                "--topology", f"{base}/system/complex.pdb",
                "--out-dir", f"{base}/prod_20ns",
                "--steps", STEPS,
                "--platform", "OpenCL"], check=True)

# The dcd is ~1.8 GB per leg and home NFS is quota-tight (a full run filled it and killed both
# 2026-10-04 legs), so the trajectory moves to scratch the moment it is written. Named by leg — the
# workspace copy is always traj.dcd — and the fetch on the Mac side reads it back from here.
dcd = os.path.join(base, "prod_20ns", "traj.dcd")
if os.path.exists(dcd):
    os.makedirs(DCD_DIR, exist_ok=True)
    shutil.move(dcd, os.path.join(DCD_DIR, f"{LEG}.dcd"))