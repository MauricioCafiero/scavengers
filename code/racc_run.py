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
import subprocess

LEG = os.environ.get("LEG", "shuffle_control_dock7")
WORKROOT = os.path.expanduser(os.environ.get("WORKROOT", "~/remote_work/openmm"))
STEPS = os.environ.get("STEPS", "10000000")
OMD = os.path.expanduser("~/.conda/envs/openmm-md/bin/omd")
base = os.path.join(WORKROOT, "runs/octinoxate/md", LEG)

# Same invocation as run_dock_pose_md.sh's dynamics stage, platform OpenCL: the conda-forge OpenMM
# build has no CUDA platform (RACC.md), and the platform string matches the OpenCL legs everywhere
# else in the project.
subprocess.run([OMD, "run",
                "--system", f"{base}/system/system.xml",
                "--topology", f"{base}/system/complex.pdb",
                "--out-dir", f"{base}/prod_20ns",
                "--steps", STEPS,
                "--platform", "OpenCL"], check=True)