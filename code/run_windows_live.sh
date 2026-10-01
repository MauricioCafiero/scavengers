#!/bin/zsh
# Compute the leading MM/GBSA windows while the dynamics is still running.
#
# There is no reason to wait for 20 ns to read the 5 ns number. The dynamics holds one CPU core on
# OpenCL; a window is analyze + prmtop build + MMPBSA.py over a few hundred frames of a ~400-atom
# solute, which this project measures at 4-16 minutes of one core. So the windows are nearly free on
# compute. What is not free is memory, since the machine is already under the run, so each window waits
# its turn and slices by streaming rather than loading a growing trajectory.
#
# The convergence series is the point of the windows, and for the docked pose it is the whole
# experiment: if Vina's pose is not a real minimum, the series should walk away from its first window
# the way orig_f12's did (-21.52 -> -13.71 over 5-20 ns), and there is no reason to learn that at the
# end rather than as it happens.
#
# run_dock_pose_md.sh skips any window whose FINAL_RESULTS_MMPBSA.dat already exists, so whatever this
# finishes early is simply not redone.
#
#   code/run_windows_live.sh <structure> [total_ns]
set -u
REPO=${0:A:h:h}
cd $REPO
OMD_ENV=${OMD_ENV:-$HOME/miniforge3/envs/openmm-md}   # override if the conda env is elsewhere
OMD=$OMD_ENV/bin/omd
PY=$OMD_ENV/bin/python
ROOT=$REPO

n=${1:?usage: run_windows_live.sh <structure> [total_ns]}
TOTAL=${2:-20}
M=runs/octinoxate/md/$n
P=$M/prod_${TOTAL}ns

for ns in 5 10 15; do
    W=$M/first_${ns}ns
    if [[ -f $W/mmgbsa/FINAL_RESULTS_MMPBSA.dat ]]; then
        echo "[$ns ns] already done"
        continue
    fi
    # Wait until the trajectory is past this window, with a margin so the slice is not taken at the
    # very head of a file being appended to.
    echo "[$ns ns] waiting for the run to pass $((ns + 1)) ns"
    while :; do
        [[ -f $P/energy.csv ]] || { sleep 60; continue }
        step=$(tail -1 $P/energy.csv | cut -d, -f1)
        [[ $step -gt $(( (ns + 1) * 500000 )) ]] && break
        sleep 120
    done
    mkdir -p $W
    echo "[$ns ns] starting $(date), step $step, $(memory_pressure | tail -1)"
    $PY code/md_window_live.py $P/traj.dcd $M/system/complex.pdb $W/traj.dcd $ns $P/energy.csv || continue
    $OMD analyze --traj $W/traj.dcd --topology $M/system/complex.pdb --out-dir $W
    # MMPBSA.py writes scratch into the WORKING directory, not --out-dir: one window left a 6 GB
    # reference.frc in a repository root. Running from inside the window's own directory contains it.
    ( cd $W && $OMD mmgbsa --protein $ROOT/$M/protein_fixed.pdb \
                           --ligand $ROOT/$M/ligand_prepped.sdf \
                           --traj traj_wrapped.xtc --topology traj_wrapped.pdb \
                           --out-dir mmgbsa --no-auto-cofactors --run )
    echo "[$ns ns] done $(date), $(memory_pressure | tail -1)"
    grep -B2 -A3 "DELTA TOTAL" $W/mmgbsa/FINAL_RESULTS_MMPBSA.dat 2>/dev/null | head -8
done
echo "WINDOWS_LIVE_DONE[$n]"
