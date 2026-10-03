#!/bin/zsh
# Convergence windows for a Modal leg, from its strided solute-only trajectory.
#
# Modal runs (code/modal_md.py) return ONLY traj_wrapped.xtc -- stride 10, ~2000
# frames over 20 ns -- with the full traj.dcd left on the Volume. So the local
# run_windows_live.sh path (slice traj.dcd, re-run `omd analyze` to wrap) does not
# apply. This slices the already-wrapped solute trajectory to each leading window
# and runs MM/GBSA directly on the slice. The subsample is uniform so the window
# means are unbiased; standard error is ~0.07 (vs ~0.04 unstrided) on 500-1500
# frames, which is fine for a convergence shape.
#
#   code/md_window_modal.sh <structure> [total_ns=20] [leg=prod_L1_modal]
set -u
REPO=${0:A:h:h}; cd $REPO; ROOT=$REPO
OMD=${OMD_ENV:-$HOME/miniforge3/envs/openmm-md}/bin/omd
PY=${OMD_ENV:-$HOME/miniforge3/envs/openmm-md}/bin/python
N=${1:?usage: md_window_modal.sh <structure> [total_ns] [leg]}
TOTAL=${2:-20}; LEG=${3:-prod_L1_modal}
M=runs/octinoxate/md/$N; P=$M/$LEG
for ns in 5 10 15; do
  W=$M/first_${ns}ns
  [[ -f $W/mmgbsa/FINAL_RESULTS_MMPBSA.dat ]] && { echo "[$ns ns] already done"; continue; }
  mkdir -p $W
  cp $P/traj_wrapped.pdb $W/traj_wrapped.pdb
  $PY -c "
import mdtraj as md
t = md.load('$P/traj_wrapped.xtc', top='$P/traj_wrapped.pdb')
k = t[:max(1, int(round(len(t)*$ns/$TOTAL)))]
k.save_xtc('$W/traj_wrapped.xtc')
print(f'[$ns ns] kept {len(k)}/{len(t)} frames')
"
  ( cd $W && $OMD mmgbsa --protein $ROOT/$M/protein_fixed.pdb --ligand $ROOT/$M/ligand_prepped.sdf \
                         --traj traj_wrapped.xtc --topology traj_wrapped.pdb \
                         --out-dir mmgbsa --no-auto-cofactors --run >/dev/null 2>&1 )
  rm -f $W/reference.frc; rm -rf $W/_MMPBSA_*   # scatter, per CLAUDE.md
  echo -n "[$ns ns] "; grep "DELTA TOTAL" $W/mmgbsa/FINAL_RESULTS_MMPBSA.dat 2>/dev/null || echo "MMGBSA FAILED"
done
echo "WINDOW_MODAL_DONE[$N]"
