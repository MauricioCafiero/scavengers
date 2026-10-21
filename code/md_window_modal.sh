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
#
# MODE/WINDOWS (2026-10-21) so this also cuts TRAILING windows, which nothing in the project did.
# Every window here and in run_dock_pose_md.sh is a LEADING cumulative mean, so a whole-run value
# below the 10 ns value implies a stronger back half rather than stating one -- and on
# shuffle_control_esm0_rf3cof_r2 the implication was large (10 ns -11.33, 20 ns -15.03, so the last
# 10 ns averages about -18.7 by subtraction). A trailing window measures that half directly, which is
# the only way to tell a seed effect from a leg that simply had not settled. MODE=trailing names the
# directory last_<ns>ns; the leading names are untouched, so no existing cell moves.
#   MODE=trailing WINDOWS=10 code/md_window_modal.sh <structure> 20 prod_20ns
set -u
REPO=${0:A:h:h}; cd $REPO; ROOT=$REPO
OMD=${OMD_ENV:-$HOME/miniforge3/envs/openmm-md}/bin/omd
PY=${OMD_ENV:-$HOME/miniforge3/envs/openmm-md}/bin/python
N=${1:?usage: md_window_modal.sh <structure> [total_ns] [leg]}
TOTAL=${2:-20}; LEG=${3:-prod_L1_modal}
# SYSNAME parameterizes the run dir, which was hardcoded to octinoxate. A racc leg returns the same
# wrapped-solute pair this slices, so the Modal path applies unchanged to oxybenzone's legs.
SYSNAME=${SYSNAME:-octinoxate}
M=runs/${SYSNAME}/md/$N; P=$M/$LEG
MODE=${MODE:-leading}
WINDOWS=${WINDOWS:-"5 10 15"}
for ns in ${=WINDOWS}; do
  if [[ $MODE == trailing ]]; then W=$M/last_${ns}ns; else W=$M/first_${ns}ns; fi
  [[ -f $W/mmgbsa/FINAL_RESULTS_MMPBSA.dat ]] && { echo "[$ns ns] already done"; continue; }
  mkdir -p $W
  cp $P/traj_wrapped.pdb $W/traj_wrapped.pdb
  $PY -c "
import mdtraj as md
t = md.load('$P/traj_wrapped.xtc', top='$P/traj_wrapped.pdb')
n = max(1, int(round(len(t)*$ns/$TOTAL)))
k = t[-n:] if '$MODE' == 'trailing' else t[:n]
k.save_xtc('$W/traj_wrapped.xtc')
print(f'[$ns ns] kept {len(k)}/{len(t)} frames ($MODE)')
"
  ( cd $W && $OMD mmgbsa --protein $ROOT/$M/protein_fixed.pdb --ligand $ROOT/$M/ligand_prepped.sdf \
                         --traj traj_wrapped.xtc --topology traj_wrapped.pdb \
                         --out-dir mmgbsa --no-auto-cofactors --run >/dev/null 2>&1 )
  # scatter, per CLAUDE.md -- traj.nc belongs in this list and was missing from it until 2026-10-21
  # (89 of them held 6.0 GB across the two campaigns before the same omission was fixed elsewhere);
  # find -delete rather than a multi-path rm, whose zsh globs abort the whole command on a miss.
  find $W \( -name 'reference.frc' -o -name '_MMPBSA_*' -o -name 'traj.nc' \) -delete 2>/dev/null
  echo -n "[$ns ns] "; grep "DELTA TOTAL" $W/mmgbsa/FINAL_RESULTS_MMPBSA.dat 2>/dev/null || echo "MMGBSA FAILED"
done
echo "WINDOW_MODAL_DONE[$N]"
