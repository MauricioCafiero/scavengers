#!/bin/zsh
# Recompute whole-run MM/GBSA after a warm-restart extension (md_warm_restart.py).
#
# Every ΔG in this project is a WHOLE-RUN average (prod_20ns/mmgbsa over all frames),
# so an extended run must be rescored the same way to stay comparable -- not as a
# trailing segment. This joins the base production trajectory and the extension into
# one trajectory of the new total length and runs the standard analyze + MM/GBSA over
# all of it, writing prod_<TOTAL>ns/ beside the original.
#
#   code/md_mmgbsa_whole_run.sh <structure> <total_ns> [base_leg=prod_20ns]
#
# Produced orig_f12_dock5's 25 ns whole-run point (2026-10-03):
#   code/md_mmgbsa_whole_run.sh orig_f12_dock5 25
set -u
REPO=${0:A:h:h}
cd $REPO
ROOT=$REPO
OMD=${OMD_ENV:-$HOME/miniforge3/envs/openmm-md}/bin/omd
PY=${OMD_ENV:-$HOME/miniforge3/envs/openmm-md}/bin/python

N=${1:?usage: md_mmgbsa_whole_run.sh <structure> <total_ns> [base_leg]}
TOTAL=${2:?need total_ns}
BASE=${3:-prod_20ns}
M=runs/octinoxate/md/$N
P=$M/$BASE
W=$M/prod_${TOTAL}ns
mkdir -p $W
echo "### whole-run MM/GBSA over $TOTAL ns for $N (base $BASE + extension)"; date

if [[ ! -f $W/traj.dcd ]]; then
    $PY -c "
import mdtraj as md
t1 = md.load('$P/traj.dcd', top='$M/system/complex.pdb')
t2 = md.load('$P/traj_ext.dcd', top='$M/system/complex.pdb')
full = t1.join(t2)
full.save_dcd('$W/traj.dcd')
print(f'joined {len(t1)} + {len(t2)} = {len(full)} frames -> $W/traj.dcd')
"
fi
[[ -f $W/traj_wrapped.xtc ]] || $OMD analyze --traj $W/traj.dcd --topology $M/system/complex.pdb --out-dir $W
if [[ ! -f $W/mmgbsa/FINAL_RESULTS_MMPBSA.dat ]]; then
    # from inside $W: MMPBSA.py scatters reference.frc + _MMPBSA_* into cwd (CLAUDE.md)
    ( cd $W && $OMD mmgbsa --protein $ROOT/$M/protein_fixed.pdb --ligand $ROOT/$M/ligand_prepped.sdf \
                           --traj traj_wrapped.xtc --topology traj_wrapped.pdb \
                           --out-dir mmgbsa --no-auto-cofactors --run )
fi
echo "$TOTAL ns whole-run DELTA:"; grep "DELTA TOTAL" $W/mmgbsa/FINAL_RESULTS_MMPBSA.dat
echo "WHOLE_RUN_MMGBSA_DONE[$N $TOTAL ns]"; date
