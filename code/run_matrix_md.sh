#!/bin/zsh
# The two-scorer pose matrix (NEXT_STEPS item 5): five 20 ns legs on Modal, MM/GBSA local.
#
# Sequential by design, one cell at a time, so a problem is visible before it has been paid for four
# more times. Any failing stage stops the chain rather than letting it cascade -- a bad system or a
# Modal error should not quietly consume the budget.
#
# Measured before launching rather than interpolated: a probe on the largest system (s3_esm2_f4_dock1,
# 14,105 particles) read 725 ns/day on A10G, 0.66 h and $0.73 for 20 ns. The five together project to
# about $3.10 against a $5 budget. Interpolating the docstring's 21k-atom rate had predicted $4.71, so
# the probe was worth its few cents -- the A10G scales better than linearly at these sizes.
#
# MM/GBSA runs locally because only dynamics needs the GPU, and because keeping it here leaves the
# convergence windows available after the fact. It runs from inside the leg directory: MMPBSA.py
# scatters scratch into the working directory and reference.frc alone reached 1.9 GB when this was run
# from the repository root.
#
# Every stage is guarded by its own output, so re-running resumes rather than restarts.
set -u
# Paths are derived, not hardcoded, so the repo can be cloned anywhere. Override the externals with
# OMD_ENV and BOLTZGEN_ROOT if they live elsewhere on this machine.
REPO=${0:A:h:h}
cd $REPO
OMD_ENV=${OMD_ENV:-$HOME/miniforge3/envs/openmm-md}
BOLTZGEN_ROOT=${BOLTZGEN_ROOT:-$BOLTZGEN_ROOT}
for need in $OMD_ENV/bin/omd $OMD_ENV/bin/python; do
    [[ -x $need ]] || { echo "missing $need -- set OMD_ENV to the openmm-md conda env"; exit 1; }
done
ROOT=$REPO
OMD=$OMD_ENV/bin/omd
STEPS=10000000            # 20 ns at 2 fs, same as every run it is compared with

caffeinate -w $$ &        # dies with this script
date
echo "### pose matrix: five legs, sequential"

# in order of information gained, not by structure (see NEXT_STEPS item 5)
CELLS=(s3_esm2_f4_dock1 bg33_4_dock1 s3_esm2_f4_dock6 s3_orig_f12_dock9 bg33_4_dock6)

for N in $CELLS; do
    M=runs/octinoxate/md/$N
    P=$M/prod_L1_modal
    echo
    echo "=== $N ==="
    date

    if [[ -f $P/mmgbsa/FINAL_RESULTS_MMPBSA.dat ]]; then
        echo "  done already, skipping"; continue
    fi

    if [[ ! -f $P/traj_wrapped.xtc ]]; then
        modal run code/modal_md.py::push --structure $N || { echo "PUSH_FAILED $N"; exit 1; }
        modal run code/modal_md.py::produce --structure $N --steps $STEPS \
            || { echo "PRODUCE_FAILED $N"; exit 1; }
        date
    else
        echo "  trajectory already here, skipping the GPU leg"
    fi

    [[ -f $P/traj_wrapped.xtc ]] || { echo "NO_TRAJECTORY $N"; exit 1; }

    ( cd $P && $OMD mmgbsa --protein $ROOT/$M/protein_fixed.pdb \
                           --ligand $ROOT/$M/ligand_prepped.sdf \
                           --traj traj_wrapped.xtc --topology traj_wrapped.pdb \
                           --out-dir mmgbsa --no-auto-cofactors --run ) \
        || { echo "MMGBSA_FAILED $N"; exit 1; }
    date

    echo "--- $N result ---"
    awk '/DELTAS|Differences/{d=1} d&&/DELTA TOTAL/{print; exit}' $P/mmgbsa/FINAL_RESULTS_MMPBSA.dat
done

echo
echo "MATRIX_DONE"
date
