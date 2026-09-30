#!/bin/zsh
# 20 ns on the DOCKED pose of s3_orig_f12, to compare against the predicted pose already run.
#
# s3_orig_f12 is the best binder measured here: MM/GBSA -24.33 +- 0.08 from the Boltz pose, 100%
# residence, zero releases in 2000 frames. AutoDock Vina's best-scoring pose for the same fold sits
# 4.09 A away -- a 2.57 A translation and a 167 degree rigid-body rotation, though only 64 degrees of
# head/tail change, so the ligand is turned but not end for end. The static metrics cannot say which
# pose is better: the docked one is equally wrapped (1.00) and slightly less enclosed (0.92 against
# 0.965), and Vina selected it for packing so its packing is not evidence. Dynamics can.
#
# Same protocol as the predicted-pose run in every respect that could confound the comparison: same
# protein_fixed.pdb (byte-identical receptor), same omd code, same 2 fs timestep, same dodecahedral
# box, same 20 ns, same MM/GBSA windows. The one difference is the platform -- this runs locally on
# OpenCL while s3_orig_f12's prod_L1_modal ran on CUDA -- which README's dynamics section covers:
# code/modal_md.py imports openmm_md.dynamics.run rather than reimplementing it, so the integrator and
# equilibration are the same code on both.
#
# 5,615 particles, measured at 1.44 ms/step on OpenCL: 12 min/ns, 120 ns/day, 4.0 h for the 20 ns.
# Do not estimate this from another system's rate. bg33_3 at 5,088 particles ran 1.11 ms/step and
# bg33_4 at 7,231 ran 1.71, and interpolating those for 5,615 predicts 1.26 -- 13% optimistic against
# the 1.44 actually measured here. The rate is not linear in particle count, so measure it from
# energy.csv once the run is going rather than borrowing a number from a different box.
#
# Every stage is guarded by its own output, so re-running resumes rather than restarts. `omd run`
# itself cannot resume -- dynamics.py writes a checkpoint but nothing loads it -- so an interrupted
# trajectory restarts from zero, and that is the one stage worth not interrupting.
set -u
cd /Users/cafierom/python_mac/peptidebuilder
OMD=~/miniforge3/envs/openmm-md/bin/omd
PY=~/miniforge3/envs/openmm-md/bin/python
ROOT=$PWD
N=s3_orig_f12_dock1
M=runs/octinoxate/md/$N
P=$M/prod_20ns
STEPS=10000000            # 20 ns at 2 fs

caffeinate -w $$ &        # dies with this script; holds off idle sleep for the whole run
date
echo "### $N: 20 ns on the docked pose"

if [[ ! -f $P/traj.dcd ]]; then
    $OMD run --system $M/system/system.xml --topology $M/system/complex.pdb \
             --out-dir $P --steps $STEPS --platform OpenCL
    echo "RUN_EXIT=$?"
    date
else
    echo "$P/traj.dcd exists, skipping the dynamics"
fi

if [[ ! -f $P/traj_wrapped.xtc ]]; then
    $OMD analyze --traj $P/traj.dcd --topology $M/system/complex.pdb --out-dir $P
    echo "ANALYZE_EXIT=$?"
fi

if [[ ! -f $P/mmgbsa/FINAL_RESULTS_MMPBSA.dat ]]; then
    # from inside the output dir: MMPBSA.py scatters scratch into the working directory, and one of
    # its files (reference.frc) reached 1.9 GB when this was run from the repository root
    ( cd $P && $OMD mmgbsa --protein $ROOT/$M/protein_fixed.pdb \
                           --ligand $ROOT/$M/ligand_prepped.sdf \
                           --traj traj_wrapped.xtc --topology traj_wrapped.pdb \
                           --out-dir mmgbsa --no-auto-cofactors --run )
    echo "MMGBSA_EXIT=$?"
    date
fi

# Leading windows, because every short window on orig_f12 read too negative -- it was still measuring
# the predicted pose rather than the ensemble. Here the question is sharper: if the docked pose is the
# wrong one, the series should walk away from its first window the way orig_f12's did.
for ns in 5 10 15; do
    W=$M/first_${ns}ns
    [[ -f $W/mmgbsa/FINAL_RESULTS_MMPBSA.dat ]] && continue
    mkdir -p $W
    if [[ ! -f $W/traj.dcd ]]; then
        $PY -c "
import mdtraj as md
t = md.load('$P/traj.dcd', top='$M/system/complex.pdb')
keep = t[: int(len(t) * $ns / 20)]
keep.save_dcd('$W/traj.dcd')
print(f'kept {len(keep)} of {len(t)} frames (leading $ns ns)')
"
    fi
    $OMD analyze --traj $W/traj.dcd --topology $M/system/complex.pdb --out-dir $W
    ( cd $W && $OMD mmgbsa --protein $ROOT/$M/protein_fixed.pdb \
                           --ligand $ROOT/$M/ligand_prepped.sdf \
                           --traj traj_wrapped.xtc --topology traj_wrapped.pdb \
                           --out-dir mmgbsa --no-auto-cofactors --run )
    echo "WINDOW_EXIT[${ns}ns]=$?"
    date
done

# What the ligand did over the run: residence and release, which one dG cannot show, and the measure
# that separated every design from every null.
$PY $ROOT/code/md_contacts.py $P --csv $M/md_contacts.csv || true
$PY $ROOT/code/md_frames.py $P || true
echo "DONE $N"
date
