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
REPO=${0:A:h:h}
cd $REPO
OMD_ENV=${OMD_ENV:-$HOME/miniforge3/envs/openmm-md}   # override if the conda env is elsewhere
OMD=$OMD_ENV/bin/omd
PY=$OMD_ENV/bin/python
ROOT=$REPO
# Parameters, env-overridable. Defaults reproduce the original s3_orig_f12 pose-1 behaviour, so an
# unparameterised call is unchanged. STRUCT is the cofolded structure whose docked pose this is; POSE
# is the 1-based pose rank dock_pose_to_sdf.py pulls from runs/<sys>/dock/<STRUCT>/. BOX defaults to
# dodecahedron now -- the co-folded baselines are dodecahedral and omd build defaults to cube, so an
# explicit shape is the one thing that must not be left to the default (CLAUDE.md). Box shape does not
# affect MM/GBSA, which strips to the solute, so a dodecahedral docked leg stays comparable to a
# cubic cofolded one.
STRUCT=${STRUCT:-s3_orig_f12}
POSE=${POSE:-1}
SYSNAME=${SYSNAME:-octinoxate}
BOX=${BOX:-dodecahedron}
BUILD_ONLY=${BUILD_ONLY:-0}   # stop after the build, to read the particle count before the 20 ns
LEGNAME=${LEGNAME:-}   # override N outright: a replicate leg (<leg>_r2) is its own staged dir with a
                       # byte-identical system/ copied from the original, so the build guard skips
                       # straight to the run; STRUCT/POSE stay set for the tail's SDF naming. See
                       # NEXT_STEPS.md, "Replicate every MM/GBSA leg".
# Honor the documented invocation `run_dock_pose_md.sh BUILD_ONLY=1` as a positional argument too:
# only reading the env var let a bare positional word through and the "build-only" launch ran a full
# 20 ns (2026-10-04, two accidental runs). `set -u` is set; loop over "$@" explicitly.
for a in "$@"; do [[ $a == BUILD_ONLY=1 ]] && BUILD_ONLY=1; done
N=${LEGNAME:-${STRUCT}_dock${POSE}}
M=runs/${SYSNAME}/md/$N                   # SYSNAME parameterizes the run dir (was hardcoded octinoxate)
SRC=${SRC:-runs/${SYSNAME}/md/$STRUCT}   # the cofolded leg: supplies the already-prepped, frame-matched receptor
                                         # boltzgen cofolds live under boltzgen_local/md/<STRUCT>; pass
                                         # SRC=/Users/cafierom/python_mac/boltzgen_local/md/<STRUCT>
P=$M/prod_20ns
STEPS=${STEPS:-10000000}         # 20 ns at 2 fs

caffeinate -w $$ &        # dies with this script; holds off idle sleep for the whole run
date
echo "### $N: 20 ns on pose $POSE of $STRUCT (box=$BOX)"
mkdir -p $M

# Build stage. Docked systems are not prebuilt, unlike the cofolded ones. Guarded by system.xml so a
# resume skips it. The receptor is copied from the cofolded leg rather than re-prepped: it is the same
# peptide, in the same frame the poses were docked into, so re-prepping could only drift it.
if [[ ! -f $M/system/system.xml ]]; then
    echo "--- building $N ($BOX) ---"
    [[ -f $M/protein_fixed.pdb ]] || cp $SRC/protein_fixed.pdb $M/protein_fixed.pdb
    [[ -f $M/${N}_ligand.sdf ]] || \
        $PY $ROOT/code/dock_pose_to_sdf.py $STRUCT --pose $POSE --system $SYSNAME \
            --out $M/${N}_ligand.sdf
    [[ -f $M/ligand_prepped.sdf ]] || \
        $OMD prep-ligand --sdf $M/${N}_ligand.sdf --out $M/ligand_prepped.sdf
    $OMD build --protein $M/protein_fixed.pdb --ligand $M/ligand_prepped.sdf \
               --out-dir $M/system --no-auto-cofactors --box-shape $BOX
    echo "BUILD_EXIT=$?"
    grep -c "^ATOM\|^HETATM" $M/system/complex.pdb | xargs echo "particles:"
    grep -m1 "^CRYST1" $M/system/complex.pdb
    date
fi
[[ $BUILD_ONLY == 1 ]] && { echo "BUILD_ONLY set, stopping after build"; exit 0; }

# Skip the run if EITHER the full dcd (local/old legs) or the wrapped solute (racc legs stripped on
# the cluster, whose full dcd never comes to this Mac) is already here -- both mean dynamics is done.
# Without the traj_wrapped.xtc clause a stripped leg would launch a 20 ns run locally on the laptop.
if [[ ! -f $P/traj.dcd && ! -f $P/traj_wrapped.xtc ]]; then
    $OMD run --system $M/system/system.xml --topology $M/system/complex.pdb \
             --out-dir $P --steps $STEPS --platform OpenCL
    echo "RUN_EXIT=$?"
    date
else
    echo "trajectory present in $P (traj.dcd or traj_wrapped.xtc), skipping the dynamics"
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
    # Two sources for the window's wrapped solute, depending on the lane. Local/Modal legs have the
    # full-system $P/traj.dcd here, so slice it and wrap as before. Racc legs stripped on the cluster
    # (racc_run.py) fetched only $P/traj_wrapped.xtc -- the full box never comes to this 8 GB Mac -- so
    # slice the leading fraction of the already-wrapped solute directly, no full-box load, no OOM.
    if [[ ! -f $W/traj_wrapped.xtc ]]; then
        if [[ -f $P/traj.dcd ]]; then
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
        else
            $PY -c "
import mdtraj as md
t = md.load('$P/traj_wrapped.xtc', top='$P/traj_wrapped.pdb')
keep = t[: int(len(t) * $ns / 20)]
keep.save_xtc('$W/traj_wrapped.xtc')
keep[0].save_pdb('$W/traj_wrapped.pdb')
print(f'kept {len(keep)} of {len(t)} frames (leading $ns ns, wrapped solute)')
"
        fi
    fi
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
