#!/bin/zsh
# Build the five MD systems of the two-scorer pose matrix (NEXT_STEPS item 5).
#
# Each cell is one docked pose of a structure that already has a predicted-pose trajectory. The
# receptor is reused byte-identically from that run -- the comparison is between ligand placements, so
# any difference in the protein preparation would be indistinguishable from the result. s3_orig_f12's
# dock1 run confirmed the reuse is exact (same md5 as runs/octinoxate/md/s3_orig_f12/protein_fixed.pdb).
#
# bg33_4's receptor comes from boltzgen_local, where its predicted-pose run lives; everything else is
# prepared here so the protocol is this repository's throughout.
#
# Local and free: dock_pose_to_sdf is RDKit, prep-ligand and build are the OpenFF/GAFF2 stack on CPU.
# Only the production leg goes to Modal. Every stage is guarded by its own output, so re-running
# resumes rather than restarts.
# --box-shape dodecahedron is NOT optional: `omd build` defaults to cube, and every co-folded baseline
# these legs are compared against was built dodecahedral. The seven legs run on 2026-10-01 inherited the
# cube default -- same physics, but ~30% more waters and about $1.30 of avoidable GPU across them, worst
# on an elongated solute (s2_esm2_control: 33,057 atoms cubic against 23,142 dodecahedral).
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
PY=.venv/bin/python

caffeinate -w $$ &        # dies with this script
date

# name                 structure      pose  receptor source
CELLS=(
  "s3_esm2_f4_dock1    s3_esm2_f4     1     runs/octinoxate/md/s3_esm2_f4/protein_fixed.pdb"
  "bg33_4_dock1        bg33_4         1     $BOLTZGEN_ROOT/md/bg33_4/protein_fixed.pdb"
  "s3_esm2_f4_dock6    s3_esm2_f4     6     runs/octinoxate/md/s3_esm2_f4/protein_fixed.pdb"
  "s3_orig_f12_dock9   s3_orig_f12    9     runs/octinoxate/md/s3_orig_f12/protein_fixed.pdb"
  "bg33_4_dock6        bg33_4         6     $BOLTZGEN_ROOT/md/bg33_4/protein_fixed.pdb"
)

for cell in $CELLS; do
    set -- ${=cell}
    N=$1 STRUCT=$2 POSE=$3 RECEPTOR=$4
    M=runs/octinoxate/md/$N
    echo "=== $N  ($STRUCT pose $POSE) ==="

    if [[ -f $M/system/system.xml ]]; then
        echo "  system built, skipping"; continue
    fi
    mkdir -p $M
    [[ -f $M/protein_fixed.pdb ]] || cp $RECEPTOR $M/protein_fixed.pdb

    if [[ ! -f $M/${N}_ligand.sdf ]]; then
        $PY code/dock_pose_to_sdf.py $STRUCT --pose $POSE --out $M/${N}_ligand.sdf || exit 1
    fi
    if [[ ! -f $M/ligand_prepped.sdf ]]; then
        $OMD prep-ligand --sdf $M/${N}_ligand.sdf --out $M/ligand_prepped.sdf || exit 1
    fi
    $OMD build --protein $M/protein_fixed.pdb --ligand $M/ligand_prepped.sdf \
               --out-dir $M/system --no-auto-cofactors --box-shape dodecahedron || exit 1
    date
done

echo "PREP_DONE"
date
