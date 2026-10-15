#!/bin/zsh
# BUILD_ONLY pass over the remaining fold-campaign cofold legs, so each Modal launch is
# push-only when its turn comes (leg 1 is producing while this runs; builds are local).
# orig_f12_o3cof is already built and is skipped by the driver's system.xml guard anyway.
set -u
REPO=${0:A:h:h}
cd $REPO
caffeinate -w $$ &
# A leg-list file ("<mol> <leg>" per line, # comments) overrides the built-in list -- used for the
# remaining cofold legs (8 OF3 BoltzGen + 21 RF3), runs/folds_md/cofold_legs.txt.
if (( $# )); then
    specs=("${(@f)$(grep -v '^#' $1 | grep -v '^$')}")
else
    specs=("octinoxate orig_f12_o3cof" "octinoxate s3_orig_f12_o3cof" "octinoxate s2_esm2_control_o3cof"
           "octinoxate s3_esm2_f4_o3cof" "oxybenzone ox1_orig_f12_o3cof" "oxybenzone ox2_orig_f8_o3cof"
           "oxybenzone ox3_orig_f8_o3cof" "oxybenzone ox1_esm1_f8_o3cof" "oxybenzone ox2_esm1_f8_o3cof")
fi
OMD=${OMD_ENV:-$HOME/miniforge3/envs/openmm-md}/bin/omd
for spec in $specs; do
    SYSNAME=${spec%% *}; STRUCT=${spec##* }
    echo "=== build $SYSNAME/$STRUCT ($(date '+%H:%M:%S'))"
    # A leg never staged by a campaign build needs its receptor prepped from the docking grid's
    # protonated protein (FOLDS.md, cofold staging chain); --out takes the FILE path.
    M=runs/$SYSNAME/md/$STRUCT
    if [[ ! -f $M/protein_fixed.pdb ]]; then
        mkdir -p $M
        $OMD prep-protein --pdb runs/$SYSNAME/dock/$STRUCT/${STRUCT}_protein.pdb --out $M/protein_fixed.pdb
    fi
    STRUCT=$STRUCT POSE=cofold SYSNAME=$SYSNAME SRC=$M code/run_dock_pose_md.sh BUILD_ONLY=1
    echo "--- exit $? ($(date '+%H:%M:%S'))"
done
echo "ALL COFOLD BUILDS DONE"