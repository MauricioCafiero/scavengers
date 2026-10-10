#!/bin/zsh
# BUILD_ONLY pass over the remaining fold-campaign cofold legs, so each Modal launch is
# push-only when its turn comes (leg 1 is producing while this runs; builds are local).
# orig_f12_o3cof is already built and is skipped by the driver's system.xml guard anyway.
set -u
REPO=${0:A:h:h}
cd $REPO
caffeinate -w $$ &
for spec in \
    "octinoxate orig_f12_o3cof" \
    "octinoxate s3_orig_f12_o3cof" \
    "octinoxate s2_esm2_control_o3cof" \
    "octinoxate s3_esm2_f4_o3cof" \
    "oxybenzone ox1_orig_f12_o3cof" \
    "oxybenzone ox2_orig_f8_o3cof" \
    "oxybenzone ox3_orig_f8_o3cof" \
    "oxybenzone ox1_esm1_f8_o3cof" \
    "oxybenzone ox2_esm1_f8_o3cof"; do
    SYSNAME=${spec%% *}; STRUCT=${spec##* }
    echo "=== build $SYSNAME/$STRUCT ($(date '+%H:%M:%S'))"
    STRUCT=$STRUCT POSE=cofold SYSNAME=$SYSNAME code/run_dock_pose_md.sh BUILD_ONLY=1
    echo "--- exit $? ($(date '+%H:%M:%S'))"
done
echo "ALL COFOLD BUILDS DONE"