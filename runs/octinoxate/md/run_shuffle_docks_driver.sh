#!/bin/zsh
# Orchestration only -- every MD step lives in code/run_dock_pose_md.sh (unchanged).
# Runs the two missing docked legs of shuffle_control sequentially (one local GPU), each paired
# with code/run_windows_live.sh so the leading MM/GBSA windows compute while the dynamics runs.
# Leg 2's window script just polls until its leg's dynamics starts writing energy.csv.
cd /Users/cafierom/python_mac/peptidebuilder

STRUCT=shuffle_control POSE=1 zsh code/run_dock_pose_md.sh \
    > runs/octinoxate/md/shuffle_control_dock1_run.log 2>&1 &
M1=$!
code/run_windows_live.sh shuffle_control_dock1 \
    > runs/octinoxate/md/shuffle_control_dock1_windows.log 2>&1 &
W1=$!
wait $M1

STRUCT=shuffle_control POSE=9 zsh code/run_dock_pose_md.sh \
    > runs/octinoxate/md/shuffle_control_dock9_run.log 2>&1 &
M2=$!
code/run_windows_live.sh shuffle_control_dock9 \
    > runs/octinoxate/md/shuffle_control_dock9_windows.log 2>&1 &
W2=$!
wait $M2

echo "ALL_LEG_SCRIPTS_DONE $(date)"