#!/bin/zsh
# Feed the four rep3 legfiles sequentially in BATCH=1 mode (drives racc_feed.sh; stages + submits).
cd /Users/cafierom/python_mac/peptidebuilder
set -u
run() {
  local lf=$1 mol=$2 time=$3
  echo "[$(date '+%H:%M:%S')] feeding $lf (SYSNAME=$mol, --time=$time)"
  BATCH=1 BATCH_NAME=$lf SYSNAME=$mol SBATCH_OPTS="--time=$time --cpus-per-task=8" \
    zsh code/racc_feed.sh runs/folds_md/$lf.legs
  echo "[$(date '+%H:%M:%S')] feeder for $lf exited: $?"
}
run rep3_oxb_1 oxybenzone 07:12:00
run rep3_oxb_2 oxybenzone   04:17:00
run rep3_oct_1 octinoxate  06:58:00
run rep3_oct_2 octinoxate  03:48:00
echo "[$(date '+%H:%M:%S')] all legfiles handled"
