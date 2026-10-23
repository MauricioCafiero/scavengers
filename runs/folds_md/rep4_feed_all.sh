#!/bin/zsh
# Feed the two rep4 legfiles sequentially in BATCH=1 (stages + submits; banks behind rep3_oct_2).
cd /Users/cafierom/python_mac/peptidebuilder
set -u
run() {
  local lf=$1 mol=$2 time=$3
  echo "[$(date '+%H:%M:%S')] feeding $lf (SYSNAME=$mol, --time=$time)"
  BATCH=1 BATCH_NAME=$lf SYSNAME=$mol SBATCH_OPTS="--time=$time --cpus-per-task=8" \
    zsh code/racc_feed.sh runs/folds_md/$lf.legs
  echo "[$(date '+%H:%M:%S')] feeder for $lf exited"
}
run rep4_oxb_1 oxybenzone 05:31:00
run rep4_oct_1 octinoxate 06:27:00
echo "[$(date '+%H:%M:%S')] rep4 legfiles handled"
