#!/bin/zsh
# Submit the remaining replicate batches through racc_feed.sh BATCH=1, in order, one per call,
# waiting for free slots -- the tried-and-true queue filler (user instruction 2026-10-23 after the
# racc_chain.sh experiment left the queue empty ~2 h). Each batch's <name>.legs, SYSNAME and
# --time come from runs/folds_md/rep5.chain + rep6.chain. Exits when every batch is submitted.
set -u
REPO=${0:A:h:h:h}
cd $REPO

caffeinate -w $$ &               # dies with this script
log() { print -r -- "[$(date '+%Y-%m-%d %H:%M:%S')] $*"; }

typeset -a NAMES MOLS TIMES
for lf in runs/folds_md/rep5.chain runs/folds_md/rep6.chain; do
  while IFS= read -r line; do
    [[ -z ${line// /} || ${line[1]} == '#' ]] && continue
    p=(${=line}); NAMES+=($p[1]); MOLS+=($p[2]); TIMES+=($p[3])
  done < $lf
done

i=1
while (( i <= ${#NAMES} )); do
  name=${NAMES[$i]}; mol=${MOLS[$i]}; t=${TIMES[$i]}
  # rep5_oxy_1 is already in the queue as job 537778 with its legs baked in -- skip it, resuming
  # is the rest of the list.
  if [[ $name == rep5_oxy_1 ]]; then (( i++ )); continue; fi
  log "feed $name ($mol, $t)"
  BATCH=1 BATCH_NAME=$name SYSNAME=$mol SBATCH_OPTS="--time=$t --cpus-per-task=8" \
    zsh code/racc_feed.sh runs/folds_md/$name.legs
  (( i++ ))
  sleep 20
done
log "all remaining batches submitted"