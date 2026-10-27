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

typeset -a NAMES MOLS TIMES CHAINS
CHAINS=("$@")                    # chain files to feed, in order; default the original rep5 + rep6
(( ${#CHAINS} )) || CHAINS=(runs/folds_md/rep5.chain runs/folds_md/rep6.chain)
for lf in $CHAINS; do
  while IFS= read -r line; do
    [[ -z ${line// /} || ${line[1]} == '#' ]] && continue
    p=(${=line}); NAMES+=($p[1]); MOLS+=($p[2]); TIMES+=($p[3])
  done < $lf
done

SOCK=$HOME/.ssh/cm-racc
rs() { ssh -o ControlPath=$SOCK -o ConnectTimeout=15 -o BatchMode=yes racc.rdg.ac.uk "$@"; }
mins_of() {  # "[-]d-hh:mm:ss" or "hh:mm:ss" -> minutes
  local el=$1 tp=${1##*-} mins=0
  [[ $el == *-* ]] && mins=$(( ${el%%-*} * 1440 ))
  print $(( mins + ${tp%%:*} * 60 + ${${tp#*:}%%:*} ))
}

i=1
while (( i <= ${#NAMES} )); do
  name=${NAMES[$i]}; mol=${MOLS[$i]}; t=${TIMES[$i]}
  # Skip batches whose job is queued/running, or that already genuinely ran (>= 10 min sacct
  # elapsed). A 0-1 s COMPLETED/FAILED record is a resubmit: it is the no-op signature. A failure
  # that ran >= 10 min is NOT touched -- it needs diagnosis (a leg that NaN'd mid-run would be
  # relaunched into the same failure forever).
  sq=$(rs "squeue -u \$USER -h -n ${name}_run" 2>/dev/null)
  if [[ -n ${sq//[[:space:]]/} ]]; then
    log "$name queued/running as $(rs "squeue -u \$USER -h -n ${name}_run -o %i" 2>/dev/null | tr -s ' ') -- skip"
    (( i++ )); continue
  fi
  rec=$(rs "sacct --name=${name}_run -X -n -o State,Elapsed" 2>/dev/null | head -1 | tr -s ' \t')
  rec=${rec##[ ]#}
  st=${rec%% *}
  if [[ $st == COMPLETED ]]; then
    (( $(mins_of ${rec##* }) >= 10 )) && { log "$name completed after ${rec##* } -- skip"; (( i++ )); continue; }
    log "$name sacct says COMPLETED after ${rec##* } -- resubmitting (no-op signature)"
  elif [[ $st == FAILED || $st == PREEMPT* || $st == CANCEL* ]]; then
    if (( $(mins_of ${rec##* }) >= 10 )); then
      log "$name $st after ${rec##* } -- NOT resubmitting; needs diagnosis"
      (( i++ )); continue
    fi
    log "$name ended instantly as $st -- resubmitting"
  elif [[ -n $st ]]; then
    log "$name state: $st -- skip"
    (( i++ )); continue
  fi
  log "feed $name ($mol, $t)"
  BATCH=1 BATCH_NAME=$name SYSNAME=$mol SBATCH_OPTS="--time=$t --cpus-per-task=8" \
    zsh code/racc_feed.sh runs/folds_md/$name.legs
  (( i++ ))
  sleep 20
done
log "all remaining batches submitted"