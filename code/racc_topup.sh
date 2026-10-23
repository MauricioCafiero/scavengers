#!/bin/zsh
# Keep the racc queue topped up to 3 jobs (1 running + 2 pending) from already-generated batches.
#
#   code/racc_topup.sh runs/folds_md/rep5.chain runs/folds_md/rep6.chain
#
# Standing policy (user ruling 2026-10-23, "it should ALWAYS queue a third job, like we have been
# doing since the start of this project" -- racc_feed.sh's CEILING=3): the second pending job is
# the standby that survives a preemption, since a preempted job's slurm hook line never fires and
# a one-pending chain would lose its buffer exactly then. This feeder sbatches each batch's
# pre-generated <name>_run.slurm (staged + hooked by racc_chain.sh) in list order whenever the
# queue drops below CEILING, and skips any batch that is queued, running or already completed, so
# it coexists with the slurm hooks without double-submitting: hooks are squeue-name-guarded and
# this feeder's skip does the rest. Chain file lines, same format racc_chain.sh reads:
#     <batch name> <SYSNAME-unused> <time>
# A lost socket costs nothing but a retry INTERVAL -- the 1-pending fallback keeps the chain
# advancing regardless; this only restores the second-pending margin. Detach via the usual
# start_new_session Popen one-liner; exits once every listed batch has been submitted.

set -u
REPO=${0:A:h:h}
SOCK=$HOME/.ssh/cm-racc
REMOTE=racc.rdg.ac.uk
CEILING=${CEILING:-3}
INTERVAL=${INTERVAL:-600}

caffeinate -w $$ &               # dies with this script

rs() { ssh -o ControlPath=$SOCK -o ConnectTimeout=15 -o BatchMode=yes $REMOTE "$@"; }
log() { print -r -- "[$(date '+%Y-%m-%d %H:%M:%S')] $*"; }

typeset -a NAMES TIMES
for lf in "$@"; do
  while IFS= read -r line; do
    [[ -z ${line// /} || ${line[1]} == '#' ]] && continue
    p=(${=line}); NAMES+=($p[1]); TIMES+=($p[3])
  done < $lf
done
n=${#NAMES}
log "topup: $n batches, ceiling $CEILING: ${NAMES[*]}"

i=1
while (( i <= n )); do
  qcnt=$(rs 'squeue -u $USER -h | wc -l' 2>/dev/null | tr -d ' ')
  if [[ -z $qcnt ]]; then
    log "socket unavailable -- retry in $INTERVAL s"
    sleep $INTERVAL
    continue
  fi
  # Skip batches already queued/running (squeue) or finished (sacct, for batches that completed
  # between polls). A batch we have not yet submitted shows neither.
  seen=$(rs "squeue -u \$USER -h -n ${NAMES[$i]}_run 2>/dev/null; sacct --name=${NAMES[$i]}_run -X -n -o State 2>/dev/null | head -1" | tr -d ' ')
  if [[ -n $seen ]]; then
    log "${NAMES[$i]} already ${seen} -- skip"
    (( i++ ))
    continue
  fi
  if (( qcnt >= CEILING )); then
    sleep $INTERVAL
    continue
  fi
  out=$(rs "sbatch --time=${TIMES[$i]} --cpus-per-task=8 \$HOME/remote_work/openmm/${NAMES[$i]}_run.slurm" 2>&1 | tail -1)
  log "submitted ${NAMES[$i]} (${TIMES[$i]}): $out"
  (( i++ ))
  sleep 20
done
log "all $n batches submitted -- hooks carry the rest"