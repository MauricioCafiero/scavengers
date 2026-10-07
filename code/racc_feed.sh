#!/bin/zsh
# Keep the racc queue topped up to the per-user ceiling from an explicit list of legs.
#
#   code/racc_feed.sh <legfile>
#
# <legfile> holds one leg per line, in the order to submit, as `<local dir>:<leg name>`; blank lines
# and `#` comments are skipped. A leg already present in squeue or already holding a finished
# energy.csv on the cluster is skipped, so re-running the feeder after a restart resumes rather than
# double-submitting.
#
# This stages and submits only. The completion side -- fetch the wrapped solute, then run the tail --
# is racc_drain.sh, which calls run_dock_pose_md.sh and run_md20.sh rather than reimplementing their
# stages, so the per-step knowledge stays in the tools that already hold it.
#
# RACC allows 3 submitted jobs per user and runs 1 at a time (QOSMaxSubmitJobPerUserLimit), so the
# whole list cannot be queued at once -- hence a feeder rather than one batch submission.
set -u

LEGFILE=${1:?usage: racc_feed.sh <legfile>}
# BATCH=1 submits every leg in the file as ONE job (racc_run.py's LEGS), waiting only for a single
# free slot. Use it for unattended stretches: gpuscavenger allows 3 submitted jobs and runs 1, so
# one-leg-per-slot needs the ssh socket re-touched every ~50 minutes, and the socket needs a password
# and TFA only the user can give -- a dropped network ends the campaign until they are back. A batch
# banks the whole remainder in one slot and needs nothing afterwards. BATCH_NAME sets the job name.
BATCH=${BATCH:-0}
BATCH_NAME=${BATCH_NAME:-oxyb_batch}
CEILING=${CEILING:-3}
INTERVAL=${INTERVAL:-300}        # seconds between queue checks
SOCK=$HOME/.ssh/cm-racc
WORKREL=runs/oxybenzone/md
REMOTE=racc.rdg.ac.uk

caffeinate -w $$ &               # dies with this script

r() { ssh -o ControlPath=$SOCK $REMOTE "$@"; }

log() { print -r -- "[$(date '+%Y-%m-%d %H:%M:%S')] $*"; }

# Read the list into parallel arrays, skipping comments and blanks.
typeset -a SRCS LEGS
while IFS= read -r line; do
    [[ -z ${line// /} || ${line[1]} == '#' ]] && continue
    SRCS+=("${line%%:*}")
    LEGS+=("${line##*:}")
done < $LEGFILE

log "feeder started: ${#LEGS} legs queued, ceiling $CEILING, interval ${INTERVAL}s, batch=$BATCH"

if [[ $BATCH == 1 ]]; then
    # Stage every leg first (idempotent -- the tar just overwrites), then wait for one slot.
    for i in {1..${#LEGS}}; do
        src=${SRCS[$i]}; leg=${LEGS[$i]}
        if [[ ! -f $src/system/system.xml ]]; then
            log "ABORT: $leg is not built ($src/system/system.xml missing)"
            exit 1
        fi
        log "staging $leg"
        tar czf - -C $src . | r "mkdir -p ~/remote_work/openmm/$WORKREL/$leg && tar xzf - -C ~/remote_work/openmm/$WORKREL/$leg" \
            || { log "ABORT: staging failed for $leg"; exit 1; }
    done

    legcsv=${(j:,:)LEGS}
    log "waiting for a free slot to submit one batch job of ${#LEGS} legs: $legcsv"
    while :; do
        njobs=$(r 'squeue -u $USER -h | wc -l' 2>/dev/null | tr -d ' ')
        [[ -n $njobs ]] && (( njobs < CEILING )) && break
        sleep $INTERVAL
    done
    out=$(r "bash -lc 'module load anaconda >/dev/null 2>&1; source activate openmm-md; cd ~/remote_work/openmm; cp racc_run.py ${BATCH_NAME}_run.py; export LEGS=$legcsv RUNREL=$WORKREL; ~/bin/submit_gpu.sh ${BATCH_NAME}_run.py'" 2>&1 | tail -3)
    log "batch submitted: $(print -r -- $out | tr '\n' ' ')"
    log "feeder done (batch mode): nothing further needs the ssh socket"
    exit 0
fi

i=1
while (( i <= ${#LEGS} )); do
    if ! r true 2>/dev/null; then
        log "ssh socket down -- waiting ${INTERVAL}s (the user reopens it; nothing is lost)"
        sleep $INTERVAL
        continue
    fi

    njobs=$(r 'squeue -u $USER -h | wc -l' 2>/dev/null | tr -d ' ')
    if [[ -z $njobs ]]; then
        log "could not read queue depth -- retrying in ${INTERVAL}s"
        sleep $INTERVAL
        continue
    fi

    if (( njobs >= CEILING )); then
        sleep $INTERVAL
        continue
    fi

    src=${SRCS[$i]}
    leg=${LEGS[$i]}

    if [[ ! -d $src ]]; then
        log "SKIP $leg: no local directory $src"
        i=$((i+1))
        continue
    fi
    if [[ ! -f $src/system/system.xml ]]; then
        log "SKIP $leg: not built ($src/system/system.xml missing) -- build it with BUILD_ONLY=1 first"
        i=$((i+1))
        continue
    fi

    log "staging $leg from $src ($(du -sh $src | cut -f1))"
    if ! tar czf - -C $src . | r "mkdir -p ~/remote_work/openmm/$WORKREL/$leg && tar xzf - -C ~/remote_work/openmm/$WORKREL/$leg"; then
        log "STAGE FAILED $leg -- will retry in ${INTERVAL}s"
        sleep $INTERVAL
        continue
    fi

    log "submitting $leg"
    out=$(r "bash -lc 'module load anaconda >/dev/null 2>&1; source activate openmm-md; cd ~/remote_work/openmm; cp racc_run.py ${leg}_run.py; export LEG=$leg RUNREL=$WORKREL; ~/bin/submit_gpu.sh ${leg}_run.py'" 2>&1 | tail -3)
    log "  $(print -r -- $out | tr '\n' ' ')"
    i=$((i+1))
    sleep 20          # let squeue reflect the new job before counting again
done

log "feeder done: every leg in $LEGFILE has been submitted"
