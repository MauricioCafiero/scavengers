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
SYSNAME=${SYSNAME:-oxybenzone}  # the molecule; one batch carries one molecule's legs
# Extra sbatch options, e.g. "--time=6:00:00 --cpus-per-task=8". submit_gpu.sh takes none, so they go
# in through a shadowing sbatch() function exported to it. The SBATCH_* env vars are no substitute:
# gpu.slurm's own --cpus-per-task=1 directive beats SBATCH_TRES_PER_TASK (tested 2026-10-14). 8 CPUs
# lifts the job's JOBSIZE priority; racc's AGE weight is 0, so waiting longer earns nothing.
SBATCH_OPTS=${SBATCH_OPTS:-}
# --time and --cpus-per-task are REQUIRED, not options to remember; the house settings are
# **estimate + 1 h 30 min**, and **8 cores**. Without --time a job takes the partition DEFAULTTIME
# (1-00:00:00 on gpuscavenger, against a 5-day cap), which is 3-5x any real batch here: a long
# request scores worse on JOBSIZE priority, cannot backfill into a short gap, and sits exposed to
# preemption for hours it will never use. Estimate a batch as the sum of
# 8 min + 39.3 min x (N/5838)^0.4 over its legs, N read off each build, then add 1.5 h.
# Five batches went out on 2026-10-18 with neither setting; the fix after the fact was
# `scontrol update JobId=<id> TimeLimit=<hh:mm:ss>` plus `NumCPUs=8`, and NumCPUs took only on the
# PENDING ones -- a running job keeps the cores it started with, so getting this right at submission
# is the only reliable route.
for required in --time --cpus-per-task; do
    if [[ $SBATCH_OPTS != *$required* ]]; then
        print -u2 "racc_feed.sh: SBATCH_OPTS must set $required. The house settings are"
        print -u2 "  SBATCH_OPTS=\"--time=<est+1h30> --cpus-per-task=8\" BATCH=1 SYSNAME=<mol> \\"
        print -u2 "      code/racc_feed.sh <legfile>"
        print -u2 "(--time: estimate + 1.5 h. Omitted, it is the partition's 24 h DEFAULTTIME --"
        print -u2 " worse JOBSIZE priority, no backfill into a short gap, hours of needless"
        print -u2 " preemption exposure. --cpus-per-task=8 lifts JOBSIZE priority; gpu.slurm's own"
        print -u2 " --cpus-per-task=1 directive beats SBATCH_TRES_PER_TASK, so it must come through"
        print -u2 " here. Each runs/folds_md/*.legs header carries that batch's own values.)"
        exit 2
    fi
done
SUBMIT="sbatch() { command sbatch $SBATCH_OPTS \"\$@\"; }; export -f sbatch; ~/bin/submit_gpu.sh"
WORKREL=runs/$SYSNAME/md
REMOTE=racc.rdg.ac.uk

caffeinate -w $$ &               # dies with this script

r() { ssh -n -o ControlPath=$SOCK -o ConnectTimeout=15 -o BatchMode=yes $REMOTE "$@"; }
# Staging pipes a tar stream INTO ssh, so it must not carry -n (which ties ssh's stdin to /dev/null
# and drops the stream -- "gzip: unexpected end of file", 2026-10-12).
rs() { ssh -o ControlPath=$SOCK -o ConnectTimeout=15 -o BatchMode=yes $REMOTE "$@"; }

log() { print -r -- "[$(date '+%Y-%m-%d %H:%M:%S')] $*"; }

# Read the list into parallel arrays, skipping comments and blanks.
typeset -a SRCS LEGS
while IFS= read -r line; do
    [[ -z ${line// /} || ${line[1]} == '#' ]] && continue
    SRCS+=("${line%%:*}")
    LEGS+=("${line##*:}")
done < $LEGFILE

log "feeder started: ${#LEGS} legs queued, ceiling $CEILING, interval ${INTERVAL}s, batch=$BATCH"
# Log the sbatch options, so what a job actually requested is recoverable from the log rather
# than inferred: macOS will not let you read a running process's environment, so after the fact
# there is otherwise no record of what SBATCH_OPTS held.
log "sbatch opts: $SBATCH_OPTS"

if [[ $BATCH == 1 ]]; then
    # Stage every leg first (idempotent -- the tar just overwrites), then wait for one slot.
    for i in {1..${#LEGS}}; do
        src=${SRCS[$i]}; leg=${LEGS[$i]}
        if [[ ! -f $src/system/system.xml ]]; then
            log "ABORT: $leg is not built ($src/system/system.xml missing)"
            exit 1
        fi
        log "staging $leg"
        tar czf - -C $src . | rs "mkdir -p ~/remote_work/openmm/$WORKREL/$leg && tar xzf - -C ~/remote_work/openmm/$WORKREL/$leg" \
            || { log "ABORT: staging failed for $leg"; exit 1; }
    done

    legcsv=${(j:,:)LEGS}
    log "waiting for a free slot to submit one batch job of ${#LEGS} legs: $legcsv"
    while :; do
        njobs=$(r 'squeue -u $USER -h | wc -l' 2>/dev/null | tr -d ' ')
        [[ -n $njobs ]] && (( njobs < CEILING )) && break
        sleep $INTERVAL
    done
    out=$(r "bash -lc 'module load anaconda >/dev/null 2>&1; source activate openmm-md; cd ~/remote_work/openmm; cp racc_run.py ${BATCH_NAME}_run.py; export LEGS=$legcsv RUNREL=$WORKREL; $SUBMIT ${BATCH_NAME}_run.py'" 2>&1 | tail -3)
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
    if ! tar czf - -C $src . | rs "mkdir -p ~/remote_work/openmm/$WORKREL/$leg && tar xzf - -C ~/remote_work/openmm/$WORKREL/$leg"; then
        log "STAGE FAILED $leg -- will retry in ${INTERVAL}s"
        sleep $INTERVAL
        continue
    fi

    log "submitting $leg"
    out=$(r "bash -lc 'module load anaconda >/dev/null 2>&1; source activate openmm-md; cd ~/remote_work/openmm; cp racc_run.py ${leg}_run.py; export LEG=$leg RUNREL=$WORKREL; $SUBMIT ${leg}_run.py'" 2>&1 | tail -3)
    log "  $(print -r -- $out | tr '\n' ' ')"
    i=$((i+1))
    sleep 20          # let squeue reflect the new job before counting again
done

log "feeder done: every leg in $LEGFILE has been submitted"
