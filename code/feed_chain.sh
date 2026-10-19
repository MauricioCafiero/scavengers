#!/bin/zsh
# Submit a sequence of racc batches, one after another, and start each one's drain as it goes.
#
#   code/feed_chain.sh <legfile base>:<molecule>:<--time> ...
#   code/feed_chain.sh o3pc_oct_1:octinoxate:06:35:00 o3pc_oxb_1:oxybenzone:05:30:00
#
# Each argument names a runs/folds_md/<base>.legs file, the molecule whose runs/<mol>/md the legs
# live under, and that batch's --time (estimate + 1 h 30 min, per the legfile's own header).
# Colons in the time are fine -- the time is everything after the second field.
#
# This orchestrates and reimplements nothing: every batch goes out through code/racc_feed.sh in
# BATCH=1 mode (one job per legfile, racc_run.py's LEGS) and is drained by code/racc_drain.sh. The
# chain exists because gpuscavenger is MaxSubmitJobsPU=3 / MaxJobsPU=1, so a campaign of more than
# three batches cannot simply be queued at once; racc_feed.sh already blocks until a slot is free,
# so calling it in sequence is all the throttling needed.
#
# A version of this lived in runs/folds_md/ for the dk* campaign and was never committed, so it was
# gone when the pocket-guided lane needed it (2026-10-19). Hence this file, in code/ and tracked.
#
# Launch it detached, or a harness reap takes the whole chain with it:
#   ~/miniforge3/bin/python -c "import subprocess; subprocess.Popen(['zsh','code/feed_chain.sh',
#     'o3pc_oct_1:octinoxate:06:35:00', ...], stdout=open('runs/folds_md/feed_chain_o3pc.log','ab'),
#     stderr=subprocess.STDOUT, stdin=subprocess.DEVNULL, start_new_session=True, cwd='$PWD')"
set -u
REPO=${0:A:h:h}
cd $REPO
caffeinate -w $$ &                 # dies with this script

log() { print -r -- "[$(date '+%Y-%m-%d %H:%M:%S')] feed_chain: $*"; }

for spec in "$@"; do
    base=${spec%%:*}; rest=${spec#*:}
    mol=${rest%%:*}; tlimit=${rest#*:}
    legfile=runs/folds_md/$base.legs
    if [[ ! -f $legfile ]]; then log "no such legfile: $legfile -- skipping"; continue; fi

    log "submitting $base ($mol, --time=$tlimit)"
    # Blocks until a queue slot is free, submits the whole legfile as one job, then returns.
    BATCH=1 BATCH_NAME=$base SYSNAME=$mol \
        SBATCH_OPTS="--time=$tlimit --cpus-per-task=8" \
        zsh $REPO/code/racc_feed.sh $legfile
    rc=$?
    log "$base feeder exit $rc"
    if (( rc != 0 )); then log "$base NOT submitted -- not starting its drain"; continue; fi

    # One drain per legfile, detached, so the chain is never blocked on a tail. Each drain holds
    # the global tail lock for its own tails, so concurrent drains cannot OOM this 8 GB Mac.
    SYSNAME=$mol nohup zsh $REPO/code/racc_drain.sh $legfile \
        >> runs/folds_md/drain_$base.log 2>&1 &
    disown
    log "$base drain started (pid $!)"
done
log "chain done: every batch submitted"
