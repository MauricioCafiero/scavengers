#!/bin/zsh
# Feed the remaining fold-cofold batches to racc in order, one racc_feed.sh BATCH=1 call each.
# racc_feed.sh waits for a free slot itself (3 submitted per user), so running them in sequence keeps
# the queue topped up as jobs finish. Each submitted batch gets its own detached racc_drain.sh.
# Time limits: 1.3 x (8 min + 39.3 min x (N/5838)^0.4 per leg), rounded up to 30 min -- the rule that
# gave oct_1 5:30 and oxb_1 6:00.
# ${0:A:h:h} -- the repo root, two levels up from code/<this>. It was :h:h:h until 2026-10-20,
# which is correct for a script at runs/folds_md/feed_chain.sh (the old handoff path it was written
# at) but lands in ~/python_mac from code/. The chain then cannot create runs/folds_md/feed_<b>.log,
# the >> redirect fails so racc_feed.sh never runs, and the submission check correctly reports
# "did not submit" and stops. That silently killed the two GNINA-pick batches at 20:35 on 10-19.
cd ${0:A:h:h}
caffeinate -w $$ &
# LEGDIR/BATCH_PREFIX so a campaign outside runs/folds_md can use this chain instead of forking it:
# the replicate campaign keeps its legfiles in runs/replicates and names its jobs rep2_* rather
# than fold_*. Defaults reproduce the fold campaign's behaviour exactly.
D=${LEGDIR:-runs/folds_md}
PREFIX=${BATCH_PREFIX:-fold}
# Specs may be given on the command line to resume a partial chain -- a transient
# "Unable to contact slurm controller" killed the 2026-10-14 run at oxb_3 with two batches left,
# and re-running the whole list would resubmit batches that already ran. No args = the full list.
SPECS=("$@")
if (( ${#SPECS} == 0 )); then
    SPECS=(oct_2:octinoxate:5:00:00 oxb_2:oxybenzone:5:00:00 oct_3:octinoxate:6:00:00 \
           oxb_3:oxybenzone:6:00:00 oct_4:octinoxate:2:30:00)
fi
for spec in $SPECS; do
    b=${spec%%:*}; rest=${spec#*:}; mol=${rest%%:*}; t=${rest#*:}
    print -r -- "[$(date '+%F %T')] chain: feeding $b ($mol, --time=$t)"
    # Check only THIS attempt's output: the log is appended to, so a retry of a batch whose earlier
    # attempt succeeded would otherwise match that old line and skip a genuine failure.
    off=0; [[ -f $D/feed_$b.log ]] && off=$(wc -c < $D/feed_$b.log)
    BATCH=1 BATCH_NAME=${PREFIX}_$b SYSNAME=$mol SBATCH_OPTS="--time=$t --cpus-per-task=8" \
        zsh code/racc_feed.sh $D/$b.legs >> $D/feed_$b.log 2>&1
    if ! tail -c +$((off+1)) $D/feed_$b.log | grep -q 'Submitted batch job'; then
        print -r -- "[$(date '+%F %T')] chain: $b did not submit -- stopping (see $D/feed_$b.log)"
        exit 1
    fi
    ~/miniforge3/bin/python -c "import subprocess,os; subprocess.Popen(['zsh','code/racc_drain.sh','$D/$b.legs'], stdout=open('$D/drain_$b.log','ab'), stderr=subprocess.STDOUT, stdin=subprocess.DEVNULL, start_new_session=True, env={**os.environ,'SYSNAME':'$mol'})"
    print -r -- "[$(date '+%F %T')] chain: $b submitted, drain launched"
done
print -r -- "[$(date '+%F %T')] chain: all batches submitted"
