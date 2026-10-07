#!/bin/zsh
# Watch racc for finished legs, fetch each one, and run its tail locally as it lands.
#
#   code/racc_drain.sh <legfile>
#
# <legfile> is the same file racc_feed.sh reads: one `<local dir>:<leg name>` per line. The two run
# side by side -- the feeder keeps the queue full, this drains completions -- so a leg's analysis
# starts as soon as it finishes instead of backlogging behind the whole campaign.
#
# This orchestrates; it does not reimplement. The tail is whichever existing script owns that leg:
#   <leg>_dock<N>   ->  code/run_dock_pose_md.sh   (STRUCT/POSE/SYSNAME/SRC; guards skip the run)
#   <leg>           ->  boltzgen_local/md/run_md20.sh <leg>
# Both are racc-aware: they skip the dynamics when only traj_wrapped.xtc is present, and slice their
# windows from the wrapped solute rather than loading a full box this 8 GB Mac cannot hold.
#
# Only the wrapped solute is fetched (~100 MB), never the full dcd (1.8-5 GB). The dcd stays on
# cluster scratch, which is also what the archive keeps.
set -u

LEGFILE=${1:?usage: racc_drain.sh <legfile>}
INTERVAL=${INTERVAL:-300}
SOCK=$HOME/.ssh/cm-racc
REMOTE=racc.rdg.ac.uk
DCD=/scratch5/gaussian/io927423/dcd
BG=$HOME/python_mac/boltzgen_local
REPO=${0:A:h:h}
SYSNAME=${SYSNAME:-oxybenzone}

caffeinate -w $$ &

# -n so a command never eats this script's stdin; ConnectTimeout/ServerAlive* so a *wedged* master --
# the socket file still there, the network underneath it gone -- fails in seconds instead of blocking.
# Measured 2026-10-07: without these, one poll hung 40 minutes against a dead master during a train
# journey, and the loop made no progress for that whole time rather than logging and retrying.
r() { ssh -n -o ControlPath=$SOCK -o ConnectTimeout=15 -o ServerAliveInterval=10 \
          -o ServerAliveCountMax=3 -o BatchMode=yes $REMOTE "$@"; }
log() { print -r -- "[$(date '+%Y-%m-%d %H:%M:%S')] $*"; }

typeset -a SRCS LEGS
while IFS= read -r line; do
    [[ -z ${line// /} || ${line[1]} == '#' ]] && continue
    SRCS+=("${line%%:*}")
    LEGS+=("${line##*:}")
done < $LEGFILE

log "drain started: watching ${#LEGS} legs, interval ${INTERVAL}s"

remaining=${#LEGS}
while (( remaining > 0 )); do
    remaining=0
    for i in {1..${#LEGS}}; do
        src=${SRCS[$i]}
        leg=${LEGS[$i]}
        [[ -f $src/.drained || -f $src/.drain_failed ]] && continue

        if ! r true 2>/dev/null; then
            log "ssh socket down -- holding (nothing is lost; the run continues on the cluster)"
            remaining=$((remaining+1))
            break
        fi

        # Completion is judged by the leg's OWN output on scratch, not by a job state: a leg run
        # inside a multi-leg batch job (racc_run.py's LEGS) has no `<leg>_run` job for sacct to find,
        # so a state query would wait on it forever. traj_wrapped.xtc is the signal either way --
        # racc_run.py writes it as the last thing it does for a leg.
        if r "test -f $DCD/${leg}_prod/traj_wrapped.xtc" 2>/dev/null; then
            state=COMPLETED
        else
            # Not finished. Only consult sacct to tell "still going / not yet started" from a leg
            # whose own single-leg job died, which is worth flagging rather than waiting on.
            state=$(r "sacct --name=${leg}_run --format=State%20 -n -X 2>/dev/null | tail -1" | tr -d ' ')
            [[ -z $state ]] && state=PENDING    # not submitted yet, or it lives inside a batch job
        fi

        case $state in
            COMPLETED)
                log "=== $leg COMPLETED -- fetching wrapped solute"
                mkdir -p $src/prod_20ns
                if r "cd $DCD/${leg}_prod 2>/dev/null && tar czf - traj_wrapped.xtc traj_wrapped.pdb energy.csv rmsd.csv analysis.png final.pdb 2>/dev/null" | tar xzf - -C $src/prod_20ns 2>/dev/null; then
                    if [[ ! -f $src/prod_20ns/traj_wrapped.xtc ]]; then
                        log "  FETCH INCOMPLETE for $leg: no traj_wrapped.xtc -- did the in-job analyze run? leaving for a retry"
                        remaining=$((remaining+1))
                        continue
                    fi
                    log "  fetched $(du -sh $src/prod_20ns | cut -f1)"
                else
                    log "  FETCH FAILED for $leg -- retrying next round"
                    remaining=$((remaining+1))
                    continue
                fi

                log "  running tail for $leg"
                if [[ $leg == *_dock<-> ]]; then
                    struct=${leg%_dock*}
                    pose=${leg##*_dock}
                    STRUCT=$struct POSE=$pose SYSNAME=$SYSNAME SRC=$BG/md/$struct \
                        zsh $REPO/code/run_dock_pose_md.sh >> $src/tail.log 2>&1
                    rc=$?
                else
                    zsh $BG/md/run_md20.sh $leg >> $src/tail.log 2>&1
                    rc=$?
                fi
                log "  tail exit $rc for $leg (log: $src/tail.log)"

                # MM/GBSA scatter: run_dock_pose_md.sh and run_md20.sh both leave reference.frc and
                # _MMPBSA_* behind, ~5-7 GB a leg if left (CLAUDE.md, MM/GBSA housekeeping). find
                # -delete rather than a multi-path rm, whose zsh globs abort the whole command when
                # one pattern matches nothing. Window xtc slices are not kept either: they are two
                # seconds of work off prod_20ns, which is kept.
                find $src \( -name 'reference.frc' -o -name '_MMPBSA_*' \) -delete 2>/dev/null
                find $src -path '*/first_*ns/traj_wrapped.xtc' -delete 2>/dev/null
                find $src -path '*/first_*ns/traj.dcd' -delete 2>/dev/null
                if [[ -f $src/prod_20ns/mmgbsa/FINAL_RESULTS_MMPBSA.dat ]]; then
                    log "  MM/GBSA present; leg now $(du -sh $src | cut -f1)"
                    touch $src/.drained
                else
                    log "  NO FINAL_RESULTS_MMPBSA.dat for $leg -- tail did not complete, leaving undrained"
                    remaining=$((remaining+1))
                fi
                ;;
            RUNNING|PENDING|REQUEUED|RESIZING|SUSPENDED)
                remaining=$((remaining+1))
                ;;
            *)
                # FAILED, CANCELLED, PREEMPTED, TIMEOUT, OUT_OF_MEMORY. A preempted leg restarts from
                # zero on gpuscavenger -- nothing to salvage -- so flag it rather than tailing a
                # partial trajectory, and let the feeder's list be re-run for it deliberately.
                log "!!! $leg ended $state -- not tailing. Re-stage and resubmit it deliberately."
                touch $src/.drain_failed
                ;;
        esac
    done
    (( remaining > 0 )) && sleep $INTERVAL
done

log "drain done: every leg is drained or flagged"
