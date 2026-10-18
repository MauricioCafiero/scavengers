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
#   <pep>_o3cof / <pep>_rf3cof  ->  code/run_dock_pose_md.sh POSE=cofold (fold-lane cofold legs)
#   <leg> under $REPO  ->  code/run_dock_pose_md.sh POSE=cofold (bare-named Boltz co-folds, e.g. *_control)
#   <leg>           ->  boltzgen_local/md/run_md20.sh <leg>
# Both are racc-aware: they skip the dynamics when only traj_wrapped.xtc is present, and slice their
# windows from the wrapped solute rather than loading a full box this 8 GB Mac cannot hold.
#
# Only the wrapped solute is fetched (~100 MB), never the full dcd (1.8-5 GB). The dcd stays on
# cluster scratch, which is also what the archive keeps.
#
# Design notes, 2026-10-16:
# * Fetch is decoupled from tail. The old drain ran one tail INLINE, which blocked the whole poll
#   loop for ~20 minutes; a finished leg behind a queued tail sat on cluster scratch all that time,
#   and a wedged poll made it sit much longer (2026-10-16: three legs through an afternoon). Now
#   every poll first fetches every completed leg that is not yet local -- fetching costs seconds --
#   and only then launches new tails, detached.
# * Tails stream through the tail lock one at a time (concurrent tails OOM this 8 GB Mac); so do the
#   per-leg record launches. Cluster production is ~1 leg / 50 min against a ~20 min tail, so the
#   local stream keeps up with the legs as they land; the queue only grows while the laptop sleeps,
#   which is exactly what capturing-at-fetch is for.
# * Bookkeeping is local and offline-safe: a tail counts as drained when its process has exited
#   leaving prod_20ns/mmgbsa/FINAL_RESULTS_MMPBSA.dat, judged per poll WITHOUT ssh; its MM/GBSA
#   scatter strip and record then happen locally too. A dead connection costs only the next fetch,
#   never the analysis of what is already here.
# * Every ssh call runs under a watchdog cap, not only keepalives: a wedged ControlMaster swallows
#   the client-side ServerAlive probes (2026-10-16: two drains hung for hours in one poll despite
#   ConnectTimeout/ServerAlive*), so the kill must be ours.
set -u

LEGFILE=${1:?usage: racc_drain.sh <legfile>}
INTERVAL=${INTERVAL:-300}
SOCK=$HOME/.ssh/cm-racc
REMOTE=racc.rdg.ac.uk
DCD=/scratch5/gaussian/io927423/dcd
BG=$HOME/python_mac/boltzgen_local
REPO=${0:A:h:h}
MDPY=$HOME/miniforge3/envs/openmm-md/bin/python
SYSNAME=${SYSNAME:-oxybenzone}

caffeinate -w $$ &

# r <seconds> <remote command>: ssh bounded by a hard watchdog. ConnectTimeout/ServerAlive* cover a
# dead network; the wedged-master case (socket path alive, master stuck) is what the cap is for.
r() { cap=$1; shift
      ssh -n -o ControlPath=$SOCK -o ConnectTimeout=15 -o ServerAliveInterval=10 \
          -o ServerAliveCountMax=3 -o BatchMode=yes $REMOTE "$@" & pid=$!
      { sleep $cap; kill -9 $pid 2>/dev/null } & wd=$!
      wait $pid; rc=$?
      kill $wd 2>/dev/null
      return $rc }
log() { print -r -- "[$(date '+%Y-%m-%d %H:%M:%S')] $*"; }

typeset -a SRCS LEGS
while IFS= read -r line; do
    [[ -z ${line// /} || ${line[1]} == '#' ]] && continue
    SRCS+=("${line%%:*}")
    LEGS+=("${line##*:}")
done < $LEGFILE

log "drain started: watching ${#LEGS} legs, interval ${INTERVAL}s"

launch_tail() {    # launch_tail <src> <leg> -- detached, queued on the tail lock by the driver
    local src=$1 leg=$2 base struct pose
    touch $src/.tail_launched
    log "  launching tail for $leg (detached; queues on the tail lock)"
    if [[ $leg == *_r2 ]]; then
        # A replicate leg (<leg>_r2) carries the word dock in the middle of its base name, so
        # STRUCT/POSE must be parsed from the base with _r2 stripped, and LEGNAME passes the real
        # leg name through (run_dock_pose_md.sh derives N=<STRUCT>_dock<POSE> otherwise). The
        # staged _r2 dir already holds a copied system/, so the driver's build guard skips and SRC
        # is never read -- it is set to the base's co-fold dir only to keep the args complete.
        base=${leg%_r2}; struct=$base; pose=1
        [[ $base == *_dock<-> ]] && { struct=${base%_dock*}; pose=${base##*_dock}; }
        ( export STRUCT=$struct POSE=$pose LEGNAME=$leg SYSNAME=$SYSNAME SRC=$BG/md/$struct
          exec zsh $REPO/code/run_dock_pose_md.sh ) >> $src/tail.log 2>&1 &
        echo $! > $src/.tail_pid
    elif [[ $leg == *_o3cof || $leg == *_rf3cof ]]; then
        # A fold-lane cofold leg (FOLDS.md): the ligand is the cofolder's own placement, the
        # receptor is the leg's own prepped protein, so SRC is the leg dir itself. NOTE: branch
        # order matters -- <pep>_o3cof_dock<N> matches here BEFORE the *_dock<-> branch below.
        ( export STRUCT=$leg POSE=cofold SYSNAME=$SYSNAME SRC=$src
          exec zsh $REPO/code/run_dock_pose_md.sh ) >> $src/tail.log 2>&1 &
        echo $! > $src/.tail_pid
    elif [[ $leg == *_dock<-> ]]; then
        struct=${leg%_dock*}
        pose=${leg##*_dock}
        ( export STRUCT=$struct POSE=$pose SYSNAME=$SYSNAME SRC=$BG/md/$struct
          exec zsh $REPO/code/run_dock_pose_md.sh ) >> $src/tail.log 2>&1 &
        echo $! > $src/.tail_pid
    elif [[ $src == $REPO/* ]]; then
        # A bare-named Boltz co-fold leg living in this repo (the unforced *_control folds,
        # NEXT_STEPS.md "Forced vs unforced co-folds"): run_md20.sh only knows boltzgen_local legs,
        # so tail it like a fold-lane cofold leg -- its ligand SDF came from cif_to_md.py.
        ( export STRUCT=$leg POSE=cofold SYSNAME=$SYSNAME SRC=$src
          exec zsh $REPO/code/run_dock_pose_md.sh ) >> $src/tail.log 2>&1 &
        echo $! > $src/.tail_pid
    else
        ( exec zsh $BG/md/run_md20.sh $leg ) >> $src/tail.log 2>&1 &
        echo $! > $src/.tail_pid
    fi
    disown
}

local_phase() {    # offline-safe: resolve recently finished tails, queue their records
    local src leg i pid rc
    src=""; leg=""; pid=""; rc=0
    for i in {1..${#LEGS}}; do
        src=${SRCS[$i]}; leg=${LEGS[$i]}
        [[ -f $src/.drained || -f $src/.drain_failed ]] && continue
        [[ -f $src/.tail_launched ]] || continue
        pid=""
        [[ -f $src/.tail_pid ]] && pid=$(< $src/.tail_pid)
        if [[ -n $pid ]] && kill -0 $pid 2>/dev/null; then
            remaining=$((remaining+1))
            continue
        fi
        if [[ -f $src/prod_20ns/mmgbsa/FINAL_RESULTS_MMPBSA.dat ]]; then
            # Same strip the inline flow did: MMPBSA scatter and window slices are regenerable
            # (traj.nc is MMPBSA.py's NetCDF conversion of kept traj_wrapped.xtc; reference.frc and
            # _MMPBSA_* pile to 5-7 GB a leg if left -- CLAUDE.md, MM/GBSA housekeeping).
            find $src \( -name 'reference.frc' -o -name '_MMPBSA_*' -o -name 'traj.nc' \) -delete 2>/dev/null
            find $src -path '*/first_*ns/traj_wrapped.xtc' -delete 2>/dev/null
            find $src -path '*/first_*ns/traj.dcd' -delete 2>/dev/null
            log "  tail done; MM/GBSA present; leg now $(du -sh $src | cut -f1)"
            touch $src/.drained
            # The leg's tables + FOLDS row, under the same tail lock so it never overlaps a tail.
            ( exec lockf -k /tmp/peptidebuilder_tail.lock $MDPY $REPO/code/record_leg.py $SYSNAME $leg ) \
                >> $REPO/runs/folds_md/record_$leg.log 2>&1 &
            disown
        else
            log "  tail exited without FINAL_RESULTS_MMPBSA.dat for $leg -- relaunching"
            launch_tail $src $leg
            remaining=$((remaining+1))
        fi
    done
}

remote_phase() {    # bounded ssh calls: status + eager fetch of every completed leg not yet local
    online=$1; src=""; leg=""
    local i state
    for i in {1..${#LEGS}}; do
        src=${SRCS[$i]}; leg=${LEGS[$i]}
        [[ -f $src/.drained || -f $src/.drain_failed ]] && continue
        if [[ $online == 0 ]]; then
            remaining=$((remaining+1)); continue
        fi

        if [[ -f $src/prod_20ns/traj_wrapped.xtc ]]; then
            # Already local, so the leg is finished: racc_run.py writes traj_wrapped.xtc as the
            # LAST thing it does on scratch. Never re-fetch it and never ask sacct about it --
            # its tail is either queued/running (re-launches only if phase A saw it die) or done.
            if [[ ! -f $src/.tail_launched ]]; then
                launch_tail $src $leg
            fi
            remaining=$((remaining+1))
        elif r 30 "test -f $DCD/${leg}_prod/traj_wrapped.xtc" 2>/dev/null; then
            state=COMPLETED
        else
            # Not finished. Only consult sacct to tell "still going / not yet started" from a leg
            # whose single-leg job died -- worth flagging rather than waiting on.
            state=$(r 60 "sacct --name=${leg}_run --format=State%20 -n -X 2>/dev/null | tail -1" | tr -d ' ')
            [[ -z $state ]] && state=PENDING    # not submitted yet, or inside a batch job
            case $state in
                RUNNING|PENDING|REQUEUED|RESIZING|SUSPENDED) remaining=$((remaining+1)) ;;
                *)
                    # FAILED, CANCELLED, PREEMPTED, TIMEOUT, OUT_OF_MEMORY. A preempted leg restarts
                    # from zero on gpuscavenger -- nothing to salvage -- so flag it rather than
                    # tailing a partial trajectory, and let the feeder's list be re-run for it
                    # deliberately.
                    log "!!! $leg ended $state -- not tailing. Re-stage and resubmit it deliberately."
                    touch $src/.drain_failed ;;
            esac
            continue
        fi

        if [[ $state == COMPLETED ]]; then
            log "=== $leg COMPLETED -- fetching wrapped solute"
            mkdir -p $src/prod_20ns
            if r 600 "cd $DCD/${leg}_prod 2>/dev/null && tar czf - traj_wrapped.xtc traj_wrapped.pdb energy.csv rmsd.csv analysis.png final.pdb 2>/dev/null" | tar xzf - -C $src/prod_20ns 2>/dev/null; then
                if [[ ! -f $src/prod_20ns/traj_wrapped.xtc ]]; then
                    log "  FETCH INCOMPLETE for $leg: no traj_wrapped.xtc -- did the in-job analyze run? retrying next round"
                    remaining=$((remaining+1)); continue
                fi
                log "  fetched $(du -sh $src/prod_20ns | cut -f1)"
            else
                log "  FETCH FAILED for $leg -- retrying next round"
                remaining=$((remaining+1)); continue
            fi
            launch_tail $src $leg
            remaining=$((remaining+1))
        fi
    done
}

remaining=${#LEGS}
while (( remaining > 0 )); do
    remaining=0
    local_phase
    online=1
    r 30 true 2>/dev/null || online=0
    if [[ $online == 0 ]]; then
        log "ssh socket down -- holding (nothing is lost; the run continues on the cluster)"
    fi
    remote_phase $online
    (( remaining > 0 )) && sleep $INTERVAL
done

log "drain done: every leg is drained or flagged"