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
#   ConnectTimeout/ServerAlive*), so the kill must be ours. Every call also forks transient zsh
#   children (the ssh job, its watchdog) that carry this script's ARGV until they exec -- in ps they
#   look exactly like fresh duplicate drains, parented to the real one, and vanish seconds later.
#   2026-10-24: twice mistaken for a runaway spawner and "killed" (harmlessly); ps-etime two samples
#   apart to prove it. The real fleet is one drain per legfile, nothing else.
set -u

LEGFILE=${1:?usage: racc_drain.sh <legfile>}
INTERVAL=${INTERVAL:-300}
SOCK=${DR_SOCK-$HOME/.ssh/cm-racc}     # empty DR_SOCK => plain ssh, no ControlMaster (Young). NOT :- (2026-10-27: :- made an empty DR_SOCK fall back to cm-racc, so every Young drain polled racc and silently saw nothing)
REMOTE=${DR_REMOTE:-racc.rdg.ac.uk}
DCD=${DR_DCD:-/scratch5/gaussian/io927423/dcd}
# How the leg's prod dir is built from the leg name ON THE CLUSTER: racc keeps it on scratch as
# "<leg>_prod"; Young keeps it in the tree as "runs/octinoxate/md/<leg>/prod_20ns" (DCD then the
# workspace root). ${(e)...} expands the inner ${leg} per leg. (Young: DR_DCD='$HOME/openmm'
# DR_LEGDIR='runs/octinoxate/md/${leg}/prod_20ns' DR_SOCK= DR_REMOTE=young-ng SYSNAME=octinoxate)
LEGDIR=${DR_LEGDIR:-'${leg}_prod'}
BG=$HOME/python_mac/boltzgen_local
REPO=${0:A:h:h}
MDPY=$HOME/miniforge3/envs/openmm-md/bin/python
SYSNAME=${SYSNAME:-oxybenzone}

caffeinate -w $$ &

# r <seconds> <remote command>: ssh bounded by a hard watchdog. ConnectTimeout/ServerAlive* cover a
# dead network; the wedged-master case (socket path alive, master stuck) is what the cap is for.
r() { cap=$1; shift
      typeset -a ssho
      ssho=()
      [[ -n $SOCK ]] && ssho=(-o ControlPath=$SOCK)
      ssh -n "${ssho[@]}" -o ConnectTimeout=15 -o ServerAliveInterval=10 \
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
    elif [[ $leg == *_o3cof || $leg == *_rf3cof || $leg == *_o3pcof ]]; then
        # A fold-lane cofold leg (FOLDS.md): the ligand is the cofolder's own placement, the
        # receptor is the leg's own prepped protein, so SRC is the leg dir itself. `_o3pcof` is the
        # pocket-guided OF3 lane (NEXT_STEPS.md Part 2); the driver picks its stage dir off the
        # same suffix.
        #
        # These patterns are END-anchored, so a DOCKED cofold leg (<pep>_o3cof_dock<N>) does NOT
        # land here -- it falls to the *_dock<-> branch below and is handed SRC=$BG/md/<struct>,
        # which for a fold-lane structure does not exist. That has been inert through 66 docked
        # legs because run_dock_pose_md.sh reads SRC only inside its `system.xml` build guard, and
        # every campaign leg is built locally first (code/build_cofolds.sh) -- so the tail never
        # looks at it. A comment here used to claim the opposite ("branch order matters -- matches
        # here BEFORE the *_dock<-> branch"); tested false 2026-10-19. If a docked leg is ever
        # drained UNBUILT, that branch needs the repo-local SRC before it will build.
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
            # The leg's tables + FOLDS row, under its OWN lock: records serialise among themselves (shared
            # tables) but never wait on a tail, which writes only inside its leg dir. Until 2026-10-27 they
            # shared the tail lock and starved behind every queued tail (one waited 2 h for a seconds-long job).
            ( exec lockf -k /tmp/peptidebuilder_record.lock $MDPY $REPO/code/record_leg.py $SYSNAME $leg ) \
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
        elif r 30 "test -f $DCD/${(e)LEGDIR}/traj_wrapped.xtc" 2>/dev/null; then
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
            # Fetch into a staging dir and move traj_wrapped.xtc into place LAST: its presence is what
            # marks the leg finished (line above), so a fetch that dies mid-stream must never leave
            # one behind (2026-10-27: a 0-byte xtc from a dropped fetch was tailed in a loop).
            stage=$src/prod_20ns/.fetch
            rm -rf $stage; mkdir -p $stage
            if r 600 "cd $DCD/${(e)LEGDIR} 2>/dev/null && tar czf - traj_wrapped.xtc traj_wrapped.pdb energy.csv rmsd.csv analysis.png final.pdb 2>/dev/null" | tar xzf - -C $stage 2>/dev/null \
               && [[ -s $stage/traj_wrapped.xtc && -s $stage/traj_wrapped.pdb && -s $stage/energy.csv ]]; then
                for f in $stage/*(N); do [[ ${f:t} == traj_wrapped.xtc ]] || mv -f $f $src/prod_20ns/; done
                mv -f $stage/traj_wrapped.xtc $src/prod_20ns/
                rm -rf $stage
                log "  fetched $(du -sh $src/prod_20ns | cut -f1)"
            else
                rm -rf $stage
                log "  FETCH FAILED or INCOMPLETE for $leg (stream died, or no xtc/pdb/energy) -- retrying next round"
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