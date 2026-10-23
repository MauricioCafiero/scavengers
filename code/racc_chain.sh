#!/bin/zsh
# *** SHELVED 2026-10-23 — DO NOT USE. Queue filling belongs to racc_feed.sh (see CLAUDE.md). ***
# This script's first real run burned all 17 replicate batches as empty no-op jobs in 40 s and
# emptied the queue for ~2 h (no LEGS/RUNREL injected, no env/cwd, hook-masked exit codes). The
# bugs are fixed in the version below and it is kept only for the record; the queue filler is
# racc_feed.sh BATCH=1 per batch, driven in order by runs/folds_md/feed_rest.sh.
#
# Chain-submit a sequence of BATCH=1 replicate batches so the hand-offs never need the ssh socket.
#
#   code/racc_chain.sh runs/folds_md/rep5.chain
#
# The chain file holds one batch per line, in run order:
#     <batch name> <SYSNAME> <time>
# with each <batch name>.legs the same file racc_feed.sh reads.
#
# What it does: stages every leg of every batch through the socket (the only socket-dependent
# work), then creates each batch's <name>_run.py and <name>_run.slurm on the cluster, appending to
# each slurm one hook line that sbatches the NEXT batch's slurm from inside the finishing job --
# cluster-side, socket-free -- and finally submits only the FIRST batch. From then on the chain
# advances itself; a lost socket costs fetch delays only.
#
# Why not racc_feed.sh + manual resubmission: gpuscavenger runs 1 job and accepts 3, and every
# hand-off otherwise costs a submission from this Mac at the exact moment the socket may be down;
# a chain of <=6 h batches keeps preemption loss bounded the way racc_feed.sh does.
#
# Preemption: a batch that dies mid-way simply stops the chain there; resubmit its slurm (sbatch
# <name>_run.slurm) and the chain resumes to its end. racc_run.py's per-leg guard skips finished
# legs either way. The squeue guard in the hook keeps a requeued batch from double-submitting its
# successor.
#
# The staging tar duplicates racc_feed.sh's BATCH=1 loop and chain_gen.sh (sent to the cluster at
# startup) duplicates submit_gpu.sh's template transform, minus the submit-at-creation: neither
# existing script can stage or create without submitting, which is exactly what a chain needs.
#
# The generated slurm is SELF-CONTAINED (fixed 2026-10-23 after the chain burned all 17 batches in
# 40 s doing nothing): racc_run.py reads LEGS and RUNREL from the environment, which neither the
# head submission nor a hook carries, so chain_gen sed-injects both into <name>_run.py; the slurm
# cd's to ~/remote_work/openmm and activates the conda env before the python line (submit-time
# cwd/env do not survive a bare-ssh submission -- the first chain version's jobs ran in $HOME and
# each produced an empty out-file); and the hook is followed by `exit $rc` so a failed batch no
# longer reads COMPLETED because its hook's sbatch succeeded.

set -u
REPO=${0:A:h:h}
SOCK=$HOME/.ssh/cm-racc
REMOTE=racc.rdg.ac.uk
CEILING=${CEILING:-3}
INTERVAL=${INTERVAL:-300}

caffeinate -w $$ &               # dies with this script

rs() { ssh -o ControlPath=$SOCK -o ConnectTimeout=15 -o BatchMode=yes $REMOTE "$@"; }
log() { print -r -- "[$(date '+%Y-%m-%d %H:%M:%S')] $*"; }

# cluster-side helper: build one batch's {name}_run.py + {name}_run.slurm from the same template
# submit_gpu.sh uses, plus the hand-off hook when a successor exists. Idempotent.
rs "cat > ~/remote_work/openmm/chain_gen.sh <<'GENEOF'
#!/bin/bash
set -e
cd ~/remote_work/openmm
name=\$1; ntime=\$2; nname=\$3; nntime=\$4; legs=\$5; sysname=\$6
cp racc_run.py \${name}_run.py
~/.conda/envs/openmm-md/bin/python - \"\$legs\" \"\$sysname\" \${name}_run.py <<'INJ'
import sys
legs, sysname, path = sys.argv[1:4]
s = open(path).read()
old_legs = 'os.environ.get(\"LEG\", \"shuffle_control_dock7\")'
old_rel = 'os.environ.get(\"RUNREL\", \"runs/octinoxate/md\")'
assert old_legs in s and old_rel in s, 'LEGS/RUNREL sentinels not found in racc_run.py'
s = s.replace(old_legs, 'os.environ.get(\"LEG\", \"' + legs + '\")')
s = s.replace(old_rel, 'os.environ.get(\"RUNREL\", \"runs/' + sysname + '/md\")')
open(path, 'w').write(s)
print('injected', legs, 'RUNREL runs/' + sysname + '/md')
INJ
cp ~/scripts/gpu.slurm \${name}_run.slurm
sed -i \"s/job-name=/job-name=\${name}_run/; s/mem=/mem=10G/\" \${name}_run.slurm
echo \"module load anaconda 2>/dev/null\" >> \${name}_run.slurm
echo \"source activate openmm-md 2>/dev/null\" >> \${name}_run.slurm
echo \"cd \\\$HOME/remote_work/openmm\" >> \${name}_run.slurm
echo \"python3.12 \${name}_run.py > \${name}_run.out\" >> \${name}_run.slurm
echo \"rc=\\\$?\" >> \${name}_run.slurm
if [ -n \"\$nname\" ]; then
  echo \"squeue -u \\\$USER -h -n \${nname}_run 2>/dev/null | grep -q . || sbatch --time=\$nntime --cpus-per-task=8 \\\$HOME/remote_work/openmm/\${nname}_run.slurm\" >> \${name}_run.slurm
fi
echo \"exit \\\$rc\" >> \${name}_run.slurm
echo \"gen \${name} ok\"
GENEOF
chmod +x ~/remote_work/openmm/chain_gen.sh" >/dev/null || { log "ABORT: could not send chain_gen.sh"; exit 1; }

typeset -a NAMES MOLS TIMES
while IFS= read -r line; do
  [[ -z ${line// /} || ${line[1]} == '#' ]] && continue
  parts=(${=line})
  NAMES+=($parts[1]); MOLS+=($parts[2]); TIMES+=($parts[3])
done < "$1"

n=${#NAMES}
log "chain of $n batches: ${NAMES[*]}"

# Parse every batch's legfile once. Needed in BOTH modes: staging uses the source paths, and the
# LEGS injection needs the leg names even when staging is skipped.
typeset -a LEGSETS
for i in {1..$n}; do
  name=${NAMES[$i]}; mol=${MOLS[$i]}
  lf=$REPO/runs/folds_md/$name.legs
  typeset -a SRCS LEGNAMES
  SRCS=(); LEGNAMES=()
  while IFS= read -r line; do
    [[ -z ${line// /} || ${line[1]} == '#' ]] && continue
    SRCS+=(${line%%:*}); LEGNAMES+=(${line##*:})
  done < $lf
  LEGSETS[$i]="${(j:,:)LEGNAMES}"
  if [[ ${SKIP_STAGE:-0} == 1 ]]; then
    log "SKIP_STAGE: not re-staging $name"
  else
    for j in {1..${#LEGNAMES}}; do
      leg=${LEGNAMES[$j]}; src=${SRCS[$j]}
      if [[ ! -f $src/system/system.xml ]]; then
        log "ABORT: $leg is not built ($src/system/system.xml missing)"; exit 1
      fi
      log "staging $leg from $src"
      tar czf - -C $src . | rs "mkdir -p ~/remote_work/openmm/runs/$mol/md/$leg && tar xzf - -C ~/remote_work/openmm/runs/$mol/md/$leg" \
        || { log "ABORT: staging failed for $leg"; exit 1; }
    done
  fi
done
log "legfile parse + staging complete"

# Build each batch's run.py + slurm, hooked to its successor. chain_gen.sh takes
# <name> <time(unused, kept for the log)> <next name> <next time>. HOOK_TAIL=<name>:<time> hooks
# this chain's LAST batch to the named successor chain head (a queued tail extending a chain).
typeset -a HOOKP
if [[ -n ${HOOK_TAIL:-} ]]; then
  HOOKP=(${HOOK_TAIL%%:*} ${HOOK_TAIL#*:})   # split on the FIRST colon only: times are h:mm:ss
fi
for i in {1..$n}; do
  name=${NAMES[$i]}
  if (( i < n )); then
    out=$(rs "bash ~/remote_work/openmm/chain_gen.sh ${name} ${TIMES[$i]} ${NAMES[$((i+1))]} ${TIMES[$((i+1))]} ${LEGSETS[$i]} ${MOLS[$i]}" 2>&1 | tail -1)
  elif [[ -n ${HOOK_TAIL:-} ]]; then
    out=$(rs "bash ~/remote_work/openmm/chain_gen.sh ${name} ${TIMES[$i]} ${HOOKP[1]} ${HOOKP[2]} ${LEGSETS[$i]} ${MOLS[$i]}" 2>&1 | tail -1)
  else
    out=$(rs "bash ~/remote_work/openmm/chain_gen.sh ${name} ${TIMES[$i]} '' '' ${LEGSETS[$i]} ${MOLS[$i]}" 2>&1 | tail -1)
  fi
  log "   $out"
  [[ $out == *"gen ${name} ok"* ]] || { log "ABORT: gen failed for $name: $out"; exit 1; }
done

# Submit the head batch; wait for a free submit-slot like racc_feed.sh does.
# STAGE_ONLY=1 skips this: a chained tail (e.g. rep6) arms its slurm files and hooks without
# submitting, because submitting a head now would occupy the third gpuscavenger slot and the
# chain's cluster-side hand-off sbatches regardless of job count (MaxSubmitJobsPU=3 means a
# fourth pending job's sbatch fails and silently stalls the chain).
if [[ ${STAGE_ONLY:-0} == 1 ]]; then
  log "STAGE_ONLY: staged + generated $n batches; head not submitted"
  exit 0
fi
# Submit the head batch; wait for a free submit-slot like racc_feed.sh does. The bash -lc wrap
# (workspace cwd + activated env) is load-bearing: sbatch exports the submitter's environment and
# cwd into the job, so a bare-ssh submission starts the job in $HOME without python -- exactly the
# 2026-10-23 failure mode where every chain job produced an empty out-file.
log "waiting for a free slot to submit the head (max $CEILING)"
while :; do
  cnt=$(rs 'squeue -u $USER -h | wc -l' 2>/dev/null | tr -d ' ')
  [[ -n $cnt && $cnt -lt $CEILING ]] && break
  sleep $INTERVAL
done
out=$(rs "bash -lc 'module load anaconda; source activate openmm-md; cd ~/remote_work/openmm && sbatch --time=${TIMES[1]} --cpus-per-task=8 ${NAMES[1]}_run.slurm'" 2>&1 | tail -1)
log "head submitted: $out"
log "chain armed: hand-offs run on the cluster; nothing further needs the ssh socket"