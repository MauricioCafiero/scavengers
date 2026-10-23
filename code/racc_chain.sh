#!/bin/zsh
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
name=\$1; ntime=\$2; nname=\$3; nntime=\$4
cp racc_run.py \${name}_run.py
cp ~/scripts/gpu.slurm \${name}_run.slurm
sed -i \"s/job-name=/job-name=\${name}_run/; s/mem=/mem=10G/\" \${name}_run.slurm
echo \"python3.12 \${name}_run.py > \${name}_run.out\" >> \${name}_run.slurm
if [ -n \"\$nname\" ]; then
  echo \"squeue -u \\\$USER -h -n \${nname}_run 2>/dev/null | grep -q . || sbatch --time=\$nntime --cpus-per-task=8 \\\$HOME/remote_work/openmm/\${nname}_run.slurm\" >> \${name}_run.slurm
fi
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

# Stage every leg of every batch -- the same tar racc_feed.sh's BATCH=1 path uses (idempotent).
# One molecule per batch; the SYSNAME field decides which runs/<mol>/md/ subtree it lands in.
for i in {1..$n}; do
  name=${NAMES[$i]}; mol=${MOLS[$i]}
  lf=$REPO/runs/folds_md/$name.legs
  typeset -a SRCS LEGNAMES
  while IFS= read -r line; do
    [[ -z ${line// /} || ${line[1]} == '#' ]] && continue
    SRCS+=(${line%%:*}); LEGNAMES+=(${line##*:})
  done < $lf
  for j in {1..${#LEGNAMES}}; do
    leg=${LEGNAMES[$j]}; src=${SRCS[$j]}
    if [[ ! -f $src/system/system.xml ]]; then
      log "ABORT: $leg is not built ($src/system/system.xml missing)"; exit 1
    fi
    log "staging $leg from $src"
    tar czf - -C $src . | rs "mkdir -p ~/remote_work/openmm/runs/$mol/md/$leg && tar xzf - -C ~/remote_work/openmm/runs/$mol/md/$leg" \
      || { log "ABORT: staging failed for $leg"; exit 1; }
  done
done
log "staging complete"

# Build each batch's run.py + slurm, hooked to its successor. chain_gen.sh takes
# <name> <time(unused, kept for the log)> <next name> <next time>.
for i in {1..$n}; do
  name=${NAMES[$i]}
  if (( i < n )); then
    out=$(rs "bash ~/remote_work/openmm/chain_gen.sh ${name} ${TIMES[$i]} ${NAMES[$((i+1))]} ${TIMES[$((i+1))]}" 2>&1 | tail -1)
  else
    out=$(rs "bash ~/remote_work/openmm/chain_gen.sh ${name} ${TIMES[$i]} '' ''" 2>&1 | tail -1)
  fi
  log "   $out"
done

# Submit the head batch; wait for a free submit-slot like racc_feed.sh does.
log "waiting for a free slot to submit the head (max $CEILING)"
while :; do
  cnt=$(rs 'squeue -u $USER -h | wc -l' 2>/dev/null | tr -d ' ')
  [[ -n $cnt && $cnt -lt $CEILING ]] && break
  sleep $INTERVAL
done
out=$(rs "sbatch --time=${TIMES[1]} --cpus-per-task=8 ~/remote_work/openmm/${NAMES[1]}_run.slurm" 2>&1 | tail -1)
log "head submitted: $out"
log "chain armed: hand-offs run on the cluster; nothing further needs the ssh socket"