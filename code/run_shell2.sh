#!/bin/zsh
# Shell 2 end to end, unattended and resumable.
#
# Every stage is guarded by its own output file, so re-running this after an interruption picks up
# where it stopped rather than repeating hours of UMA relaxation. Stage 5 is the expensive one --
# twelve complexes at 30-70 minutes each -- and `binding_energy.py` writes a partial_<name>.json per
# structure as it finishes, so even a kill part way through stage 5 keeps what it has paid for.
#
# Names are deliberately distinct from shell 1's throughout (condense_shell2.csv, design_shell2.json,
# s2_* folds, fold_check_shell2.csv, binding_shell2.csv). Shell 1's results are committed and must
# not be touched.
#
# Usage:  zsh code/run_shell2.sh
set -u
cd /Users/cafierom/python_mac/peptidebuilder

# Holds off idle and system sleep for exactly as long as this script lives, independently of whatever
# session-level assertion happens to be armed. Dies with the script.
caffeinate -is -w $$ &

PY=.venv/bin/python
RUN=runs/octinoxate
LOGS=$RUN/logs
UMALOGS=$RUN/uma_logs
mkdir -p $LOGS $UMALOGS

SHELL2="phenylalanine:2,arginine:3,lysine:10,aspartic:21,glutamic:15,isoleucine:8,leucine:12,serine:32,tryptophan:5,glutamic:10,leucine:13,serine:2"

say() { print -r -- "\n=== $1 : $(date '+%H:%M:%S') ===" }

# ---------------------------------------------------------------------------------------------
say "stage 1: wait for the reachability sweep"
# The sweep was launched separately and is already running. Wait rather than restart it.
while pgrep -f 'condense.py sweep' > /dev/null; do sleep 60; done
pairs=$(cut -d, -f1-4 $RUN/condense_shell2.csv | tail -n +2 | sort -u | wc -l | tr -d ' ')
print -r -- "sweep finished with $pairs/132 ordered pairs"
if [ "$pairs" -lt 132 ]; then
    print -r -- "STOPPING: the sweep did not cover all 132 pairs, so the assignment would be"
    print -r -- "searching an incomplete graph and could miss the best ordering."
    exit 1
fi
# The sweep log is pure geometry -- no UMA runs in it -- so it belongs with its siblings in logs/.
[ -f runs/condense_shell2.log ] && mv runs/condense_shell2.log $LOGS/condense_shell2.log

# ---------------------------------------------------------------------------------------------
say "stage 2: ordering and spacer counts"
if [ -f $RUN/design_shell2.json ]; then
    print -r -- "design_shell2.json exists, skipping"
else
    $PY code/assign.py $RUN --poses "$SHELL2" --pairs-csv condense_shell2.csv \
        --save $RUN/design_shell2.json 2>&1 | tee $LOGS/assign_shell2.log
    [ -f $RUN/design_shell2.json ] || { print -r -- "STOPPING: assign.py wrote no design"; exit 1; }
fi
S2SEQ=$($PY -c "import json;print(json.load(open('$RUN/design_shell2.json'))['sequence'])")
print -r -- "shell 2 design: $S2SEQ"

# ---------------------------------------------------------------------------------------------
say "stage 3: two ESM2 linker variants"
# fill() is the in-process API design_test.py uses. Calling it directly keeps sequences.csv
# untouched -- fill_linkers.py's own main() rewrites that file, and it is committed shell-1 data.
# Sampling, not argmax: with this many masks the 35M checkpoint is nearly flat and the argmax
# collapses every linker to leucine.
if [ -f $RUN/variants_shell2.txt ]; then
    print -r -- "variants_shell2.txt exists, skipping"
else
    $PY -c "
import sys, json
sys.path.insert(0, 'code')
from fill_linkers import fill
seq = json.load(open('$RUN/design_shell2.json'))['sequence']
out = []
for item in fill([seq], variants=2):
    out.append(item['filled'] if isinstance(item, dict) else item)
for s in out:
    print(s)
open('$RUN/variants_shell2.txt','w').write('\n'.join(out) + '\n')
" 2>&1 | tee $UMALOGS/fill_shell2.log
fi
ESM1=$(sed -n 1p $RUN/variants_shell2.txt 2>/dev/null)
ESM2=$(sed -n 2p $RUN/variants_shell2.txt 2>/dev/null)
print -r -- "esm1: ${ESM1:-<none>}"
print -r -- "esm2: ${ESM2:-<none>}"

# ---------------------------------------------------------------------------------------------
say "stage 4: twelve folds"
# force AND --use_potentials are both required. `force: false` is discarded by the featurizer, and
# `force: true` alone is only conditioning -- the guidance that acts on the structure during
# sampling is gated behind --use_potentials, which is off by default.
fold_one() {
    local name=$1 seq=$2 top=$3
    local cif=$RUN/boltz/boltz_results_$name/predictions/$name/${name}_model_0.cif
    if [ -f "$cif" ]; then print -r -- "  $name: already folded, skipping"; return 0; fi
    local cmd=($PY code/boltz_hints.py $RUN design_shell2.json --name $name --run)
    if [ "$top" = "control" ]; then
        cmd+=(--no-hints)
    else
        # --boltz-args=--use_potentials must be ONE argv entry. Written space-separated on the night
        # of 2026-09-26, which argparse rejects with "expected one argument" because the value starts
        # with a dash -- all nine forced folds died at parse time before this was caught.
        cmd+=(--top $top --force --boltz-args=--use_potentials)
    fi
    [ -n "$seq" ] && cmd+=(--sequence $seq)
    print -r -- "  folding $name"
    "${cmd[@]}" >> $LOGS/folds_shell2.log 2>&1 || print -r -- "  $name: FAILED, see folds_shell2.log"
}

# suffix:top -- the suffix has to read f4/f8/f12 because that is the ladder make_figures.py sorts by
for set_tag seq in orig "" esm1 "$ESM1" esm2 "$ESM2"; do
    [ "$set_tag" != "orig" ] && [ -z "$seq" ] && continue
    for pair in control:control f4:4 f8:8 f12:12; do
        fold_one s2_${set_tag}_${pair%%:*} "$seq" ${pair##*:}
    done
done

FOLDED=()
for name in s2_orig_control s2_orig_f4 s2_orig_f8 s2_orig_f12 \
            s2_esm1_control s2_esm1_f4 s2_esm1_f8 s2_esm1_f12 \
            s2_esm2_control s2_esm2_f4 s2_esm2_f8 s2_esm2_f12; do
    [ -f $RUN/boltz/boltz_results_$name/predictions/$name/${name}_model_0.cif ] && FOLDED+=($name)
done
print -r -- "${#FOLDED} of 12 folds present"
[ ${#FOLDED} -eq 0 ] && { print -r -- "STOPPING: nothing folded"; exit 1; }

# ---------------------------------------------------------------------------------------------
say "stage 5: geometry check (free, and the primary measure)"
# --out keeps this off shell 1's committed fold_check.csv. check_fold.py scans every fold under
# boltz/, so this file carries both shells; that is a superset, not a loss.
$PY code/check_fold.py $RUN --out fold_check_shell2.csv 2>&1 | tee $LOGS/check_fold_shell2.log

# ---------------------------------------------------------------------------------------------
say "stage 6: score every fold -- THE BOTTLENECK, 30-70 min each on CPU"
# Everything must be scored, not a subset. strain_peptide is left out deliberately: a free peptide
# collapses in vacuum, so that term measures collapse rather than strain.
# This is also the first scoring run since the relaxed-geometry saving was added, so boltz/structures/
# should appear and fill as it goes.
$PY code/binding_energy.py $RUN --terms interaction,strain_ligand \
    --out binding_shell2.csv --structures "${(j:,:)FOLDED}" 2>&1 | tee -a $UMALOGS/score_shell2.log

say "stage 7: figures, into their own directory"
# --fig-dir and --match keep shell 1's figures/ and its manifest intact.
$PY code/make_figures.py $RUN --fig-dir figures_shell2 --match s2_ \
    --design-sequence "$S2SEQ" 2>&1 | tee $LOGS/make_figures_shell2.log

say "shell 2 complete"
