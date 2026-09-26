#!/bin/zsh
# The nine forced folds of shell 2, plus the downstream stages, with the argument bug fixed.
#
# Why this exists: `run_shell2.sh` passes `--boltz-args '--use_potentials'`, which reaches argparse
# as two separate argv entries, so `--use_potentials` is read as an unknown flag and every forced
# fold dies at parse time. The working form is `--boltz-args=--use_potentials`, as a single entry.
# That script could not be corrected in place because zsh reads a script incrementally as it runs.
#
# Both halves are required and neither is optional: `force: false` is discarded by Boltz's
# featurizer rather than softened, and `force: true` alone is only conditioning -- the guidance that
# acts on the structure during sampling is gated behind --use_potentials, which is off by default.
#
# Safe to run alongside or after run_shell2.sh: every stage is guarded by its own output, so a fold
# that already exists is skipped rather than repeated, and nothing of shell 1's is written.
#
# Usage:  zsh code/run_shell2_forced.sh
set -u
cd /Users/cafierom/python_mac/peptidebuilder

caffeinate -is -w $$ &

PY=.venv/bin/python
RUN=runs/octinoxate
LOGS=$RUN/logs
UMALOGS=$RUN/uma_logs
mkdir -p $LOGS $UMALOGS

say() { print -r -- ""; print -r -- "=== $1 : $(date '+%H:%M:%S') ===" }

[ -f $RUN/design_shell2.json ] || { print -r -- "no design_shell2.json yet -- let run_shell2.sh reach stage 2 first"; exit 1; }
S2SEQ=$($PY -c "import json;print(json.load(open('$RUN/design_shell2.json'))['sequence'])")
ESM1=$(sed -n 1p $RUN/variants_shell2.txt 2>/dev/null)
ESM2=$(sed -n 2p $RUN/variants_shell2.txt 2>/dev/null)
print -r -- "design: $S2SEQ"
print -r -- "esm1:   ${ESM1:-<none>}"
print -r -- "esm2:   ${ESM2:-<none>}"

fold_one() {
    local name=$1 seq=$2 top=$3
    local cif=$RUN/boltz/boltz_results_$name/predictions/$name/${name}_model_0.cif
    if [ -f "$cif" ]; then print -r -- "  $name: already folded, skipping"; return 0; fi
    # --boltz-args=--use_potentials as ONE argv entry: the space-separated form fails to parse.
    local cmd=($PY code/boltz_hints.py $RUN design_shell2.json --name $name --run
               --top $top --force --boltz-args=--use_potentials)
    [ -n "$seq" ] && cmd+=(--sequence $seq)
    print -r -- "  folding $name"
    "${cmd[@]}" >> $LOGS/folds_shell2.log 2>&1 || print -r -- "  $name: FAILED, see folds_shell2.log"
}

say "the nine forced folds"
for set_tag seq in orig "" esm1 "$ESM1" esm2 "$ESM2"; do
    [ "$set_tag" != "orig" ] && [ -z "$seq" ] && continue
    for pair in f4:4 f8:8 f12:12; do
        fold_one s2_${set_tag}_${pair%%:*} "$seq" ${pair##*:}
    done
done

# Score only the forced folds. binding_energy.py has no "already scored" check -- it writes a
# partial_<name>.json per structure but never reads one back to skip work -- so handing it the three
# controls that run_shell2.sh has already scored would repeat 2-3 hours of relaxation for nothing.
# The controls land in binding_shell2.csv, these nine in binding_shell2_forced.csv, and
# make_figures.py globs binding_*.csv so the manifest sees all twelve.
FOLDED=()
for name in s2_orig_f4 s2_orig_f8 s2_orig_f12 \
            s2_esm1_f4 s2_esm1_f8 s2_esm1_f12 \
            s2_esm2_f4 s2_esm2_f8 s2_esm2_f12; do
    [ -f $RUN/boltz/boltz_results_$name/predictions/$name/${name}_model_0.cif ] && FOLDED+=($name)
done
ALL=(); for name in $RUN/boltz/boltz_results_s2_*(N); do ALL+=(${name:t}); done
print -r -- "${#FOLDED} of 9 forced folds present; ${#ALL} s2_ fold directories in total"
[ ${#FOLDED} -eq 0 ] && { print -r -- "STOPPING: no forced folds to score"; exit 1; }

say "geometry check"
$PY code/check_fold.py $RUN --out fold_check_shell2.csv 2>&1 | tee $LOGS/check_fold_shell2.log

say "score every fold -- 30-70 min each on CPU"
# binding_energy.py writes a partial_<name>.json per structure as it finishes, and skips work it
# already has, so this is safe to re-run after an interruption.
# Its own --out, so this cannot clobber whatever run_shell2.sh already wrote to binding_shell2.csv.
# make_figures.py globs binding_*.csv, so both files are read into one manifest.
$PY code/binding_energy.py $RUN --terms interaction,strain_ligand \
    --out binding_shell2_forced.csv --structures "${(j:,:)FOLDED}" 2>&1 \
    | tee -a $UMALOGS/score_shell2_forced.log

say "figures, into their own directory"
$PY code/make_figures.py $RUN --fig-dir figures_shell2 --match s2_ \
    --design-sequence "$S2SEQ" 2>&1 | tee $LOGS/make_figures_shell2.log

say "shell 2 complete"
