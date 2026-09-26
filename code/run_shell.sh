#!/bin/zsh
# One shell, end to end, unattended and resumable. Parameterised so a new shell needs no new script.
#
#   zsh code/run_shell.sh <prefix> <shell-number> "<pose spec>"
#
#   zsh code/run_shell.sh s3 3 "tryptophan:2,tyrosine:0,phenylalanine:4,lysine:1,aspartic:21,glutamic:8,leucine:5,serine:4,tryptophan:5,lysine:8,leucine:13,serine:32"
#
# This exists because the shell-2 run used two hand-copied scripts, and copying the Boltz command line
# by hand is how nine folds were lost: `--boltz-args '--use_potentials'` reaches argparse as two argv
# entries, and a value beginning with a dash is rejected with "expected one argument". The working form
# is `--boltz-args=--use_potentials`, as a single entry. It is written once here, and **stage 0
# preflights it** -- one fold's yaml is written without folding and checked for `force: true`, which
# takes seconds and fails loudly rather than after hours of wasted folding.
#
# Both halves are required and neither is optional: `force: false` is discarded by Boltz's featurizer
# rather than softened, and `force: true` alone is only conditioning -- the guidance that acts on the
# structure during sampling is gated behind --use_potentials, which is off by default.
#
# Every stage is guarded by its own output, so re-running after an interruption resumes rather than
# repeats. `binding_energy.py` also writes a partial_<name>.json per structure as it finishes, so even
# a kill part way through scoring keeps what it has paid for. It has no "already scored" check of its
# own, so stage 6 passes it only the structures with no partial yet.
#
# Nothing of another shell's is written: every output is named from <prefix> or <shell-number>.
#
# Ligand strain comes from stage 5c, not from the scoring stage: one reference conformer per run
# (ligand_reference.json) and one cheap pass over the folds (strain_shell<N>.csv). The scoring stage
# computes the interaction energy alone.
set -u

if [ $# -lt 3 ]; then
    print -r -- "usage: zsh code/run_shell.sh <prefix> <shell-number> \"<pose spec>\""
    print -r -- "   eg: zsh code/run_shell.sh s3 3 \"tryptophan:2,tyrosine:0,...\""
    exit 2
fi
PREFIX=$1
NUM=$2
POSES=$3

cd "${0:A:h:h}"          # repo root, from this script's own location -- no absolute path baked in

# Holds off idle and system sleep for exactly as long as this script lives. Dies with the script.
caffeinate -is -w $$ &

PY=.venv/bin/python
RUN=runs/octinoxate
LOGS=$RUN/logs
UMALOGS=$RUN/uma_logs
mkdir -p $LOGS $UMALOGS

SWEEP=condense_shell$NUM.csv
DESIGN=design_shell$NUM.json
VARIANTS=$RUN/variants_shell$NUM.txt
FOLDCHECK=fold_check_shell$NUM.csv
BINDING=binding_shell$NUM.csv
FIGDIR=figures_shell$NUM

say() { print -r -- ""; print -r -- "=== $1 : $(date '+%H:%M:%S') ===" }

say "shell $NUM, prefix $PREFIX"
print -r -- "poses: $POSES"

# ---------------------------------------------------------------------------------------------
say "stage 0: preflight the forced-contact flags"
# Seconds, and it is the check whose absence cost nine folds. Needs a design, so it runs again after
# stage 2 if there is none yet; here it only fires when resuming.
preflight() {
    local probe=${PREFIX}_preflight_probe
    $PY code/boltz_hints.py $RUN $DESIGN --name $probe --top 4 --force \
        --boltz-args=--use_potentials > $LOGS/preflight_shell$NUM.log 2>&1
    local rc=$?
    local yaml=$RUN/boltz/$probe.yaml
    if [ $rc -ne 0 ] || [ ! -f "$yaml" ]; then
        print -r -- "PREFLIGHT FAILED (exit $rc) -- see $LOGS/preflight_shell$NUM.log"
        rm -f "$yaml"; return 1
    fi
    local n=$(grep -c 'force: true' "$yaml")
    rm -f "$yaml"
    if [ "$n" -lt 1 ]; then
        print -r -- "PREFLIGHT FAILED: no 'force: true' in the yaml, so contacts would not be enforced"
        return 1
    fi
    print -r -- "preflight ok: flags parse and the yaml carries $n forced contacts"
    return 0
}
if [ -f $RUN/$DESIGN ]; then
    preflight || exit 1
else
    print -r -- "no $DESIGN yet, preflight deferred to after stage 2"
fi

# ---------------------------------------------------------------------------------------------
say "stage 1: reachability sweep"
if [ -f $RUN/$SWEEP ] && [ "$(cut -d, -f1-4 $RUN/$SWEEP | tail -n +2 | sort -u | wc -l | tr -d ' ')" -ge 1 ]; then
    print -r -- "$SWEEP exists with $(cut -d, -f1-4 $RUN/$SWEEP | tail -n +2 | sort -u | wc -l | tr -d ' ') pairs, skipping"
else
    # Pure geometry, no UMA. k=0 is skipped deliberately: settled, 0 successes in 810 attempts.
    $PY code/condense.py sweep $RUN --poses "$POSES" --kmin 1 --kmax 4 --restarts 24 \
        --out $SWEEP > $LOGS/condense_shell$NUM.log 2>&1
fi
NPOSE=$(print -r -- "$POSES" | tr ',' '\n' | wc -l | tr -d ' ')
WANT=$(( NPOSE * (NPOSE - 1) ))
GOT=$(cut -d, -f1-4 $RUN/$SWEEP | tail -n +2 | sort -u | wc -l | tr -d ' ')
print -r -- "sweep: $GOT of $WANT ordered pairs"
if [ "$GOT" -lt "$WANT" ]; then
    print -r -- "STOPPING: incomplete sweep, so the assignment would search a partial graph"
    exit 1
fi

# ---------------------------------------------------------------------------------------------
say "stage 2: ordering and spacer counts"
if [ -f $RUN/$DESIGN ]; then
    print -r -- "$DESIGN exists, skipping"
else
    $PY code/assign.py $RUN --poses "$POSES" --pairs-csv $SWEEP --save $RUN/$DESIGN \
        2>&1 | tee $LOGS/assign_shell$NUM.log
    [ -f $RUN/$DESIGN ] || { print -r -- "STOPPING: assign.py wrote no design"; exit 1; }
    preflight || exit 1
fi
SEQ=$($PY -c "import json;print(json.load(open('$RUN/$DESIGN'))['sequence'])")
print -r -- "design: $SEQ"

# The designed shell, which overlay.py and make_figures.py both need. assign.py does not write it.
# The command word is quoted because zsh rejects `$PY - <<HEREDOC` with "redirection with no command".
"$PY" - <<PYEOF
import sys, json
sys.path.insert(0, 'code')
from design_test import write_design_shell
write_design_shell('$RUN', json.load(open('$RUN/$DESIGN')), 'octinoxate')
PYEOF

# ---------------------------------------------------------------------------------------------
say "stage 3: two ESM2 linker variants"
# fill() in process, so sequences.csv is left alone -- fill_linkers.py's own main() rewrites it, and
# that file is committed data. Sampling, not argmax: with this many masks the 35M checkpoint is nearly
# flat and the argmax collapses every linker to leucine. Expect a large net-charge shift and treat any
# variant whose charge moves far from the design's as suspect.
if [ -f $VARIANTS ]; then
    print -r -- "$(basename $VARIANTS) exists, skipping"
else
    "$PY" - > $UMALOGS/fill_shell$NUM.log 2>&1 <<PYEOF
import sys, json
sys.path.insert(0, 'code')
from fill_linkers import fill
seq = json.load(open('$RUN/$DESIGN'))['sequence']
out = [i['filled'] if isinstance(i, dict) else i for i in fill([seq], variants=2)]
open('$VARIANTS', 'w').write('\n'.join(out) + '\n')
PYEOF
fi
ESM1=$(sed -n 1p $VARIANTS 2>/dev/null)
ESM2=$(sed -n 2p $VARIANTS 2>/dev/null)
chg() { print -r -- "$1" | awk '{n=gsub(/[RK]/,"");m=gsub(/[DE]/,"");printf "%+d", n-m}' }
print -r -- "design charge $(chg $SEQ)   esm1 $(chg ${ESM1:-}) ${ESM1:-<none>}   esm2 $(chg ${ESM2:-}) ${ESM2:-<none>}"

# ---------------------------------------------------------------------------------------------
say "stage 4: twelve folds"
fold_one() {
    local name=$1 seq=$2 top=$3
    local cif=$RUN/boltz/boltz_results_$name/predictions/$name/${name}_model_0.cif
    if [ -f "$cif" ]; then print -r -- "  $name: already folded, skipping"; return 0; fi
    local cmd=($PY code/boltz_hints.py $RUN $DESIGN --name $name --run)
    if [ "$top" = "control" ]; then
        cmd+=(--no-hints)
    else
        cmd+=(--top $top --force --boltz-args=--use_potentials)
    fi
    [ -n "$seq" ] && cmd+=(--sequence $seq)
    print -r -- "  folding $name"
    "${cmd[@]}" >> $LOGS/folds_shell$NUM.log 2>&1 || print -r -- "  $name: FAILED, see folds_shell$NUM.log"
}
# The f4/f8/f12 suffixes matter: make_figures.py sorts the ladder by them.
for set_tag seq in orig "" esm1 "$ESM1" esm2 "$ESM2"; do
    [ "$set_tag" != "orig" ] && [ -z "$seq" ] && continue
    for pair in control:control f4:4 f8:8 f12:12; do
        fold_one ${PREFIX}_${set_tag}_${pair%%:*} "$seq" ${pair##*:}
    done
done

ALL=(); TODO=()
for set_tag in orig esm1 esm2; do
    for tag in control f4 f8 f12; do
        name=${PREFIX}_${set_tag}_${tag}
        [ -f $RUN/boltz/boltz_results_$name/predictions/$name/${name}_model_0.cif ] || continue
        ALL+=($name)
        [ -f $RUN/boltz/partial_$name.json ] || TODO+=($name)
    done
done
print -r -- "${#ALL} folds present, ${#TODO} still to score"
[ ${#ALL} -eq 0 ] && { print -r -- "STOPPING: nothing folded"; exit 1; }

# ---------------------------------------------------------------------------------------------
say "stage 5: geometry check -- free, and the primary measure"
$PY code/check_fold.py $RUN --out $FOLDCHECK 2>&1 | tee $LOGS/check_fold_shell$NUM.log

say "stage 5b: designed positions reproduced -- also free"
# overlay.py matches a fold to sequences/<foldname>.xyz, so the design's shell is copied in under each
# fold's name. The poses are identical across the twelve; only sequence and constraint level differ.
for name in $ALL; do
    [ -f $RUN/sequences/$name.xyz ] || cp $RUN/sequences/$SEQ.xyz $RUN/sequences/$name.xyz
done
$PY code/overlay.py $RUN --match ${PREFIX}_ --out overlay_shell$NUM.csv \
    2>&1 | tee $LOGS/overlay_shell$NUM.log

# ---------------------------------------------------------------------------------------------
say "stage 6: score -- THE BOTTLENECK, 40-50 min each on CPU"
# strain_peptide is left out deliberately: a free peptide collapses in vacuum, so that term measures
# collapse rather than strain. If a CUDA device is around, binding_energy.py finds it and this takes
# minutes per structure instead -- see code/modal_score.py for renting one.
#
# SKIP_SCORING=1 stops here and writes the structure list, for scoring on a GPU elsewhere. Everything
# before this stage is cheap; this stage is the only expensive one. Default is unchanged: score here.
if [ "${SKIP_SCORING:-0}" = "1" ]; then
    mkdir -p $RUN/modal
    print -rl -- $TODO > $RUN/modal/shell$NUM.txt
    print -r -- "SKIP_SCORING=1: ${#TODO} structures left unscored, listed in $RUN/modal/shell$NUM.txt"
    print -r -- ""
    print -r -- "  modal run code/modal_score.py --structures-file $RUN/modal/shell$NUM.txt \\"
    print -r -- "      --out $BINDING --gpu L4"
    print -r -- ""
    print -r -- "then re-run this script to pick up from stage 5 with the scores in place."
    exit 0
elif [ ${#TODO} -eq 0 ]; then
    print -r -- "every fold already scored, skipping"
else
    # interaction only. strain_ligand is stage 5c's job now, against the shared reference, and
    # leaving it out here also removes the per-structure free-ligand relaxation from the slow stage.
    $PY code/binding_energy.py $RUN --terms interaction \
        --out $BINDING --structures "${(j:,:)TODO}" 2>&1 | tee -a $UMALOGS/score_shell$NUM.log
fi

say "stage 6b: ligand strain against one shared reference -- cheap"
# Runs after scoring because the bound state is the ligand as it exists in the complex,
# hydrogens included, and binding_energy.py is what writes that geometry
# (structures/<name>_ligand_bound.xyz). Re-placing those hydrogens on the isolated ligand
# would change the initial state to one the ligand never occupies and relax away part of
# the energy being measured.
#
# The reference is a property of the ligand, not of any fold, so it is computed once per run and
# reused. Relaxing each complex's own bound pose instead measures every structure against a different
# local minimum, which is what made the old strain numbers depend on the force cutoff, the step cap
# and pdbfixer's non-deterministic hydrogens. See HANDOFF.md section 10.
if [ -f $RUN/ligand_reference.json ]; then
    print -r -- "ligand_reference.json exists, reusing it"
else
    $PY code/ligand_reference.py $RUN 2>&1 | tee $UMALOGS/ligand_reference.log
    [ -f $RUN/ligand_reference.json ] || { print -r -- "STOPPING: no reference written"; exit 1; }
fi
# Needs only the Boltz heavy atoms plus RDKit hydrogens, so this runs without any complex relaxation.
$PY code/strain_global.py $RUN --match ${PREFIX}_ --out strain_shell$NUM.csv --save-structures \
    2>&1 | tee $UMALOGS/strain_shell$NUM.log

say "stage 7: figures, into their own directory"
$PY code/make_figures.py $RUN --fig-dir $FIGDIR --match ${PREFIX}_ --design-sequence "$SEQ" \
    2>&1 | tee $LOGS/make_figures_shell$NUM.log

say "shell $NUM complete"
print -r -- "next: correlate.py $RUN --match ${PREFIX}_ --group '${PREFIX}_(orig|esm1|esm2)_'"
