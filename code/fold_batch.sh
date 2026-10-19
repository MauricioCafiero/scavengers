#!/bin/zsh
# Fold-comparison batch driver: fold every peptide in runs/<mol>/folds/inputs/ with one model and
# write the structure beside that molecule's manifest. See FOLDS.md for what this campaign is.
#
#   LANE=esmfold code/fold_batch.sh                 # ESMFold: ligand-free fold, sequence only
#   LANE=openfold3 code/fold_batch.sh               # OpenFold3: protein+ligand cofold on the A10G
#   LANE=of3apo code/fold_batch.sh                  # OpenFold3: protein-only (apo) query
#   LANE=rf3 code/fold_batch.sh                     # RosettaFold3 cofold, second opinion
#   LANE=boltzapo code/fold_batch.sh                # Boltz-2 apo, local MPS (no Modal)
#   LANE=openfold3p code/fold_batch.sh              # OpenFold3 pocket-guided cofold (inputs_pocket/)
#   MOL=oxybenzone LANE=esmfold code/fold_batch.sh  # restrict to one molecule's peptide set
#   ONLY="s3_orig_f12" LANE=openfold3p code/fold_batch.sh   # restrict to named entries (smoke runs)
#
# Resumable: an entry whose output already exists is skipped, so re-running after an interruption
# continues where it stopped. Runs are sequential because one Modal run at a time is the house rule;
# each invocation carries its own caffeinate -i -s guard (-s since a dark sleep dropped a -i-only
# Modal client on 2026-10-10, CLAUDE.md "Modal") -- applied here directly rather than via
# code/modal_run.sh, because the fold repo's apps must run under the fold repo's `uv run modal` for
# its package to be importable (modal_run.sh hardcodes peptidebuilder's modal_md.py invocation).
#
# ESMFold writes pdb directly (--out-path). The cofold lanes run with --input-file (the SMILES
# rides along) -- OpenFold3 drops outputs/<job>/ under the fold repo and names one model.cif per
# sample; RF3 drops outputs/rf3_<job>/ and promotes its best model. The driver converts the
# best-ranked cif (fold.results.find_best_cif) to pdb in the fold lane dir and leaves the raw cifs
# in the fold repo; a later re-run of the metrics script reads the converted pdb.
# of3apo submits the SEQUENCE: line alone (no SMILES) as a single-chain apo query and names its
# job <pep>_apo so it cannot collide with the cofold lane's outputs. Requires the fold repo's
# guards relaxed to allow one-chain queries (done 2026-10-09). rf3apo is the same for RF3 -- the
# app itself prefixes outputs with rf3_, so its job dir is outputs/rf3_<pep>_apo.
# boltzapo folds the SEQUENCE: line alone with the same local Boltz-2 that made the co-folds
# (boltz_env.run_boltz, msa: empty as in runs/<mol>/boltz/*.yaml, no ligand, no hints). It runs on
# this Mac, so no Modal; the yaml, log and Boltz's raw output stay under boltzapo/raw/.
# openfold3p is the OpenFold3 cofold with a pocket_constraint: it reads runs/<mol>/folds/inputs_pocket/
# (written by code/pocket_inputs.py from each forced Boltz yaml -- same SEQUENCE/SMILES as inputs/,
# plus POCKET/POCKET_MAX_DISTANCE) and names its job <pep>_pocket so it cannot collide with the
# unguided lane's outputs/<pep>. Pocket-GUIDED, not forced: see NEXT_STEPS.md "Part 2".
set -u
REPO=${0:A:h:h}
cd $REPO
FOLD=${FOLD_ROOT:-$HOME/python_mac/fold}
LANE=${LANE:-esmfold}

case $LANE in
    esmfold|openfold3|openfold3p|rf3|of3apo|rf3apo|boltzapo) ;;
    *) print -u2 "LANE must be esmfold, openfold3, openfold3p, rf3, of3apo, rf3apo or boltzapo"; exit 1 ;;
esac
INDIR=inputs; [[ $LANE == openfold3p ]] && INDIR=inputs_pocket
MOL=${MOL:-*}
ONLY=${ONLY:-}

for f in runs/${~MOL}/folds/$INDIR/*.txt(N); do
    mol=${f:h:h:h:t}
    name=${f:t:r}
    [[ -n $ONLY && " $ONLY " != *" $name "* ]] && continue
    log=$REPO/runs/$mol/folds/$LANE.batch.log
    echo "=== $LANE $mol/$name ($(date '+%H:%M:%S'))" >> $log
    if [[ $LANE == esmfold ]]; then
        out=$REPO/runs/$mol/folds/esmfold/$name.pdb
        [[ -e $out ]] && { echo "skip $mol/$name (exists)" >> $log; continue; }
        seq=$(awk '/^SEQUENCE:/ {print $2}' $f)
        ( cd $FOLD && caffeinate -i -s uv run modal run src/fold/esmfold_app.py::main \
              --sequence $seq --out-path $out ) >> $log 2>&1
        rc=$?
        echo "--- $LANE $mol/$name exit $rc ($(date '+%H:%M:%S'))" >> $log
        [[ $rc -eq 0 && -s $out ]] || echo "FAILED $mol/$name (exit $rc) -- see $log" | tee -a $log
    elif [[ $LANE == boltzapo ]]; then
        out=$REPO/runs/$mol/folds/boltzapo/$name.pdb
        [[ -e $out ]] && { echo "skip $mol/$name (exists)" >> $log; continue; }
        raw=$REPO/runs/$mol/folds/boltzapo/raw
        mkdir -p $raw
        seq=$(awk '/^SEQUENCE:/ {print $2}' $f)
        printf 'version: 1\nsequences:\n  - protein:\n      id: A\n      sequence: %s\n      msa: empty\n' \
            $seq > $raw/$name.yaml
        caffeinate -i .venv/bin/python -c "
import sys; sys.path.insert(0, 'code')
from boltz_env import run_boltz
sys.exit(run_boltz('$raw/$name.yaml', '$raw', '$raw/$name.log'))" >> $log 2>&1
        rc=$?
        cif=$raw/boltz_results_$name/predictions/$name/${name}_model_0.cif
        if [[ -f $cif ]]; then
            ( cd $FOLD && uv run python -m fold.analyze cif-to-pdb $cif --out-path $out ) >> $log 2>&1
            echo "--- $LANE $mol/$name exit $rc, converted ($(date '+%H:%M:%S'))" >> $log
            [[ -s $out ]] || echo "CONVERT FAILED $mol/$name -- cif found but pdb empty, see $log" \
                | tee -a $log
        else
            echo "NO CIF for $mol/$name (boltz exit $rc) -- see $raw/$name.log" | tee -a $log
        fi
    elif [[ $LANE == of3apo || $LANE == rf3apo ]]; then
        out=$REPO/runs/$mol/folds/$LANE/$name.pdb
        [[ -e $out ]] && { echo "skip $mol/$name (exists)" >> $log; continue; }
        [[ $LANE == rf3apo ]] && APP=src/fold/rf3_app.py::main || APP=src/fold/app.py::main
        [[ $LANE == rf3apo ]] && dirjob=rf3_${name}_apo || dirjob=${name}_apo
        seq=$(awk '/^SEQUENCE:/ {print $2}' $f)
        ( cd $FOLD && caffeinate -i -s uv run modal run $APP \
              --sequence $seq --job-name ${name}_apo ) >> $log 2>&1
        rc=$?
        cif=$( cd $FOLD && uv run python -c "
from fold.results import find_best_cif
print(find_best_cif('outputs/$dirjob'))" 2>/dev/null | tail -1 )
        if [[ -n $cif && -f $FOLD/$cif ]]; then
            mkdir -p runs/$mol/folds/$LANE
            ( cd $FOLD && uv run python -m fold.analyze cif-to-pdb $cif --out-path $out ) >> $log 2>&1
            echo "--- $LANE $mol/$name exit $rc, converted ($(date '+%H:%M:%S'))" >> $log
            [[ -s $out ]] || echo "CONVERT FAILED $mol/$name -- cif found but pdb empty, see $log" \
                | tee -a $log
        else
            echo "NO CIF for $mol/$name (modal exit $rc) -- see $log" | tee -a $log
        fi
    else
        out=$REPO/runs/$mol/folds/$LANE/$name.pdb
        [[ -e $out ]] && { echo "skip $mol/$name (exists)" >> $log; continue; }
        APP=src/fold/app.py::main; prefix=""; job=$name
        [[ $LANE == rf3 ]] && { APP=src/fold/rf3_app.py::main; prefix=rf3_; }
        [[ $LANE == openfold3p ]] && job=${name}_pocket
        ( cd $FOLD && caffeinate -i -s uv run modal run $APP --input-file $REPO/$f --job-name $job ) \
            >> $log 2>&1
        rc=$?
        # best-ranked cif -> pdb beside the manifest; raw cifs stay in the fold repo's outputs/
        cif=$( cd $FOLD && uv run python -c "
from fold.results import find_best_cif
print(find_best_cif('outputs/$prefix$job'))" 2>/dev/null | tail -1 )
        if [[ -n $cif && -f $FOLD/$cif ]]; then
            mkdir -p runs/$mol/folds/$LANE
            ( cd $FOLD && uv run python -m fold.analyze cif-to-pdb $cif --out-path $out ) >> $log 2>&1
            echo "--- $LANE $mol/$name exit $rc, converted ($(date '+%H:%M:%S'))" >> $log
            [[ -s $out ]] || echo "CONVERT FAILED $mol/$name -- cif found but pdb empty, see $log" \
                | tee -a $log
        else
            echo "NO CIF for $mol/$name (modal exit $rc) -- see $log" | tee -a $log
        fi
    fi
    sleep 5
done
echo "batch $LANE complete: $(ls runs/*/folds/*/*.pdb(N) | wc -l) pdb files total"