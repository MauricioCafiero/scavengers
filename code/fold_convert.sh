#!/bin/zsh
# Convert the best-ranked cofold cif for every peptide that has one but is missing its pdb.
# Complements fold_batch.sh's conversion step: fold_batch.sh passed the out-path as a positional
# argument, which typer rejects, so every OF3 conversion failed while the (expensive, complete)
# modal run and its cifs were fine -- this script does the conversion pass afterwards instead of
# re-running modal. Resumable: skips entries whose pdb already exists.
#
#   LANE=openfold3 code/fold_convert.sh   # openfold3 (outputs/<job>/) or rf3 (outputs/rf3_<job>/)
set -u
REPO=${0:A:h:h}
cd $REPO
FOLD=${FOLD_ROOT:-$HOME/python_mac/fold}
LANE=${LANE:-openfold3}
[[ $LANE == openfold3 || $LANE == rf3 ]] || { print -u2 "LANE must be openfold3 or rf3"; exit 1; }
prefix=""; [[ $LANE == rf3 ]] && prefix=rf3_

for f in runs/*/folds/inputs/*.txt(N); do
    mol=${f:h:h:h:t}
    name=${f:t:r}
    out=$REPO/runs/$mol/folds/$LANE/$name.pdb
    [[ -e $out && -s $out ]] && { echo "skip $mol/$name"; continue; }
    cif=$( cd $FOLD && uv run python -c "
from fold.results import find_best_cif
print(find_best_cif('outputs/$prefix$name'))" 2>/dev/null | tail -1 )
    if [[ -n $cif && -f $FOLD/$cif ]]; then
        mkdir -p runs/$mol/folds/$LANE
        ( cd $FOLD && uv run python -m fold.analyze cif-to-pdb $cif --out-path $out ) >/dev/null 2>&1
        [[ -s $out ]] && echo "converted $mol/$name" || echo "FAILED $mol/$name (cif $cif)"
    else
        echo "NO CIF for $mol/$name"
    fi
done