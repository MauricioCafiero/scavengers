#!/bin/bash
# Rescore ten poses per structure with GNINA. Run this on the Linux machine, from inside this folder.
#
#   ./run_gnina.sh              # heavy-atom receptor, the one Vina scored against
#   ./run_gnina.sh receptor_h   # PDBFixer receptor with hydrogens, the one MD used
#   ./run_gnina.sh receptor minimize   # local minimisation before scoring
#
# --score_only evaluates the pose exactly as given, with no search and no minimisation, which is what
# makes it a rescoring of Vina's output rather than a fresh docking. That matters here: the comparison
# is only meaningful if GNINA is shown the same coordinates Vina produced.
#
# --minimize is offered as a second pass because CNN scores are sensitive to small clashes, and these
# poses came out of a rigid-receptor search that never relaxed them. If the two passes disagree about
# the ranking, the clash sensitivity is the likely reason and worth knowing. It moves the coordinates,
# so its output is written separately and must not be mixed with the score_only numbers.
set -u
REC=${1:-receptor}
MODE=${2:-score_only}
OUT=results_${REC}_${MODE}
mkdir -p "$OUT"

command -v gnina >/dev/null || { echo "gnina is not on PATH"; exit 1; }
gnina --version | head -2 | tee "$OUT/gnina_version.txt"

if [[ "$MODE" == "minimize" ]]; then FLAG="--minimize"; else FLAG="--score_only"; fi

for d in */; do
    n="${d%/}"
    [[ -f "$n/poses_with_reference.sdf" ]] || continue
    [[ -f "$n/${REC}.pdb" ]] || { echo "SKIP $n: no ${REC}.pdb"; continue; }
    echo "################ $n"
    # -o writes an SDF carrying the scores as tags (CNNscore, CNNaffinity, minimizedAffinity), which
    # collect_gnina.py reads. Parsing stdout instead would be fragile across GNINA versions.
    gnina $FLAG -r "$n/${REC}.pdb" -l "$n/poses_with_reference.sdf" \
          -o "$OUT/${n}_scored.sdf" --seed 42 2>&1 | tee "$OUT/${n}.log"
    echo "EXIT[$n]=$?"
done
echo "GNINA_DONE -> $OUT"
echo
echo "Bring back the whole $OUT directory. Then, on either machine:"
echo "    python collect_gnina.py $OUT"
