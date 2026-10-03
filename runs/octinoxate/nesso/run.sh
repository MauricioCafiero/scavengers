#!/bin/zsh
caffeinate -w $$ &
export HF_HUB_DISABLE_XET=1
root=/Users/cafierom/python_mac/nesso
out=/Users/cafierom/python_mac/peptidebuilder/runs/octinoxate/nesso
echo "START $(date)"
"$root/.venv/bin/nesso" predict "$out/inputs" --out_dir "$out" --cache "$root/.cache" --seed 0
echo "END $(date)"
