#!/bin/zsh
# Stage a replicate of an already-run leg: same system, different seed.
#   stage_replicate.sh <mol> <leg> [n]        n = 2 (default) for _r2, 3 for _r3, ...
#
# A replicate is the SAME system.xml run again. Nothing here sets a seed: racc_run.py builds a
# fresh integrator on the cluster, so the initial velocities differ between runs and that is the
# whole difference. Copying the inputs byte-for-byte is therefore the entire job, and the check
# below (cmp on system.xml) is what makes the run a replicate rather than a different experiment.
#
# Only the four inputs are copied -- system/, protein_fixed.pdb and the two ligand sdfs. Everything
# else in a finished leg (prod_20ns, first_*ns, md_contacts.csv, tail.log) is output and must NOT
# come along: racc_drain.sh treats the mere presence of prod_20ns/traj_wrapped.xtc as "leg
# finished, tail it", so a stray copy of the parent's trajectory would be recorded as the
# replicate's own result.
set -eu
REPO=${0:A:h:h}
mol=$1; leg=$2; n=${3:-2}
src=$REPO/runs/$mol/md/$leg
dst=$REPO/runs/$mol/md/${leg}_r$n

[[ -d $src ]] || { print -u2 "no such leg: $src"; exit 1; }
[[ -f $src/system/system.xml ]] || { print -u2 "no system.xml in $src -- cannot replicate"; exit 1; }
[[ -e $dst ]] && { print -u2 "$dst already exists; refusing to overwrite"; exit 1; }

mkdir -p $dst
cp -R $src/system $dst/system
for f in protein_fixed.pdb ligand_prepped.sdf ${leg}_ligand.sdf; do
  [[ -f $src/$f ]] && cp $src/$f $dst/$f
done

cmp -s $src/system/system.xml $dst/system/system.xml \
  || { print -u2 "system.xml differs after copy -- not a replicate"; exit 1; }
print -r -- "staged ${leg}_r$n  ($(grep -c '<Particle ' $dst/system/system.xml) particles, system.xml identical)"
