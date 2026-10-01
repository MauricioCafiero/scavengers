#!/bin/zsh
# Row 5: s2_esm2_control, both docked poses.
# Vina cell if the clock allows.
#
# bg33_3 is chosen because it is the case most likely to break the pattern rather than confirm it. Five
# of the six cells run so far end below their co-folded value, which tempts the conclusion that docking
# then running is a worse route than co-folding. bg33_3's co-folded pose is the worst-retained structure
# anywhere in this project -- 41% residence, 18 release episodes -- so what docking has to beat is
# genuinely bad. If a docked start wins anywhere, here is where. It also completes the BoltzGen pair and
# extends the range from -19.66 down to -11.86, where the three rows already done are the top three
# binders and nothing has been tested at the weak end.
#
# s2_esm2_control tests the hypothesis the fourth row produced: that whether a docked start beats
# co-folding depends on how good the co-folded pose was. Boltz won on its two fully-retained structures
# (s3_orig_f12 and bg33_4, both 100% residence, zero releases) and lost on bg33_3, whose co-folded pose
# holds the ligand only 41% of the time across 18 releases. s2_esm2_control sits between them at 59.6%
# and one release -- the lowest retention left unrun -- so it is where the hypothesis makes its sharpest
# prediction: a docked start should win here too, and by less than it did on bg33_3.
#
# GNINA's pick is pose 6 by CNNaffinity. No deadline tonight, so the cutoff is nominal.
# --box-shape dodecahedron is NOT optional: `omd build` defaults to cube, and every co-folded baseline
# these legs are compared against was built dodecahedral. The seven legs run on 2026-10-01 inherited the
# cube default -- same physics, but ~30% more waters and about $1.30 of avoidable GPU across them, worst
# on an elongated solute (s2_esm2_control: 33,057 atoms cubic against 23,142 dodecahedral).
set -u
# Paths are derived, not hardcoded, so the repo can be cloned anywhere. Override the externals with
# OMD_ENV and BOLTZGEN_ROOT if they live elsewhere on this machine.
REPO=${0:A:h:h}
cd $REPO
OMD_ENV=${OMD_ENV:-$HOME/miniforge3/envs/openmm-md}
BOLTZGEN_ROOT=${BOLTZGEN_ROOT:-$BOLTZGEN_ROOT}
for need in $OMD_ENV/bin/omd $OMD_ENV/bin/python; do
    [[ -x $need ]] || { echo "missing $need -- set OMD_ENV to the openmm-md conda env"; exit 1; }
done
ROOT=$REPO
OMD=$OMD_ENV/bin/omd
PY=.venv/bin/python
STEPS=10000000
CUTOFF=2359               # no deadline tonight

caffeinate -w $$ &
date
echo "### row 5: s2_esm2_control, Vina p1 and GNINA p6"

# name              structure        pose  receptor
CELLS=(
  "s2_esm2_control_dock1 s2_esm2_control 1 runs/octinoxate/md/s2_esm2_control/protein_fixed.pdb"
  "s2_esm2_control_dock6 s2_esm2_control 6 runs/octinoxate/md/s2_esm2_control/protein_fixed.pdb"
)

for cell in $CELLS; do
    set -- ${=cell}
    N=$1 STRUCT=$2 POSE=$3 RECEPTOR=$4
    M=runs/octinoxate/md/$N
    P=$M/prod_L1_modal
    echo; echo "=== $N  ($STRUCT pose $POSE) ==="; date

    [[ -f $P/mmgbsa/FINAL_RESULTS_MMPBSA.dat ]] && { echo "  done already"; continue; }

    NOW=$(date +%H%M)
    if (( NOW > CUTOFF )); then
        echo "  past cutoff $CUTOFF (now $NOW) -- not starting $N"
        continue
    fi

    # --- prep (local, free) ---
    if [[ ! -f $M/system/system.xml ]]; then
        mkdir -p $M
        [[ -f $M/protein_fixed.pdb ]] || cp $RECEPTOR $M/protein_fixed.pdb
        [[ -f $M/${N}_ligand.sdf ]] || $PY code/dock_pose_to_sdf.py $STRUCT --pose $POSE \
            --out $M/${N}_ligand.sdf || { echo "PREP_FAILED $N"; exit 1; }
        [[ -f $M/ligand_prepped.sdf ]] || $OMD prep-ligand --sdf $M/${N}_ligand.sdf \
            --out $M/ligand_prepped.sdf || { echo "PREP_FAILED $N"; exit 1; }
        $OMD build --protein $M/protein_fixed.pdb --ligand $M/ligand_prepped.sdf \
                   --out-dir $M/system --no-auto-cofactors --box-shape dodecahedron || { echo "BUILD_FAILED $N"; exit 1; }
    fi

    # --- dynamics on Modal ---
    if [[ ! -f $P/traj_wrapped.xtc ]]; then
        modal run code/modal_md.py::push --structure $N || { echo "PUSH_FAILED $N"; exit 1; }
        modal run code/modal_md.py::produce --structure $N --steps $STEPS \
            || { echo "PRODUCE_FAILED $N"; exit 1; }
        date
    fi
    [[ -f $P/traj_wrapped.xtc ]] || { echo "NO_TRAJECTORY $N"; exit 1; }

    # --- MM/GBSA local, from inside the leg dir (MMPBSA scatters GBs of scratch into cwd) ---
    ( cd $P && $OMD mmgbsa --protein $ROOT/$M/protein_fixed.pdb \
                           --ligand $ROOT/$M/ligand_prepped.sdf \
                           --traj traj_wrapped.xtc --topology traj_wrapped.pdb \
                           --out-dir mmgbsa --no-auto-cofactors --run ) \
        || { echo "MMGBSA_FAILED $N"; exit 1; }
    date

    # QC that the pose matrix taught us to print: how far the ligand moved before production
    $OMD_ENV/bin/python - "$M" <<'PYQC' || true
import sys, numpy as np, mdtraj as md
M=sys.argv[1]
try:
    s=md.load(f"{M}/system/complex.pdb"); p=md.load(f"{M}/prod_L1_modal/traj_wrapped.pdb")
    f=lambda t:(np.array([a.index for a in t.topology.atoms if a.residue.is_protein]),
                np.array([a.index for a in t.topology.atoms if not a.residue.is_protein
                          and a.residue.name not in ("HOH","WAT","NA","CL")]))
    ps,ls=f(s); pp,lp=f(p); S=s.xyz[0]*10; P=p.xyz[0]*10
    A,B=S[ps]-S[ps].mean(0), P[pp]-P[pp].mean(0)
    U,_,Vt=np.linalg.svd(A.T@B); R=Vt.T@np.diag([1,1,np.sign(np.linalg.det(Vt.T@U.T))])@U.T
    Sl=(S[ls]-S[ps].mean(0))@R.T; Pl=P[lp]-P[pp].mean(0)
    print(f"[qc] ligand RMSD input -> production start: {np.sqrt(((Sl-Pl)**2).sum(-1).mean()):.2f} A")
except Exception as e:
    print(f"[qc] skipped: {e}")
PYQC

    echo "--- $N result ---"
    awk '/DELTAS|Differences/{d=1} d&&/DELTA TOTAL/{print; exit}' $P/mmgbsa/FINAL_RESULTS_MMPBSA.dat
done

echo; echo "ROW5_DONE"; date
