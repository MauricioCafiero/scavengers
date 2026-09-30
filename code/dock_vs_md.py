#!/usr/bin/env python3
"""Which ligand pose does the dynamics actually prefer -- Boltz's or Vina's?

The docking in `vina_redock.py` disagrees with Boltz about where the ligand sits, by 3-11 A, mostly
by turning it end for end. That leaves the question of which pose is right, and neither the Vina
score nor the wrapping metric can answer it: Vina's poses were *selected* to pack tightly against a
rigid receptor, so their packing is not independent evidence, and wrapping saturates whichever way
round the ligand lies.

The MD trajectories can answer it, because they were started from the Boltz pose in explicit solvent
with a flexible peptide and left to run. If the ligand is still near where Boltz put it in the medoid
frame, the Boltz pose is dynamically stable and the Vina poses are artefacts of a rigid receptor and
an empirical score. If it has drifted towards one of the Vina poses, the opposite.

Unlike `dock_compare.py`, this one *does* need a superposition. The MD frame has its own coordinate
frame -- the peptide tumbles and the box is re-wrapped -- so the peptide is superimposed onto the
docking receptor first and that transform is applied to the ligand. `peptide_rmsd_A` reports how well
that superposition worked: it is also the peptide's own conformational drift over the run, and when it
is large the ligand comparison is correspondingly blurred, so it is reported rather than hidden.

Atom order is relied upon, not perceived: the MD prep renames ligand atoms (C38 -> C1x) but preserves
their order, which is checked here against the reference's elements and bond graph before anything is
compared, and the run aborts if it does not hold.

    ~/python_mac/dock_assist/dock-env/bin/python code/dock_vs_md.py --all
"""
import argparse
import csv
import os
import sys

import numpy as np
from rdkit import Chem

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.expanduser("~/python_mac/dock_assist/code"))

from dock_compare import kabsch  # noqa: E402
from vina_redock import ref_order_poses  # noqa: E402

# Which production leg holds each structure's medoid frame. orig_f12 was run locally before the
# Modal legs existed, hence the odd one out.
LEGS = ["prod_L1_modal", "prod_20ns"]


def read_pdb(path):
    """({(resseq, name): xyz} for the peptide, ligand xyz in file order, ligand elements)."""
    pep, lig_xyz, lig_el = {}, [], []
    for line in open(path):
        if not line.startswith(("ATOM", "HETATM")) or line[76:78].strip() == "H":
            continue
        xyz = [float(line[30:38]), float(line[38:46]), float(line[46:54])]
        if line[17:20].strip() == "LIG":
            lig_xyz.append(xyz)
            lig_el.append(line[76:78].strip())
        else:
            pep[(int(line[22:26]), line[12:16].strip())] = xyz
    return pep, np.array(lig_xyz), lig_el


def compare(structure, dock_root, md_roots):
    work = os.path.join(dock_root, structure)
    ref_mol = Chem.MolFromMolFile(os.path.join(work, "ligand_ref.sdf"), removeHs=True)
    ref = ref_mol.GetConformer().GetPositions()
    poses = ref_order_poses(os.path.join(work, "poses.pdbqt"),
                            os.path.join(work, "ligand.pdbqt"), ref)
    rec_pep, _, _ = read_pdb(os.path.join(work, f"{structure}_protein.pdb"))

    medoid = None
    for root in md_roots:
        for leg in LEGS:
            p = os.path.join(root, structure, leg, "frame_medoid.pdb")
            if os.path.exists(p):
                medoid, leg_used = p, leg
                break
        if medoid:
            break
    if medoid is None:
        return None

    md_pep, md_lig, md_el = read_pdb(medoid)
    if len(md_lig) != len(ref):
        raise SystemExit(f"{structure}: medoid has {len(md_lig)} ligand heavy atoms, "
                         f"reference has {len(ref)}")
    # The MD prep renames ligand atoms but keeps their order; verify before relying on it.
    if md_el != [a.GetSymbol() for a in ref_mol.GetAtoms()]:
        raise SystemExit(f"{structure}: ligand element order differs between medoid and reference")

    def adj(X):
        d = np.linalg.norm(X[:, None] - X[None], axis=2)
        return (d < 1.85) & (d > 0)
    if not (adj(md_lig) == adj(ref)).all():
        raise SystemExit(f"{structure}: ligand bond graph differs between medoid and reference; "
                         "atom order cannot be assumed")

    shared = sorted(set(rec_pep) & set(md_pep))
    if len(shared) < 20:
        raise SystemExit(f"{structure}: only {len(shared)} peptide atoms matched by name")
    A = np.array([md_pep[k] for k in shared])
    B = np.array([rec_pep[k] for k in shared])
    R, t, pep_rmsd = kabsch(A, B)
    md_in_dock_frame = (R @ md_lig.T).T + t

    autos = ref_mol.GetSubstructMatches(ref_mol, uniquify=False, useChirality=False) or \
        [tuple(range(len(ref)))]

    def rmsd(P, Q):
        return min(float(np.sqrt(((P[list(a)] - Q) ** 2).sum(axis=1).mean())) for a in autos)

    to_boltz = rmsd(md_in_dock_frame, ref)
    to_poses = [rmsd(md_in_dock_frame, p) for p in poses]
    nearest = int(np.argmin(to_poses))
    return {
        "structure": structure,
        "leg": leg_used,
        "peptide_atoms_matched": len(shared),
        "peptide_rmsd_A": round(pep_rmsd, 2),
        "md_ligand_vs_boltz_A": round(to_boltz, 2),
        "md_ligand_vs_nearest_pose_A": round(to_poses[nearest], 2),
        "nearest_pose": nearest + 1,
        "md_ligand_vs_pose1_A": round(to_poses[0], 2),
        "verdict": ("stayed at Boltz pose" if to_boltz < to_poses[nearest]
                    else f"moved towards Vina pose {nearest + 1}"),
    }


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("structures", nargs="*")
    ap.add_argument("--system", default="octinoxate")
    ap.add_argument("--all", action="store_true")
    ap.add_argument("--md-root", action="append", default=None, metavar="DIR",
                    help="extra directory of MD runs to look for medoid frames in, repeatable")
    a = ap.parse_args(argv)

    dock_root = os.path.join(REPO, "runs", a.system, "dock")
    md_roots = [os.path.join(REPO, "runs", a.system, "md")]
    md_roots += [os.path.expanduser(r) for r in (a.md_root or [])]
    names = list(a.structures)
    if a.all:
        names += [d for d in sorted(os.listdir(dock_root))
                  if os.path.exists(os.path.join(dock_root, d, "poses.pdbqt"))]
        names = sorted(set(names))
    if not names:
        ap.error("name at least one structure, or pass --all")

    rows = [r for r in (compare(n, dock_root, md_roots) for n in names) if r]
    out = os.path.join(dock_root, "md_vs_dock.csv")
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print(f"wrote {out}\n")
    print(f"{'structure':<22}{'pepRMSD':>9}{'MD vs Boltz':>13}{'MD vs best pose':>17}"
          f"{'(which)':>9}   verdict")
    for r in rows:
        print(f"{r['structure']:<22}{r['peptide_rmsd_A']:>9.2f}{r['md_ligand_vs_boltz_A']:>13.2f}"
              f"{r['md_ligand_vs_nearest_pose_A']:>17.2f}{r['nearest_pose']:>9}   {r['verdict']}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
