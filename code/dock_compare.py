#!/usr/bin/env python3
"""How does a docked pose actually differ from the Boltz pose it was meant to reproduce?

An RMSD says the two disagree but not how, and 8 A of RMSD can mean the ligand slid along a groove,
flipped end for end in the same cavity, or uncoiled its alkyl tail without moving at all. Those have
different consequences for the design, so this separates them.

No superposition of the receptors is involved or needed: a docked pose and its reference already sit
in one frame, because the receptor PDB Vina was given is the same file the reference ligand came
from, coordinates untouched. The whole difference is therefore the ligand's own, and it splits into
three parts that are measured separately:

* **rigid-body** -- the Kabsch transform taking the reference onto the pose. `translation_A` is how
  far the centroid moved, `rotation_deg` how far it turned about that transform's axis. A large
  rotation with a small translation is a reorientation in place; the reverse is a slide.
* **internal** -- `aligned_rmsd_A`, the RMSD left *after* that transform is removed. This is the
  ligand changing its own shape, which for this molecule means the 9 rotatable torsions, mostly in
  the 2-ethylhexyl tail. Compare it against `inplace_rmsd_A` (the total): if the two are close, the
  pose difference is conformational and the ligand barely moved; if `aligned_rmsd_A` is small while
  the total is large, the ligand moved as a rigid body and kept its shape.
* **orientation of the two ends** -- octinoxate is amphiphilic, a methoxyphenyl head and a branched
  alkyl tail. `headtail_deg` is the angle between the head->tail vectors of the two poses, so ~180
  means the ligand is sitting the other way round in the pocket. That is the difference that matters
  chemically and the one an RMSD hides most thoroughly.

Then, because a pose can score well in a different place entirely, `contact_jaccard` compares which
*residues* line the ligand in each pose (heavy-atom contact within --contact A). A high Jaccard with
a high RMSD is one cavity holding the ligand several ways; a low Jaccard is a different site.

Reads what code/vina_redock.py wrote. Writes pose_differences.csv beside it, and a PyMOL script per
structure so the poses can also just be looked at.

    ~/python_mac/dock_assist/dock-env/bin/python code/dock_compare.py
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
os.environ["PATH"] = (os.path.expanduser("~/python_mac/dock_assist/dock-env/bin")
                      + os.pathsep + os.environ["PATH"])

from vina_redock import ref_order_poses  # noqa: E402


def kabsch(P, Q):
    """Rotation R and translation t with R @ p + t ~= q, both centred first.

    Returns (R, t, rmsd_after). R is a proper rotation: the sign correction on the smallest singular
    value keeps det(R) = +1, so a reflection is never reported as a large rotation.
    """
    pc, qc = P.mean(axis=0), Q.mean(axis=0)
    H = (P - pc).T @ (Q - qc)
    U, S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1.0, 1.0, d]) @ U.T
    t = qc - R @ pc
    resid = (R @ P.T).T + t - Q
    return R, t, float(np.sqrt((resid ** 2).sum(axis=1).mean()))


def head_tail_axis(mol, coords):
    """Vector from the aromatic ring centroid to the chain atom furthest from it, through bonds.

    Graph distance, not spatial: a folded tail can curl back near the ring, and a spatial choice
    would then pick an atom that is not the end of the chain at all.
    """
    ring_idx = [a.GetIdx() for a in mol.GetAtoms() if a.GetIsAromatic()]
    if not ring_idx:
        return None
    dmat = Chem.GetDistanceMatrix(mol)
    far = max(range(mol.GetNumAtoms()), key=lambda i: min(dmat[i][r] for r in ring_idx))
    head = coords[ring_idx].mean(axis=0)
    return coords[far] - head


def residue_contacts(pdb_path, lig_coords, cutoff):
    """{residue number} whose heavy atoms come within `cutoff` of any ligand heavy atom."""
    res_atoms = {}
    with open(pdb_path) as fh:
        for line in fh:
            if not line.startswith(("ATOM", "HETATM")) or line[76:78].strip() == "H":
                continue
            key = (line[21], int(line[22:26]), line[17:20].strip())
            res_atoms.setdefault(key, []).append(
                [float(line[30:38]), float(line[38:46]), float(line[46:54])])
    hit = set()
    for key, xyz in res_atoms.items():
        if np.linalg.norm(np.array(xyz)[:, None, :] - lig_coords[None, :, :], axis=2).min() < cutoff:
            hit.add(key)
    return hit


def pymol_script(work, structure, n_poses):
    """A session that loads the receptor, the reference ligand and every pose, reference in view."""
    lines = [
        f"# {structure}: Boltz reference (green) vs Vina poses",
        "# No alignment commands: every object is already in one frame.",
        f"load {structure}_protein.pdb, receptor",
        "load ligand_ref.sdf, boltz_ref",
        "load poses.sdf, vina_poses",
        "hide everything",
        "show cartoon, receptor",
        "color grey80, receptor",
        "show sticks, boltz_ref",
        "color green, boltz_ref",
        "show sticks, vina_poses",
        "color cyan, vina_poses",
        "set all_states, on, vina_poses",
        "set stick_radius, 0.12, vina_poses",
        "set stick_radius, 0.20, boltz_ref",
        "orient boltz_ref",
        "set ray_opaque_background, 0",
    ]
    path = os.path.join(work, "compare.pml")
    with open(path, "w") as fh:
        fh.write("\n".join(lines) + "\n")
    return path


def compare(structure, dock_root, contact):
    work = os.path.join(dock_root, structure)
    ref_sdf = os.path.join(work, "ligand_ref.sdf")
    rec_pdb = os.path.join(work, f"{structure}_protein.pdb")
    ref_mol = Chem.MolFromMolFile(ref_sdf, removeHs=True)
    ref = ref_mol.GetConformer().GetPositions()
    poses = ref_order_poses(os.path.join(work, "poses.pdbqt"),
                            os.path.join(work, "ligand.pdbqt"), ref)

    scores = []
    import re
    row_re = re.compile(r"^\s*\d+\s+(-?\d+\.\d+)\s+")
    for line in open(os.path.join(work, "vina.log")):
        m = row_re.match(line)
        if m:
            scores.append(float(m.group(1)))

    autos = ref_mol.GetSubstructMatches(ref_mol, uniquify=False, useChirality=False) or \
        [tuple(range(ref_mol.GetNumAtoms()))]
    ref_axis = head_tail_axis(ref_mol, ref)
    ref_contacts = residue_contacts(rec_pdb, ref, contact)

    rows = []
    for i, pose in enumerate(poses):
        # Take the automorphism that makes the two poses correspond best, so a symmetry-equivalent
        # ring flip is not reported as a real rotation.
        best = min(autos, key=lambda a: ((pose[list(a)] - ref) ** 2).sum())
        P = pose[list(best)]
        inplace = float(np.sqrt(((P - ref) ** 2).sum(axis=1).mean()))
        R, t, aligned = kabsch(ref, P)
        angle = float(np.degrees(np.arccos(np.clip((np.trace(R) - 1) / 2, -1, 1))))
        axis = head_tail_axis(ref_mol, P)
        ht = float(np.degrees(np.arccos(np.clip(
            ref_axis @ axis / (np.linalg.norm(ref_axis) * np.linalg.norm(axis)), -1, 1))))
        con = residue_contacts(rec_pdb, pose, contact)
        union = ref_contacts | con
        rows.append({
            "structure": structure,
            "pose": i + 1,
            "vina_score": scores[i] if i < len(scores) else None,
            "inplace_rmsd_A": round(inplace, 3),
            "aligned_rmsd_A": round(aligned, 3),
            "translation_A": round(float(np.linalg.norm(P.mean(axis=0) - ref.mean(axis=0))), 3),
            "rotation_deg": round(angle, 1),
            "headtail_deg": round(ht, 1),
            "n_contacts_ref": len(ref_contacts),
            "n_contacts_pose": len(con),
            "n_shared": len(ref_contacts & con),
            "contact_jaccard": round(len(ref_contacts & con) / len(union), 3) if union else "",
            "lost_residues": " ".join(f"{r[2]}{r[1]}" for r in sorted(ref_contacts - con,
                                                                     key=lambda r: r[1])),
            "gained_residues": " ".join(f"{r[2]}{r[1]}" for r in sorted(con - ref_contacts,
                                                                       key=lambda r: r[1])),
        })
    pymol_script(work, structure, len(poses))
    return rows


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("structures", nargs="*", help="structure names under runs/<system>/dock")
    ap.add_argument("--system", default="octinoxate")
    ap.add_argument("--all", action="store_true", help="every structure with a docking run")
    ap.add_argument("--contact", type=float, default=4.0,
                    help="heavy-atom separation counting as a residue contact, A (default 4.0)")
    a = ap.parse_args(argv)

    dock_root = os.path.join(REPO, "runs", a.system, "dock")
    names = list(a.structures)
    if a.all:
        names += [d for d in sorted(os.listdir(dock_root))
                  if os.path.exists(os.path.join(dock_root, d, "poses.pdbqt"))]
        names = sorted(set(names))
    if not names:
        ap.error("name at least one structure, or pass --all")

    rows = []
    for n in names:
        rows.extend(compare(n, dock_root, a.contact))

    out = os.path.join(dock_root, "pose_differences.csv")
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print(f"wrote {out}  ({len(rows)} poses)")
    print(f"wrote a compare.pml beside each structure's poses")

    print(f"\n{'structure':<22}{'pose':>5}{'score':>7}{'rmsd':>7}{'aligned':>8}"
          f"{'trans':>7}{'rot':>6}{'head/tail':>10}{'jaccard':>9}")
    for r in rows:
        print(f"{r['structure']:<22}{r['pose']:>5}{r['vina_score']:>7.1f}"
              f"{r['inplace_rmsd_A']:>7.2f}{r['aligned_rmsd_A']:>8.2f}"
              f"{r['translation_A']:>7.2f}{r['rotation_deg']:>6.0f}"
              f"{r['headtail_deg']:>10.0f}{r['contact_jaccard']:>9.2f}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
