#!/usr/bin/env python3
"""Assemble a self-contained folder for rescoring the docked poses with GNINA elsewhere.

GNINA's CNN was trained to recognise near-native poses, which is the one thing AutoDock Vina
demonstrably failed at here: across eight structures Vina ranked the pose closest to the predicted one
6th to 9th of nine every single time, including the one genuine 1.30 A hit in `bg33_4`. So the question
this bundle is built to answer is not "what is the affinity" but "does a CNN rescorer put the predicted
pose on top when Vina would not". Both flip-degenerate orientations and the predicted pose itself go in
as scorable ligands so that ranking can be read directly.

Each structure gets one multi-model SDF holding ten poses in a fixed order:

    pose 0   the predicted pose -- what Boltz or BoltzGen placed, and what MD was run from
    pose 1-9 Vina's output, in Vina's own ranking order

The poses are NOT copied out of `poses.sdf`. Open Babel's pdbqt reader loses valence, so that file's
molecules come back with bracketed `[C]` atoms and no implicit hydrogens, which would corrupt the atom
typing every CNN score depends on. Instead the reference SDF supplies the topology and only the
coordinates are replaced, reindexed by `vina_redock.ref_order_poses` -- the same approach
`dock_pose_to_sdf.py` uses for the MD inputs.

Two receptors are shipped per structure because the right choice is not obvious and a round trip to
find out is expensive: `receptor.pdb` is the heavy-atom file Vina actually scored against, so GNINA
numbers computed on it are directly comparable to the Vina column; `receptor_h.pdb` is the PDBFixer
output with hydrogens, which is what the MD used. Score against the first for comparability, and
against the second if the CNN turns out to be sensitive to it.

    python code/make_gnina_bundle.py --all
"""
import argparse
import csv
import os
import shutil
import sys

import numpy as np
from rdkit import Chem
from rdkit.Geometry import Point3D

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.expanduser(os.environ.get("DOCK_ASSIST", "~/python_mac/dock_assist") + "/code"))

from vina_redock import ref_order_poses  # noqa: E402

# Everything measured for these structures already, so a GNINA number can be correlated the moment it
# comes back rather than after another round of joining files. MM/GBSA dG and its standard error,
# residence within 10 A, and release episodes, from runs/octinoxate/md/mmgbsa_summary.csv and the
# release analysis in README.md; the BoltzGen pair from boltzgen_local/README.md.
REFERENCE = {
    "s3_orig_f12":          dict(dg=-24.33, dg_err=0.08, residence=100.0, releases=0,  source="peptidebuilder"),
    "s3_esm2_f4":           dict(dg=-21.08, dg_err=0.08, residence=99.6,  releases=0,  source="peptidebuilder"),
    "bg33_4":               dict(dg=-19.66, dg_err=0.02, residence=100.0, releases=0,  source="boltzgen"),
    "s2_esm2_control":      dict(dg=-16.25, dg_err=0.08, residence=59.6,  releases=1,  source="peptidebuilder"),
    # flat docked-leg fields ride along with the pose rows (REFERENCE values would otherwise
    # describe only the co-folded start); measured 2026-10-04, from mmgbsa_summary.csv + md_contacts
    "shuffle_control":      dict(dg=-15.13, dg_err=0.13, residence=75.9,  releases=5,  source="peptidebuilder (null)",
                                 vina_p1_dg=-16.15, vina_p1_residence=57.8, vina_p1_releases=19,
                                 gnina_p9_dg=-11.15, gnina_p9_residence=46.4, gnina_p9_releases=10,
                                 cnnscore_p7_dg=-8.40, cnnscore_p7_residence=22.4, cnnscore_p7_releases=24),
    "shuffle_control_esm0": dict(dg=-14.32, dg_err=0.10, residence=80.2,  releases=5,  source="peptidebuilder (null)",
                                 vina_p1_dg=-15.50, vina_p1_residence=76.1, vina_p1_releases=6,
                                 gnina_p9_dg=-14.76, gnina_p9_residence=90.5, gnina_p9_releases=8),
    "orig_f12":             dict(dg=-13.71, dg_err=0.04, residence=77.3,  releases=0,  source="peptidebuilder"),
    "bg33_3":               dict(dg=-11.86, dg_err=0.02, residence=41.0,  releases=18, source="boltzgen"),
    "bg33_2":               dict(dg=-8.02,  dg_err=0.04, residence=50.9,  releases=9,  source="boltzgen",
                                 gnina_p4_dg=-10.44, gnina_p4_residence=57.2, gnina_p4_releases=2),
    "bg33_1":               dict(dg=-11.31, dg_err=0.03, residence=14.5,  releases=16, source="boltzgen",
                                 gnina_p2_dg=-12.55, gnina_p2_residence=42.5, gnina_p2_releases=11),
}


def build(structure, dock_root, out_root):
    work = os.path.join(dock_root, structure)
    ref_h = Chem.MolFromMolFile(os.path.join(work, "ligand_ref.sdf"), removeHs=False)
    ref = Chem.MolFromMolFile(os.path.join(work, "ligand_ref.sdf"), removeHs=True)
    ref_xyz = ref.GetConformer().GetPositions()
    poses = ref_order_poses(os.path.join(work, "poses.pdbqt"),
                            os.path.join(work, "ligand.pdbqt"), ref_xyz)

    dest = os.path.join(out_root, structure)
    os.makedirs(dest, exist_ok=True)
    shutil.copy2(os.path.join(work, f"{structure}_protein.pdb"), os.path.join(dest, "receptor.pdb"))
    for md_root in (os.path.join(REPO, "runs", "octinoxate", "md"),
                    os.path.expanduser("~/python_mac/boltzgen_local/md")):
        h = os.path.join(md_root, structure, "protein_fixed.pdb")
        if os.path.exists(h):
            shutil.copy2(h, os.path.join(dest, "receptor_h.pdb"))
            break

    autos = ref.GetSubstructMatches(ref, uniquify=False, useChirality=False) or \
        [tuple(range(ref.GetNumAtoms()))]

    def rmsd(P):
        return min(float(np.sqrt(((P[list(a)] - ref_xyz) ** 2).sum(axis=1).mean())) for a in autos)

    scores = []
    import re
    row_re = re.compile(r"^\s*\d+\s+(-?\d+\.\d+)\s+")
    for line in open(os.path.join(work, "vina.log")):
        m = row_re.match(line)
        if m:
            scores.append(float(m.group(1)))

    rows = []
    # The predicted pose on its own as well as in the bundle: --autobox_ligand needs a single molecule
    # to build the search box from, and it is also the reference any docked pose is scored against.
    writer = Chem.SDWriter(os.path.join(dest, "poses_with_reference.sdf"))
    for i, xyz in enumerate([ref_xyz] + poses):
        mol = Chem.RWMol(ref_h)
        conf = mol.GetConformer()
        for j, p in enumerate(xyz):
            conf.SetAtomPosition(j, Point3D(*[float(v) for v in p]))
        out = Chem.RemoveHs(mol.GetMol())
        Chem.AssignStereochemistryFrom3D(out)
        name = f"{structure}_pose{i}" + ("_reference" if i == 0 else "")
        out.SetProp("_Name", name)
        out.SetProp("structure", structure)
        out.SetProp("pose", str(i))
        out.SetProp("is_reference", "1" if i == 0 else "0")
        r = 0.0 if i == 0 else rmsd(xyz)
        out.SetProp("rmsd_to_reference_A", f"{r:.3f}")
        vina = "" if i == 0 else f"{scores[i - 1]:.1f}"
        out.SetProp("vina_score", vina)
        writer.write(out)
        if i == 0:
            w0 = Chem.SDWriter(os.path.join(dest, "ref_pose.sdf"))
            w0.write(out)
            w0.close()
        rows.append(dict(structure=structure, pose=i, is_reference=int(i == 0),
                         rmsd_to_reference_A=round(r, 3), vina_score=vina,
                         **REFERENCE.get(structure, {})))
    writer.close()
    return rows


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("structures", nargs="*")
    ap.add_argument("--system", default="octinoxate")
    ap.add_argument("--all", action="store_true")
    a = ap.parse_args(argv)

    dock_root = os.path.join(REPO, "runs", a.system, "dock")
    out_root = os.path.join(REPO, "runs", a.system, "gnina")
    os.makedirs(out_root, exist_ok=True)
    names = list(a.structures)
    if a.all:
        names += [d for d in sorted(os.listdir(dock_root))
                  if os.path.exists(os.path.join(dock_root, d, "poses.pdbqt"))]
        names = sorted(set(names))
    if not names:
        ap.error("name at least one structure, or pass --all")

    rows = []
    for n in names:
        rows.extend(build(n, dock_root, out_root))
        print(f"  {n}: receptor + 10 poses + ref_pose.sdf")

    ref_csv = os.path.join(out_root, "reference.csv")
    with open(ref_csv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print(f"wrote {ref_csv} ({len(rows)} poses across {len(names)} structures)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
