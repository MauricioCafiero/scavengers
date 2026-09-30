#!/usr/bin/env python3
"""Write one docked pose as an SDF the MD prep will accept.

`omd prep-ligand` wants a molecule with correct bond orders and a stereocentre it can perceive, and
the poses Vina writes have neither: `poses.sdf` comes back from Open Babel's pdbqt reader with
bracketed `[C]` atoms and no valence, so parameterising it would either fail or silently produce a
different molecule than the reference.

So the pose is not read as a molecule at all. The reference SDF supplies the topology — it is the same
molecule, bond orders, stereocentre and all — and only its coordinates are replaced, with the pose's,
reindexed into the reference's atom order by `vina_redock.ref_order_poses`. The result is the reference
molecule sitting where Vina put it, which is exactly what should go into dynamics.

    python code/dock_pose_to_sdf.py s3_orig_f12 --pose 1
"""
import argparse
import os
import sys

from rdkit import Chem
from rdkit.Geometry import Point3D

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.expanduser("~/python_mac/dock_assist/code"))

from vina_redock import ref_order_poses  # noqa: E402


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("structure")
    ap.add_argument("--pose", type=int, default=1, help="1-based pose number (default 1, Vina's best)")
    ap.add_argument("--system", default="octinoxate")
    ap.add_argument("--out", default=None)
    a = ap.parse_args(argv)

    work = os.path.join(REPO, "runs", a.system, "dock", a.structure)
    ref = Chem.MolFromMolFile(os.path.join(work, "ligand_ref.sdf"), removeHs=False)
    if ref is None:
        ap.error(f"could not read the reference ligand in {work}")
    heavy = Chem.MolFromMolFile(os.path.join(work, "ligand_ref.sdf"), removeHs=True)
    poses = ref_order_poses(os.path.join(work, "poses.pdbqt"),
                            os.path.join(work, "ligand.pdbqt"),
                            heavy.GetConformer().GetPositions())
    if not 1 <= a.pose <= len(poses):
        ap.error(f"--pose must be 1..{len(poses)}")
    pose = poses[a.pose - 1]

    # The reference is read twice: with hydrogens for the molecule that gets written, without them for
    # the atom mapping, since the pdbqt is united-atom. Heavy atoms come first in this SDF, so the
    # mapping is positional; the hydrogens are then dropped, because prep-ligand adds its own and
    # leaving the reference's would put them on the old pose's geometry.
    out = Chem.RWMol(ref)
    conf = out.GetConformer()
    for i, xyz in enumerate(pose):
        assert out.GetAtomWithIdx(i).GetSymbol() != "H", "heavy atoms are not first in the reference"
        conf.SetAtomPosition(i, Point3D(*[float(v) for v in xyz]))
    mol = Chem.RemoveHs(out.GetMol())
    Chem.AssignStereochemistryFrom3D(mol)

    dest = a.out or os.path.join(work, f"pose{a.pose}_ligand.sdf")
    mol.SetProp("_Name", f"{a.structure}_dockpose{a.pose}")
    w = Chem.SDWriter(dest)
    w.write(mol)
    w.close()
    print(f"wrote {dest}")
    print(f"  SMILES    {Chem.MolToSmiles(mol)}")
    print(f"  reference {Chem.MolToSmiles(heavy)}")
    print(f"  match     {Chem.MolToSmiles(mol) == Chem.MolToSmiles(heavy)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
