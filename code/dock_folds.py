#!/usr/bin/env python3
"""Stage every folded structure (cofolds and apo folds) for Vina and dock it.

The fold campaign's structures (runs/<mol>/folds/<lane>/<pep>.pdb) contain no
Vina inputs, so before any of them can become an MD leg they need the input
pair vina_redock.redock() wants: <name>_protein.pdb (the peptide's heavy atoms,
ligand stripped) and <name>_ligand.sdf (the reference ligand with bond orders).
This script writes those and then calls redock() unchanged -- same 22 A box,
exhaustiveness 16, seed 42, same cutoffs, pose 0 = the reference placement.

The reference SDF is built as:
- cofold lanes (openfold3, rf3): the fold contains the ligand. Its atoms are
  perceived by obabel (bonds from geometry) and the molecule's reference
  template -- copied with all its bond orders and stereocentres from the
  co-folding leg's own <pep>_ligand.sdf (these folds and the Boltz cofolds
  dock the same molecule for a given system) -- is rigidly aligned onto the
  perceived ligand. The reference then sits exactly where the cofolder put the
  ligand, which is what pose 0 measures against.
- apo lanes (esmfold, of3apo, rf3apo): the fold has no ligand to perceive. The
  peptide's OF3-cofold reference SDF is superposed onto the apo fold's C-alpha
  (resseq-matched Kabsch) and the transformed coordinates become the reference
  pose -- so a ligand-free fold still gets its search box centred on the part
  of the chain where its cofold packs the ligand. This is a site definition,
  not a prediction: the box says "look here", and nothing else.

Work lands in runs/<mol>/dock/<name>/ (names carry the lane suffix, so they
cannot collide with the committed Boltz-cofold dock folders); the rows go to
runs/<mol>/folds/metrics/dock_folds_summary.csv and _poses.csv. Appending to
the committed dock_summary.csv / dock_poses.csv is a separate step, done
additions-only after the pass. Resumable: a name already in the summary is
skipped.

Run with dock-env's python (rdkit + vina):
    ~/python_mac/dock_assist/dock-env/bin/python code/dock_folds.py
"""
import csv
import os
import sys
import time
from collections import defaultdict

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolAlign
from rdkit.Geometry import Point3D


def apply_transform(mol, mat):
    """Rigidly place a molecule's every atom by a 4x4 transform. (TransformMol is
    absent from some RDKit builds, so the affine product is done directly.)"""
    conf = mol.GetConformer()
    for i in range(mol.GetNumAtoms()):
        x = np.array(conf.GetAtomPosition(i))
        conf.SetAtomPosition(i, Point3D(*(x @ mat[:3, :3].T + mat[:3, 3])))

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
BG_MD = os.path.expanduser("~/python_mac/boltzgen_local/md")
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(os.path.expanduser("~/python_mac/dock_assist"), "code"))

from vina_redock import DockError, redock, convert  # noqa: E402

# lane -> (suffix, has_ligand). Cofold lanes first: the apo lanes derive their
# reference placement from the openfold3 cofold lane's staged SDF.
LANES = [("openfold3", "o3cof", True),
         ("rf3", "rf3cof", True),
         ("of3apo", "o3apo", False),
         ("rf3apo", "rf3apo", False),
         ("esmfold", "esmapo", False)]

WATERS = {"HOH", "WAT", "NA", "CL"}


def fold_lines(path):
    """Peptide and ligand pdb line groups from a converted fold structure.

    Peptide = the chain holding the most distinct residue numbers (OF3 ligand
    is chain Z, 1 residue; RF3's is chain B; apo folds are single-chain).
    Ligand = every other chain's non-water HETATM residues. Returns
    ([peptide ATOM/HETATM lines], [ligand lines]).
    """
    chains = defaultdict(list)
    for line in open(path):
        if not line.startswith(("ATOM", "HETATM")):
            continue
        resname = line[17:20].strip()
        if resname in WATERS:
            continue
        chains[line[21]].append(line)
    if len(chains) == 1:
        return list(chains.values())[0], []
    pep_key = max(chains, key=lambda c: len({x[22:26] for x in chains[c]}))
    pep = [x for x in chains[pep_key] if not is_hydrogen(x)]
    lig = [x for c, lines in chains.items() if c != pep_key for x in lines
           if not is_hydrogen(x)]
    return pep, lig


def is_hydrogen(line):
    elem = line[76:78].strip()
    if elem == "H":
        return True
    name = line[12:16].strip()
    return elem == "" and name.startswith("H") and len(name) == 1


def write_pdb(lines, out):
    with open(out, "w") as fh:
        fh.writelines(lines)
        fh.write("END\n")


def template_sdf_path(mol, pep):
    """The co-folding leg's own reference ligand: bond orders + stereo, Boltz pose."""
    for root in (os.path.join(REPO, "runs", mol, "md"), BG_MD):
        p = os.path.join(root, pep, f"{pep}_ligand.sdf")
        if os.path.exists(p):
            return p
    raise DockError(f"no template ligand SDF for {mol}/{pep}")


def cofold_ref_sdf(pep_lines, lig_lines, template_path, out):
    """Reference molecule placed where the cofold put the ligand.

    The fold's ligand bonds are perceived from geometry (RDKit proximity
    bonding; obabel fallback), then the template -- read for real, with bond
    orders, stereocentres and hydrogens -- has its bond orders assigned back and
    is rigidly aligned onto the perceived ligand via GetBestAlignmentTransform.
    Geometry perception misreads bonds on strained cofold conformers (obabel
    loses the alkene, RDKit closes an occasional ring bond), but those are
    exactly what AssignBondOrdersFromTemplate repairs: the template's
    connectivity is imposed wherever the perceived graph admits it, so the
    substructure match that follows is over chemically valid mappings only.
    Only the template's coordinates change -- the molecule, its orders and its
    atom order are untouched. Returns the alignment rmsd."""
    t_heavy = Chem.MolFromMolFile(template_path, removeHs=True)
    t_full = Chem.MolFromMolFile(template_path, removeHs=False)
    tmp = out + ".ligand.pdb"
    write_pdb(lig_lines, tmp)
    probe = Chem.MolFromPDBFile(tmp, removeHs=True, proximityBonding=True)
    if probe is None:
        perceived = out + ".perceived.sdf"
        convert(["obabel", "-ipdb", tmp, "-osdf", "-O", perceived],
                "fold ligand pdb->SDF perception")
        probe = Chem.MolFromMolFile(perceived, removeHs=True)
        if probe is None:
            raise DockError(f"neither RDKit nor obabel read the fold's ligand: {tmp}")
    if probe.GetNumAtoms() != t_heavy.GetNumAtoms():
        raise DockError(f"fold ligand {probe.GetNumAtoms()} heavy atoms != "
                        f"template {t_heavy.GetNumAtoms()}")
    probe = AllChem.AssignBondOrdersFromTemplate(t_heavy, probe)
    rms, mat, _ = rdMolAlign.GetBestAlignmentTransform(t_heavy, probe)
    apply_transform(t_full, mat)
    w = Chem.SDWriter(out)
    w.write(t_full)
    w.close()
    return float(rms)


def cas(lines):
    out = {}
    for line in lines:
        if not line.startswith("ATOM") or line[12:16].strip() != "CA":
            continue
        if is_hydrogen(line):
            continue
        out[int(line[22:26])] = np.array(
            [float(line[30:38]), float(line[38:46]), float(line[46:54])])
    return out


def kabsch(P, Q):
    """(R, t) with P @ R + t superposed onto Q (both (n,3), same order)."""
    Pc, Qc = P - P.mean(0), Q - Q.mean(0)
    U, S, Vt = np.linalg.svd(Pc.T @ Qc)
    d = np.sign(np.linalg.det(U @ Vt))
    R = U @ np.diag([1, 1, d]) @ Vt
    return R, Q.mean(0) - P.mean(0) @ R


def apo_ref_sdf(ref_sdf, cofold_cas, apo_cas, out):
    """Reference pose carried onto an apo fold by C-alpha superposition."""
    keys = sorted(set(cofold_cas) & set(apo_cas))
    if len(keys) < 10:
        raise DockError(f"only {len(keys)} matched CA residues -- refusing the superposition")
    R, t = kabsch(np.array([apo_cas[r] for r in keys]),
                  np.array([cofold_cas[r] for r in keys]))
    mol = Chem.MolFromMolFile(ref_sdf, removeHs=False)
    conf = mol.GetConformer()
    for i in range(mol.GetNumAtoms()):
        x = np.array(conf.GetAtomPosition(i))
        conf.SetAtomPosition(i, Point3D(*(x @ R + t)))
    w = Chem.SDWriter(out)
    w.write(mol)
    w.close()


def cofold_pose_sdf(lig_lines, template_path, out):
    """The cofold's own ligand conformation, as an SDF the MD prep will accept.

    The reference SDF built by cofold_ref_sdf carries the template's conformer
    rigidly placed at the fold ligand (~1 A off). For an MD leg that starts
    FROM the cofold pose, that approximation is unnecessary: the perceived
    ligand with template bond orders (AssignBondOrdersFromTemplate's output)
    IS the cofold's ligand to its own geometry. This writes it in template atom
    order (the mapping GetBestAlignmentTransform returns), hydrogens removed --
    the same shape dock_pose_to_sdf.py produces for a docked pose, and
    prep-ligand adds hydrogens itself."""
    t_heavy = Chem.MolFromMolFile(template_path, removeHs=True)
    t_full = Chem.MolFromMolFile(template_path, removeHs=False)
    tmp = out + ".ligand.pdb"
    write_pdb(lig_lines, tmp)
    probe = Chem.MolFromPDBFile(tmp, removeHs=True, proximityBonding=True)
    if probe is None:
        raise DockError(f"RDKit could not read the fold's ligand: {tmp}")
    if probe.GetNumAtoms() != t_heavy.GetNumAtoms():
        raise DockError(f"fold ligand {probe.GetNumAtoms()} heavy atoms != "
                        f"template {t_heavy.GetNumAtoms()}")
    probe = AllChem.AssignBondOrdersFromTemplate(t_heavy, probe)
    _, _, amap = rdMolAlign.GetBestAlignmentTransform(t_heavy, probe)
    pc = probe.GetConformer().GetPositions()
    rw = Chem.RWMol(t_full)
    conf = rw.GetConformer()
    mapped = 0
    for prb_i, ref_i in amap:
        conf.SetAtomPosition(int(prb_i), tuple(pc[int(ref_i)]))
        mapped += 1
    if mapped != t_heavy.GetNumAtoms():
        raise DockError(f"template substructure map covers {mapped}/{t_heavy.GetNumAtoms()} "
                        f"heavy atoms -- cofold pose unmappable")
    w = Chem.SDWriter(out)
    w.write(Chem.RemoveHs(rw.GetMol()))
    w.close()


def main():
    t0 = time.time()
    done = 0
    for mol in ("octinoxate", "oxybenzone"):
        man = os.path.join(REPO, "runs", mol, "folds", "manifest.csv")
        if not os.path.exists(man):
            continue
        peps = [row["peptide"] for row in csv.DictReader(open(man))]
        stage = os.path.join(REPO, "runs", mol, "folds", "dock_stage")
        out_root = os.path.join(REPO, "runs", mol, "dock")
        sum_csv = os.path.join(REPO, "runs", mol, "folds", "metrics", "dock_folds_summary.csv")
        pose_csv = os.path.join(REPO, "runs", mol, "folds", "metrics", "dock_folds_poses.csv")
        summary_rows, pose_rows = [], []
        if os.path.exists(sum_csv):
            summary_rows = list(csv.DictReader(open(sum_csv)))
            done_names = {r["structure"] for r in summary_rows}
        else:
            done_names = set()
        os.makedirs(stage, exist_ok=True)

        # Leg-ready cofold poses for the o3cof lane, independent of the docking
        # resume guard (which skips docked names without staging these).
        for pep in peps:
            fold = os.path.join(REPO, "runs", mol, "folds", "openfold3", f"{pep}.pdb")
            pose_out = os.path.join(stage, "openfold3", f"{pep}_o3cof_ligand_cofold.sdf")
            if not (os.path.exists(fold) and not os.path.exists(pose_out)):
                continue
            _, lig_lines = fold_lines(fold)
            if lig_lines:
                cofold_pose_sdf(lig_lines, template_sdf_path(mol, pep), pose_out)
                print(f"cofold pose staged: {pep}_o3cof")

        # Cofold lanes first: the apo lanes need their staged reference SDFs.
        for pep in peps:
            for lane, suffix, has_lig in LANES:
                fold = os.path.join(REPO, "runs", mol, "folds", lane, f"{pep}.pdb")
                if not os.path.exists(fold):
                    print(f"no fold {mol}/{pep}/{lane}, skipping")
                    continue
                name = f"{pep}_{suffix}"
                if name in done_names:
                    continue
                started = time.time()
                try:
                    pep_lines, lig_lines = fold_lines(fold)
                    md_dir = os.path.join(stage, lane)
                    os.makedirs(md_dir, exist_ok=True)
                    write_pdb(pep_lines, os.path.join(md_dir, f"{name}_protein.pdb"))
                    ref_sdf = os.path.join(md_dir, f"{name}_ligand.sdf")
                    if has_lig:
                        if not lig_lines:
                            raise DockError("cofold fold without a ligand")
                        rms = cofold_ref_sdf(pep_lines, lig_lines,
                                             template_sdf_path(mol, pep), ref_sdf)
                        extra = f"align_rmsd {rms:.2f}"
                        if lane == "openfold3":
                            # the apo lanes re-use this SDF transformed
                            openfold3_ref[pep] = ref_sdf
                    else:
                        src_ref = openfold3_ref.get(pep)
                        if src_ref is None or not os.path.exists(src_ref):
                            raise DockError(f"no staged OF3-cofold reference to build {name} from")
                        if pep not in cofold_cas_cache:
                            p, _ = fold_lines(os.path.join(
                                REPO, "runs", mol, "folds", "openfold3", f"{pep}.pdb"))
                            cofold_cas_cache[pep] = cas(p)
                        apo_ref_sdf(src_ref, cofold_cas_cache[pep], cas(pep_lines), ref_sdf)
                        extra = "superposed from o3cof"
                    s, p = redock(name, md_dir, out_root)
                    summary_rows.append(s)
                    pose_rows.extend(p)
                    el = time.time() - started
                    print(f"{mol[:3]} {name:24s} best {s['vina_best_score']:>6} "
                          f"pose1_rmsd {s['pose1_rmsd_to_boltz']:>6} | {extra} | {el:5.1f}s")
                except Exception as e:  # noqa: E722 -- keep the pass going, report at the end
                    print(f"FAILED {mol}/{name}: {e}")
                done += 1
        # write the scratch tables after each peptide so a kill keeps the rows
        if summary_rows:
            for path, rows in ((sum_csv, summary_rows), (pose_csv, pose_rows)):
                if not rows:
                    continue
                os.makedirs(os.path.dirname(path), exist_ok=True)
                with open(path, "w") as f:
                    w = csv.DictWriter(f, list(rows[0]))
                    w.writeheader()
                    w.writerows(rows)
    print(f"ALL DONE: {done} structures docked, {(time.time() - t0) / 60:.0f} min total")


if __name__ == "__main__":
    openfold3_ref = {}
    cofold_cas_cache = {}
    main()