"""Overlay a Boltz-2 co-folded complex on the fragment-based design it came from.

Both structures contain the same ligand, so the ligand is the common reference: the Boltz
complex is rigidly superposed onto the design by its ligand atoms, and then the question is
whether the folded peptide puts its side chains where the fragment search put them.

Reports, per structure:
  - ligand heavy-atom RMSD after superposition (how much the ligand conformer differs)
  - for each placed fragment, the distance to the nearest side chain of the same residue
    type in the folded peptide
  - how many fragments are matched within a cutoff

Usage:
    python code/overlay.py runs/octinoxate [--cutoff 3.0] [--save-overlays]
"""
import argparse
import csv
import glob
import os
import sys

import numpy as np

BACKBONE = {"N", "CA", "C", "O", "OXT"}

THREE_TO_ONE = {"ARG": "R", "LYS": "K", "ASP": "D", "GLU": "E", "ILE": "I", "LEU": "L",
                "SER": "S", "TRP": "W", "TYR": "Y", "PHE": "F", "GLY": "G"}


def kabsch(mobile, target):
    """Rotation + translation that best maps `mobile` onto `target` (both (N,3))."""
    mc, tc = mobile.mean(axis=0), target.mean(axis=0)
    u, _, vt = np.linalg.svd((mobile - mc).T @ (target - tc))
    d = np.sign(np.linalg.det(vt.T @ u.T))
    rot = vt.T @ np.diag([1.0, 1.0, d]) @ u.T
    return rot, tc - rot @ mc


# the rigid part of the ligand: methoxyphenyl + vinyl + ester. The hexyl tail is flexible,
# so including it would fold conformational differences into the alignment frame.
CORE_SMARTS = "COc1ccc(cc1)C=CC(=O)O"


def ligand_mol(symbols, xyz, charge=0):
    """RDKit mol with 3D coordinates from element symbols and positions."""
    from rdkit import Chem
    from rdkit.Chem import rdDetermineBonds

    block = f"{len(symbols)}\nligand\n" + "\n".join(
        f"{s} {c[0]} {c[1]} {c[2]}" for s, c in zip(symbols, xyz))
    mol = Chem.MolFromXYZBlock(block)
    rdDetermineBonds.DetermineBonds(mol, charge=charge)
    return mol


def ligand_correspondence(design_symbols, design_xyz, boltz_lig, smiles, core_only=True):
    """Corresponding ligand atom positions (boltz, design), matched through the ligand's
    rigid core and resolved over symmetry-equivalent mappings by best fit."""
    from rdkit import Chem
    from rdkit.Chem import AllChem

    design_mol = Chem.RemoveHs(ligand_mol(design_symbols, design_xyz))
    design_heavy_xyz = np.array([design_xyz[i] for i, s in enumerate(design_symbols) if s != "H"])

    boltz_heavy = [a for a in boltz_lig if a["element"] != "H"]
    boltz_xyz = np.array([a["xyz"] for a in boltz_heavy])
    lig_block = "\n".join(
        f"HETATM{i + 1:5d} {a['name']:<4s} LIG A   1    "
        f"{a['xyz'][0]:8.3f}{a['xyz'][1]:8.3f}{a['xyz'][2]:8.3f}  1.00  0.00          {a['element']:>2s}"
        for i, a in enumerate(boltz_heavy))
    boltz_mol = Chem.MolFromPDBBlock(lig_block + "\nEND\n", removeHs=False, sanitize=False)
    boltz_mol = AllChem.AssignBondOrdersFromTemplate(Chem.MolFromSmiles(smiles), boltz_mol)
    Chem.SanitizeMol(boltz_mol)

    query = Chem.MolFromSmarts(CORE_SMARTS) if core_only else Chem.MolFromSmiles(smiles)
    design_matches = design_mol.GetSubstructMatches(query, uniquify=False)
    boltz_matches = boltz_mol.GetSubstructMatches(query, uniquify=False)
    if not design_matches or not boltz_matches:
        raise RuntimeError("the ligand core could not be matched in both structures")

    # try every pairing of symmetry images and keep the one that fits best
    best = None
    for dm in design_matches:
        target = design_heavy_xyz[list(dm)]
        for bm in boltz_matches:
            mobile = boltz_xyz[list(bm)]
            rot, trans = kabsch(mobile, target)
            rmsd = float(np.sqrt(((mobile @ rot.T + trans - target) ** 2).sum(axis=1).mean()))
            if best is None or rmsd < best[0]:
                best = (rmsd, mobile, target)
    return best[1], best[2]


def read_xyz(path):
    lines = open(path).read().splitlines()
    n = int(lines[0])
    symbols, xyz = [], []
    for line in lines[2:2 + n]:
        f = line.split()
        symbols.append(f[0])
        xyz.append([float(f[1]), float(f[2]), float(f[3])])
    return symbols, np.array(xyz), lines[1]


def design_fragment_centres(outdir, sequence_file, n_ligand_atoms):
    """Centroid and residue type of each fragment placed in a design structure."""
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    import peptide_builder as pb

    symbols, xyz, comment = read_xyz(sequence_file)
    names = comment.replace("placed:", "").split()
    frags = {f["name"]: f for f in pb.define_fragments()}
    centres = []
    offset = n_ligand_atoms
    for name in names:
        n = frags[name]["num_atoms"]
        block = [j for j in range(offset, offset + n) if symbols[j] != "H"]
        centres.append((pb.ONE_LETTER[name], xyz[block].mean(axis=0)))
        offset += n
    return centres, symbols[:n_ligand_atoms], xyz[:n_ligand_atoms]


def sidechain_centres(peptide_atoms):
    """Centroid and one-letter type of each residue side chain (glycine skipped)."""
    residues = {}
    for a in peptide_atoms:
        residues.setdefault((a["resseq"], a["resname"]), []).append(a)
    out = []
    for (resseq, resname), atoms in sorted(residues.items(), key=lambda kv: int(kv[0][0])):
        side = [a for a in atoms if a["name"] not in BACKBONE]
        if not side or resname == "GLY":
            continue
        out.append((THREE_TO_ONE.get(resname, "?"), int(resseq),
                    np.array([a["xyz"] for a in side]).mean(axis=0)))
    return out


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("outdir", help="a peptide_builder run directory, e.g. runs/octinoxate")
    parser.add_argument("--cutoff", type=float, default=3.0, help="a fragment counts as matched within this, A (default: 3.0)")
    parser.add_argument("--save-overlays", action="store_true", help="write the superposed complex next to each design")
    parser.add_argument("--out", default="overlay.csv",
                        help="output csv name inside boltz/ (default overlay.csv). Give a second "
                             "shell its own, or it overwrites the first shell's results")
    parser.add_argument("--match",
                        help="only folds whose name contains this, e.g. s2_ for a second shell")
    args = parser.parse_args(argv)

    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from boltz_check import ligand_smiles
    from uma_binding import parse_cif
    smiles = ligand_smiles(args.outdir)

    boltz_dir = os.path.join(args.outdir, "boltz")
    seq_dir = os.path.join(args.outdir, "sequences")
    overlay_dir = os.path.join(boltz_dir, "overlays")
    if args.save_overlays:
        os.makedirs(overlay_dir, exist_ok=True)

    cifs = sorted(glob.glob(os.path.join(boltz_dir, "boltz_results_*", "predictions", "*", "*_model_0.cif")))
    if args.match:
        cifs = [c for c in cifs if args.match in os.path.basename(c)]
    if not cifs:
        sys.exit(f"no Boltz structures found under {boltz_dir}"
                 + (f" matching {args.match!r}" if args.match else ""))

    rows = []
    for path in cifs:
        name = os.path.basename(path).replace("_model_0.cif", "")
        design = os.path.join(seq_dir, f"{name}.xyz")
        if not os.path.exists(design):
            print(f"{name}: no design structure in {seq_dir}, skipping")
            continue

        pep_atoms, lig_atoms = parse_cif(path)
        # the design file starts with the ligand, so its atom count comes from the run's ligand.xyz
        n_lig = int(open(os.path.join(args.outdir, "ligand.xyz")).readline())
        centres, lig_symbols, lig_xyz = design_fragment_centres(args.outdir, design, n_lig)

        try:
            boltz_pairs, design_pairs = ligand_correspondence(lig_symbols, lig_xyz, lig_atoms, smiles)
        except Exception as exc:
            print(f"{name}: ligand matching failed ({type(exc).__name__}: {exc})")
            continue

        rot, trans = kabsch(boltz_pairs, design_pairs)
        lig_rmsd = float(np.sqrt(((boltz_pairs @ rot.T + trans - design_pairs) ** 2).sum(axis=1).mean()))

        moved = [(t, r, rot @ c + trans) for t, r, c in sidechain_centres(pep_atoms)]
        matched, distances, any_type = 0, [], []
        detail = []
        for one_letter, centre in centres:
            if moved:
                any_type.append(min(float(np.linalg.norm(m[2] - centre)) for m in moved))
            same = [m for m in moved if m[0] == one_letter]
            if not same:
                detail.append(f"{one_letter}:none")
                continue
            d = min(float(np.linalg.norm(m[2] - centre)) for m in same)
            distances.append(d)
            detail.append(f"{one_letter}:{d:.1f}")
            if d <= args.cutoff:
                matched += 1
        median = float(np.median(distances)) if distances else float("nan")
        median_any = float(np.median(any_type)) if any_type else float("nan")
        print(f"{name:<28} ligand RMSD {lig_rmsd:5.2f} A   {matched}/{len(centres)} fragments within "
              f"{args.cutoff} A   median {median:5.2f} A (any type {median_any:5.2f} A)", flush=True)
        print(f"{'':28} per fragment: {' '.join(detail)}")
        rows.append({"name": name, "ligand_rmsd_A": f"{lig_rmsd:.3f}", "n_fragments": len(centres),
                     "matched_within_cutoff": matched, "median_distance_A": f"{median:.3f}",
                     "median_distance_any_type_A": f"{median_any:.3f}"})

        if args.save_overlays:
            pep_moved = np.array([a["xyz"] for a in pep_atoms]) @ rot.T + trans
            lig_moved = np.array([a["xyz"] for a in lig_atoms]) @ rot.T + trans
            symbols = [a["element"] for a in pep_atoms] + [a["element"] for a in lig_atoms]
            coords = np.vstack([pep_moved, lig_moved])
            import peptide_builder as pb
            pb.write_xyz(os.path.join(overlay_dir, f"{name}_boltz_on_design.xyz"), symbols, coords,
                         f"Boltz complex superposed on {name} design by ligand")

    if rows:
        csv_path = os.path.join(boltz_dir, args.out)
        with open(csv_path, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
            writer.writeheader()
            writer.writerows(rows)
        print(f"\nResults in {csv_path}")


if __name__ == "__main__":
    main()
