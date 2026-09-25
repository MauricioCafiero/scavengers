"""How much do the folded complexes actually differ from each other?

Two different questions, and they do not answer each other:

**Superposed on the ligand.** The ligand is the one thing every complex has in common, so aligning
on it and then measuring how far the peptide's CA atoms sit apart says how differently each peptide
is *arranged around the ligand*. That is the quantity the design cares about -- two folds could be
the same shape and still present opposite faces to the ligand.

**Superposed on the peptide's own CA atoms.** This ignores the ligand and asks whether the folds are
the same shape at all. A low number here with a high number above means the same fold docked
differently.

All the structures compared must be the same length, which the designs are (33 residues), so CA
positions correspond one-to-one even where the sequences differ.

Usage:
    python code/compare_folds.py runs/octinoxate
    python code/compare_folds.py runs/octinoxate --reference orig_f4
"""
import argparse
import csv
import glob
import itertools
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from uma_binding import parse_cif  # noqa: E402


def kabsch(mobile, target):
    """Rotation and translation that best fits `mobile` onto `target`."""
    mc, tc = mobile.mean(axis=0), target.mean(axis=0)
    v, _, w = np.linalg.svd((mobile - mc).T @ (target - tc))
    d = np.sign(np.linalg.det(v @ w))
    rot = v @ np.diag([1.0, 1.0, d]) @ w
    return rot, mc, tc


def apply_fit(coords, rot, mc, tc):
    return (coords - mc) @ rot + tc


def rmsd(a, b):
    return float(np.sqrt(np.mean(np.sum((a - b) ** 2, axis=1))))


def load(cif):
    """(name, CA coordinates by residue, ligand heavy-atom coordinates)."""
    pep, lig = parse_cif(cif)
    ca = np.array([a["xyz"] for a in pep if a["name"] == "CA"])
    heavy = np.array([a["xyz"] for a in lig if a["element"] != "H"])
    return ca, heavy


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("outdir")
    parser.add_argument("--out", default="fold_similarity.csv")
    args = parser.parse_args(argv)

    boltz_dir = os.path.join(args.outdir, "boltz")
    cifs = sorted(glob.glob(os.path.join(boltz_dir, "boltz_results_*", "predictions", "*",
                                         "*_model_0.cif")))
    if len(cifs) < 2:
        sys.exit(f"need at least two folded complexes under {boltz_dir}")

    data = {}
    for cif in cifs:
        name = os.path.basename(cif).replace("_model_0.cif", "")
        ca, lig = load(cif)
        data[name] = (ca, lig)
    names = sorted(data)
    lengths = {len(ca) for ca, _ in data.values()}
    if len(lengths) > 1:
        print(f"WARNING: differing residue counts {sorted(lengths)}; CA positions do not "
              f"correspond one-to-one, so these numbers are not comparable")
    print(f"{len(names)} complexes, {sorted(lengths)[0]} residues each\n")

    on_ligand, on_self = {}, {}
    for a, b in itertools.combinations(names, 2):
        ca_a, lig_a = data[a]
        ca_b, lig_b = data[b]
        # align on the ligand, then compare where the peptides sit
        rot, mc, tc = kabsch(lig_b, lig_a)
        on_ligand[(a, b)] = rmsd(apply_fit(ca_b, rot, mc, tc), ca_a)
        # align on the peptide itself: is it even the same shape?
        rot, mc, tc = kabsch(ca_b, ca_a)
        on_self[(a, b)] = rmsd(apply_fit(ca_b, rot, mc, tc), ca_a)

    def matrix(title, table):
        print(f"\n{title}")
        short = [n.replace("RGGDGGKGGGGLGGGGKGIGEGWGGDGSGGEGS", "orig")[:14] for n in names]
        print("               " + "".join(f"{s[:7]:>8}" for s in short))
        for i, a in enumerate(names):
            row = f"{short[i]:<14} "
            for j, b in enumerate(names):
                if i == j:
                    row += f"{'-':>8}"
                else:
                    key = (a, b) if (a, b) in table else (b, a)
                    row += f"{table[key]:>8.1f}"
            print(row)

    matrix("CA RMSD after superposing on the LIGAND (A) -- how differently each peptide "
           "is arranged around it:", on_ligand)
    matrix("CA RMSD after superposing on the PEPTIDE (A) -- whether the folds are the "
           "same shape at all:", on_self)

    lv = np.array(list(on_ligand.values()))
    sv = np.array(list(on_self.values()))
    print(f"\non the ligand:  median {np.median(lv):.1f} A, range {lv.min():.1f}-{lv.max():.1f}")
    print(f"on the peptide: median {np.median(sv):.1f} A, range {sv.min():.1f}-{sv.max():.1f}")
    print("\nA large ligand-superposed number with a small peptide-superposed one means the same"
          "\nfold presenting itself to the ligand differently.")

    out = os.path.join(boltz_dir, args.out)
    with open(out, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["structure_a", "structure_b", "ca_rmsd_on_ligand", "ca_rmsd_on_peptide"])
        for a, b in itertools.combinations(names, 2):
            w.writerow([a, b, f"{on_ligand[(a, b)]:.3f}", f"{on_self[(a, b)]:.3f}"])
    print(f"\nresults in {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
