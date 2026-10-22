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


def replicate_report(outdir, out_csv):
    """Each unseeded re-fold against its own original. See code/fold_replicates.py.

    Three numbers per draw, because they answer different questions and a single "RMSD" would hide
    the interesting case:

    `ca_on_peptide`  -- superpose the two peptides' CAs and compare them. Is it even the same fold?
    `ca_on_ligand`   -- superpose on the LIGAND, then compare CAs. How differently is the peptide
                        arranged around the ligand? This is the quantity the design cares about.
    `lig_on_peptide` -- superpose on the peptide's CAs, then compare ligand heavy atoms. Did the
                        ligand land in the same pocket of its own peptide? Same construction as
                        ligand_slide.py, so the number is read on the same scale as a leg's slide.

    A small `ca_on_peptide` with a large `ca_on_ligand` is the same fold presenting a different face
    to the ligand -- which would matter more for this project than a different shape would.
    """
    reps = sorted(glob.glob(os.path.join(outdir, "folds", "foldrep", "*", "rep*")))
    if not reps:
        sys.exit(f"no replicate folds under {outdir}/folds/foldrep/ -- run code/fold_replicates.py")
    rows = []
    for rep in reps:
        name = os.path.basename(os.path.dirname(rep))
        k = os.path.basename(rep)
        rep_cif = os.path.join(rep, f"boltz_results_{name}", "predictions", name,
                               f"{name}_model_0.cif")
        orig_cif = os.path.join(outdir, "boltz", f"boltz_results_{name}", "predictions", name,
                                f"{name}_model_0.cif")
        if not (os.path.exists(rep_cif) and os.path.exists(orig_cif)):
            print(f"skipping {name} {k}: missing {'replicate' if not os.path.exists(rep_cif) else 'original'} cif")
            continue
        ca_o, lig_o = load(orig_cif)
        ca_r, lig_r = load(rep_cif)
        if len(ca_o) != len(ca_r) or len(lig_o) != len(lig_r):
            print(f"skipping {name} {k}: atom counts differ from the original")
            continue
        rot, mc, tc = kabsch(ca_r, ca_o)
        ca_self = rmsd(apply_fit(ca_r, rot, mc, tc), ca_o)
        lig_self = rmsd(apply_fit(lig_r, rot, mc, tc), lig_o)
        rot, mc, tc = kabsch(lig_r, lig_o)
        ca_lig = rmsd(apply_fit(ca_r, rot, mc, tc), ca_o)
        rows.append({"structure": name, "draw": k, "residues": len(ca_o),
                     "ca_on_peptide": round(ca_self, 2), "ca_on_ligand": round(ca_lig, 2),
                     "lig_on_peptide": round(lig_self, 2)})
        print(f"{name:34s} {k:6s} ca_on_peptide {ca_self:6.2f}   ca_on_ligand {ca_lig:6.2f}   "
              f"lig_on_peptide {lig_self:6.2f}")
    if not rows:
        return 1
    path = os.path.join(outdir, "folds", out_csv)
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    for key in ("ca_on_peptide", "ca_on_ligand", "lig_on_peptide"):
        v = np.array([r[key] for r in rows])
        print(f"\n{key:15s} median {np.median(v):.2f} A   range {v.min():.2f}-{v.max():.2f}")
    print(f"\n{len(rows)} draws over {len({r['structure'] for r in rows})} structures -> {path}")
    print("These are independent unseeded draws, NOT failed reproductions: the spread is the "
          "measurement.\nIt gives the resolution at which single-draw comparisons between fold arms "
          "can be read, and\ndoes not retract any forced-vs-unforced finding.")
    return 0


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("outdir")
    parser.add_argument("--out", default="fold_similarity.csv")
    parser.add_argument("--replicates", action="store_true",
                        help="compare each unseeded re-fold under folds/foldrep/ against its own "
                             "original, instead of all pairs of different peptides")
    args = parser.parse_args(argv)
    if args.replicates:
        out = args.out if args.out != "fold_similarity.csv" else "fold_replicate_rms.csv"
        return replicate_report(args.outdir, out)

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
    skipped = 0
    for a, b in itertools.combinations(names, 2):
        ca_a, lig_a = data[a]
        ca_b, lig_b = data[b]
        # CA positions only correspond one-to-one between peptides of the same length. The warning
        # above anticipated mixed lengths but nothing acted on it, so the subtraction raised
        # instead -- oxybenzone's set holds 25-, 31- and 33-residue peptides, which made the
        # all-pairs mode unusable for that molecule (found 2026-10-22). Skip the pair and say so.
        if len(ca_a) != len(ca_b) or len(lig_a) != len(lig_b):
            skipped += 1
            continue
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
                key = (a, b) if (a, b) in table else (b, a)
                if i == j:
                    row += f"{'-':>8}"
                elif key not in table:
                    row += f"{'.':>8}"      # skipped: lengths differ, no correspondence
                else:
                    row += f"{table[key]:>8.1f}"
            print(row)

    matrix("CA RMSD after superposing on the LIGAND (A) -- how differently each peptide "
           "is arranged around it:", on_ligand)
    matrix("CA RMSD after superposing on the PEPTIDE (A) -- whether the folds are the "
           "same shape at all:", on_self)

    if not on_ligand:
        sys.exit(f"no comparable pairs: all {len(names)} complexes differ in length from each other")
    if skipped:
        print(f"\n({skipped} pairs skipped: different residue or ligand atom counts, so CA "
              f"positions do not correspond)")
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
            if (a, b) not in on_ligand:     # skipped above: lengths differ
                continue
            w.writerow([a, b, f"{on_ligand[(a, b)]:.3f}", f"{on_self[(a, b)]:.3f}"])
    print(f"\nresults in {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
