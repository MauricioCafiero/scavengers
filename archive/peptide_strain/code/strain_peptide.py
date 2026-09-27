"""Peptide strain: what the peptide gives up conformationally in order to bind.

    peptide strain = E(peptide as it sits in the complex) - E(that sequence's apo fold)

Both states are Boltz heavy atoms with hydrogens relaxed at the same cutoff, so the difference between
them is the fold change and nothing else.

This replaces a term that was disabled for good reason. Relaxing a free peptide in vacuum collapses it
into a compact hydrogen-bonded ball -- 237.8 kcal/mol on a 33-mer, not converged -- so the old reference
measured collapse rather than strain. The apo fold is the sequence's predicted structure in the ligand's
absence, which is the state binding actually competes with.

**One reference per sequence, shared across its four constraint levels.** The apo state depends on the
sequence alone, so `s3_orig_control`, `_f4`, `_f8` and `_f12` all measure against the same apo fold.
That is what makes a ladder comparable: the differences along it are what forcing contacts costs the
peptide, with the sequence held fixed.

`E(peptide in the complex)` is the peptide portion of `structures/<name>_complex_relaxed_h.xyz`, which
is the peptide with the hydrogens the complex relaxation gave it, in the ligand's field. It is not
re-relaxed and its hydrogens are not re-placed: strain is the energy released going from the bound state
to the apo one, and the bound state is the peptide as it exists in the complex. The peptide comes first
in that file, so the first `n_peptide_atoms` rows are it -- that count is in the binding csv.

Usage:
    python code/strain_peptide.py runs/octinoxate --match s3_ --out strain_peptide_shell3.csv
"""
import argparse
import csv
import glob
import json
import os
import sys

import numpy as np


def read_xyz(path):
    lines = open(path).read().splitlines()
    n = int(lines[0])
    sym, xyz = [], []
    for ln in lines[2:2 + n]:
        f = ln.split()
        sym.append(f[0])
        xyz.append([float(f[1]), float(f[2]), float(f[3])])
    return sym, np.array(xyz)


def apo_key(name):
    """Which apo reference a fold uses: its sequence set, which its name carries.

    `s3_esm1_f8` -> `apo_s3_esm1`. The constraint level is dropped deliberately, because the apo state
    depends on the sequence and not on how many contacts were forced.
    """
    parts = name.split("_")
    if len(parts) < 3 or not parts[0].startswith("s"):
        return None
    return f"apo_{parts[0]}_{parts[1]}"


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0],
                               formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("outdir")
    p.add_argument("--match", help="only folds whose name contains this, e.g. s3_")
    p.add_argument("--reference", default="peptide_reference.json")
    p.add_argument("--out", default="strain_peptide.csv")
    p.add_argument("--model", default="uma-s-1p2p1")
    args = p.parse_args(argv)

    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from binding_energy import energy, EV_TO_KCAL
    from fairchem.core import pretrained_mlip, FAIRChemCalculator

    ref_path = os.path.join(args.outdir, args.reference)
    if not os.path.exists(ref_path):
        sys.exit(f"no {ref_path}\nrun: python code/peptide_apo.py fold {args.outdir} --shell N\n"
                 f"then: modal run code/modal_peptide_ref.py --shell N   (or peptide_apo.py energy)")
    refs = json.load(open(ref_path))

    boltz_dir = os.path.join(args.outdir, "boltz")
    # n_peptide_atoms tells us where the peptide ends in the combined file.
    natoms, charge = {}, {}
    for f in glob.glob(os.path.join(boltz_dir, "binding_*.csv")):
        for r in csv.DictReader(open(f)):
            if r.get("n_peptide_atoms"):
                natoms[r["name"]] = int(r["n_peptide_atoms"])
                charge[r["name"]] = int(r.get("peptide_charge") or 0)

    struct = os.path.join(boltz_dir, "structures")
    names = sorted(os.path.basename(p_)[:-len("_complex_relaxed_h.xyz")]
                   for p_ in glob.glob(os.path.join(struct, "*_complex_relaxed_h.xyz")))
    if args.match:
        names = [n for n in names if args.match in n]
    if not names:
        sys.exit("no complex geometries found; score the folds first")

    calc = FAIRChemCalculator(pretrained_mlip.get_predict_unit(args.model, device="cpu"),
                             task_name="omol")
    rows, skipped = [], []
    for i, name in enumerate(names, 1):
        key = apo_key(name)
        if key not in refs:
            skipped.append(f"{name} (no {key})")
            continue
        if name not in natoms:
            skipped.append(f"{name} (no n_peptide_atoms in any binding csv)")
            continue
        ref = refs[key]
        sym, xyz = read_xyz(os.path.join(struct, f"{name}_complex_relaxed_h.xyz"))
        n = natoms[name]
        # Sanity: the apo reference must be the same molecule, or the difference is meaningless.
        if ref["n_atoms"] != n:
            skipped.append(f"{name} ({n} peptide atoms vs {ref['n_atoms']} in {key})")
            continue
        e_bound = energy(sym[:n], xyz[:n], charge.get(name, 0), calc) * EV_TO_KCAL
        strain = e_bound - float(ref["e_apo_kcal"])
        rows.append({"name": name, "apo_reference": key,
                     "e_peptide_bound_kcal": f"{e_bound:.3f}",
                     "e_apo_kcal": f"{ref['e_apo_kcal']:.3f}",
                     "strain_peptide_kcal": f"{strain:.3f}"})
        print(f"[{i}/{len(names)}] {name:<20} strain {strain:>9.3f}   [{key}]", flush=True)

    if skipped:
        print(f"\nskipped {len(skipped)}:")
        for s in skipped:
            print(f"  {s}")
    if not rows:
        sys.exit("nothing to write")

    out = os.path.join(boltz_dir, args.out)
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    vals = [float(r["strain_peptide_kcal"]) for r in rows]
    print(f"\nstrain range {min(vals):.2f} to {max(vals):.2f} kcal/mol")
    if min(vals) < -1.0:
        print("note: a negative strain means the bound fold is lower in energy than the apo fold, "
              "which is possible here -- the apo fold is one Boltz prediction, not a minimum, and the "
              "bound fold was predicted with the ligand's constraints. Read it as the apo prediction "
              "being the worse structure, not as negative strain.")
    print(f"wrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
