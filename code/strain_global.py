"""Recompute ligand strain for every fold against the one global reference.

`strain_i = E(bound_i) - E_ref`, where `E_ref` comes from `ligand_reference.py` and is the same number
for every structure in the run. That replaces the old definition, which relaxed each complex's own
bound pose and so measured every structure against a different local minimum.

Two things make this reproducible where the old one was not, and cheap where the old one was not:

* **`E(bound_i)` needs only the ligand's heavy atoms**, which come from the Boltz CIF, plus hydrogens
  added by RDKit's `AddHs(addCoords=True)` -- geometric, deterministic -- and then relaxed with the
  heavy atoms held fixed. No complex relaxation, so this runs on folds whose complex geometries were
  never saved, and takes seconds per structure instead of tens of minutes.
* **The hydrogens are optimised on the isolated ligand, not inside the complex.** Previously they were
  relaxed in the peptide's field, starting from a `pdbfixer` protonation that is not deterministic, and
  that was worth 3.6 kcal/mol of apparent strain on one structure tested twice. Optimising them alone
  also separates conformational strain from interaction, which is the point of reporting them apart.

Results go to a csv of their own; nothing existing is modified.

Usage:
    python code/strain_global.py runs/octinoxate
    python code/strain_global.py runs/octinoxate --match s2_ --out strain_global_shell2.csv
"""
import argparse
import csv
import glob
import json
import os
import sys

import numpy as np


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0],
                               formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("outdir")
    p.add_argument("--match", help="only folds whose name contains this")
    p.add_argument("--reference", default="ligand_reference.json")
    p.add_argument("--out", default="strain_global.csv")
    p.add_argument("--fmax", type=float, default=0.01, help="hydrogen relaxation cutoff (default 0.01)")
    p.add_argument("--steps", type=int, default=300, help="hydrogen relaxation cap (default 300)")
    p.add_argument("--model", default="uma-s-1p2p1")
    p.add_argument("--save-structures", action="store_true",
                   help="write each bound ligand with its relaxed hydrogens")
    args = p.parse_args(argv)

    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from boltz_check import ligand_smiles
    from uma_binding import parse_cif, protonate_ligand
    from binding_energy import energy, relax_hydrogens, EV_TO_KCAL
    from peptide_builder import write_xyz
    from fairchem.core import pretrained_mlip, FAIRChemCalculator

    ref_path = os.path.join(args.outdir, args.reference)
    if not os.path.exists(ref_path):
        sys.exit(f"no {ref_path}; run code/ligand_reference.py first")
    ref = json.load(open(ref_path))
    e_ref = float(ref["reference_kcal"])
    print(f"reference: {e_ref:.3f} kcal/mol-equivalent "
          f"(conformer {ref['reference_conformer']} of {ref['n_conformers']}, "
          f"spread {ref['spread_kcal']} across {ref['keep']} relaxed)")

    # The reference is a property of the ligand and the potential, so it is computed once and reused
    # for every shell of the same ligand -- which is the point. That also makes it the kind of file
    # that outlives the assumptions behind it, so check them: a reference from a different molecule or
    # a different model is not a smaller error than a missing one, it is a silently wrong strain for
    # every structure.
    smiles_now = ligand_smiles(args.outdir)
    if ref.get("smiles") != smiles_now:
        sys.exit(f"reference was computed for a different ligand:\n"
                 f"  reference: {ref.get('smiles')}\n"
                 f"  this run:  {smiles_now}\n"
                 f"delete {os.path.basename(ref_path)} and re-run ligand_reference.py")
    if ref.get("model") != args.model:
        sys.exit(f"reference was computed with model {ref.get('model')!r} but this run uses "
                 f"{args.model!r}. Energies from different potentials are not comparable; re-run "
                 f"ligand_reference.py with --model {args.model}")

    boltz_dir = os.path.join(args.outdir, "boltz")
    cifs = sorted(glob.glob(os.path.join(boltz_dir, "boltz_results_*", "predictions", "*",
                                         "*_model_0.cif")))
    if args.match:
        cifs = [c for c in cifs if args.match in os.path.basename(c)]
    if not cifs:
        sys.exit("no folded complexes found")
    print(f"{len(cifs)} folds\n")

    smiles = ligand_smiles(args.outdir)
    calc = FAIRChemCalculator(pretrained_mlip.get_predict_unit(args.model, device="cpu"),
                             task_name="omol")
    struct_dir = os.path.join(boltz_dir, "structures")
    if args.save_structures:
        os.makedirs(struct_dir, exist_ok=True)

    rows = []
    for i, cif in enumerate(cifs, 1):
        name = os.path.basename(cif).replace("_model_0.cif", "")
        _pep, lig = parse_cif(cif)
        sym, xyz, chg = protonate_ligand(lig, smiles)
        # Hydrogens relaxed with the heavy atoms fixed: Boltz's coordinates are the prediction and
        # must not move, but RDKit's geometric hydrogens carry an artificial strain of their own.
        xyz = relax_hydrogens(sym, xyz, chg, calc, fmax=args.fmax, steps=args.steps,
                              label=f"{name} ligand")
        e_bound = energy(sym, xyz, chg, calc) * EV_TO_KCAL
        strain = e_bound - e_ref
        if args.save_structures:
            write_xyz(os.path.join(struct_dir, f"{name}_ligand_bound_hrelaxed.xyz"), sym, xyz,
                      f"{name}: Boltz heavy atoms, hydrogens relaxed alone at fmax {args.fmax}")
        rows.append({"name": name, "n_atoms": len(sym), "ligand_charge": chg,
                     "e_bound_kcal": f"{e_bound:.3f}",
                     "e_reference_kcal": f"{e_ref:.3f}",
                     "strain_global_kcal": f"{strain:.3f}"})
        print(f"[{i}/{len(cifs)}] {name:<38} strain {strain:>8.3f}", flush=True)

    out = os.path.join(boltz_dir, args.out)
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    vals = [float(r["strain_global_kcal"]) for r in rows]
    print(f"\nstrain range {min(vals):.2f} to {max(vals):.2f} kcal/mol")
    if min(vals) < -0.5:
        print(f"WARNING: a negative strain means a bound pose below the reference, so the conformer "
              f"search missed the global minimum. Re-run ligand_reference.py with more conformers.")
    print(f"wrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
