"""Null control: do the designed peptides bind better than random ones of the same length?

Generates random sequences of a given length, co-folds each with the run's ligand in
Boltz-2, computes the UMA interaction energy of the predicted complex, and compares the
result with the designed sequences of that same length.

The comparison that matters is the UMA interaction energy on the predicted geometry,
together with the minimum peptide-ligand contact distance, since a positive energy or a
contact under about 2.6 A means the predicted structure is sterically broken rather than
weakly bound.

Usage:
    python code/random_control.py runs/octinoxate --length 19 --count 5 [--seed 1]
    python code/random_control.py runs/octinoxate --length 19 --shuffle RYGLSGIKWDKSFEGDGGE
"""
import argparse
import csv
import glob
import os
import random
import sys

import numpy as np

# the 20 standard amino acids
AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("outdir", help="a peptide_builder run directory, e.g. runs/octinoxate")
    parser.add_argument("--length", type=int, required=True, help="peptide length to test")
    parser.add_argument("--count", type=int, default=5, help="how many random peptides (default: 5)")
    parser.add_argument("--seed", type=int, default=1, help="random seed (default: 1)")
    parser.add_argument("--shuffle", help="instead of random sequences, shuffle this one (keeps composition)")
    parser.add_argument("--model", default="uma-s-1p2p1", help="fairchem model (default: uma-s-1p2p1)")
    parser.add_argument("--boltz-repo", default=os.path.expanduser("~/python_mac/boltz_local"))
    parser.add_argument("--fixer-venv", default=os.path.expanduser("~/python_mac/pocket_assist/venv"))
    args = parser.parse_args(argv)

    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from boltz_check import ligand_smiles, run_boltz
    import uma_binding as ub

    random.seed(args.seed)
    if args.shuffle:
        pool = list(args.shuffle)
        sequences = []
        while len(sequences) < args.count:
            random.shuffle(pool)
            candidate = "".join(pool)
            if candidate != args.shuffle and candidate not in sequences:
                sequences.append(candidate)
        tag = "shuf"
    else:
        sequences = ["".join(random.choice(AMINO_ACIDS) for _ in range(args.length))
                     for _ in range(args.count)]
        tag = "rand"

    smiles = ligand_smiles(args.outdir)
    boltz_dir = os.path.join(args.outdir, "boltz")
    os.makedirs(boltz_dir, exist_ok=True)
    print(f"ligand SMILES: {smiles}")
    print(f"{len(sequences)} {'shuffled' if args.shuffle else 'random'} sequences of length "
          f"{len(sequences[0])}, seed {args.seed}\n")

    import torch
    from fairchem.core import FAIRChemCalculator, pretrained_mlip
    from ase import Atoms

    device = "cuda" if torch.cuda.is_available() else "cpu"
    predictor = pretrained_mlip.get_predict_unit(args.model, device=device)
    calculator = FAIRChemCalculator(predictor, task_name="omol")

    def energy(symbols, positions, charge):
        atoms = Atoms(symbols, positions=positions)
        atoms.info.update({"charge": int(charge), "spin": 1})
        atoms.calc = calculator
        return atoms.get_potential_energy()

    rows = []
    for i, sequence in enumerate(sequences, 1):
        name = f"{tag}{args.seed}_{i}_{sequence}"
        print(f"[{i}/{len(sequences)}] {sequence}", flush=True)
        out = run_boltz(name, sequence, smiles, False, boltz_dir, args.boltz_repo)
        if out is None:
            continue
        cif = glob.glob(os.path.join(boltz_dir, f"boltz_results_{name}", "predictions", "*", "*_model_0.cif"))
        if not cif:
            print("  no structure written")
            continue
        pep_atoms, lig_atoms = ub.parse_cif(cif[0])
        try:
            pep_symbols, pep_xyz = ub.protonate_peptide(pep_atoms, args.fixer_venv)
            lig_symbols, lig_xyz, lig_charge = ub.protonate_ligand(lig_atoms, smiles)
        except Exception as exc:
            print(f"  could not add hydrogens ({type(exc).__name__})")
            continue
        pep_charge, n_res = ub.peptide_charge(pep_atoms)

        symbols = pep_symbols + lig_symbols
        positions = np.vstack([pep_xyz, lig_xyz])
        ie = 23.06035 * (energy(symbols, positions, pep_charge + lig_charge)
                         - energy(pep_symbols, pep_xyz, pep_charge)
                         - energy(lig_symbols, lig_xyz, lig_charge))
        P = np.array([a["xyz"] for a in pep_atoms])
        L = np.array([a["xyz"] for a in lig_atoms])
        contact = float(np.linalg.norm(P[:, None, :] - L[None, :, :], axis=2).min())

        print(f"  UMA IE {ie:8.2f} kcal/mol   min contact {contact:.2f} A   "
              f"Boltz dG {out['dG_kcal_mol']:.2f}   binder P {out['binder_probability']:.2f}", flush=True)
        rows.append({"sequence": sequence, "kind": "shuffled" if args.shuffle else "random",
                     "n_residues": n_res, "uma_ie_kcal_mol": f"{ie:.3f}",
                     "min_contact_A": f"{contact:.3f}",
                     "boltz_dG_kcal_mol": f"{out['dG_kcal_mol']:.3f}",
                     "binder_probability": f"{out['binder_probability']:.3f}"})

        csv_path = os.path.join(boltz_dir, "controls.csv")
        existing = []
        if os.path.exists(csv_path):
            with open(csv_path) as f:
                existing = [r for r in csv.DictReader(f) if r["sequence"] not in {x["sequence"] for x in rows}]
        with open(csv_path, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
            writer.writeheader()
            writer.writerows(existing + rows)

    # compare with the designed sequences of the same length
    uma_path = os.path.join(boltz_dir, "uma_binding.csv")
    if rows and os.path.exists(uma_path):
        with open(uma_path) as f:
            designed = [r for r in csv.DictReader(f)
                        if int(r["n_residues"]) == args.length and not r["name"].startswith(("rand", "shuf"))]
        ours = np.array([float(r["uma_ie_kcal_mol"]) for r in designed])
        theirs = np.array([float(r["uma_ie_kcal_mol"]) for r in rows])
        print(f"\nlength {args.length}: {len(designed)} designed vs {len(theirs)} "
              f"{'shuffled' if args.shuffle else 'random'}")
        if len(ours):
            print(f"  designed UMA IE: {', '.join(f'{v:.1f}' for v in sorted(ours))}"
                  f"   (mean {ours.mean():.2f}, best {ours.min():.2f})")
        print(f"  control  UMA IE: {', '.join(f'{v:.1f}' for v in sorted(theirs))}"
              f"   (mean {theirs.mean():.2f}, best {theirs.min():.2f})")
        if len(ours):
            print(f"  designed better than {(theirs > ours.mean()).sum()}/{len(theirs)} controls "
                  f"on the mean; best design beats {(theirs > ours.min()).sum()}/{len(theirs)}")
    print(f"\nResults in {os.path.join(boltz_dir, 'controls.csv')}")


if __name__ == "__main__":
    main()
