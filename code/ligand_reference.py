"""Find the ligand's global minimum once, as the single reference for every strain calculation.

Ligand strain is `E(ligand in the complex) - E(ligand free)`, and the free term was being obtained by
relaxing each complex's own bound pose. That makes the reference different for every structure: which
local basin the relaxation falls into depends on where it started and how long it ran, so the strains
were not on a common footing. Measured on shell 2, tightening the force cutoff from 0.1 to 0.01 moved
individual strains by +1.5 to +7.7 kcal/mol, and the largest shifts came with up to 1.0 A of
heavy-atom drift -- a rotor flipping, not a relaxation.

The free ligand is the same molecule for every structure in a run, so its energy should be computed
once and reused. That is what this does:

1. embed 20 conformers with ETKDG from the run's own SMILES, perceived from `ligand.xyz`
2. optimise all of them with MMFF94, which is cheap, and rank
3. take the five lowest and relax each with UMA to a tight cutoff
4. the lowest UMA energy is the reference, written to `ligand_reference.json`

Strain then becomes `E(bound_i) - E_ref` with one shared `E_ref`, so differences between structures
reflect only their bound geometries. It is a strain against the ligand's own best conformer, which is
the quantity that matters for binding, and it is reproducible because the reference no longer depends
on any individual complex.

MMFF ranks before UMA because ranking 20 conformers with UMA would cost twenty relaxations to find
the same five; MMFF is approximate but easily good enough to discard the obviously bad ones, and the
five survivors are then judged by the potential that everything else here uses.

Usage:
    python code/ligand_reference.py runs/octinoxate
    python code/ligand_reference.py runs/octinoxate --n-conformers 40 --keep 8
"""
import argparse
import json
import os
import sys

import numpy as np


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0],
                               formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("outdir", help="a run directory, e.g. runs/octinoxate")
    p.add_argument("--n-conformers", type=int, default=20, help="ETKDG embeddings (default 20)")
    p.add_argument("--keep", type=int, default=5,
                   help="how many of the MMFF-lowest to relax with UMA (default 5)")
    p.add_argument("--fmax", type=float, default=0.01,
                   help="UMA force cutoff, eV/A (default 0.01). Tight on purpose: this is the "
                        "reference, and it should be a real minimum")
    p.add_argument("--steps", type=int, default=1000, help="UMA step cap (default 1000)")
    p.add_argument("--seed", type=int, default=0xC0FFEE, help="ETKDG seed, for reproducibility")
    p.add_argument("--model", default="uma-s-1p2p1")
    p.add_argument("--out", default="ligand_reference.json")
    args = p.parse_args(argv)

    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from rdkit import Chem
    from rdkit.Chem import AllChem
    from boltz_check import ligand_smiles
    from binding_energy import energy, relaxed_energy, EV_TO_KCAL
    from peptide_builder import write_xyz
    from fairchem.core import pretrained_mlip, FAIRChemCalculator

    smiles = ligand_smiles(args.outdir)
    print(f"ligand SMILES (perceived from ligand.xyz): {smiles}")

    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    params = AllChem.ETKDGv3()
    params.randomSeed = args.seed
    cids = AllChem.EmbedMultipleConfs(mol, numConfs=args.n_conformers, params=params)
    print(f"embedded {len(cids)} conformers")

    mmff = AllChem.MMFFOptimizeMoleculeConfs(mol, maxIters=2000)
    ranked = sorted(((e, cid) for cid, (rc, e) in zip(cids, mmff)), key=lambda t: t[0])
    print(f"\nMMFF94 energies, kcal/mol (lowest {args.keep} go to UMA):")
    for i, (e, cid) in enumerate(ranked):
        print(f"  {'->' if i < args.keep else '  '} conf {cid:>3}  {e:>10.3f}")

    symbols = [a.GetSymbol() for a in mol.GetAtoms()]
    n_heavy = sum(1 for s in symbols if s != "H")
    print(f"\n{len(symbols)} atoms ({n_heavy} heavy)")

    calc = FAIRChemCalculator(pretrained_mlip.get_predict_unit(args.model, device="cpu"),
                             task_name="omol")

    results = []
    print(f"\nUMA relaxation of the {args.keep} lowest, fmax {args.fmax}:")
    for e_mmff, cid in ranked[:args.keep]:
        xyz = np.array(mol.GetConformer(cid).GetPositions())
        e0 = energy(symbols, xyz, 0, calc) * EV_TO_KCAL
        e_rel, x_rel, steps, conv = relaxed_energy(symbols, xyz, 0, calc,
                                                   fmax=args.fmax, steps=args.steps, label=None)
        e_rel_kcal = e_rel * EV_TO_KCAL
        results.append({"conformer": int(cid), "mmff_kcal": round(float(e_mmff), 3),
                        "uma_start_kcal": round(float(e0), 3),
                        "uma_relaxed_kcal": round(float(e_rel_kcal), 3),
                        "steps": int(steps), "converged": bool(conv),
                        "xyz": x_rel.tolist()})
        print(f"  conf {cid:>3}  MMFF {e_mmff:>9.3f}   UMA start {e0:>13.3f}   "
              f"UMA relaxed {e_rel_kcal:>13.3f}   ({steps} steps, conv={bool(conv)})")

    best = min(results, key=lambda r: r["uma_relaxed_kcal"])
    spread = max(r["uma_relaxed_kcal"] for r in results) - best["uma_relaxed_kcal"]
    print(f"\nreference: conformer {best['conformer']} at {best['uma_relaxed_kcal']:.3f} kcal/mol")
    print(f"spread across the {len(results)} relaxed conformers: {spread:.3f} kcal/mol")
    if not all(r["converged"] for r in results):
        print("WARNING: not every conformer converged; raise --steps")
    # If MMFF's ranking disagrees with UMA's, the five kept may not contain the true best.
    if results.index(best) != 0:
        print(f"note: MMFF's lowest was not UMA's lowest (UMA prefers rank "
              f"{results.index(best) + 1} of {len(results)}), so MMFF's ranking is only a filter. "
              f"Raise --keep if the spread above is small enough to worry about.")

    struct_dir = os.path.join(args.outdir, "boltz", "structures")
    os.makedirs(struct_dir, exist_ok=True)
    ref_xyz = os.path.join(struct_dir, "ligand_global_reference.xyz")
    write_xyz(ref_xyz, symbols, np.array(best["xyz"]),
              f"ligand global minimum: conformer {best['conformer']} of {args.n_conformers}, "
              f"MMFF-ranked, UMA-relaxed at fmax {args.fmax}; "
              f"{best['uma_relaxed_kcal']:.3f} kcal/mol-equivalent")

    payload = {
        "smiles": smiles, "n_atoms": len(symbols), "n_heavy": n_heavy,
        "n_conformers": args.n_conformers, "keep": args.keep,
        "fmax": args.fmax, "steps": args.steps, "seed": args.seed, "model": args.model,
        "reference_kcal": best["uma_relaxed_kcal"],
        "reference_conformer": best["conformer"],
        "spread_kcal": round(float(spread), 3),
        "candidates": [{k: v for k, v in r.items() if k != "xyz"} for r in results],
        "reference_xyz": os.path.relpath(ref_xyz, args.outdir),
    }
    out = os.path.join(args.outdir, args.out)
    with open(out, "w") as fh:
        json.dump(payload, fh, indent=2)
    print(f"\nwrote {out}\nwrote {ref_xyz}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
