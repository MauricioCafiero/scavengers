"""Fold each sequence on its own, as the reference state for peptide strain.

Peptide strain was disabled because its reference was meaningless: relaxing a free peptide in vacuum
collapses it into a compact hydrogen-bonded ball, giving 237.8 kcal/mol on a 33-mer and not converging.
That measures collapse, not strain.

The reference used here is the peptide's **apo fold** -- Boltz's prediction for the sequence with no
ligand present -- with hydrogens added and relaxed while the heavy atoms stay fixed. Then

    peptide strain = E(peptide as it sits in the complex) - E(peptide folded alone)

which asks what the peptide gives up conformationally in order to bind.

Two things make this well defined where the vacuum reference was not:

* **Only hydrogens are relaxed, heavy atoms fixed.** Relaxing the whole apo structure with UMA would
  collapse it toward the same vacuum ball and reintroduce the original problem. Holding Boltz's apo
  heavy atoms also makes the reference symmetric with the bound state, which is likewise Boltz heavy
  atoms with hydrogens relaxed, so the difference between them is the fold change and nothing else.
* **One reference per sequence, not per fold.** The apo state depends on the sequence alone, so a
  sequence's four constraint levels all share it -- which is what makes the comparison across a ladder
  meaningful: it isolates what forcing contacts costs the peptide.

Two caveats worth keeping in view. This is strain against Boltz's *predicted* apo fold, a
model-internal quantity rather than an experimental free state; it is well defined but it is not the
same kind of reference as the ligand's, which is a genuine conformer-search minimum (a search is not
available for a 34-mer). And Boltz returns one structure where a real free peptide is an ensemble, so
the reference understates the free state's flexibility. Read the apo confidence alongside it,
especially for the glycine designs, which are 3-6% helical and may fold apo with little conviction.

Split by cost: `fold` runs Boltz, which is fast on this laptop's GPU. `energy` runs UMA, which is not.

**`energy` runs in either place, with no second implementation.** It picks CUDA when it can see a
device and CPU otherwise, so the same code path serves `code/modal_peptide_ref.py` on a rented GPU and
this laptop when there is no credit to spend. On CPU expect roughly what a scoring pass costs per
structure -- tens of minutes for a 400-600 atom peptide -- against a few minutes on an L4. There are
only three apo folds per shell, so the local route is a long coffee rather than an overnight run.

Usage:
    python code/peptide_apo.py fold runs/octinoxate --shell 3          # Boltz, locally
    modal run code/modal_peptide_ref.py --shell 3                      # UMA on a GPU
    python code/peptide_apo.py energy runs/octinoxate --shell 3        # UMA here instead
    python code/peptide_apo.py energy runs/octinoxate --shell 3 --device cpu   # force it
"""
import argparse
import glob
import json
import os
import sys

import numpy as np

APO_DIR = "apo"


def sequences_for_shell(outdir, shell):
    """The design sequence and its ESM2 variants, tagged the way the folds are named."""
    design = os.path.join(outdir, f"design_shell{shell}.json")
    if not os.path.exists(design):
        sys.exit(f"no {design}")
    out = [("orig", json.load(open(design))["sequence"])]
    variants = os.path.join(outdir, f"variants_shell{shell}.txt")
    if os.path.exists(variants):
        for i, ln in enumerate(open(variants).read().split(), 1):
            out.append((f"esm{i}", ln.strip()))
    return out


def do_fold(args):
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    import boltz_env

    apo = os.path.join(args.outdir, APO_DIR)
    os.makedirs(apo, exist_ok=True)
    for tag, seq in sequences_for_shell(args.outdir, args.shell):
        name = f"apo_s{args.shell}_{tag}"
        cif = os.path.join(apo, f"boltz_results_{name}", "predictions", name, f"{name}_model_0.cif")
        if os.path.exists(cif):
            print(f"{name}: already folded, skipping")
            continue
        yaml_path = os.path.join(apo, f"{name}.yaml")
        # Protein only -- no ligand block and no constraints. That is the whole point: this is the
        # sequence's preferred fold in the ligand's absence.
        with open(yaml_path, "w") as fh:
            fh.write("\n".join(["version: 1", "sequences:", "  - protein:", "      id: A",
                                f"      sequence: {seq}", "      msa: empty", ""]))
        log = os.path.join(apo, f"{name}.log")
        print(f"folding {name} ({len(seq)} residues) -> {log}", flush=True)
        code = boltz_env.run_boltz(yaml_path, apo, log)
        print(f"  exit {code}; {'wrote ' + cif if os.path.exists(cif) else 'NO STRUCTURE, see ' + log}")
    return 0


def do_energy(args):
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from uma_binding import parse_cif, protonate_peptide, peptide_charge
    from binding_energy import energy, relax_hydrogens, EV_TO_KCAL
    from peptide_builder import write_xyz
    from fairchem.core import pretrained_mlip, FAIRChemCalculator

    apo = os.path.join(args.outdir, APO_DIR)
    cifs = sorted(glob.glob(os.path.join(apo, "boltz_results_apo_s*", "predictions", "*",
                                         "*_model_0.cif")))
    if args.shell:
        cifs = [c for c in cifs if f"apo_s{args.shell}_" in os.path.basename(c)]
    if not cifs:
        sys.exit(f"no apo structures under {apo}; run `fold` first")

    calc = FAIRChemCalculator(pretrained_mlip.get_predict_unit(args.model, device=args.device),
                             task_name="omol")
    out_path = os.path.join(args.outdir, args.out)
    store = json.load(open(out_path)) if os.path.exists(out_path) else {}

    for i, cif in enumerate(cifs, 1):
        name = os.path.basename(cif).replace("_model_0.cif", "")
        pep, _lig = parse_cif(cif)
        sym, xyz = protonate_peptide(pep)
        # peptide_charge returns (charge, n_residues), and the cyclic flag follows the name, exactly as
        # binding_energy.py does it -- so the apo charge matches the bound charge for the same sequence.
        chg, n_res = peptide_charge(pep, name.startswith("cyclo_"))
        # Hydrogens only, heavy atoms fixed, at the same cutoff the complex uses -- so the apo and
        # bound states are treated identically and their difference is the fold change alone.
        xyz = relax_hydrogens(sym, xyz, chg, calc, fmax=args.fmax, steps=args.steps,
                              label=f"{name} apo")
        e = energy(sym, xyz, chg, calc) * EV_TO_KCAL
        write_xyz(os.path.join(apo, f"{name}_relaxed_h.xyz"), sym, xyz,
                  f"{name}: Boltz apo fold, hydrogens relaxed (fmax {args.fmax}, "
                  f"max {args.steps} steps)")
        store[name] = {"e_apo_kcal": round(float(e), 3), "n_atoms": len(sym),
                       "peptide_charge": int(chg), "fmax": args.fmax, "steps": args.steps,
                       "model": args.model, "n_residues": n_res}
        print(f"[{i}/{len(cifs)}] {name:<22} E_apo {e:>14.3f}  ({len(sym)} atoms, charge {chg:+d})",
              flush=True)

    with open(out_path, "w") as fh:
        json.dump(store, fh, indent=2)
    print(f"\nwrote {out_path}")
    return 0


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0],
                               formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="action", required=True)

    f = sub.add_parser("fold", help="Boltz-fold each sequence with no ligand (fast; run locally)")
    f.add_argument("outdir")
    f.add_argument("--shell", required=True, help="shell number, e.g. 3")

    e = sub.add_parser("energy", help="protonate, relax hydrogens, single point (slow; run on a GPU)")
    e.add_argument("outdir")
    e.add_argument("--shell", help="only this shell's apo folds")
    e.add_argument("--out", default="peptide_reference.json")
    e.add_argument("--fmax", type=float, default=0.10,
                   help="hydrogen relaxation cutoff; must match the complex's (default 0.10)")
    e.add_argument("--steps", type=int, default=75,
                   help="step cap; must match the complex's (default 75)")
    e.add_argument("--model", default="uma-s-1p2p1")
    e.add_argument("--device", default=None,
                   help="cpu or cuda; default is cuda when available")

    args = p.parse_args(argv)
    if args.action == "fold":
        return do_fold(args)
    if args.device is None:
        import torch
        args.device = "cuda" if torch.cuda.is_available() else "cpu"
    return do_energy(args)


if __name__ == "__main__":
    raise SystemExit(main())
