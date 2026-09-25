"""Binding energy for a folded complex: interaction, desolvation and strain.

`uma_binding.py` reports the rigid gas-phase interaction energy, and the traps list in
`NEXT_STEPS.md` is blunt about what that leaves out: no desolvation, no entropy, and a size bias
(Spearman -0.52 against peptide length). This adds the two terms that can be computed the same way
UMADock computes them, with the explicit-water machinery copied into `solvate.py`.

**Why the peptide term is the one that matters here.** UMADock desolvates the ligand, which is
right for its usual problem -- many ligands against one protein, where the protein's desolvation is
common to every candidate and cancels. Our problem is the transpose: many peptides against one
ligand. The ligand's desolvation is then near-constant across designs and cancels instead, and the
term carrying the variation is the **peptide's**. So both are computed here.

**Why not simply add the two full desolvations.** Stripping every water off the peptide and off the
ligand and adding the costs would badly overcount: in the complex only the buried interface is
desolvated, and for a peptide wrapped round a small ligand that is a small patch of a large
surface. The thermodynamic cycle handles it properly. With

    dE_solv(X) = E(X + waters) - E(X) - E(waters)          (negative; how much solvation helps X)

the desolvation penalty for forming the complex is

    desolvation = dE_solv(complex) - dE_solv(peptide) - dE_solv(ligand)      (positive)

which counts only the solvation that binding destroys. Strain is separate, because the interaction
energy is evaluated with both partners at their complex geometry and so contains none of it:

    strain(X) = E(X at its geometry in the complex) - E(X relaxed on its own)

and the total is

    binding = interaction + desolvation + strain(peptide) + strain(ligand)

**Two honest caveats.** The waters are placed at random, so every desolvation number is noisy;
`--replicates` runs the placement several times and reports the spread, which is the only way to
see how much of a difference between designs is real. And an explicit-water shell round a
200-atom peptide is a large system for an MLIP on CPU -- this is expensive, so the terms are
individually selectable.

Usage:
    python code/binding_energy.py runs/octinoxate --terms interaction,strain
    python code/binding_energy.py runs/octinoxate --terms all --replicates 3 --limit 2
"""
import argparse
import csv
import glob
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from peptide_builder import write_xyz  # noqa: E402
from solvate import add_waters, waters_for  # noqa: E402
from uma_binding import (FIXER_VENV, parse_cif, peptide_charge,  # noqa: E402
                         protonate_ligand, protonate_peptide)

EV_TO_KCAL = 23.06035


def relax_hydrogens(symbols, positions, charge, calculator, fmax=0.10, steps=75, spin=1,
                    label="", log=True):
    """Relax only the hydrogens, heavy atoms fixed. Returns the new positions.

    pdbfixer places rotatable hydrogens -- hydroxyls, amines -- wherever it likes, and it is not
    deterministic: the same folded structure scored twice gave interaction energies of -14.47 and
    -13.31 kcal/mol, a 1.2 kcal/mol spread from hydrogen placement alone. That is the same size as
    the difference the design-versus-random comparison rested on, so it is not a rounding detail.

    Relaxing the hydrogens with the heavy atoms held converges them to a common minimum and removes
    the arbitrariness, while leaving the predicted heavy-atom geometry exactly as Boltz produced it
    -- which is the thing being scored.

    Loose by default, on the same budget as the strain relaxations: the aim is to get the hydrogens
    off pdbfixer's arbitrary starting positions, not to minimise them precisely, and precision here
    would exceed what a gas-phase single point on a predicted structure can support anyway.
    """
    from ase import Atoms
    from ase.constraints import FixAtoms
    from ase.optimize import BFGS
    atoms = Atoms(symbols, positions=positions)
    atoms.info.update({"charge": int(charge), "spin": int(spin)})
    atoms.calc = calculator
    heavy = [i for i, sym in enumerate(symbols) if sym != "H"]
    atoms.set_constraint(FixAtoms(indices=heavy))
    if log:
        print(f"    relaxing {len(symbols) - len(heavy)} hydrogens of {label or 'system'} "
              f"({len(heavy)} heavy atoms fixed)", flush=True)
    opt = BFGS(atoms, logfile="-" if log else None)
    opt.run(fmax=fmax, steps=steps)
    if log:
        print(f"    -> {opt.get_number_of_steps()} steps", flush=True)
    return atoms.get_positions()


def energy(symbols, positions, charge, calculator, spin=1):
    from ase import Atoms
    atoms = Atoms(symbols, positions=positions)
    atoms.info.update({"charge": int(charge), "spin": int(spin)})
    atoms.calc = calculator
    return atoms.get_potential_energy()


def relaxed_energy(symbols, positions, charge, calculator, fmax=0.10, steps=75, spin=1,
                   label="", log=True):
    """BFGS relaxation at a deliberately loose convergence.

    0.05 eV/A over 400 steps is what the rest of the repo uses, and on a 200-atom glycine-rich
    peptide it does not finish -- 31 minutes produced no completed complex. At 0.10 over 75 steps
    the free reference is not a true minimum, so the strain it implies is an **underestimate**; that
    is acceptable while designs are being compared against each other on the same budget, and not
    acceptable if a strain value is ever quoted on its own.
    """

    from ase import Atoms
    from ase.optimize import BFGS
    atoms = Atoms(symbols, positions=positions)
    atoms.info.update({"charge": int(charge), "spin": int(spin)})
    atoms.calc = calculator
    # logfile='-' sends the step table to stdout. Passing None, as this did before, throws the
    # optimiser's progress away: a 31-minute relaxation then looks identical to a stuck one, with
    # no step count and no energy trace to tell them apart.
    if log:
        print(f"    relaxing {label or 'system'}: {len(atoms)} atoms, "
              f"fmax {fmax}, max {steps} steps", flush=True)
    opt = BFGS(atoms, logfile="-" if log else None)
    converged = opt.run(fmax=fmax, steps=steps)
    taken = opt.get_number_of_steps()
    if log:
        print(f"    -> {taken} steps, converged={bool(converged)}", flush=True)
    return atoms.get_potential_energy(), atoms.get_positions(), taken, bool(converged)


def solvation_energy(symbols, positions, charge, calculator, n_waters=None, replicates=1,
                     fmax=0.30, steps=200, seed=0, min_contact=2.4, verbose=True):
    """dE_solv = E(solute + waters) - E(solute) - E(waters), averaged over water placements.

    Everything is relaxed at a loose fmax, as UMADock does: the point is to let the shell settle,
    not to find a minimum. Returns (mean, spread, per-replicate list) in kcal/mol.
    """
    from ase import Atoms
    solute = Atoms(symbols, positions=positions)
    if n_waters is None:
        n_waters = waters_for(solute)
    values = []
    for rep in range(replicates):
        if verbose:
            print(f"    solvation replicate {rep + 1}/{replicates}, target {n_waters} waters")
        solvated, n_added = add_waters(solute, n_waters, min_contact=min_contact,
                                      seed=seed + rep, verbose=verbose)
        if not n_added:
            continue
        n_solute = len(solute)
        e_solvated, xyz, _, _ = relaxed_energy(solvated.get_chemical_symbols(),
                                               solvated.get_positions(), charge, calculator,
                                               fmax=fmax, steps=steps,
                                               label="solvated shell")
        # both parts evaluated at their geometry in the relaxed solvated system
        e_solute = energy(symbols, xyz[:n_solute], charge, calculator)
        water_symbols = solvated.get_chemical_symbols()[n_solute:]
        e_waters = energy(water_symbols, xyz[n_solute:], 0, calculator)
        values.append((e_solvated - e_solute - e_waters) * EV_TO_KCAL)
        if verbose:
            print(f"      dE_solv = {values[-1]:.2f} kcal/mol ({n_added} waters)")
    if not values:
        return float("nan"), float("nan"), []
    return float(np.mean(values)), float(np.std(values)), values


def cavity_desolvation(pep_symbols, pep_xyz, lig_xyz, pep_charge, calculator, replicates=1,
                       max_waters=40, region_radius=2.4, min_contact=2.4, fmax=0.30, steps=200,
                       seed=0, verbose=True):
    """The cost of emptying the binding cavity of water, in kcal/mol (positive = a penalty).

    Waters go **only where the ligand will sit**, not over the whole peptide. Water on the
    peptide's outer surface is there before and after binding, so it cancels; including it would
    mean extracting a small cavity term as the difference between two large numbers whose random
    placements do not cancel, and the noise would swamp the answer. So the ligand's own coordinates
    define the region, the peptide is present, and the ligand is not.

        penalty = E(peptide) + E(waters) - E(peptide + waters)

    evaluated with everything at the geometry of the relaxed water-filled cavity. What this leaves
    out is the bulk solvation those waters regain once displaced, which is roughly constant per
    water -- so compare designs at similar water counts, and treat `n_waters` as part of the result
    rather than a detail.
    """
    from ase import Atoms
    peptide = Atoms(pep_symbols, positions=pep_xyz)
    values, counts = [], []
    for rep in range(replicates):
        solvated, n_added = add_waters(peptide, max_waters, region=lig_xyz,
                                       region_radius=region_radius, min_contact=min_contact,
                                       seed=seed + rep, verbose=False, stop_on_streak=False)
        if not n_added:
            continue
        n_pep = len(peptide)
        _, xyz, _, _ = relaxed_energy(solvated.get_chemical_symbols(), solvated.get_positions(),
                                      pep_charge, calculator, fmax=fmax, steps=steps,
                                      label="peptide + cavity waters")
        e_solvated = energy(solvated.get_chemical_symbols(), xyz, pep_charge, calculator)
        e_pep = energy(pep_symbols, xyz[:n_pep], pep_charge, calculator)
        water_symbols = solvated.get_chemical_symbols()[n_pep:]
        e_waters = energy(water_symbols, xyz[n_pep:], 0, calculator)
        values.append((e_pep + e_waters - e_solvated) * EV_TO_KCAL)
        counts.append(n_added)
        if verbose:
            print(f"    cavity replicate {rep + 1}: {n_added} waters, "
                  f"penalty {values[-1]:.2f} kcal/mol")
    if not values:
        return float("nan"), float("nan"), 0
    return float(np.mean(values)), float(np.std(values)), int(round(float(np.mean(counts))))


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("outdir")
    parser.add_argument("--terms", default="interaction,strain_ligand,cavity",
                        help="comma-separated: interaction, strain_ligand, strain_peptide, cavity, "
                             "desolvation, or all. `strain_peptide` is off by default because in "
                             "vacuum it measures the peptide collapsing, not strain. "
                             "`cavity` empties only the binding site, which is what varies between "
                             "designs; `desolvation` runs the full cycle over the whole peptide and "
                             "is far more expensive and much noisier "
                             "(default: interaction,strain,cavity)")
    parser.add_argument("--no-relax-h", action="store_true",
                        help="skip the hydrogen-only relaxation (not advised: pdbfixer's hydrogen "
                             "placement alone moves the interaction energy by ~1 kcal/mol)")
    parser.add_argument("--fmax", type=float, default=0.10,
                        help="force convergence for the strain relaxations, eV/A (default 0.10)")
    parser.add_argument("--steps", type=int, default=75,
                        help="step cap for the strain relaxations (default 75)")
    parser.add_argument("--replicates", type=int, default=1,
                        help="water placements per solvation energy (default 1)")
    parser.add_argument("--waters", type=int, help="fixed water count (default: scaled by size)")
    parser.add_argument("--region-radius", type=float, default=2.4,
                        help="how far from a ligand atom a cavity water may sit, A (default 2.4)")
    parser.add_argument("--min-contact", type=float, default=2.4,
                        help="closest a water may be placed, A (UMADock used 1.7)")
    parser.add_argument("--limit", type=int, help="only the first N complexes")
    parser.add_argument("--structures", help="comma-separated structure names to use")
    parser.add_argument("--fixer-venv", default=FIXER_VENV)
    parser.add_argument("--model", default="uma-s-1p2p1")
    parser.add_argument("--out", default="binding_energy.csv")
    args = parser.parse_args(argv)

    terms = {t.strip() for t in args.terms.split(",")}
    if "all" in terms:
        terms = {"interaction", "strain_ligand", "cavity"}
    if "strain" in terms:                      # the old name meant both
        terms |= {"strain_peptide", "strain_ligand"}

    boltz_dir = os.path.join(args.outdir, "boltz")
    cifs = sorted(glob.glob(os.path.join(boltz_dir, "boltz_results_*", "predictions", "*",
                                         "*_model_0.cif")))
    if args.structures:
        wanted = {s.strip() for s in args.structures.split(",")}
        cifs = [c for c in cifs
                if os.path.basename(c).replace("_model_0.cif", "") in wanted]
    if args.limit:
        cifs = cifs[:args.limit]
    if not cifs:
        sys.exit(f"no folded complexes under {boltz_dir}")
    print(f"{len(cifs)} complexes; terms: {sorted(terms)}\n")

    import torch
    from fairchem.core import FAIRChemCalculator, pretrained_mlip
    device = "cuda" if torch.cuda.is_available() else "cpu"
    calculator = FAIRChemCalculator(pretrained_mlip.get_predict_unit(args.model, device=device),
                                    task_name="omol")

    smiles_path = os.path.join(args.outdir, "ligand.xyz")
    from boltz_check import ligand_smiles
    smiles = ligand_smiles(args.outdir)

    fields = ["name", "n_peptide_atoms", "n_ligand_atoms", "peptide_charge", "ligand_charge",
              "interaction_kcal", "strain_peptide_kcal", "strain_ligand_kcal",
              "desolv_complex_kcal", "desolv_peptide_kcal", "desolv_ligand_kcal",
              "desolvation_kcal", "desolv_spread_kcal",
              "cavity_kcal", "cavity_spread_kcal", "cavity_waters",
              "strain_peptide_steps", "strain_peptide_converged",
              "strain_ligand_steps", "strain_ligand_converged", "binding_kcal"]
    out_path = os.path.join(boltz_dir, args.out)
    with open(out_path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields)
        writer.writeheader()
        for path in cifs:
            name = os.path.basename(path).replace("_model_0.cif", "")
            print(f"{name}")
            pep_atoms, lig_atoms = parse_cif(path)
            try:
                pep_symbols, pep_xyz = protonate_peptide(pep_atoms, args.fixer_venv)
                lig_symbols, lig_xyz, lig_charge = protonate_ligand(lig_atoms, smiles)
            except Exception as exc:
                print(f"  could not add hydrogens ({type(exc).__name__}: {exc})")
                continue
            pep_charge, _ = peptide_charge(pep_atoms, name.startswith("cyclo_"))
            row = {"name": name, "n_peptide_atoms": len(pep_symbols),
                   "n_ligand_atoms": len(lig_symbols), "peptide_charge": pep_charge,
                   "ligand_charge": lig_charge}

            # A partial result is worth keeping: when the cavity term was cancelled mid-run the
            # interaction and strain values computed before it were lost from the csv, because the
            # row was only written once every requested term had finished. Each term now lands in
            # a sidecar as it is computed, so nothing already paid for is thrown away.
            partial = os.path.join(boltz_dir, f"partial_{name}.json")

            def save_partial():
                import json
                with open(partial, "w") as ph:
                    json.dump(row, ph, indent=2)

            all_symbols = pep_symbols + lig_symbols
            all_xyz = np.vstack([pep_xyz, lig_xyz])
            total_charge = pep_charge + lig_charge

            # Every relaxation here is paid for in hours, so the geometry it produces is written
            # out, not just the number derived from it. `condense_strain.py` has saved its built,
            # free and held structures from the start; this did not, and a full scoring pass was
            # run without keeping a single relaxed structure.
            struct_dir = os.path.join(boltz_dir, "structures")
            os.makedirs(struct_dir, exist_ok=True)

            if not args.no_relax_h:
                all_xyz = relax_hydrogens(all_symbols, all_xyz, total_charge, calculator,
                                          fmax=args.fmax, steps=args.steps, label="complex")
                pep_xyz, lig_xyz = all_xyz[:len(pep_symbols)], all_xyz[len(pep_symbols):]
                write_xyz(os.path.join(struct_dir, f"{name}_complex_relaxed_h.xyz"),
                          all_symbols, all_xyz,
                          f"{name}: heavy atoms as Boltz predicted, hydrogens relaxed "
                          f"(fmax {args.fmax}, max {args.steps} steps)")
                write_xyz(os.path.join(struct_dir, f"{name}_ligand_bound.xyz"),
                          lig_symbols, lig_xyz,
                          f"{name}: ligand in its bound pose, the reference for ligand strain")

            if "interaction" in terms:
                e_c = energy(all_symbols, all_xyz, total_charge, calculator)
                e_p = energy(pep_symbols, pep_xyz, pep_charge, calculator)
                e_l = energy(lig_symbols, lig_xyz, lig_charge, calculator)
                row["interaction_kcal"] = f"{(e_c - e_p - e_l) * EV_TO_KCAL:.3f}"
                print(f"  interaction {float(row['interaction_kcal']):9.2f} kcal/mol")
                save_partial()

            if "strain_peptide" in terms:
                # Kept available but off by default. On a 33-residue peptide this came out at
                # 237.8 kcal/mol *without converging* in 75 steps: relaxed in vacuum a peptide
                # collapses into a compact hydrogen-bonded ball, so the number measures how far
                # that collapse got on the step budget, not strain. It needs solvent to mean
                # anything.
                e_p_bound = energy(pep_symbols, pep_xyz, pep_charge, calculator)
                e_p_free, pep_free_xyz, st_p, cv_p = relaxed_energy(
                    pep_symbols, pep_xyz, pep_charge, calculator,
                    fmax=args.fmax, steps=args.steps, label="free peptide")
                write_xyz(os.path.join(struct_dir, f"{name}_peptide_relaxed.xyz"),
                          pep_symbols, pep_free_xyz,
                          f"{name}: peptide relaxed free, {st_p} steps, converged={bool(cv_p)}")
                row["strain_peptide_kcal"] = f"{(e_p_bound - e_p_free) * EV_TO_KCAL:.3f}"
                row["strain_peptide_steps"] = f"{st_p}"
                row["strain_peptide_converged"] = int(cv_p)
                print(f"  strain peptide {float(row['strain_peptide_kcal']):8.2f} kcal/mol "
                      f"({st_p} steps, converged={bool(cv_p)})")

            if "strain_ligand" in terms:
                e_l_bound = energy(lig_symbols, lig_xyz, lig_charge, calculator)
                e_l_free, lig_free_xyz, st_l, cv_l = relaxed_energy(
                    lig_symbols, lig_xyz, lig_charge, calculator,
                    fmax=args.fmax, steps=args.steps, label="free ligand")
                write_xyz(os.path.join(struct_dir, f"{name}_ligand_relaxed.xyz"),
                          lig_symbols, lig_free_xyz,
                          f"{name}: ligand relaxed free, {st_l} steps, converged={bool(cv_l)}; "
                          f"against _ligand_bound.xyz this is the strain")
                row["strain_ligand_kcal"] = f"{(e_l_bound - e_l_free) * EV_TO_KCAL:.3f}"
                row["strain_ligand_steps"] = f"{st_l}"
                row["strain_ligand_converged"] = int(cv_l)
                print(f"  strain ligand  {float(row['strain_ligand_kcal']):8.2f} kcal/mol "
                      f"({st_l} steps, converged={bool(cv_l)})")
                save_partial()

            if "desolvation" in terms:
                kw = dict(calculator=calculator, replicates=args.replicates,
                          n_waters=args.waters, min_contact=args.min_contact)
                print("  complex:")
                d_c, s_c, _ = solvation_energy(all_symbols, all_xyz, total_charge, **kw)
                print("  peptide:")
                d_p, s_p, _ = solvation_energy(pep_symbols, pep_xyz, pep_charge, **kw)
                print("  ligand:")
                d_l, s_l, _ = solvation_energy(lig_symbols, lig_xyz, lig_charge, **kw)
                desolv = d_c - d_p - d_l
                row.update({"desolv_complex_kcal": f"{d_c:.3f}",
                            "desolv_peptide_kcal": f"{d_p:.3f}",
                            "desolv_ligand_kcal": f"{d_l:.3f}",
                            "desolvation_kcal": f"{desolv:.3f}",
                            "desolv_spread_kcal": f"{max(s_c, s_p, s_l):.3f}"})
                print(f"  desolvation {desolv:9.2f} kcal/mol "
                      f"(worst replicate spread {max(s_c, s_p, s_l):.2f})")

            if "cavity" in terms:
                cav, spread, n_w = cavity_desolvation(
                    pep_symbols, pep_xyz, lig_xyz, pep_charge, calculator,
                    replicates=args.replicates, max_waters=args.waters or 40,
                    min_contact=args.min_contact, region_radius=args.region_radius)
                row.update({"cavity_kcal": f"{cav:.3f}", "cavity_spread_kcal": f"{spread:.3f}",
                            "cavity_waters": n_w})
                print(f"  cavity desolvation {cav:8.2f} kcal/mol "
                      f"({n_w} waters, spread {spread:.2f})")
                save_partial()

            parts = [float(row[k]) for k in ("interaction_kcal", "strain_peptide_kcal",
                                             "strain_ligand_kcal", "desolvation_kcal",
                                             "cavity_kcal")
                     if row.get(k)]
            if "interaction" in terms:
                row["binding_kcal"] = f"{sum(parts):.3f}"
                print(f"  binding {sum(parts):13.2f} kcal/mol"
                      + ("" if terms & {"desolvation", "cavity"} else "   (no desolvation term)"))
            writer.writerow(row)
            fh.flush()
    print(f"\nresults in {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
