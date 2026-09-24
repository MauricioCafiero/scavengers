"""UMA binding energies for the Boltz-2 co-folded complexes.

For each predicted structure: IE = E(complex) - E(peptide) - E(ligand), all three at the
same geometry, so any strain in the predicted structure cancels and what is left is the
rigid interaction energy in kcal/mol.

Boltz writes heavy atoms only, so hydrogens are added first: the peptide with pdbfixer
(pH 7, so Arg/Lys are protonated, Asp/Glu deprotonated, termini charged), the ligand with
RDKit from the SMILES the run used. The net charge handed to UMA follows from that.

Usage:
    python code/uma_binding.py runs/octinoxate [--relax-h] [--limit N]

pdbfixer is called in a separate environment (--fixer-venv), so nothing needs installing here.
"""
import argparse
import csv
import glob
import os
import subprocess
import sys
import tempfile

import numpy as np

FIXER_VENV = os.path.expanduser("~/python_mac/pocket_assist/venv")

# side chains that carry a charge at pH 7, and the termini
CHARGED = {"ARG": 1, "LYS": 1, "ASP": -1, "GLU": -1}

FIXER_SCRIPT = """
import sys
from pdbfixer import PDBFixer
from openmm.app import PDBFile
fixer = PDBFixer(filename=sys.argv[1])
fixer.findMissingResidues()
fixer.findMissingAtoms()
fixer.addMissingAtoms()
fixer.addMissingHydrogens(float(sys.argv[3]))
PDBFile.writeFile(fixer.topology, fixer.positions, open(sys.argv[2], "w"), keepIds=True)
"""


def parse_cif(path):
    """Boltz model CIF -> (peptide atoms, ligand atoms) as lists of dicts."""
    peptide, ligand = [], []
    for line in open(path):
        if not line.startswith(("ATOM", "HETATM")):
            continue
        f = line.split()
        atom = {"element": f[2], "name": f[3], "resname": f[5], "resseq": f[7],
                "xyz": (float(f[10]), float(f[11]), float(f[12]))}
        (ligand if line.startswith("HETATM") else peptide).append(atom)
    return peptide, ligand


def write_pdb(atoms, path, hetatm=False):
    record = "HETATM" if hetatm else "ATOM  "
    with open(path, "w") as f:
        for i, a in enumerate(atoms, 1):
            x, y, z = a["xyz"]
            f.write(f"{record}{i:5d} {a['name']:<4s} {a['resname'][:3]:>3s} A{int(a['resseq']):4d}    "
                    f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {a['element']:>2s}\n")
        f.write("END\n")


def protonate_peptide(atoms, fixer_venv, ph=7.0):
    """Add hydrogens with pdbfixer; returns (symbols, positions)."""
    with tempfile.TemporaryDirectory() as tmp:
        raw, fixed = os.path.join(tmp, "in.pdb"), os.path.join(tmp, "out.pdb")
        write_pdb(atoms, raw)
        script = os.path.join(tmp, "fix.py")
        open(script, "w").write(FIXER_SCRIPT)
        out = subprocess.run([os.path.join(fixer_venv, "bin", "python"), script, raw, fixed, str(ph)],
                             capture_output=True, text=True)
        if not os.path.exists(fixed):
            raise RuntimeError(f"pdbfixer failed:\n{out.stdout}\n{out.stderr}")
        symbols, positions = [], []
        for line in open(fixed):
            if line.startswith(("ATOM", "HETATM")):
                element = line[76:78].strip() or line[12:16].strip()[0]
                symbols.append(element.capitalize())
                positions.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
    return symbols, np.array(positions)


def protonate_ligand(atoms, smiles):
    """Add hydrogens to the ligand using its known SMILES; returns (symbols, positions)."""
    from rdkit import Chem
    from rdkit.Chem import AllChem

    with tempfile.TemporaryDirectory() as tmp:
        path = os.path.join(tmp, "lig.pdb")
        write_pdb(atoms, path, hetatm=True)
        mol = Chem.MolFromPDBFile(path, removeHs=False, sanitize=False)
    template = Chem.MolFromSmiles(smiles)
    mol = AllChem.AssignBondOrdersFromTemplate(template, mol)
    Chem.SanitizeMol(mol)
    mol = Chem.AddHs(mol, addCoords=True)
    conf = mol.GetConformer()
    symbols = [a.GetSymbol() for a in mol.GetAtoms()]
    positions = np.array([list(conf.GetAtomPosition(i)) for i in range(mol.GetNumAtoms())])
    return symbols, positions, Chem.GetFormalCharge(mol)


def peptide_charge(atoms, cyclic=False):
    """Net charge at pH 7 from the residue composition, plus the termini when linear."""
    residues = {}
    for a in atoms:
        residues[a["resseq"]] = a["resname"]
    charge = sum(CHARGED.get(name, 0) for name in residues.values())
    if not cyclic:
        charge += 1 - 1  # NH3+ and COO- cancel
    return charge, len(residues)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("outdir", help="a peptide_builder run directory, e.g. runs/octinoxate")
    parser.add_argument("names", nargs="*", help="the structures to compute: a sequence name, or a path to a .cif")
    parser.add_argument("--all", action="store_true", help="every structure in the run that is not already done")
    parser.add_argument("--relax-h", action="store_true",
                        help="relax hydrogens in the complex first (heavy atoms fixed)")
    parser.add_argument("--force", action="store_true", help="recompute even if already in the csv")
    parser.add_argument("--model", default="uma-s-1p2p1", help="fairchem model (default: uma-s-1p2p1)")
    parser.add_argument("--fixer-venv", default=FIXER_VENV, help=f"venv with pdbfixer (default: {FIXER_VENV})")
    args = parser.parse_args(argv)

    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from boltz_check import ligand_smiles
    smiles = ligand_smiles(args.outdir)
    print(f"ligand SMILES: {smiles}")

    boltz_dir = os.path.join(args.outdir, "boltz")
    csv_path = os.path.join(boltz_dir, "uma_binding.csv")
    done = {}
    if os.path.exists(csv_path):
        with open(csv_path) as f:
            done = {r["name"]: r for r in csv.DictReader(f)}

    def cif_for(name):
        if name.endswith(".cif"):
            return name
        hits = glob.glob(os.path.join(boltz_dir, f"boltz_results_{name}", "predictions", "*", "*_model_0.cif"))
        return hits[0] if hits else None

    if args.names:
        cifs = []
        for name in args.names:
            path = cif_for(name)
            if path is None:
                print(f"{name}: no Boltz structure found, skipping")
            elif name in done and not args.force:
                print(f"{name}: already done (use --force to redo)")
            else:
                cifs.append(path)
    elif args.all:
        cifs = sorted(glob.glob(os.path.join(boltz_dir, "boltz_results_*", "predictions", "*", "*_model_0.cif")))
        if not args.force:
            cifs = [p for p in cifs if os.path.basename(p).replace("_model_0.cif", "") not in done]
    else:
        sys.exit("name the structures to compute, or pass --all")

    if not cifs:
        sys.exit("nothing to compute")
    print(f"{len(cifs)} structures to compute\n")

    if not cifs:
        print("nothing new to compute")
    import torch
    from fairchem.core import FAIRChemCalculator, pretrained_mlip
    from ase import Atoms
    from ase.constraints import FixAtoms
    from ase.optimize import BFGS

    device = "cuda" if torch.cuda.is_available() else "cpu"
    predictor = pretrained_mlip.get_predict_unit(args.model, device=device)
    calculator = FAIRChemCalculator(predictor, task_name="omol")

    def energy(symbols, positions, charge, spin=1):
        atoms = Atoms(symbols, positions=positions)
        atoms.info.update({"charge": int(charge), "spin": int(spin)})
        atoms.calc = calculator
        return atoms.get_potential_energy()

    results = [r for r in done.values() if r["name"] not in
               {os.path.basename(p).replace("_model_0.cif", "") for p in cifs}]
    for path in cifs:
        name = os.path.basename(path).replace("_model_0.cif", "")
        cyclic = name.startswith("cyclo_")
        pep_atoms, lig_atoms = parse_cif(path)
        try:
            pep_symbols, pep_xyz = protonate_peptide(pep_atoms, args.fixer_venv)
            lig_symbols, lig_xyz, lig_charge = protonate_ligand(lig_atoms, smiles)
        except Exception as exc:
            print(f"{name}: could not add hydrogens ({type(exc).__name__}: {exc})")
            continue
        pep_charge, n_res = peptide_charge(pep_atoms, cyclic)

        symbols = pep_symbols + lig_symbols
        positions = np.vstack([pep_xyz, lig_xyz])
        total_charge = pep_charge + lig_charge

        if args.relax_h:
            complex_atoms = Atoms(symbols, positions=positions)
            complex_atoms.info.update({"charge": int(total_charge), "spin": 1})
            complex_atoms.calc = calculator
            heavy = [i for i, s in enumerate(symbols) if s != "H"]
            complex_atoms.set_constraint(FixAtoms(indices=heavy))
            BFGS(complex_atoms, logfile=None).run(fmax=0.05)
            positions = complex_atoms.get_positions()
            pep_xyz, lig_xyz = positions[:len(pep_symbols)], positions[len(pep_symbols):]

        e_complex = energy(symbols, positions, total_charge)
        e_pep = energy(pep_symbols, pep_xyz, pep_charge)
        e_lig = energy(lig_symbols, lig_xyz, lig_charge)
        ie = 23.06035 * (e_complex - e_pep - e_lig)

        print(f"{name:<28} {n_res:>3} res, {len(symbols):>3} atoms, peptide charge {pep_charge:+d}  "
              f"IE {ie:8.2f} kcal/mol", flush=True)
        results.append({"name": name, "n_residues": n_res, "n_atoms": len(symbols),
                        "peptide_charge": pep_charge, "uma_ie_kcal_mol": f"{ie:.3f}"})

        with open(csv_path, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=list(results[0].keys()))
            writer.writeheader()
            writer.writerows(results)

    # line up with the Boltz affinities and the fragment sums
    comparison = os.path.join(boltz_dir, "comparison.csv")
    if os.path.exists(comparison) and len(results) > 2:
        with open(comparison) as f:
            boltz = {r["sequence"].replace("cyclo-", "cyclo_"): r for r in csv.DictReader(f)}
        rows = [(r, boltz[r["name"]]) for r in results if r["name"] in boltz]
        if len(rows) > 2:
            uma = np.array([float(r["uma_ie_kcal_mol"]) for r, _ in rows])
            pic50 = np.array([-float(b["pIC50"]) for _, b in rows])
            frag = np.array([float(b["total_ie_kcal_mol"]) for _, b in rows])
            rank = lambda v: np.argsort(np.argsort(v))
            print(f"\n{len(rows)} structures compared (Spearman rank correlation):")
            print(f"  UMA IE vs Boltz pIC50   {np.corrcoef(rank(uma), rank(pic50))[0, 1]:+.2f}")
            print(f"  UMA IE vs fragment sum  {np.corrcoef(rank(uma), rank(frag))[0, 1]:+.2f}")
            print(f"  fragment sum vs pIC50   {np.corrcoef(rank(frag), rank(pic50))[0, 1]:+.2f}")
    print(f"\nResults in {os.path.join(boltz_dir, 'uma_binding.csv')}")


if __name__ == "__main__":
    main()
