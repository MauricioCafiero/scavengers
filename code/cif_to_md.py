"""Split a Boltz complex into the two files the OpenMM MD pipeline wants.

Boltz writes one CIF holding the peptide and the ligand together, heavy atoms only. `omd` (the pipeline
in `~/python_mac/openmm`) wants them apart: a protein PDB for `prep-protein`, which runs PDBFixer over
it, and a ligand SDF for `prep-ligand`, which needs correct bond orders to hand to GAFF2 or OpenFF.

The ligand is the part that needs care. A PDB or CIF records no bond orders, so the SDF is built the
same way `uma_binding.protonate_ligand` does it -- read the heavy atoms, assign bond orders from the
run's own SMILES as a template, then add hydrogens with coordinates. Without the template step the
aromatic ring and the ester come out as single bonds and the ligand is parameterised as something it
is not.

Hydrogens are added to the ligand here because the SDF has to carry them for parameterisation, and
they are placed geometrically from the heavy atoms, which Boltz fixed. The protein's hydrogens are
left to PDBFixer inside `prep-protein`, which is where that pipeline expects to add them.

Usage:
    python code/cif_to_md.py runs/octinoxate --structure orig_f12 --out-dir md/orig_f12
"""
import argparse
import glob
import os
import sys


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0],
                               formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("outdir", help="a run directory, e.g. runs/octinoxate")
    p.add_argument("--structure", required=True, help="fold name, e.g. orig_f12 or s3_orig_f12")
    p.add_argument("--out-dir", help="where to write (default <outdir>/md/<structure>)")
    args = p.parse_args(argv)

    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from uma_binding import parse_cif, write_pdb, protonate_ligand
    from boltz_check import ligand_smiles
    from rdkit import Chem

    name = args.structure
    hits = glob.glob(os.path.join(args.outdir, "boltz", f"boltz_results_{name}",
                                  "predictions", name, f"{name}_model_0.cif"))
    if not hits:
        sys.exit(f"no folded complex for {name!r} under {args.outdir}/boltz")
    cif = hits[0]

    out = args.out_dir or os.path.join(args.outdir, "md", name)
    os.makedirs(out, exist_ok=True)

    pep, lig = parse_cif(cif)
    smiles = ligand_smiles(args.outdir)
    print(f"{name}: {len(pep)} peptide atoms, {len(lig)} ligand atoms")
    print(f"ligand SMILES: {smiles}")

    pdb = os.path.join(out, f"{name}_protein.pdb")
    write_pdb(pep, pdb)
    print(f"wrote {pdb}  (heavy atoms only -- prep-protein adds hydrogens)")

    # Bond orders from the SMILES template, then hydrogens with coordinates. Same route as
    # uma_binding.protonate_ligand, but kept as an RDKit mol so it can be written as an SDF.
    tmp_pdb = os.path.join(out, f"{name}_ligand_heavy.pdb")
    write_pdb(lig, tmp_pdb, hetatm=True)
    mol = Chem.MolFromPDBFile(tmp_pdb, removeHs=False, sanitize=False)
    template = Chem.MolFromSmiles(smiles)
    from rdkit.Chem import AllChem
    mol = AllChem.AssignBondOrdersFromTemplate(template, mol)
    Chem.SanitizeMol(mol)
    mol = Chem.AddHs(mol, addCoords=True)
    mol.SetProp("_Name", f"{name}_ligand")
    sdf = os.path.join(out, f"{name}_ligand.sdf")
    with Chem.SDWriter(sdf) as w:
        w.write(mol)
    os.remove(tmp_pdb)

    formula = Chem.rdMolDescriptors.CalcMolFormula(mol)
    n_arom = sum(1 for a in mol.GetAtoms() if a.GetIsAromatic())
    print(f"wrote {sdf}  ({mol.GetNumAtoms()} atoms, {formula}, {n_arom} aromatic, "
          f"charge {Chem.GetFormalCharge(mol):+d})")
    # If the formula is wrong the template match silently mis-assigned, and everything downstream
    # would be parameterising the wrong molecule.
    print("\ncheck the formula against the run's ligand before going further.")
    print(f"\nnext, in the openmm-md environment:")
    omd = "~/miniforge3/envs/openmm-md/bin/omd"
    print(f"  {omd} prep-protein --pdb {pdb} --out {out}/protein_fixed.pdb")
    print(f"  {omd} prep-ligand  --sdf {sdf} --out {out}/ligand_prepped.sdf")
    print(f"  {omd} build --protein {out}/protein_fixed.pdb --ligand {out}/ligand_prepped.sdf "
          f"--out-dir {out}/system --no-auto-cofactors")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
