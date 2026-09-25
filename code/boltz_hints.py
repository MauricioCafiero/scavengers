"""Co-fold a condensation design in Boltz-2 with the designed shell given as contact hints.

`boltz_check.py` folds a sequence and hopes the designed side-chain arrangement re-emerges. It
mostly does not. Boltz-2 will however accept **contact constraints** -- a peptide residue and a
ligand atom, with a maximum separation -- and the fragment search already knows exactly which
ligand atoms each designed side chain was placed against. So the shell can be asked for instead of
hoped for.

The two runs answer different questions and both are worth having:

    unconstrained   does the *sequence* encode the shell?      (the design question)
    constrained     can *any* backbone hold the shell?         (the feasibility question)

A constrained fold that reproduces the arrangement proves nothing about the sequence -- the answer
was supplied. What it does test is whether the geometry is reachable at all, by a route completely
independent of this repo's own backbone code. If Boltz cannot satisfy the contacts either, that is
real corroboration; if it can, there is a structure to relax and score.

Ligand atom naming is the fiddly part. Boltz names the atoms of a SMILES ligand
`element + canonical rank`, computed after `standardize()` (applied whenever affinity is
requested) and `AddHs`. Reproducing that by hand would be guesswork, so the names are taken from
Boltz's own code path, run in Boltz's own environment, and mapped back onto the run's `ligand.xyz`
by substructure match.

Usage:
    python code/assign.py runs/octinoxate --pairs-csv condense_shell12.csv \
        --save runs/octinoxate/design.json
    python code/boltz_hints.py runs/octinoxate design.json --cutoff 4.5
"""
import argparse
import json
import os
import subprocess
import sys
import tempfile

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from boltz_check import ligand_smiles  # noqa: E402
from peptide_builder import ONE_LETTER, define_fragments, load_ligand  # noqa: E402

# Boltz is resolved by boltz_env at call time. These stay only for the affinity naming path,
# which needs Boltz's own standardize(), and are unset unless the environment names a venv.
BOLTZ_VENV = os.environ.get("PEPTIDEBUILDER_BOLTZ_VENV")
BOLTZ_REPO = os.environ.get("PEPTIDEBUILDER_BOLTZ_REPO")

# Runs inside the Boltz environment: name the ligand atoms exactly as Boltz will, then map those
# names onto the atom order of the run's own ligand.xyz.
NAMER = r'''
import json, sys
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdDetermineBonds
from boltz.data.parse.schema import standardize

smiles, xyz_path, affinity = sys.argv[1], sys.argv[2], sys.argv[3] == "1"
seq = standardize(smiles) if affinity else smiles
mol = AllChem.MolFromSmiles(seq)
mol = AllChem.AddHs(mol)
order = AllChem.CanonicalRankAtoms(mol)
Chem.AssignStereochemistry(mol, force=True, cleanIt=True)
names = [a.GetSymbol().upper() + str(i + 1) for a, i in zip(mol.GetAtoms(), order)]

ref = Chem.MolFromXYZFile(xyz_path)
rdDetermineBonds.DetermineBonds(ref, charge=0)
match = mol.GetSubstructMatch(ref)          # match[i] = index in `mol` of xyz atom i
if not match:
    ref_noH, mol_noH = Chem.RemoveHs(ref), Chem.RemoveHs(mol)
    raise SystemExit(json.dumps({"error": "no substructure match between ligand.xyz and the SMILES",
                                 "xyz_atoms": ref.GetNumAtoms(), "smiles_atoms": mol.GetNumAtoms()}))
print(json.dumps({"names_in_xyz_order": [names[j] for j in match]}))
'''


def ligand_atom_names_local(outdir, smiles):
    """Boltz's atom names for a SMILES ligand, computed here with RDKit alone.

    Boltz names them `element + canonical rank`, and the only Boltz-specific step in that path is
    `standardize()` -- which it calls **only when affinity is requested**. We never request affinity
    (its head is broken against this RDKit and unreliable anyway), so no Boltz code is involved and
    there is no reason to shell out to its environment for this. Verified identical across all 44
    atoms of octinoxate against the subprocess version below.

    Returns the names in the atom order of the run's own `ligand.xyz`, matched by substructure.
    """
    from rdkit import Chem
    from rdkit.Chem import AllChem, rdDetermineBonds
    mol = AllChem.AddHs(AllChem.MolFromSmiles(smiles))
    order = AllChem.CanonicalRankAtoms(mol)
    Chem.AssignStereochemistry(mol, force=True, cleanIt=True)
    names = [a.GetSymbol().upper() + str(i + 1) for a, i in zip(mol.GetAtoms(), order)]
    ref = Chem.MolFromXYZFile(os.path.join(outdir, "ligand.xyz"))
    rdDetermineBonds.DetermineBonds(ref, charge=0)
    match = mol.GetSubstructMatch(ref)
    if not match:
        raise RuntimeError("ligand.xyz does not match the SMILES as a substructure")
    return [names[j] for j in match]


def ligand_atom_names(outdir, smiles, affinity=False, boltz_venv=BOLTZ_VENV):
    """Atom names for the run's ligand, in the atom order of ligand.xyz.

    Computed locally unless `affinity` is set, which is the one case where Boltz's own
    `standardize()` changes the naming and its environment therefore has to be used.
    """
    if not affinity:
        return ligand_atom_names_local(outdir, smiles)
    return _ligand_atom_names_via_boltz(outdir, smiles, affinity, boltz_venv)


def _ligand_atom_names_via_boltz(outdir, smiles, affinity=True, boltz_venv=BOLTZ_VENV):
    """The naming as Boltz itself computes it, run in Boltz's environment.

    Only needed for the affinity path, and kept as the reference the local version was checked
    against.
    """
    python = os.path.join(boltz_venv, "bin", "python")
    if not os.path.exists(python):
        sys.exit(f"no Boltz environment at {boltz_venv}")
    with tempfile.NamedTemporaryFile("w", suffix=".py", delete=False) as fh:
        fh.write(NAMER)
        script = fh.name
    try:
        out = subprocess.run([python, script, smiles, os.path.join(outdir, "ligand.xyz"),
                              "1" if affinity else "0"],
                             capture_output=True, text=True)
    finally:
        os.unlink(script)
    if out.returncode != 0:
        sys.exit(f"naming the ligand atoms failed:\n{out.stdout}\n{out.stderr}")
    return json.loads(out.stdout.strip())["names_in_xyz_order"]


def residue_positions(design):
    """1-based residue index in the sequence for each pose on the path.

    The path alternates a designed residue and its glycine spacers, so position i+1 follows
    position i by the number of linkers on that step plus one.
    """
    positions, at = [], 1
    for step, (name, pose) in enumerate(design["path"]):
        positions.append((at, name, pose))
        if step < len(design["linkers"]):
            at += design["linkers"][step] + 1
    return positions


def contacts(outdir, design, cutoff, names):
    """One contact per designed residue: the ligand atom its fragment pose sits closest to.

    Only the closest is used. Pinning every contact a side chain makes would over-constrain the
    fold into simply reproducing the input, and the point is to see whether a backbone can get
    there, not to dictate the whole arrangement.
    """
    import ase.io
    frags = {f["name"]: f for f in define_fragments()}
    ligand = load_ligand(design.get("ligand", "octinoxate"))
    lig_xyz = np.asarray(ligand["coords"])
    heavy = [i for i, s in enumerate(ligand["atoms"]) if s != "H"]

    out = []
    for res_index, name, pose in residue_positions(design):
        path = os.path.join(outdir, "poses", f"{ligand['name']}_w_{name}{pose}.xyz")
        coords = ase.io.read(path).get_positions()[ligand["num_atoms"]:]
        frag_heavy = [i for i, s in enumerate(frags[name]["atoms"]) if s != "H"]
        d = np.linalg.norm(coords[frag_heavy][:, None, :] - lig_xyz[heavy][None, :, :], axis=-1)
        fi, li = np.unravel_index(d.argmin(), d.shape)
        distance = float(d[fi, li])
        if distance <= cutoff:
            out.append({"residue": res_index, "letter": ONE_LETTER[name],
                        "fragment": f"{name}{pose}", "atom": names[heavy[li]],
                        "distance": round(distance, 2)})
    return out


def write_yaml(path, sequence, smiles, hints, margin, force, cyclic=False, affinity=False):
    """Boltz-2 input for one design.

    `affinity` is off by default: Boltz 2.2.1 computes the binder's molecular weight with
    `AllChem.Descriptors.MolWt`, and RDKit 2026.03.6 no longer exposes `Descriptors` through
    `AllChem`, so asking for affinity fails at parse time in this installation. Dropping it costs
    little -- the earlier work already found that head barely discriminates and is not
    geometry-consistent -- and the interaction energy is scored with UMA here anyway. It also
    changes ligand atom naming, since Boltz only calls `standardize()` when affinity is requested,
    which is why the namer takes the same flag.
    """
    lines = ["version: 1", "sequences:", "  - protein:", "      id: A",
             f"      sequence: {sequence}", "      msa: empty"]
    if cyclic:
        lines.append("      cyclic: true")
    lines += ["  - ligand:", "      id: B", f"      smiles: '{smiles}'"]
    if hints:
        lines.append("constraints:")
        for h in hints:
            lines += ["  - contact:",
                      f"      token1: [A, {h['residue']}]",
                      f"      token2: [B, {h['atom']}]",
                      f"      max_distance: {h['distance'] + margin:.1f}",
                      f"      force: {'true' if force else 'false'}"]
    if affinity:
        lines += ["properties:", "  - affinity:", "      binder: B"]
    lines.append("")
    with open(path, "w") as fh:
        fh.write("\n".join(lines))


def fold(yaml_path, boltz_dir, boltz_repo=None, extra_args=(), boltz_cmd=None, boltz_venv=None):
    """Co-fold one yaml, however Boltz is available on this machine.

    `boltz_env.resolve_boltz` works out the command: an explicit one, a named venv, boltz installed
    here, boltz on PATH, or the MPS wrapper if that layout happens to exist. `boltz_repo` is the old
    hardwired argument, honoured as a venv if a caller still passes one.

    `extra_args` goes straight to Boltz's own predict command. The one that matters here is
    `--use_potentials`: a contact constraint with `force: true` becomes a `ContactPotential`, but
    the guidance loop that applies potentials during sampling is gated on `fk_steering` and
    `physical_guidance_update`, and both are off unless that flag is passed. Without it a forced
    contact is conditioning only -- the model is told, and may decline. (With `force: false` the
    constraint is dropped by the featurizer outright, so it is not even conditioning.)
    """
    from boltz_env import run_boltz

    name = os.path.splitext(os.path.basename(yaml_path))[0]
    log_path = os.path.join(boltz_dir, f"{name}.log")
    print(f"folding {name} -> {log_path}", flush=True)
    code = run_boltz(yaml_path, boltz_dir, log_path, extra_args=extra_args,
                     boltz_cmd=boltz_cmd,
                     boltz_venv=boltz_venv or (os.path.join(boltz_repo, ".venv")
                                               if boltz_repo else None))
    if code:
        print(f"  boltz failed (exit {code}), see {log_path}")
        return None
    cif = os.path.join(boltz_dir, f"boltz_results_{name}", "predictions", name,
                       f"{name}_model_0.cif")
    print(f"  {'wrote ' + cif if os.path.exists(cif) else 'no structure written, see ' + log_path}")
    return cif if os.path.exists(cif) else None


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("outdir", help="the run directory, e.g. runs/octinoxate")
    parser.add_argument("design", help="design JSON from `assign.py --save`")
    parser.add_argument("--sequence", help="use this sequence instead of the design's "
                                           "(for an ESM2-substituted variant of the same design)")
    parser.add_argument("--cutoff", type=float, default=4.5,
                        help="only hint contacts the pose already makes within this, A "
                             "(default 4.5)")
    parser.add_argument("--margin", type=float, default=1.0,
                        help="added to each measured distance to give Boltz room (default 1.0 A)")
    parser.add_argument("--force", action="store_true",
                        help="enforce the contacts rather than hinting at them")
    parser.add_argument("--name", help="output basename (default: the sequence)")
    parser.add_argument("--no-hints", action="store_true",
                        help="write the same design with no constraints, as the control")
    parser.add_argument("--top", type=int,
                        help="hint only the N tightest designed contacts, so that forcing them "
                             "and failing points at something specific")
    parser.add_argument("--residues", help="hint only these residue numbers, comma-separated")
    parser.add_argument("--run", action="store_true", help="co-fold it straight away")
    parser.add_argument("--boltz-args", default="",
                        help="extra arguments for Boltz predict, e.g. '--use_potentials' to "
                             "actually apply forced contacts during sampling")
    parser.add_argument("--affinity", action="store_true",
                        help="also ask for the affinity head (broken in this Boltz/RDKit pair)")
    parser.add_argument("--boltz-venv", default=BOLTZ_VENV,
                        help="venv with Boltz (or $PEPTIDEBUILDER_BOLTZ_VENV); not needed if "
                             "boltz is installed here or on PATH")
    parser.add_argument("--boltz-cmd", default=os.environ.get("PEPTIDEBUILDER_BOLTZ_CMD"),
                        help="complete command that runs Boltz, e.g. 'python -m boltz'")
    args = parser.parse_args(argv)

    design_path = args.design
    if not os.path.exists(design_path):
        design_path = os.path.join(args.outdir, args.design)
    design = json.load(open(design_path))
    sequence = args.sequence or design["sequence"]
    if args.sequence and len(args.sequence) != len(design["sequence"]):
        sys.exit(f"substituted sequence is {len(args.sequence)} residues but the design is "
                 f"{len(design['sequence'])}; the residue numbering would not line up")

    smiles = ligand_smiles(args.outdir)
    print(f"ligand SMILES: {smiles}")
    names = ligand_atom_names(args.outdir, smiles, affinity=args.affinity,
                              boltz_venv=args.boltz_venv)
    print(f"{len(names)} ligand atoms named by Boltz's own code path\n")

    hints = [] if args.no_hints else contacts(args.outdir, design, args.cutoff, names)
    if hints and args.residues:
        wanted = {int(r) for r in args.residues.split(",")}
        hints = [h for h in hints if h["residue"] in wanted]
    elif hints and args.top:
        # the tightest designed contacts first: those are the ones the pose is most committed to,
        # and forcing a handful makes a failure attributable where forcing all of them does not
        hints = sorted(hints, key=lambda h: h["distance"])[:args.top]
        hints.sort(key=lambda h: h["residue"])
    print(f"sequence: {sequence}  ({len(sequence)} residues)")
    if hints:
        print(f"{len(hints)} contact hints (closest ligand atom per designed residue):")
        for h in hints:
            print(f"  residue {h['residue']:>3} {h['letter']}  ({h['fragment']:<14}) "
                  f"-> ligand {h['atom']:<4} at {h['distance']:.2f} A "
                  f"-> max_distance {h['distance'] + args.margin:.1f}")
    else:
        print("no contact hints (control run)" if args.no_hints else "no contacts within cutoff")

    boltz_dir = os.path.join(args.outdir, "boltz")
    os.makedirs(boltz_dir, exist_ok=True)
    base = args.name or (sequence + ("_nohints" if args.no_hints else "_hints"))
    out = os.path.join(boltz_dir, base + ".yaml")
    write_yaml(out, sequence, smiles, hints, args.margin, args.force, affinity=args.affinity)
    print(f"\nwrote {out}")
    if args.run:
        fold(out, boltz_dir, extra_args=tuple(args.boltz_args.split()) if args.boltz_args else ())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
