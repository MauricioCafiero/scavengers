"""Co-fold the peptides designed by peptide_builder with their ligand in Boltz-2,
and compare Boltz's affinity ranking with the fragment interaction energies.

The ligand SMILES is derived from the run's own ligand.xyz, so it always matches the
structure the fragments were sampled against. MSAs are skipped (`msa: empty`): these
are designed peptides with no natural homologs, and it keeps them off the public MSA
server.

Usage:
    python code/boltz_check.py runs/octinoxate [--limit N] [--sequences S,S] [--force]

Needs the Boltz environment from the parallel boltz_local repo (--boltz-repo).
"""
import argparse
import csv
import json
import os
import subprocess
import sys

# resolved by boltz_env at call time; unset unless the environment names one
BOLTZ_REPO = os.environ.get("PEPTIDEBUILDER_BOLTZ_REPO")

YAML = """version: 1
sequences:
  - protein:
      id: A
      sequence: {sequence}
      msa: empty{cyclic}
  - ligand:
      id: B
      smiles: '{smiles}'
properties:
  - affinity:
      binder: B
"""


def ligand_smiles(outdir: str):
    """SMILES for the ligand of a run, perceived from its saved geometry."""
    from rdkit import Chem
    from rdkit.Chem import rdDetermineBonds

    with open(os.path.join(outdir, "ligand.xyz")) as f:
        mol = Chem.MolFromXYZBlock(f.read())
    rdDetermineBonds.DetermineBonds(mol, charge=0)
    return Chem.MolToSmiles(Chem.RemoveHs(mol))


def run_boltz(name: str, sequence: str, smiles: str, cyclic: bool, boltz_dir: str,
              boltz_repo: str = None, boltz_cmd: str = None, boltz_venv: str = None):
    """Write the yaml for one peptide, co-fold it, and return its affinity results.

    `boltz_env.resolve_boltz` finds Boltz however it is installed here; `boltz_repo` is the old
    hardwired argument, honoured as a venv if a caller still passes one.
    """
    from boltz_env import run_boltz as _run

    yaml_path = os.path.join(boltz_dir, f"{name}.yaml")
    with open(yaml_path, "w") as f:
        f.write(YAML.format(sequence=sequence, smiles=smiles,
                            cyclic="\n      cyclic: true" if cyclic else ""))

    log_path = os.path.join(boltz_dir, f"{name}.log")
    code = _run(yaml_path, boltz_dir, log_path, boltz_cmd=boltz_cmd,
                boltz_venv=boltz_venv or (os.path.join(boltz_repo, ".venv")
                                          if boltz_repo else None))
    if code:
        print(f"  {name}: boltz failed (exit {code}), see {log_path}")
        return None

    path = os.path.join(boltz_dir, f"boltz_results_{name}", "predictions", name, f"affinity_{name}.json")
    if not os.path.exists(path):
        print(f"  {name}: no affinity written, see {log_path}")
        return None
    with open(path) as f:
        data = json.load(f)
    # Boltz reports affinity_pred_value as log10(IC50 in uM)
    pIC50 = 6 - data["affinity_pred_value"]
    return {"pIC50": pIC50, "dG_kcal_mol": -1.364 * pIC50,
            "binder_probability": data["affinity_probability_binary"]}


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("outdir", help="a peptide_builder run directory, e.g. runs/octinoxate")
    parser.add_argument("--limit", type=int, help="only the N best-scoring sequences")
    parser.add_argument("--sequences", help="comma-separated sequences to run instead of the csv order")
    parser.add_argument("--force", action="store_true", help="re-run sequences already in the comparison csv")
    parser.add_argument("--boltz-repo", default=BOLTZ_REPO,
                        help="a boltz checkout to use instead of however boltz is installed here")
    parser.add_argument("--boltz-venv", default=os.environ.get("PEPTIDEBUILDER_BOLTZ_VENV"),
                        help="venv with Boltz (or $PEPTIDEBUILDER_BOLTZ_VENV)")
    parser.add_argument("--boltz-cmd", default=os.environ.get("PEPTIDEBUILDER_BOLTZ_CMD"),
                        help="complete command that runs Boltz")
    args = parser.parse_args(argv)

    with open(os.path.join(args.outdir, "sequences.csv")) as f:
        rows = list(csv.DictReader(f))
    if args.sequences:
        wanted = args.sequences.split(",")
        rows = [r for r in rows if r["sequence"].replace("cyclo-", "") in wanted or r["sequence"] in wanted]
    rows.sort(key=lambda r: float(r["total_ie_kcal_mol"]))
    if args.limit:
        rows = rows[:args.limit]

    boltz_dir = os.path.join(args.outdir, "boltz")
    os.makedirs(boltz_dir, exist_ok=True)
    csv_path = os.path.join(boltz_dir, "comparison.csv")
    done = {}
    if os.path.exists(csv_path) and not args.force:
        with open(csv_path) as f:
            done = {r["sequence"]: r for r in csv.DictReader(f)}

    smiles = ligand_smiles(args.outdir)
    print(f"ligand SMILES (from {args.outdir}/ligand.xyz): {smiles}")
    print(f"{len(rows)} sequences to co-fold\n")

    results = dict(done)
    for i, row in enumerate(rows, 1):
        label = row["sequence"]
        if label in done:
            print(f"[{i}/{len(rows)}] {label}: already done, skipping")
            continue
        cyclic = label.startswith("cyclo-")
        sequence = label.replace("cyclo-", "")
        name = ("cyclo_" if cyclic else "") + sequence
        print(f"[{i}/{len(rows)}] {label} ({len(sequence)} residues){' cyclic' if cyclic else ''}", flush=True)
        out = run_boltz(name, sequence, smiles, cyclic, boltz_dir, args.boltz_repo,
                        boltz_cmd=args.boltz_cmd, boltz_venv=args.boltz_venv)
        if out is None:
            continue
        results[label] = {"sequence": label, "n_residues": len(sequence),
                          "total_ie_kcal_mol": row["total_ie_kcal_mol"],
                          "pIC50": f"{out['pIC50']:.3f}",
                          "boltz_dG_kcal_mol": f"{out['dG_kcal_mol']:.3f}",
                          "binder_probability": f"{out['binder_probability']:.3f}"}
        print(f"  pIC50 {out['pIC50']:.2f}, dG {out['dG_kcal_mol']:.2f} kcal/mol, "
              f"binder probability {out['binder_probability']:.2f}", flush=True)

        fields = ["sequence", "n_residues", "total_ie_kcal_mol", "pIC50", "boltz_dG_kcal_mol", "binder_probability"]
        with open(csv_path, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=fields)
            writer.writeheader()
            for rec in sorted(results.values(), key=lambda r: float(r["total_ie_kcal_mol"])):
                writer.writerow(rec)

    if len(results) > 2:
        import numpy as np
        ours = np.array([float(r["total_ie_kcal_mol"]) for r in results.values()])
        theirs = np.array([-float(r["pIC50"]) for r in results.values()])
        rank = lambda v: np.argsort(np.argsort(v))
        ro, rt = rank(ours), rank(theirs)
        spearman = np.corrcoef(ro, rt)[0, 1]
        print(f"\n{len(results)} sequences compared; Spearman rank correlation "
              f"(fragment IE vs Boltz pIC50): {spearman:+.2f}")
    print(f"Results in {csv_path}")


if __name__ == "__main__":
    main()
