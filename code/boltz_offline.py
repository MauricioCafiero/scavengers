"""Run Boltz-2 somewhere else: export the inputs, fold them anywhere, bring the results back.

The rest of this repo assumes Boltz is installed alongside it, which suits the machine it was written
on and anyone who sets up both repositories the same way. This gives the other route. `export` writes
every Boltz input into one self-contained folder with the command to run and the output names to
expect; you fold them wherever you like — a GPU box, Colab, a cluster queue — and `import` puts the
structures back where the scoring and analysis scripts look for them.

It is also the practical way to use a GPU. Folding is quick anywhere, but scoring a 33-residue complex
takes 30 to 70 minutes on CPU against a couple of minutes on an A100, so the same split is worth using
for the energies: fold and score elsewhere, bring back the CIFs and the CSVs.

Usage:
    python code/boltz_offline.py export runs/octinoxate --design design.json --variants 2
    #  ... copy runs/octinoxate/boltz_inputs/ elsewhere, follow its RUN.md ...
    python code/boltz_offline.py import runs/octinoxate --from ~/Downloads/boltz_out

`import` accepts whatever layout Boltz produced — its own `boltz_results_*/predictions/...` tree, a
flat directory of `.cif` files, or a zip of either.
"""
import argparse
import glob
import json
import os
import shutil
import sys
import zipfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

RUN_MD = """# Boltz-2 inputs for {run}

{n} structures to fold. Every one needs **both** of these or the contact constraints do nothing:

* `force: true` on each constraint — already set in these files. Without it Boltz's featurizer
  discards the constraint outright.
* `--use_potentials` on the command line — the guidance that acts on the structure during sampling
  is gated behind it, and it is off by default. A constraint without it is conditioning only, and
  the model may decline.

## Fold them

```bash
for y in *.yaml; do
    boltz predict "$y" --num_workers 0 --out_dir results --use_potentials
done
```

On Colab with a GPU, `pip install boltz` then the same loop. Expect a couple of minutes per
structure on an A100 against two to three on a laptop GPU.

## Bring the results back

Zip or copy the `results` directory, then on the machine holding the run:

```bash
python code/boltz_import.py {run} --from /path/to/results
```

or equivalently `python code/boltz_offline.py import {run} --from /path/to/results`.

## Expected structures

{names}

## Note on affinity

These inputs do not request the affinity head. Boltz 2.2.1 computes the binder's molecular weight
with `AllChem.Descriptors.MolWt`, which RDKit 2026.03.6 no longer exposes, so asking for affinity
fails at parse time. Little is lost: that head spans 1.65 kcal/mol across 37 sequences and has given
confident scores to complexes with the ligand 6.6 A away. Interaction energies are computed here
instead, with UMA.
"""


def cmd_export(args):
    from boltz_hints import contacts, ligand_atom_names, residue_positions, write_yaml  # noqa
    from boltz_check import ligand_smiles

    design_path = args.design
    if not os.path.exists(design_path):
        design_path = os.path.join(args.outdir, args.design)
    design = json.load(open(design_path))

    out_dir = os.path.join(args.outdir, args.into)
    os.makedirs(out_dir, exist_ok=True)

    smiles = ligand_smiles(args.outdir)
    names = ligand_atom_names(args.outdir, smiles, affinity=False,
                              boltz_venv=args.boltz_venv) if not args.no_hints else []
    print(f"ligand SMILES: {smiles}")

    sequences = [(design["sequence"], "orig")]
    if args.variants:
        try:
            from fill_linkers import fill
            for i, item in enumerate(fill([design["sequence"]], variants=args.variants), 1):
                seq = item["filled"] if isinstance(item, dict) else item
                sequences.append((seq, f"esm{i}"))
                print(f"  esm{i}: {seq}")
        except Exception as exc:
            print(f"  ESM2 unavailable ({type(exc).__name__}: {exc}); glycine design only")

    all_hints = contacts(args.outdir, design, args.cutoff, names) if names else []
    written = []
    for seq, tag in sequences:
        # one unconstrained control per sequence, then the forced ladders
        write_yaml(os.path.join(out_dir, f"{tag}_control.yaml"), seq, smiles, [], args.margin,
                   False)
        written.append(f"{tag}_control")
        for n in args.forced:
            subset = sorted(all_hints, key=lambda h: h["distance"])[:n] if n < len(all_hints) \
                else list(all_hints)
            subset.sort(key=lambda h: h["residue"])
            write_yaml(os.path.join(out_dir, f"{tag}_f{n}.yaml"), seq, smiles, subset,
                       args.margin, True)
            written.append(f"{tag}_f{n}")

    listing = "\n".join(f"* `boltz_results_{n}/predictions/{n}/{n}_model_0.cif`" for n in written)
    with open(os.path.join(out_dir, "RUN.md"), "w") as fh:
        fh.write(RUN_MD.format(run=args.outdir, n=len(written), names=listing))
    with open(os.path.join(out_dir, "expected.json"), "w") as fh:
        json.dump({"structures": written, "smiles": smiles,
                   "sequence": design["sequence"]}, fh, indent=2)

    print(f"\n{len(written)} inputs -> {out_dir}")
    print(f"read {os.path.join(out_dir, 'RUN.md')} for the command, then bring `results` back with")
    print(f"  python code/boltz_offline.py import {args.outdir} --from /path/to/results")
    return 0


def cmd_import(args):
    src = os.path.expanduser(args.source)
    boltz_dir = os.path.join(args.outdir, "boltz")
    os.makedirs(boltz_dir, exist_ok=True)

    scratch = None
    if os.path.isfile(src) and src.endswith(".zip"):
        scratch = os.path.join(boltz_dir, "_unzip")
        os.makedirs(scratch, exist_ok=True)
        with zipfile.ZipFile(src) as z:
            z.extractall(scratch)
        src = scratch
        print(f"unpacked {args.source}")

    cifs = glob.glob(os.path.join(src, "**", "*.cif"), recursive=True)
    cifs = [c for c in cifs if "_model_" in os.path.basename(c)] or cifs
    if not cifs:
        sys.exit(f"no .cif files found under {src}")

    brought = []
    for cif in sorted(cifs):
        base = os.path.basename(cif)
        name = base.replace("_model_0.cif", "").replace(".cif", "")
        # rebuild the layout the rest of the pipeline globs for
        dest_dir = os.path.join(boltz_dir, f"boltz_results_{name}", "predictions", name)
        os.makedirs(dest_dir, exist_ok=True)
        dest = os.path.join(dest_dir, f"{name}_model_0.cif")
        shutil.copy(cif, dest)
        brought.append(name)

    for extra in glob.glob(os.path.join(src, "**", "*.csv"), recursive=True):
        shutil.copy(extra, os.path.join(boltz_dir, os.path.basename(extra)))

    expected_path = os.path.join(args.outdir, args.into, "expected.json")
    if os.path.exists(expected_path):
        expected = json.load(open(expected_path))["structures"]
        missing = [e for e in expected if e not in brought]
        print(f"\n{len(brought)} structures imported, {len(expected)} expected")
        if missing:
            print("missing: " + ", ".join(missing))
    else:
        print(f"\n{len(brought)} structures imported")
    for n in sorted(brought):
        print(f"  {n}")

    if scratch:
        shutil.rmtree(scratch, ignore_errors=True)
    print("\nnext:\n  python code/check_fold.py " + args.outdir +
          "\n  python code/binding_energy.py " + args.outdir +
          " --terms interaction,strain_ligand")
    return 0


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="cmd", required=True)

    p_ex = sub.add_parser("export", help="write every Boltz input into one portable folder")
    p_ex.add_argument("outdir")
    p_ex.add_argument("--design", default="design.json")
    p_ex.add_argument("--into", default="boltz_inputs")
    p_ex.add_argument("--variants", type=int, default=2, help="ESM2 linker variants (default 2)")
    p_ex.add_argument("--forced", type=int, nargs="+", default=[4, 8, 12],
                      help="how many contacts to force, one input each (default 4 8 12)")
    p_ex.add_argument("--cutoff", type=float, default=4.5)
    p_ex.add_argument("--margin", type=float, default=1.0)
    p_ex.add_argument("--no-hints", action="store_true", help="controls only, no constraints")
    p_ex.add_argument("--boltz-venv", default=None,
                      help="only needed for ligand atom naming, which uses Boltz's own code")
    p_ex.set_defaults(func=cmd_export)

    p_im = sub.add_parser("import", help="bring folded structures back into the run directory")
    p_im.add_argument("outdir")
    p_im.add_argument("--from", dest="source", required=True,
                      help="directory or .zip of Boltz output, in any layout")
    p_im.add_argument("--into", default="boltz_inputs",
                      help="where the export wrote expected.json, for checking completeness")
    p_im.set_defaults(func=cmd_import)

    args = parser.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
