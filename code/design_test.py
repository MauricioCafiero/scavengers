"""End-to-end test of a condensation design: order it, fill it, fold it, score it.

  1. `assign.py`      choose the order through the placed side chains, and how many glycine
                      spacers each step needs, from the measured reachability sweep
  2. `condense_chain` **optional, `--chain-gate`.** Build the whole path as one chain and measure
                      the worst junction, as a gate: a path that cannot be built is not worth
                      folding. Off by default because it is expensive -- tens of minutes for a long
                      path -- and the folds themselves are cheap. Worth turning on when folds are
                      the scarce resource, or to check that a path's spacer counts really are
                      sufficient rather than only pairwise sufficient
  3. `fill_linkers`   ESM2 variants of those spacers, so the linkers are not all glycine
  4. `boltz_hints`    write a Boltz input per sequence, with the designed shell as contact hints,
                      plus one unhinted control. Hints need `force: true` *and* `--use_potentials`,
                      or they do nothing at all
  5. `boltz_hints`    co-fold each of them
  6. `check_fold`     **how much of the ligand each fold wraps**, plus whether it is bound at all
                      and whether the hints were honoured
  7. `binding_energy` interaction energy and ligand strain (cavity desolvation optional)
  8. `overlay`        how many of the designed side-chain positions the fold reproduces

**Step 6 is the one that decides things.** Wrapping is the design objective; hint satisfaction and
interaction energy are proxies for it and have both inverted the verdict -- forced contacts satisfied
almost no hints while raising ligand engagement from 75% to 100%, and the ESM2 variants matched the
original on summed energy while engaging 60% and 40% against its 75%.

The hinted fold asks whether the shell is reachable by a backbone at all; the unhinted control asks
whether the *sequence* encodes it. The first sets the ceiling for the second, and the old pipeline
only ever tried the second.

The design's own shell is written into `sequences/` so `overlay.py` has something to compare the
fold against -- the ligand plus the placed fragments, i.e. the arrangement being asked for.

Usage:
    python code/design_test.py runs/octinoxate --pairs-csv condense_shell12.csv \
        --poses arginine:3,lysine:5,... [--variants 2] [--skip-fold]
"""
import argparse
import json
import os
import subprocess
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ase.io  # noqa: E402

from condense import parse_poses  # noqa: E402
from peptide_builder import combine_poses, define_fragments, load_ligand  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
PY = sys.executable


def step(title):
    print(f"\n{'=' * 78}\n{title}\n{'=' * 78}", flush=True)


def run(cmd):
    print("$ " + " ".join(cmd), flush=True)
    return subprocess.call(cmd)


def write_design_shell(outdir, design, ligand_name):
    """The ligand plus this design's placed fragments, named after the sequence.

    `overlay.py` compares a folded complex against one of these, so a design that was not produced
    by `build_sequence` needs its shell written out before it can be scored.
    """
    frags = define_fragments()
    ligand = load_ligand(ligand_name)
    by = {f["name"]: i for i, f in enumerate(frags)}
    n_poses = {}
    for name, pose in design["path"]:
        n_poses[by[name]] = max(n_poses.get(by[name], 0), pose + 1)
    mols = [[] for _ in frags]
    for idx, count in n_poses.items():
        for j in range(count):
            path = os.path.join(outdir, "poses",
                                f"{ligand['name']}_w_{frags[idx]['name']}{j}.xyz")
            mols[idx].append(ase.io.read(path).get_positions())
    selected = [(by[name], pose) for name, pose in design["path"]]
    seq_dir = os.path.join(outdir, "sequences")
    os.makedirs(seq_dir, exist_ok=True)
    out = os.path.join(seq_dir, design["sequence"] + ".xyz")
    combine_poses(ligand, frags, mols, selected, out)
    print(f"wrote the design's own shell to {out}")
    return out


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("outdir")
    parser.add_argument("--poses", required=True, help="the shell, frag:pose,frag:pose,...")
    parser.add_argument("--pairs-csv", required=True, help="reachability sweep for that shell")
    parser.add_argument("--ligand", default="octinoxate")
    parser.add_argument("--tol", type=float, default=0.2)
    parser.add_argument("--variants", type=int, default=2, help="ESM2 linker variants (default 2)")
    parser.add_argument("--cutoff", type=float, default=4.5, help="contact hint cutoff, A")
    parser.add_argument("--margin", type=float, default=1.0, help="slack added to each hint, A")
    parser.add_argument("--chain-gate", action="store_true",
                        help="build the whole path as one chain first and refuse to fold it if the "
                             "worst junction exceeds --max-junction (expensive; off by default)")
    parser.add_argument("--max-junction", type=float, default=1.0,
                        help="largest junction miss, in A, still worth folding. A real backbone "
                             "flexes a few degrees, so a near miss is absorbable; several angstroms "
                             "is not (default 1.0)")
    parser.add_argument("--chain-restarts", type=int, default=40)
    parser.add_argument("--max-combos", type=int, default=12)
    parser.add_argument("--replicates-water", type=int, default=2,
                        help="water placements per cavity desolvation (default 2, to show noise)")
    parser.add_argument("--skip-fold", action="store_true", help="prepare everything but do not fold")
    args = parser.parse_args(argv)

    design_path = os.path.join(args.outdir, "design.json")

    step("1. order the shell and count the spacers each step needs")
    if run([PY, os.path.join(HERE, "assign.py"), args.outdir, "--poses", args.poses,
            "--pairs-csv", args.pairs_csv, "--tol", str(args.tol), "--save", design_path]):
        sys.exit("assign.py failed")
    design = json.load(open(design_path))
    print(f"\ndesign: {design['sequence']}  ({design['residues']} residues, "
          f"{len(design['path'])} side chains, total fragment IE "
          f"{design['total_ie_kcal_mol']} kcal/mol)")
    print(f"spacers per step: {design['linkers']}  (sum {design['sum_k']}, "
          f"budget {design['budget']})")
    if design["sum_k"] < design["budget"]:
        print("WARNING: below the linker budget, so this path is not expected to close as spelled")

    write_design_shell(args.outdir, design, args.ligand)

    if args.chain_gate:
        step("2. condense the whole path as one chain -- the gate before spending folds "
             "(skipped unless --chain-gate)")
        chain_json = os.path.join(args.outdir, "chain.json")
        cmd = [PY, os.path.join(HERE, "condense_chain.py"), args.outdir,
               "--path"] + [f"{name}:{pose}" for name, pose in design["path"]] + \
              ["--k"] + [str(k) for k in design["linkers"]] + \
              ["--steric", "--restarts", str(args.chain_restarts),
               "--max-combos", str(args.max_combos), "--tol", str(args.max_junction),
               "--save", chain_json]
        if run(cmd):
            sys.exit("condense_chain.py failed")
        chain = json.load(open(chain_json))
        worst = chain["worst_junction"]
        print(f"\nworst junction {worst:.3f} A over {len(chain['closures'])} junctions "
              f"(threshold {args.max_junction} A)")
        print("  per junction: " + "  ".join(f"{c:.2f}" for c in chain["closures"]))
        if worst > args.max_junction:
            print(f"\nSTOPPING: the path cannot be built as one chain within {args.max_junction} A, so "
                  f"folding it\nwould be testing a sequence whose designed geometry is not realisable. "
                  f"Raise --max-junction\nto override, or give the ordering search more linkers.")
            return 1
        print("  within tolerance -- a real backbone can absorb this, so the fold is worth running")

    step("3. ESM2 variants of the glycine spacers")
    sequences = [(design["sequence"], "hints")]
    if args.variants:
        from fill_linkers import fill
        try:
            filled = fill([design["sequence"]], variants=args.variants)
            for i, item in enumerate(filled, 1):
                seq = item["filled"] if isinstance(item, dict) else item
                print(f"  variant {i}: {seq}")
                sequences.append((seq, f"esm{i}"))
        except Exception as exc:
            print(f"  ESM2 fill unavailable ({type(exc).__name__}: {exc}); "
                  f"continuing with the glycine design only")

    step("4. write the Boltz inputs (hinted, plus one unhinted control)")
    yamls = []
    for seq, tag in sequences:
        cmd = [PY, os.path.join(HERE, "boltz_hints.py"), args.outdir, "design.json",
               "--cutoff", str(args.cutoff), "--margin", str(args.margin),
               "--name", f"{design['sequence']}_{tag}"]
        if seq != design["sequence"]:
            cmd += ["--sequence", seq]
        if run(cmd):
            sys.exit("boltz_hints.py failed")
        yamls.append(f"{design['sequence']}_{tag}")
    if run([PY, os.path.join(HERE, "boltz_hints.py"), args.outdir, "design.json", "--no-hints",
            "--name", f"{design['sequence']}_control"]):
        sys.exit("boltz_hints.py failed")
    yamls.append(f"{design['sequence']}_control")

    if args.skip_fold:
        print("\n--skip-fold: stopping before the folds")
        return 0

    step("5. co-fold each one")
    from boltz_hints import fold
    boltz_dir = os.path.join(args.outdir, "boltz")
    for name in yamls:
        fold(os.path.join(boltz_dir, name + ".yaml"), boltz_dir)

    step("6. did each fold bind the ligand, and did it honour its hints?")
    run([PY, os.path.join(HERE, "check_fold.py"), args.outdir, "--quiet"])

    step("7. binding energy on every folded complex: interaction, strain, cavity desolvation")
    run([PY, os.path.join(HERE, "binding_energy.py"), args.outdir,
         "--terms", "interaction,strain,cavity", "--replicates", str(args.replicates_water)])

    step("8. how much of the designed shell each fold reproduces")
    run([PY, os.path.join(HERE, "overlay.py"), args.outdir])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
