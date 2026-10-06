"""A null control for the whole pipeline: the same residues, arranged at random.

Every conclusion in this repository rests on designs whose side-chain *arrangement* was chosen --
greedy pose selection, then an exact Held-Karp ordering, then spacer counts from a reachability
sweep. This asks what that arrangement is worth, by holding everything else fixed:

* the same 34 residues, so identical composition and identical glycine fraction
* the same **spacer pattern**, so the same linker lengths in the same places
* only which side chain occupies which slot is randomised

That last restriction is the point. A naive shuffle of the whole string also randomises the spacing
and produces adjacent side-chain pairs with no glycine between them -- 4 to 5 of them in practice,
where the design has none -- so a naive shuffle changes two variables and cannot say which one
mattered. Keeping the spacer pattern leaves exactly one variable.

There is deliberately **no scoring**. The control has no placed poses, so it has no fragment
interaction energy, and inventing one to get it through `sequences.csv` would put a fake number in
the file that records designed shells. It is folded directly and sent to dynamics, where the
comparison is MM/GBSA against the designed folds.

Usage:
    python code/shuffle_control.py runs/octinoxate --design design_shell3.json --seed 2
    python code/shuffle_control.py runs/octinoxate --design design_shell3.json --seed 2 --fold
    python code/shuffle_control.py runs/oxybenzone --design design_shell2.json --seed 2 \
        --name ox2_shuffle --fold      # one control per shell needs one --name per shell
"""
import argparse
import csv
import json
import os
import random
import sys


def shuffled(sequence, seed, linker="G"):
    """Randomise which side chain sits in which non-linker slot, keeping the linker pattern."""
    slots = [i for i, c in enumerate(sequence) if c != linker]
    side = [sequence[i] for i in slots]
    rng = random.Random(seed)
    out = side[:]
    rng.shuffle(out)
    s = list(linker * len(sequence))
    for i, c in zip(slots, out):
        s[i] = c
    unmoved = sum(1 for a, b in zip(side, out) if a == b)
    return "".join(s), "".join(side), "".join(out), unmoved, slots


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0],
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("outdir", help="a run directory, e.g. runs/octinoxate")
    p.add_argument("--design", default="design.json",
                   help="the design json to shuffle, relative to outdir (default design.json)")
    p.add_argument("--seed", type=int, default=2, help="shuffle seed (default 2)")
    p.add_argument("--name", default="shuffle_control",
                   help="name for this control, used for its records and its CSV (default "
                        "shuffle_control). Give each shell's control its own name when a run has "
                        "more than one, e.g. ox1_shuffle, or the second overwrites the first")
    p.add_argument("--variants", type=int, default=1, help="ESM2 linker-filled variants (default 1)")
    p.add_argument("--fold", action="store_true", help="co-fold the control and its variant in Boltz")
    p.add_argument("--boltz-repo", default=os.environ.get("PEPTIDEBUILDER_BOLTZ_REPO"))
    p.add_argument("--boltz-venv", default=os.environ.get("PEPTIDEBUILDER_BOLTZ_VENV"))
    p.add_argument("--boltz-cmd", default=os.environ.get("PEPTIDEBUILDER_BOLTZ_CMD"))
    args = p.parse_args(argv)

    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    design = json.load(open(os.path.join(args.outdir, args.design)))["sequence"]
    ctrl, side_in, side_out, unmoved, slots = shuffled(design, args.seed)

    print(f"design   {design}")
    print(f"control  {ctrl}")
    print(f"  {len(design)} residues both, {design.count('G')} glycines both")
    print(f"  side chains {side_in} -> {side_out}")
    print(f"  {unmoved}/{len(side_in)} side chains left in their designed slot")
    if ctrl == design:
        sys.exit("the shuffle reproduced the design; choose another --seed")

    records = [{"name": args.name, "sequence": ctrl, "source": f"{args.design} seed {args.seed}"}]

    if args.variants:
        from fill_linkers import fill
        # Same ESM2 route the designed variants used, so the control's variant is not a different
        # kind of object from the ones it is being compared against. One record per variant, each
        # carrying the filled sequence as a string.
        for r in fill([ctrl], variants=args.variants, seed=args.seed):
            filled, j = r["filled"], r["variant"]
            print(f"  esm variant {j}: {filled}  ({len(r['steps'])} linkers filled)")
            records.append({"name": f"{args.name}_esm{j}", "sequence": filled,
                            "source": f"esm2 fill of {args.name}, seed {args.seed}"})

    out_csv = os.path.join(args.outdir, f"{args.name}.csv")
    with open(out_csv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["name", "sequence", "source"])
        w.writeheader()
        w.writerows(records)
    print(f"\nwrote {out_csv} ({len(records)} sequences); sequences.csv untouched")

    if args.fold:
        from boltz_check import ligand_smiles, run_boltz
        smiles = ligand_smiles(args.outdir)
        boltz_dir = os.path.join(args.outdir, "boltz")
        os.makedirs(boltz_dir, exist_ok=True)
        print(f"\nligand SMILES: {smiles}")
        for rec in records:
            print(f"\nfolding {rec['name']} ({len(rec['sequence'])} residues), no hints", flush=True)
            out = run_boltz(rec["name"], rec["sequence"], smiles, False, boltz_dir,
                            args.boltz_repo, boltz_cmd=args.boltz_cmd, boltz_venv=args.boltz_venv)
            if out is None:
                print(f"  {rec['name']}: Boltz produced nothing")
                continue
            # No affinity: the structure is the deliverable, and check_fold.py is what reads it.
            print(f"  folded -> {os.path.relpath(out['cif'], args.outdir)}")
    else:
        print("\nre-run with --fold to co-fold these, or fold them however you prefer")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
