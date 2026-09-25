"""Test the design the pipeline actually produced against the measured reachability criterion.

`build_sequence` chains the placed fragments by nearest centroid and converts each gap into
`round(d / 3.8) - 1` glycines. Nothing in that decides whether a backbone can run from one
fragment's attachment point to the next: the distance it uses is between centroids, and the
direction the attachment points face is not consulted at all.

This walks the design's own ordering, and for each consecutive step asks the question the sweep
answers: at the number of linkers the design assigned, can an ideal backbone make the connection?
And if not, what would it take? The answer is a per-step verdict for the design that produced the
3.6% side-chain reproduction rate.

Usage:
    python code/check_design.py runs/octinoxate --pairs-csv condense_pairs_geom.csv
"""
import argparse
import csv
import os
import sys
from collections import defaultdict

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from condense import selected_attachments  # noqa: E402
from peptide_builder import (ONE_LETTER, build_sequence, define_fragments,  # noqa: E402
                             load_ligand)


def load_sweep(path):
    """{(frag_a, pose_a, frag_b, pose_b, k): best row} over the CA candidate combinations."""
    best = {}
    with open(path) as fh:
        for row in csv.DictReader(fh):
            key = (row["frag_a"], int(row["pose_a"]), row["frag_b"], int(row["pose_b"]),
                   int(row["k"]))
            row["closure"] = float(row["closure"])
            if key not in best or row["closure"] < best[key]["closure"]:
                best[key] = row
    return best


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("outdir")
    parser.add_argument("--ligand", default="octinoxate")
    parser.add_argument("--pairs-csv", default="condense_pairs_geom.csv")
    parser.add_argument("--tol", type=float, default=0.2)
    parser.add_argument("--linker-slack", type=int, default=0)
    args = parser.parse_args(argv)

    frags = define_fragments()
    ligand = load_ligand(args.ligand)
    atts = selected_attachments(args.outdir, frags, ligand)
    centers = [np.mean(a["frag_coords"], axis=0) for a in atts]
    names = [a["name"] for a in atts]
    poses = {a["name"]: a["pose"] for a in atts}

    sequence, steps = build_sequence(centers, names, linker_slack=args.linker_slack, verbose=False)
    print(f"design sequence: {sequence}   ({len(sequence)} residues, "
          f"{sum(1 for c in sequence if c == 'G')} glycine)\n")

    path = args.pairs_csv
    if not os.path.isabs(path) and not os.path.exists(path):
        path = os.path.join(args.outdir, path)
    sweep = load_sweep(path)
    ks_available = sorted({key[4] for key in sweep})

    print(f"{'step':<28} {'centroid d':>10} {'linkers':>8} {'d_CA':>6} {'theta_a':>8} "
          f"{'theta_b':>8} {'closure':>8}  verdict")
    verdicts = defaultdict(int)
    for name_from, name_to, distance, n_linkers in steps:
        key = (name_from, poses[name_from], name_to, poses[name_to], n_linkers)
        row = sweep.get(key)
        label = f"{ONE_LETTER[name_from]}{poses[name_from]} -> {ONE_LETTER[name_to]}{poses[name_to]}"
        if row is None:
            print(f"{label:<28} {distance:>10.2f} {n_linkers:>8} "
                  f"{'':>6} {'':>8} {'':>8} {'':>8}  not swept at k={n_linkers}")
            verdicts["not swept"] += 1
            continue
        ok = row["closure"] < args.tol
        # the smallest number of linkers that would have worked
        works = [k for k in ks_available
                 if (name_from, poses[name_from], name_to, poses[name_to], k) in sweep
                 and sweep[(name_from, poses[name_from], name_to, poses[name_to],
                            k)]["closure"] < args.tol]
        if ok:
            verdict = "connects"
            verdicts["connects"] += 1
        elif works:
            verdict = f"NO - needs k={min(works)}, design gave {n_linkers}"
            verdicts["wrong linker count"] += 1
        else:
            verdict = f"NO - no k up to {max(ks_available)} connects this step"
            verdicts["unreachable"] += 1
        print(f"{label:<28} {distance:>10.2f} {n_linkers:>8} {float(row['d_ca']):>6.2f} "
              f"{float(row['theta_a']):>8.1f} {float(row['theta_b']):>8.1f} "
              f"{row['closure']:>8.3f}  {verdict}")

    total = sum(verdicts.values())
    print(f"\n{verdicts['connects']} of {total} steps in the design are connectable as drawn.")
    for label in ("wrong linker count", "unreachable", "not swept"):
        if verdicts[label]:
            print(f"  {verdicts[label]} {label}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
