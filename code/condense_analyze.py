"""Turn a condensation sweep into the reachability criterion step 3 needs.

`condense.py sweep` writes one row per (ordered pose pair, CA candidate pair, spacer count k)
with the closure error in angstroms. This reads that back and answers the question the plan
asks: given two placed fragments and k glycines between them, when can a backbone actually make
the connection?

The answer comes out as a lookup table rather than a fitted formula, because the two variables
that matter are not interchangeable. Distance is the obvious one. The other is how far *behind*
its own fragment each partner sits -- the angle at CA between the partner and that residue's own
CB -- because the backbone leaves CA pointing away from CB, so a partner at a small angle is one
the chain has to double back to reach. A pair 11 A apart can close easily while a pair 11 A apart
at a small angle cannot close at all.

Usage:
    python code/condense_analyze.py runs/octinoxate/condense_pairs_geom.csv
    python code/condense_analyze.py runs/octinoxate/condense_pairs_geom.csv --tol 0.2
"""
import argparse
import csv
from collections import defaultdict

import numpy as np

THETA_BINS = [(0, 30), (30, 60), (60, 90), (90, 180)]
D_BINS = [(0, 4), (4, 6), (6, 8), (8, 10), (10, 12), (12, 14), (14, 30)]


def load(path):
    with open(path) as fh:
        rows = list(csv.DictReader(fh))
    for r in rows:
        for key in ("d_ca", "d_cb", "closure", "cos_att", "theta_a", "theta_b",
                    "clash_frag", "clash_lig", "penalty"):
            if key in r and r[key] != "":
                r[key] = float(r[key])
        r["k"] = int(r["k"])
    return rows


def bin_of(value, bins):
    for lo, hi in bins:
        if lo <= value < hi:
            return (lo, hi)
    return bins[-1]


def summarise(rows, tol, clash_min=None):
    """Per k: how many combinations close, and how far they reach."""
    print(f"closes = closure < {tol} A"
          + (f" and no fragment/ligand contact under {clash_min} A" if clash_min else ""))
    print(f"\n{'k':>2} {'rows':>6} {'closes':>7} {'rate':>6} {'median closure':>15} "
          f"{'max d_CA that closes':>22}")
    for k in sorted({r["k"] for r in rows}):
        sub = [r for r in rows if r["k"] == k]
        good = [r for r in sub if closes(r, tol, clash_min)]
        med = np.median([r["closure"] for r in sub])
        reach = max((r["d_ca"] for r in good), default=float("nan"))
        print(f"{k:>2} {len(sub):>6} {len(good):>7} {len(good) / len(sub):>6.2f} "
              f"{med:>15.3f} {reach:>22.2f}")


def closes(row, tol, clash_min=None):
    if row["closure"] >= tol:
        return False
    if clash_min is not None:
        if min(row.get("clash_frag", 9e9), row.get("clash_lig", 9e9)) < clash_min:
            return False
    return True


def reach_table(rows, tol, clash_min=None):
    """The criterion: largest CA-CA distance that closes, by spacer count and approach angle."""
    print("\n\nlargest CA-CA distance (A) that closes, by spacer count and the smaller of the two"
          "\napproach angles (the angle at CA between the partner and that residue's own CB):\n")
    header = "  ".join(f"{lo:>3}-{hi:<3}" for lo, hi in THETA_BINS)
    print(f"{'k':>2}  {header}")
    for k in sorted({r["k"] for r in rows}):
        cells = []
        for lo, hi in THETA_BINS:
            good = [r["d_ca"] for r in rows
                    if r["k"] == k and lo <= min(r["theta_a"], r["theta_b"]) < hi
                    and closes(r, tol, clash_min)]
            cells.append(f"{max(good):>7.1f}" if good else f"{'-':>7}")
        print(f"{k:>2}  " + "  ".join(cells))


def rate_grid(rows, tol, clash_min=None):
    """Closing rate over distance and angle, pooled across k: the shape of the constraint."""
    print("\n\nfraction of candidate attachments that close, pooled over k "
          "(blank = no rows in that cell):\n")
    header = "  ".join(f"{lo:>3}-{hi:<3}" for lo, hi in THETA_BINS)
    print(f"{'d_CA':>9}  {header}")
    for dlo, dhi in D_BINS:
        cells = []
        for tlo, thi in THETA_BINS:
            sub = [r for r in rows
                   if dlo <= r["d_ca"] < dhi and tlo <= min(r["theta_a"], r["theta_b"]) < thi]
            if not sub:
                cells.append(f"{'':>7}")
            else:
                rate = sum(closes(r, tol, clash_min) for r in sub) / len(sub)
                cells.append(f"{rate:>7.2f}")
        print(f"{dlo:>4}-{dhi:<4}  " + "  ".join(cells))


def per_pair(rows, tol, clash_min=None):
    """The smallest k each ordered pair needs, and the pairs no k reaches."""
    best = defaultdict(lambda: None)
    for r in rows:
        if not closes(r, tol, clash_min):
            continue
        key = (r["frag_a"], r["pose_a"], r["frag_b"], r["pose_b"])
        if best[key] is None or r["k"] < best[key]["k"]:
            best[key] = r
    all_pairs = {(r["frag_a"], r["pose_a"], r["frag_b"], r["pose_b"]) for r in rows}
    unreachable = sorted(all_pairs - {k for k, v in best.items() if v})
    print(f"\n\n{len(all_pairs) - len(unreachable)} of {len(all_pairs)} ordered pairs can be "
          f"connected by some k; {len(unreachable)} cannot.")
    counts = defaultdict(int)
    for v in best.values():
        if v:
            counts[v["k"]] += 1
    print("smallest k that works: " + ", ".join(f"k={k}: {counts[k]}" for k in sorted(counts)))
    if unreachable:
        print("\nnot connectable at any k tried:")
        for fa, pa, fb, pb in unreachable:
            sub = [r for r in rows if (r["frag_a"], r["pose_a"], r["frag_b"], r["pose_b"])
                   == (fa, pa, fb, pb)]
            near = min(sub, key=lambda r: r["closure"])
            print(f"  {fa}{pa} -> {fb}{pb}   d_CA {near['d_ca']:5.2f}  "
                  f"theta {near['theta_a']:5.1f}/{near['theta_b']:5.1f}  "
                  f"best closure {near['closure']:.2f} A at k={near['k']}")


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("csv_path")
    parser.add_argument("--tol", type=float, default=0.2,
                        help="closure error a real backbone can absorb (default 0.2 A)")
    parser.add_argument("--clash-min", type=float, default=None,
                        help="also require no fragment/ligand contact closer than this")
    args = parser.parse_args(argv)

    rows = load(args.csv_path)
    print(f"{len(rows)} rows from {args.csv_path}\n")
    summarise(rows, args.tol, args.clash_min)
    reach_table(rows, args.tol, args.clash_min)
    rate_grid(rows, args.tol, args.clash_min)
    per_pair(rows, args.tol, args.clash_min)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
