"""Correlate the measured properties of folded complexes against each other.

Every number this pipeline produces about a fold lands in one of three CSVs -- `fold_check*.csv` for
geometry, `binding_*.csv` for energies, `overlay*.csv` for how much of the design was reproduced --
and the questions worth asking span them: does enclosure predict interaction energy, does interaction
energy just track peptide size, does forcing contacts move side chains toward their designed
positions. This joins all three on structure name and reports the correlations.

Two features exist because of specific mistakes made on this project.

**`--group` is first class.** The correlations here invert between subsets. Enclosure ordered the
interaction energy perfectly within shell 1's glycine ladder and nearly reversed it for that shell's
`esm2` variant; across shell 2's `esm2` folds the interaction energy correlates +0.93 with ligand
centroid separation and +0.18 with enclosure, because a +5 peptide binds by long-range electrostatics
rather than by wrapping. A correlation computed over a pooled set of folds that bind by different
mechanisms is close to meaningless, so `--group` splits by a regex capture and reports each subset
separately alongside the pooled figure.

**n is printed next to every coefficient, and small n is flagged.** Two patterns announced on this
project were killed by the next measurement. With four structures per ladder, a rank correlation of
1.0 happens by chance about 4% of the time, so the flag is not decorative.

Both Pearson (linear) and Spearman (rank) are given. Spearman is the honest one for small n and for
monotonic-but-curved relationships, which is most of them here; Pearson is reported because the
README quotes it in places.

Usage:
    python code/correlate.py runs/octinoxate
    python code/correlate.py runs/octinoxate --match s2_
    python code/correlate.py runs/octinoxate --match s2_ --group 's2_(esm1|esm2|orig)_'
    python code/correlate.py runs/octinoxate --pairs enclosed_fraction:interaction_kcal
    python code/correlate.py runs/octinoxate --list
"""
import argparse
import csv
import glob
import math
import os
import re
import sys

# Properties worth correlating, in the order they get reported. Anything numeric in the joined CSVs
# is available; this list is what `--all` skips having to spell out.
DEFAULT_PROPS = [
    "enclosed_fraction", "wrapped_fraction", "centroid_separation", "peptide_rg",
    "contacts_under_cutoff", "mean_ligand_distance", "engaged",
    "interaction_kcal", "strain_ligand_kcal", "sum_kcal", "interaction_per_atom",
    "n_peptide_atoms", "peptide_charge", "abs_charge",
    "hints_satisfied", "n_hints",
    "matched_within_cutoff", "median_distance_A", "median_distance_any_type_A", "ligand_rmsd_A",
]


def load(outdir):
    """Join every results CSV on structure name, and add the derived properties."""
    boltz = os.path.join(outdir, "boltz")
    rows = {}
    patterns = ["fold_check*.csv", "binding_*.csv", "overlay*.csv"]
    files = [p for pat in patterns for p in sorted(glob.glob(os.path.join(boltz, pat)))]
    if not files:
        sys.exit(f"no results CSVs under {boltz} (looked for {', '.join(patterns)})")
    for path in files:
        with open(path) as fh:
            for row in csv.DictReader(fh):
                name = row.get("name")
                if not name:
                    continue
                rows.setdefault(name, {}).update(
                    {k: v for k, v in row.items() if v not in ("", None)})

    # `binding_energy.py` writes its csv only when the whole batch finishes, but drops a
    # partial_<name>.json as each structure completes. Reading those too makes this usable while a
    # scoring run is still going -- which is when the numbers are most worth looking at, and the
    # reason the partials exist at all. They lose to the csv where both have a value, since the csv
    # is the finished record.
    import json
    for path in sorted(glob.glob(os.path.join(boltz, "partial_*.json"))):
        try:
            d = json.load(open(path))
        except (ValueError, OSError):
            continue
        name = d.get("name") or os.path.basename(path)[len("partial_"):-len(".json")]
        have = rows.setdefault(name, {})
        for k, v in d.items():
            if v not in ("", None) and k not in have:
                have[k] = v
        files.append(path)

    for name, r in rows.items():
        def num(key):
            try:
                return float(r[key])
            except (KeyError, TypeError, ValueError):
                return None
        i, s = num("interaction_kcal"), num("strain_ligand_kcal")
        if i is not None and s is not None:
            r["sum_kcal"] = i + s
        # The interaction energy grows with system size (Spearman -0.52 against peptide length), so
        # a per-atom figure is the one to rank by. Kept alongside, not instead of, the raw value.
        n_at = num("n_peptide_atoms")
        if i is not None and n_at:
            r["interaction_per_atom"] = i / n_at
        q = num("peptide_charge")
        if q is not None:
            r["abs_charge"] = abs(q)
    return rows, files


def numeric(rows, names, prop):
    """The (name, value) pairs where `prop` parses as a number."""
    out = []
    for n in names:
        v = rows[n].get(prop)
        if v is None:
            continue
        try:
            out.append((n, float(v)))
        except (TypeError, ValueError):
            continue
    return out


def ranks(xs):
    """Average ranks, so ties do not bias the rank correlation."""
    order = sorted(range(len(xs)), key=lambda i: xs[i])
    out = [0.0] * len(xs)
    i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and xs[order[j + 1]] == xs[order[i]]:
            j += 1
        avg = (i + j) / 2 + 1
        for k in range(i, j + 1):
            out[order[k]] = avg
        i = j + 1
    return out


def pearson(a, b):
    n = len(a)
    if n < 3:
        return None
    ma, mb = sum(a) / n, sum(b) / n
    va = sum((x - ma) ** 2 for x in a)
    vb = sum((y - mb) ** 2 for y in b)
    if va == 0 or vb == 0:
        return None
    return sum((x - ma) * (y - mb) for x, y in zip(a, b)) / math.sqrt(va * vb)


def correlate(rows, names, x, y):
    """Returns (n, pearson, spearman) over the structures that have both properties."""
    ax = dict(numeric(rows, names, x))
    ay = dict(numeric(rows, names, y))
    shared = sorted(set(ax) & set(ay))
    if len(shared) < 3:
        return len(shared), None, None
    xs = [ax[n] for n in shared]
    ys = [ay[n] for n in shared]
    return len(shared), pearson(xs, ys), pearson(ranks(xs), ranks(ys))


def flag(n, rho):
    """A word of warning, not a p-value: this is a diagnostic tool, not a statistics package."""
    if rho is None:
        return ""
    if n < 5:
        return "  (n<5, treat as anecdote)"
    if n < 8 and abs(rho) > 0.9:
        return "  (n<8 and near-perfect: check for a confound)"
    return ""


def report(rows, names, props, label):
    print(f"\n{'=' * 78}\n{label}  (n={len(names)} structures)\n{'=' * 78}")
    present = [p for p in props if len(numeric(rows, names, p)) >= 3]
    missing = [p for p in props if p not in present]
    if missing:
        print(f"not enough data for: {', '.join(missing)}\n")
    for i, x in enumerate(present):
        for y in present[i + 1:]:
            n, r, rho = correlate(rows, names, x, y)
            if rho is None:
                continue
            print(f"  {x:<26} {y:<26} n={n:<3} pearson {r:>+6.2f}  spearman {rho:>+6.2f}"
                  f"{flag(n, rho)}")


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0],
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("outdir", help="a run directory, e.g. runs/octinoxate")
    parser.add_argument("--match", help="only structures whose name contains this, e.g. s2_")
    parser.add_argument("--group", help="regex with one capture group; correlations are reported per "
                                       "captured value as well as pooled")
    parser.add_argument("--pairs", help="only these pairs, as x:y,x:y. Default is every pair of the "
                                        "standard properties")
    parser.add_argument("--props", help="comma-separated properties to use instead of the defaults")
    parser.add_argument("--list", action="store_true", help="list the joined properties and exit")
    args = parser.parse_args(argv)

    rows, files = load(args.outdir)
    print(f"joined {len(rows)} structures from {len(files)} csv(s): "
          f"{', '.join(os.path.basename(f) for f in files)}")

    if args.list:
        keys = sorted({k for r in rows.values() for k in r})
        print(f"\n{len(keys)} properties:")
        for k in keys:
            print(f"  {k:<30} numeric for {len(numeric(rows, list(rows), k))} structures")
        return 0

    names = sorted(rows)
    if args.match:
        names = [n for n in names if args.match in n]
        if not names:
            sys.exit(f"no structures match {args.match!r}")

    if args.pairs:
        pairs = [tuple(p.split(":", 1)) for p in args.pairs.split(",")]
        print(f"\n{'=' * 78}")
        for x, y in pairs:
            n, r, rho = correlate(rows, names, x, y)
            if rho is None:
                print(f"  {x} vs {y}: only {n} structures have both, need 3")
                continue
            print(f"  {x:<26} {y:<26} n={n:<3} pearson {r:>+6.2f}  spearman {rho:>+6.2f}"
                  f"{flag(n, rho)}")
        return 0

    props = args.props.split(",") if args.props else DEFAULT_PROPS
    report(rows, names, props, "POOLED" + (f", match {args.match!r}" if args.match else ""))

    if args.group:
        rx = re.compile(args.group)
        groups = {}
        for n in names:
            m = rx.search(n)
            if m:
                groups.setdefault(m.group(1), []).append(n)
        for key in sorted(groups):
            if len(groups[key]) >= 3:
                report(rows, groups[key], props, f"GROUP {key!r}")
            else:
                print(f"\nGROUP {key!r}: only {len(groups[key])} structures, skipped")
        print("\nWhere a pooled correlation and a per-group one disagree, the per-group one is the "
              "\nreal relationship and the pooled figure is mixing binding mechanisms.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
