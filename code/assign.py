"""Reachability-constrained assignment: choose which poses to use, and in what order (method B).

`build_sequence` chains fragment centroids by nearest neighbour and converts each gap into
`round(d / 3.8) - 1` glycines. It has no notion of whether a backbone can make the step, so it
routinely prescribes connections no peptide can follow -- see `check_design.py`.

This picks the path instead. Nodes are placed poses, each worth its interaction energy. An edge
a -> b exists only if the sweep found a backbone that connects them, and it costs the number of
residues that connection takes: k glycine linkers plus residue b itself. The result is an ordering
that is known to be traversable before a single second is spent folding it, which is the whole
point of doing the geometry first.

The search is exact. With ten poses the subset-and-order space is 2^10 * 10, so a Held-Karp
dynamic program over (visited set, current node) finds the true optimum rather than a greedy walk.

Usage:
    python code/assign.py runs/octinoxate --pairs-csv condense_pairs_geom.csv
    python code/assign.py runs/octinoxate --lam 0.5      # charge 0.5 kcal/mol per residue
"""
import argparse
import csv
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from condense import parse_poses, selected_attachments  # noqa: E402
from peptide_builder import ONE_LETTER, define_fragments, load_ligand  # noqa: E402


def load_edges(path, tol, max_penalty=None):
    """Cheapest connection for each ordered pair: {(a, b): (k, closure, penalty)}.

    `a` and `b` are (fragment name, pose). The cheapest connection is the one with the fewest
    linkers, since every linker is a residue the peptide has to pay for.
    """
    edges = {}
    with open(path) as fh:
        for row in csv.DictReader(fh):
            closure = float(row["closure"])
            if closure >= tol:
                continue
            penalty = float(row.get("penalty") or 0.0)
            if max_penalty is not None and penalty > max_penalty:
                continue
            a = (row["frag_a"], int(row["pose_a"]))
            b = (row["frag_b"], int(row["pose_b"]))
            k = int(row["k"])
            if (a, b) not in edges or k < edges[(a, b)][0]:
                edges[(a, b)] = (k, closure, penalty)
    return edges


def best_paths(nodes, values, edges, lam):
    """Held-Karp over (visited subset, last node), minimising sum(IE) + lam * residues.

    Returns {number of poses used: (score, total_ie, residues, path)} -- the best path for each
    number of poses, so the trade-off between covering the shell and paying for chain length is
    visible rather than hidden in one number.
    """
    n = len(nodes)
    index = {node: i for i, node in enumerate(nodes)}
    adjacency = [[] for _ in range(n)]
    for (a, b), (k, _, _) in edges.items():
        if a in index and b in index:
            adjacency[index[a]].append((index[b], k))

    # state: (mask, last) -> (score, total_ie, residues, predecessor state)
    start = {}
    for i in range(n):
        start[(1 << i, i)] = (values[i] + lam, values[i], 1, None)
    dp = dict(start)
    order = sorted(dp)
    frontier = list(order)
    while frontier:
        nxt = []
        for state in frontier:
            mask, last = state
            score, total_ie, residues, _ = dp[state]
            for j, k in adjacency[last]:
                if mask & (1 << j):
                    continue
                new_mask = mask | (1 << j)
                added = k + 1                       # k linkers plus residue j
                cand = (score + values[j] + lam * added, total_ie + values[j],
                        residues + added, state)
                key = (new_mask, j)
                if key not in dp or cand[0] < dp[key][0]:
                    dp[key] = cand
                    nxt.append(key)
        frontier = nxt

    best = {}
    for (mask, last), (score, total_ie, residues, _) in dp.items():
        count = bin(mask).count("1")
        if count not in best or score < best[count][0]:
            best[count] = (score, total_ie, residues, (mask, last))
    out = {}
    for count, (score, total_ie, residues, state) in best.items():
        path = []
        while state is not None:
            path.append(nodes[state[1]])
            state = dp[state][3]
        out[count] = (score, total_ie, residues, list(reversed(path)))
    return out


def linker_budget(n_poses):
    """Glycines a path of this many poses needs before it can generically close.

    Counting the whole path at once: m = n_poses - 1 segments give 1 + m + 2*sum(k) free dihedrals
    against 4m closure conditions, so sum(k) >= (3m - 1) / 2. Surplus freedom in one segment pays
    for a deficit in the next, which is why this is a budget over the path rather than a test on
    each step. Necessary, not sufficient -- sterics and the actual geometry still decide.
    """
    m = max(n_poses - 1, 0)
    return int(np.ceil((3 * m - 1) / 2)) if m else 0


def spell(path, edges):
    """The peptide the path prescribes: residues with glycine linkers between them."""
    letters = ONE_LETTER[path[0][0]]
    for a, b in zip(path, path[1:]):
        k = edges[(a, b)][0]
        letters += "G" * k + ONE_LETTER[b[0]]
    return letters


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("outdir")
    parser.add_argument("--ligand", default="octinoxate")
    parser.add_argument("--pairs-csv", default="condense_pairs_geom.csv")
    parser.add_argument("--tol", type=float, default=0.2,
                        help="closure that counts as connected (A)")
    parser.add_argument("--max-penalty", type=float, default=None,
                        help="reject connections whose steric penalty exceeds this (A^2)")
    parser.add_argument("--poses", help="explicit shell as frag:pose,frag:pose,... "
                        "(default: the poses flagged selected in energies.csv)")
    parser.add_argument("--save", help="write the chosen design to this JSON file")
    parser.add_argument("--poses-used", type=int,
                        help="how many poses the saved design should use (default: the most)")
    parser.add_argument("--lam", type=float, default=0.0,
                        help="kcal/mol charged per residue, to discourage long linkers")
    args = parser.parse_args(argv)

    frags = define_fragments()
    ligand = load_ligand(args.ligand)
    atts = selected_attachments(args.outdir, frags, ligand, parse_poses(args.poses))
    nodes = [(a["name"], a["pose"]) for a in atts]
    values = [a["ie"] for a in atts]

    path_csv = args.pairs_csv
    if not os.path.isabs(path_csv) and not os.path.exists(path_csv):
        path_csv = os.path.join(args.outdir, path_csv)
    edges = load_edges(path_csv, args.tol, args.max_penalty)
    possible = len(nodes) * (len(nodes) - 1)
    print(f"{len(edges)} of {possible} ordered pairs are connectable "
          f"(closure < {args.tol} A)\n")

    results = best_paths(nodes, values, edges, args.lam)
    print(f"{'poses':>5} {'residues':>9} {'sum k':>6} {'budget':>7} {'pad':>4} "
          f"{'total IE':>10} {'IE/residue':>11}  sequence")
    for count in sorted(results):
        score, total_ie, residues, path = results[count]
        sum_k = residues - count
        budget = linker_budget(count)
        pad = max(budget - sum_k, 0)
        print(f"{count:>5} {residues:>9} {sum_k:>6} {budget:>7} {pad:>4} "
              f"{total_ie:>10.2f} {total_ie / residues:>11.2f}  {spell(path, edges)}")
    print("\n'budget' is the glycine count the whole-path parameter count demands; 'pad' is how "
          "many\nmore are needed than the per-edge minimums supply. A path short of budget cannot "
          "close\nas spelled, however connectable each step looks on its own.")

    if args.save:
        count = args.poses_used or max(results)
        if count not in results:
            sys.exit(f"no path uses {count} poses; available: {sorted(results)}")
        score, total_ie, residues, path = results[count]
        design = {
            "ligand": args.ligand,
            "path": [[name, pose] for name, pose in path],
            "linkers": [edges[(a, b)][0] for a, b in zip(path, path[1:])],
            "sequence": spell(path, edges),
            "total_ie_kcal_mol": round(total_ie, 4),
            "residues": residues,
            "sum_k": residues - count,
            "budget": linker_budget(count),
            "closures": [round(edges[(a, b)][1], 4) for a, b in zip(path, path[1:])],
        }
        with open(args.save, "w") as fh:
            json.dump(design, fh, indent=2)
        print(f"\nsaved {count}-pose design to {args.save}: {design['sequence']}")

    if results:
        best_count = max(results)
        _, total_ie, residues, path = results[best_count]
        sum_k = residues - best_count
        budget = linker_budget(best_count)
        print(f"\nlongest traversable path uses {best_count} of {len(nodes)} poses "
              f"(sum k = {sum_k}, budget = {budget}"
              + (f", SHORT BY {budget - sum_k}" if sum_k < budget else ", budget met") + "):")
        for a, b in zip(path, path[1:]):
            k, closure, penalty = edges[(a, b)]
            print(f"  {a[0]}{a[1]} -> {b[0]}{b[1]}   {k} linker(s), "
                  f"closure {closure:.3f} A, steric penalty {penalty:.2f}")
        unused = [n for n in nodes if n not in path]
        if unused:
            print("  unused poses: " + ", ".join(f"{n[0]}{n[1]}" for n in unused))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
