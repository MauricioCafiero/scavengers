"""What the construction actually delivers: how many designed side-chain PAIRS reach the ligand.

The design pipeline does not deliver a structure. Boltz rarely reproduces the designed arrangement --
`s3_orig_f12`, the best binder measured, reproduces **none** of its twelve designed positions -- so the
question of what the construction is worth cannot be answered by asking how closely the fold matches it.

What it delivers is a set of positions each chosen to contact the ligand, and at n positions that is
n(n-1)/2 pairwise combinations. Chain connectivity hands you the n-1 adjacent pairs for free, so the
informative band is **[n-1, n(n-1)/2]** and the informative content is the non-adjacent pairs: those
require the fold to bring sequence-distant slots onto the ligand together. This measures, over a
trajectory, how much of that band a fold realises.

The measure separates a design from its own shuffle where mean contact count does not. Shell 3's design
and a shuffle of it -- identical residues, identical spacer pattern, identical length and charge --
realise 53/66 and 18/66 pairs; the design reaches every gap out to 11, both termini on the ligand at
once, while the shuffle reaches nothing beyond gap 8 and sustains only 8 pairs, below the free floor
of 11. It also explains ligand release without extra assumptions: 3.82 side chains engaged on average
means losing one does not detach the ligand, where 1.75 means it often might.

Designed slots come from the parent design's non-glycine positions. Glycine is never in the fragment
library, so in a glycine-spaced design every non-G position is a placed fragment, and an ESM2 variant
keeps those slots because `fill_linkers` only rewrites the linkers.

Usage:
    ~/miniforge3/envs/openmm-md/bin/python code/pair_contacts.py \
        runs/octinoxate/md/s3_orig_f12/prod_L1_modal --design design_shell3.json
    ... several legs at once; --design is guessed from each structure name if omitted
"""
import argparse
import csv
import itertools
import json
import os
import re

AA3 = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E",
       "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F",
       "PRO": "P", "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V"}


def design_for(structure):
    """Which design a structure's slots come from, by name. Shuffles inherit their parent's."""
    if structure.startswith("s2_"):
        return "design_shell2.json"
    if structure.startswith("s3_") or structure.startswith("shuffle_control"):
        return "design_shell3.json"
    return "design.json"


def slots_of(run_dir, design_file, linker="G"):
    seq = json.load(open(os.path.join(run_dir, design_file)))["sequence"]
    return [i for i, c in enumerate(seq) if c != linker], seq


def analyse(leg_dir, slots, cutoff=4.0, persist=0.05):
    import mdtraj as md
    import numpy as np

    xtc = os.path.join(leg_dir, "traj_wrapped.xtc")
    pdb = os.path.join(leg_dir, "traj_wrapped.pdb")
    if not (os.path.exists(xtc) and os.path.exists(pdb)):
        return None
    t = md.load(xtc, top=pdb)
    top = t.topology
    lig_res = {"LIG", "UNK", "UNL"}
    lig = [a.index for a in top.atoms if a.residue.name in lig_res and a.element.symbol != "H"]
    res = [r for r in top.residues if r.name not in lig_res]

    keys, con = [], []
    for k in slots:
        if k >= len(res):
            continue
        r = res[k]
        # side chain only: the backbone contacts the ligand for reasons the design did not choose
        sc = [a.index for a in r.atoms if a.element.symbol != "H"
              and a.name not in ("N", "CA", "C", "O", "OXT")]
        if not sc:                       # glycine has none, and is never a designed slot anyway
            continue
        keys.append((k, AA3.get(r.name, "?")))
        pairs = np.array([(a, b) for a in sc for b in lig])
        con.append((md.compute_distances(t, pairs) * 10 < cutoff).any(axis=1))
    con = np.array(con)
    n = len(keys)
    if n < 2:
        return None

    n_sim = con.sum(0)
    occ = con.mean(1)
    by_gap, realised, held = {}, 0, 0
    for i, j in itertools.combinations(range(n), 2):
        gap = j - i
        both = con[i] & con[j]
        d = by_gap.setdefault(gap, [0, 0, 0])
        d[0] += 1
        if both.any():
            d[1] += 1
            realised += 1
        if both.mean() >= persist:
            d[2] += 1
            held += 1
    ceiling = n * (n - 1) // 2
    floor = n - 1
    adj = by_gap.get(1, [0, 0, 0])
    nonadj_poss = sum(v[0] for g, v in by_gap.items() if g > 1)
    nonadj_real = sum(v[1] for g, v in by_gap.items() if g > 1)
    longr_poss = sum(v[0] for g, v in by_gap.items() if g >= 7)
    longr_real = sum(v[1] for g, v in by_gap.items() if g >= 7)
    return {"leg": leg_dir, "frames": int(t.n_frames), "slots": n,
            "floor": floor, "ceiling": ceiling,
            "simultaneous_mean": round(float(n_sim.mean()), 2),
            "simultaneous_max": int(n_sim.max()),
            "two_or_more_pct": round(100 * float((n_sim >= 2).mean()), 0),
            "pairs_realised": realised,
            "pairs_realised_pct_of_ceiling": round(100 * realised / ceiling, 0),
            "band_position_pct": round(100 * (realised - floor) / (ceiling - floor), 0),
            "pairs_held": held,
            "adjacent": f"{adj[1]}/{adj[0]}",
            "nonadjacent": f"{nonadj_real}/{nonadj_poss}",
            "nonadjacent_pct": round(100 * nonadj_real / nonadj_poss, 0) if nonadj_poss else 0,
            "longrange": f"{longr_real}/{longr_poss}",
            "dead_slots": sum(1 for x in occ if x == 0),
            "labels": " ".join(f"{c}{k}" for k, c in keys),
            "occupancy": " ".join(f"{occ[i]:.2f}" for i in range(n)),
            "by_gap": by_gap}


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0],
                               formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("legs", nargs="+", help="production directories with traj_wrapped.xtc/.pdb")
    p.add_argument("--run-dir", default="runs/octinoxate", help="where the design jsons live")
    p.add_argument("--design", help="design json for every leg, instead of guessing per name")
    p.add_argument("--cutoff", type=float, default=4.0, help="contact cutoff, A (default 4.0)")
    p.add_argument("--persist", type=float, default=0.05,
                   help="a pair counts as held if simultaneous in this fraction (default 0.05)")
    p.add_argument("--csv", help="write the summary here")
    p.add_argument("--gaps", action="store_true", help="also print the per-gap breakdown")
    args = p.parse_args(argv)

    rows = []
    for leg in args.legs:
        structure = os.path.basename(os.path.dirname(leg.rstrip("/")))
        design = args.design or design_for(structure)
        slots, seq = slots_of(args.run_dir, design)
        r = analyse(leg, slots, args.cutoff, args.persist)
        if r is None:
            print(f"skipping {leg}: no wrapped trajectory, or too few slots")
            continue
        r["structure"], r["design"] = structure, design
        rows.append(r)
        print(f"\n=== {structure}  (slots from {design})")
        print(f"  {r['slots']} slots: {r['labels']}")
        print(f"  simultaneous side chains on the ligand: mean {r['simultaneous_mean']}, "
              f"max {r['simultaneous_max']}, >=2 in {r['two_or_more_pct']:.0f}% of frames")
        print(f"  pairs realised {r['pairs_realised']}/{r['ceiling']} "
              f"({r['pairs_realised_pct_of_ceiling']:.0f}% of ceiling, "
              f"{r['band_position_pct']:.0f}% of the way up the [{r['floor']}, {r['ceiling']}] band)")
        print(f"  held >={args.persist:.0%} of the run: {r['pairs_held']}   "
              f"adjacent {r['adjacent']}   non-adjacent {r['nonadjacent']} "
              f"({r['nonadjacent_pct']:.0f}%)   long-range gap>=7 {r['longrange']}")
        print(f"  slots never in contact: {r['dead_slots']}")
        if args.gaps:
            print(f"    {'gap':>5}{'poss':>6}{'real':>6}{'held':>6}")
            for g in sorted(r["by_gap"]):
                v = r["by_gap"][g]
                print(f"    {g:>5}{v[0]:>6}{v[1]:>6}{v[2]:>6}")

    if len(rows) > 1:
        print(f"\n{'structure':22}{'slots':>6}{'sim':>6}{'>=2%':>6}{'pairs':>9}{'band%':>7}"
              f"{'held':>6}{'adjacent':>10}{'non-adj':>9}{'long':>7}{'dead':>6}")
        for r in rows:
            pairs = f"{r['pairs_realised']}/{r['ceiling']}"
            print(f"{r['structure']:22}{r['slots']:>6}{r['simultaneous_mean']:>6.2f}"
                  f"{r['two_or_more_pct']:>6.0f}{pairs:>9}"
                  f"{r['band_position_pct']:>7.0f}{r['pairs_held']:>6}{r['adjacent']:>10}"
                  f"{r['nonadjacent']:>9}{r['longrange']:>7}{r['dead_slots']:>6}")

    if args.csv and rows:
        fields = [k for k in rows[0] if k != "by_gap"]
        with open(args.csv, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=fields, extrasaction="ignore")
            w.writeheader()
            w.writerows(rows)
        print(f"\nwrote {args.csv}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
