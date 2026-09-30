#!/usr/bin/env python3
"""Join GNINA's scored SDFs to everything already known about these poses.

Reads the score tags GNINA writes into its output SDF -- CNNscore, CNNaffinity, minimizedAffinity --
and merges them with reference.csv, which carries each pose's RMSD to the predicted pose, its Vina
score, and the structure's MM/GBSA dG, residence and release count.

Then it answers the question the bundle was built for. Vina ranked the pose nearest the predicted one
6th to 9th of nine in all eight structures, so `reference_rank` is the headline: where does GNINA place
pose 0, the pose that was actually predicted and simulated, among the ten it was shown? A rank of 1
means the CNN succeeds where Vina's scoring function failed.

    python collect_gnina.py results_receptor_score_only
"""
import argparse
import csv
import glob
import os
import sys

TAGS = ["CNNscore", "CNNaffinity", "CNNvariance", "minimizedAffinity", "Affinity"]


def read_sdf_tags(path):
    """[(name, {tag: value})] per molecule in an SDF.

    Deliberately not RDKit: this has to run on whatever machine GNINA ran on, and no chemistry is
    needed to read the score tags. Records are split on the '$$$$' terminator, the name is the first
    line of each record, and a '> <tag>' line is followed by its value on the next non-blank line.
    """
    out = []
    for chunk in open(path).read().split("$$$$"):
        lines = chunk.split("\n")
        while lines and lines[0].strip() == "":
            lines.pop(0)
        if not lines or len(lines) < 4:
            continue
        name = lines[0].strip()
        props, i = {}, 0
        while i < len(lines):
            line = lines[i].strip()
            if line.startswith("> <") or line.startswith(">  <"):
                tag = line[line.index("<") + 1:line.rindex(">")]
                j = i + 1
                while j < len(lines) and lines[j].strip() == "":
                    j += 1
                if j < len(lines):
                    props[tag] = lines[j].strip()
                i = j
            i += 1
        out.append((name, props))
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("results_dir")
    ap.add_argument("--reference", default="reference.csv")
    a = ap.parse_args(argv)

    ref = {}
    with open(a.reference) as fh:
        for r in csv.DictReader(fh):
            ref[(r["structure"], int(r["pose"]))] = r

    rows = []
    for path in sorted(glob.glob(os.path.join(a.results_dir, "*_scored.sdf"))):
        for name, props in read_sdf_tags(path):
            st = props.get("structure")
            pose = props.get("pose")
            if st is None or pose is None:
                continue
            base = ref.get((st, int(pose)), {})
            row = {"structure": st, "pose": int(pose),
                   "is_reference": base.get("is_reference", ""),
                   "rmsd_to_reference_A": base.get("rmsd_to_reference_A", ""),
                   "vina_score": base.get("vina_score", ""),
                   "dg": base.get("dg", ""), "residence": base.get("residence", ""),
                   "releases": base.get("releases", "")}
            for t in TAGS:
                row[t] = props.get(t, "")
            rows.append(row)

    if not rows:
        sys.exit(f"no *_scored.sdf with score tags found in {a.results_dir}")
    rows.sort(key=lambda r: (r["structure"], r["pose"]))
    out = os.path.join(a.results_dir, "gnina_scores.csv")
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print(f"wrote {out} ({len(rows)} poses)\n")

    # Where GNINA puts the predicted pose among the ten, by each score it reported.
    by_struct = {}
    for r in rows:
        by_struct.setdefault(r["structure"], []).append(r)
    print(f"{'structure':<22}{'n':>3}  " + "".join(f"{t[:14]:>16}" for t in
          ("CNNscore rank", "CNNaffinity rk", "vina rank")))
    for st, group in sorted(by_struct.items()):
        line = f"{st:<22}{len(group):>3}  "
        for tag, better_is_high in (("CNNscore", True), ("CNNaffinity", True), ("vina_score", False)):
            vals = [(r["pose"], r.get(tag, "")) for r in group if str(r.get(tag, "")).strip() != ""]
            try:
                vals = [(p, float(v)) for p, v in vals]
            except ValueError:
                line += f"{'-':>16}"
                continue
            vals.sort(key=lambda t: -t[1] if better_is_high else t[1])
            rank = next((i + 1 for i, (p, _) in enumerate(vals) if p == 0), None)
            line += f"{(str(rank) + '/' + str(len(vals))) if rank else '-':>16}"
        print(line)
    print("\nrank 1 = the predicted pose scored best of the ten. Vina never managed better than 6th.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
