#!/usr/bin/env python3
"""Append any MM/GBSA leg/window results not yet in mmgbsa_summary.csv.

Scans runs/<ligand>/md/<structure>/<leg>/mmgbsa/FINAL_RESULTS_MMPBSA.dat for the
named structures, parses the DELTA (Differences) block, and appends a row for each
(structure, leg) pair missing from the summary -- matching the column order
structure,leg,vdw,eel,egb,esurf,dg_gas,dg_solv,dg_bind,std_dev,std_err. The
presence check is line-anchored on "<structure>,<leg>," so near-names like
s3_orig_f12 vs orig_f12 do not false-match. Idempotent: re-running adds only what
is new. Used to backfill the convergence windows computed by md_window_modal.sh.

    python code/collect_mmgbsa_rows.py s3_orig_f12 s3_esm2_f4 ...
    python code/collect_mmgbsa_rows.py --all        # every md/<structure>
"""
import argparse, re, sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parent
TERMS = ("VDWAALS", "EEL", "EGB", "ESURF", "DELTA G gas", "DELTA G solv", "DELTA TOTAL")


def parse(fp):
    blk = False; d = {}
    for line in open(fp):
        if line.startswith("Differences"):
            blk = True; continue
        if not blk:
            continue
        m = re.match(r"^(VDWAALS|EEL|EGB|ESURF|DELTA G gas|DELTA G solv|DELTA TOTAL)\s+"
                     r"(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)", line)
        if m:
            d[m.group(1)] = (float(m.group(2)), float(m.group(3)), float(m.group(4)))
    return d if all(t in d for t in TERMS) else None


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("structures", nargs="*")
    ap.add_argument("--all", action="store_true")
    ap.add_argument("--outdir", default="runs/octinoxate")
    ap.add_argument("--md-root", action="append", default=None, metavar="DIR",
                    help="extra directory of prepared MD runs to look in, repeatable. Co-folded "
                         "BoltzGen legs live outside this repo (~/python_mac/boltzgen_local/md), but "
                         "their rows belong in the same mmgbsa_summary.csv -- octinoxate's bg33_1 and "
                         "bg33_2 co-folds are already there. runs/<system>/md is searched first.")
    a = ap.parse_args(argv)
    md = REPO / a.outdir / "md"
    summary = md / "mmgbsa_summary.csv"
    # Roots to resolve a named structure against, in order. The summary stays where it is: one file
    # per ligand, whichever lane a leg ran in.
    roots = [md] + [Path(r).expanduser() for r in (a.md_root or [])]
    existing = summary.read_text().splitlines()
    present = {(l.split(",")[0], l.split(",")[1]) for l in existing if l and "," in l}
    structs = ([p.name for p in md.iterdir() if p.is_dir()] if a.all else a.structures)
    added = []
    for s in structs:
        sdir = next((r / s for r in roots if (r / s).is_dir()), md / s)
        if not sdir.is_dir():
            print(f"skip {s}: no dir in {', '.join(str(r) for r in roots)}"); continue
        for legdir in sorted(sdir.iterdir()):
            fp = legdir / "mmgbsa" / "FINAL_RESULTS_MMPBSA.dat"
            if not fp.exists():
                continue
            leg = legdir.name
            if (s, leg) in present:
                continue
            d = parse(fp)
            if d is None:
                print(f"skip {s}/{leg}: unparsed"); continue
            row = (f"{s},{leg},{d['VDWAALS'][0]:.4f},{d['EEL'][0]:.4f},{d['EGB'][0]:.4f},"
                   f"{d['ESURF'][0]:.4f},{d['DELTA G gas'][0]:.4f},{d['DELTA G solv'][0]:.4f},"
                   f"{d['DELTA TOTAL'][0]:.4f},{d['DELTA TOTAL'][1]:.4f},{d['DELTA TOTAL'][2]:.4f}")
            added.append(row); present.add((s, leg))
    if added:
        with open(summary, "a") as f:
            for r in added:
                f.write(r + "\n")
    print(f"appended {len(added)} rows")
    for r in added:
        print("  " + ",".join(r.split(",")[:2]) + " -> " + r.split(",")[8])
    return 0


if __name__ == "__main__":
    sys.exit(main())
