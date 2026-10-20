#!/usr/bin/env python3
"""Wait for the next leg to finish, print its result, and exit.

The one waiter for every campaign. Two throwaway versions of this lived in a scratchpad for
weeks -- `watch_leg.sh` (keyed on the record log) and `watch_rep_tails.sh` (keyed on MM/GBSA) --
and because neither was reusable, the second was written from scratch rather than extended and
shipped with a path bug that silently swallowed three legs' reports on 2026-10-20.

Two properties it must have, both learned the hard way:

**It exits on the first completion.** The harness re-invokes on a background job's *exit*, not on
its output, so a waiter that reports and keeps looping is mute however much it prints. Report,
then re-arm it.

**Its trigger is verified against a leg that has already landed** (`--selftest`), because the two
bugs that hid longest were a wrong path and a signal that fires early -- neither of which shows up
by reading the code.

Usage:
    code/watch_legs.py octinoxate/s3_esm2_f4_dock1_r2 oxybenzone/ox2_orig_f8_dock4_r2
    code/watch_legs.py --legfile runs/replicates/rep_oxb_1.legs
    code/watch_legs.py --any                 # any leg that records, campaign-wide
    code/watch_legs.py --legfile … --selftest # resolve every trigger, print, exit
"""

from __future__ import annotations

import argparse
import csv
import glob
import os
import sys
import time
from datetime import datetime

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# `record` is the true done-signal: record_leg.py creates the log empty and writes its `audit:`
# line last, so keying on the file appearing fires ~16 min early (seen 2026-10-19). `mmgbsa` is
# the earlier, cruder signal -- the whole-run free energy is readable while the window slices are
# still computing, which is worth having when the number is the thing being waited on.
SIGNALS = ("record", "mmgbsa")


def leg_dir(mol, leg):
    return os.path.join(REPO, "runs", mol, "md", leg)


def resolve(spec):
    """`mol/leg` as given, or a bare leg name looked up under the molecules that exist."""
    if "/" in spec:
        mol, leg = spec.split("/", 1)
        return mol, leg
    mols = sorted(d for d in os.listdir(os.path.join(REPO, "runs"))
                  if os.path.isdir(os.path.join(REPO, "runs", d, "md")))
    hits = [m for m in mols if os.path.isdir(leg_dir(m, spec))]
    if len(hits) != 1:
        sys.exit(f"{spec}: matches {hits or 'no molecule'} under runs/*/md/ -- qualify it as mol/leg")
    return hits[0], spec


def trigger(mol, leg, signal):
    if signal == "mmgbsa":
        return os.path.join(leg_dir(mol, leg), "prod_20ns", "mmgbsa", "FINAL_RESULTS_MMPBSA.dat")
    return os.path.join(REPO, "runs", "folds_md", f"record_{leg}.log")


def fired(path, signal):
    if not os.path.exists(path):
        return False
    if signal == "mmgbsa":
        return True
    # the record log exists from the start of the tail; the audit line is written last
    with open(path) as fh:
        return any(line.startswith("audit:") for line in fh)


def read_legfile(path):
    """`<local dir>:<leg name>` lines, the same format racc_feed.sh and racc_drain.sh take."""
    specs = []
    for line in open(path):
        line = line.split("#")[0].strip()
        if not line or ":" not in line:
            continue
        d, leg = line.rsplit(":", 1)
        parts = os.path.normpath(d).split(os.sep)
        mol = parts[parts.index("runs") + 1] if "runs" in parts else None
        specs.append(f"{mol}/{leg}" if mol else leg)
    return specs


def dg(path):
    """The Differences block, not the complex section -- a naive VDWAALS grep returns -hundreds."""
    body = open(path).read().split("Differences")[1]
    line = next(x for x in body.split("\n") if x.startswith("DELTA TOTAL"))
    f = line.split()
    return float(f[2]), float(f[4])


def report(mol, leg):
    d = leg_dir(mol, leg)
    out = [f"=== {leg}  ({mol})  {datetime.now():%H:%M:%S}"]

    windows = []
    for n in (5, 10, 15):
        hits = glob.glob(f"{d}/first_{n}ns/**/FINAL_RESULTS_MMPBSA.dat", recursive=True)
        windows.append(dg(hits[0])[0] if hits else None)
    whole = [f for f in glob.glob(f"{d}/**/FINAL_RESULTS_MMPBSA.dat", recursive=True)
             if "/first_" not in f]
    if not whole:
        out.append("  no whole-run MM/GBSA yet")
        return "\n".join(out)
    value, se = dg(whole[0])
    windows.append(value)
    shown = " / ".join(f"{w:.2f}" if w is not None else "--" for w in windows)
    out.append(f"  dG {value:.2f}  se {se:.2f}   windows {shown}")

    # A replicate is scored against its own run 1. Neither run is the reference; the pair measures
    # how far the seed moves the dynamics, so the sign is kept and nothing is called a failure.
    if leg.endswith("_r2"):
        summary = os.path.join(REPO, "runs", mol, "md", "mmgbsa_summary.csv")
        rows = [r for r in csv.DictReader(open(summary))
                if r["structure"] == leg[:-3] and not r["leg"].startswith("first_")]
        if rows:
            first = float(rows[0]["dg_bind"])
            out.append(f"  run 1 {first:.2f}   run 2 {value:.2f}   signed ddG {value - first:+.2f}")

    # md_contacts.csv is BINNED -- 10 rows of 2 ns, whose `contacts` column is the bin mean. So a
    # mean over it is fair, but counting its zero rows is not a release measure: a 2 ns average is
    # essentially never exactly 0 even on a leg that releases constantly. The real `released_pct`
    # (fraction of all 2000 frames with zero atom-atom contacts under 4 A, md_contacts.py:74) is
    # per-frame, is never written to this csv, and reaches us only via the record log below. Report
    # it from there or say it is pending -- do not substitute anything derived from the bins.
    contacts = os.path.join(d, "md_contacts.csv")
    if os.path.exists(contacts):
        c = [float(r["contacts"]) for r in csv.DictReader(open(contacts))]
        out.append(f"  mean contacts {sum(c) / len(c):.1f}  (binned: {len(c)} x 2 ns means)")

    log = os.path.join(REPO, "runs", "folds_md", f"record_{leg}.log")
    lines = [l.rstrip() for l in open(log)] if os.path.exists(log) else []
    shown_audit = False
    for line in lines:
        if line.startswith(("audit:", f"{leg}:")):
            out.append("  " + line)
            shown_audit = True
    if not shown_audit:
        out.append("  release/residence: pending (record_leg.py has not written its line yet)")
    return "\n".join(out)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("legs", nargs="*", help="mol/leg (or a bare leg name)")
    ap.add_argument("--legfile", action="append", default=[],
                    help="a racc legfile; its legs are watched (repeatable)")
    ap.add_argument("--any", action="store_true",
                    help="any leg that records, not a named set")
    ap.add_argument("--on", default="record", choices=SIGNALS,
                    help="record = the audit line (true done-signal); "
                         "mmgbsa = the whole-run dG, ~16 min earlier")
    ap.add_argument("--max-minutes", type=float, default=120)
    ap.add_argument("--interval", type=float, default=30)
    ap.add_argument("--selftest", action="store_true",
                    help="resolve every trigger path, say which already exist, and exit")
    args = ap.parse_args()

    specs = list(args.legs) + [s for f in args.legfile for s in read_legfile(f)]
    watch = [resolve(s) for s in specs]
    if not watch and not args.any:
        sys.exit("name some legs, pass --legfile, or pass --any")

    if args.selftest:
        for mol, leg in watch:
            path = trigger(mol, leg, args.on)
            print(f"{'FIRED  ' if fired(path, args.on) else 'waiting'}  "
                  f"{os.path.relpath(path, REPO)}")
        return 0

    def pending():
        if args.any:
            return {os.path.basename(f)[len("record_"):-len(".log")]
                    for f in glob.glob(f"{REPO}/runs/folds_md/record_*.log")
                    if fired(f, "record")}
        return {leg for mol, leg in watch if fired(trigger(mol, leg, args.on), args.on)}

    seen = pending()
    deadline = time.time() + args.max_minutes * 60
    while time.time() < deadline:
        time.sleep(args.interval)
        new = pending() - seen
        if new:
            for leg in sorted(new):
                mol = next((m for m, l in watch if l == leg), None) or resolve(leg)[0]
                print(report(mol, leg), flush=True)
            return 0            # exit so the harness re-invokes; the caller re-arms
    print(f"no leg finished within {args.max_minutes:g} min "
          f"(as of {datetime.now():%H:%M:%S})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
