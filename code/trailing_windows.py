"""Trailing 10 ns MM/GBSA for every replicate pair, run between tails.

Why this exists: every convergence window in this project is a LEADING cumulative mean (5/10/15 ns
from the start), so a leg whose binding tightens late has its weak early half averaged into the
whole-run value that gets quoted. Against a converged partner that inflates the pair's ddG, and
unlike noise the bias has a direction. The back half, measured directly, separates "the seed moved
the answer" from "20 ns was not enough for this system".

Measured on the first pair (shuffle_control_esm0_rf3cof, 2026-10-21): run 1 flat at -21.58 against
its -21.53 whole run, run 2 at -18.72 against its -15.03 -- headline ddG +6.50, back-half ddG +2.86.

It does not reimplement the slicing. `md_window_modal.sh MODE=trailing WINDOWS=<ns>` cuts the slice
and runs MM/GBSA; this only decides which legs need one, and serialises them.

**One lock acquisition per window, not one for the whole run.** Each window is ~5 min of CPU and the
machine is simultaneously tailing finished legs, which take 10-18 min and are the campaign's critical
path. Holding /tmp/peptidebuilder_tail.lock across all 24 would park every tail behind two hours of
diagnostics; taking it per window lets a waiting tail in at each boundary, which is what "between
tails" means.

Resumable: a leg whose last_<ns>ns/mmgbsa/FINAL_RESULTS_MMPBSA.dat exists is skipped, so a run
interrupted by a reap or a sleep continues where it stopped.

    python3 code/trailing_windows.py --list      # the plan, nothing computed
    python3 code/trailing_windows.py             # compute, skipping what is done
"""

import argparse
import os
import subprocess
import sys
import time

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BG = os.path.expanduser("~/python_mac/boltzgen_local/md")
LOCK = "/tmp/peptidebuilder_tail.lock"
sys.path.insert(0, os.path.join(REPO, "code"))
os.chdir(REPO)


def locate(mol, leg):
    """(leg dir, production dir name) -- legs sit in two roots and two production layouts.

    Modal legs keep prod_L1_modal and a stride-10 trajectory (2,000 frames), so their trailing
    window is the last 1,000 frames: a uniform subsample, unbiased in the mean, with ~0.07 standard
    error against ~0.03 on an unstrided racc leg. Fine for a convergence shape, which is all this is.
    """
    for root in (f"runs/{mol}/md", BG):
        d = f"{root}/{leg}"
        for prod in ("prod_20ns", "prod_L1_modal"):
            if os.path.exists(f"{d}/{prod}/traj_wrapped.xtc"):
                return d, prod
    return None, None


def plan(ns, only=()):
    """`only` names structures (the PAIRS cell, not the _r2 leg) to restrict the plan to.

    The back-fill of 2026-10-21 showed the back half changes no pair's verdict -- median paired
    shift +0.80 -- so the whole plan is rarely worth 1.75 h of CPU. What is worth it is one pair
    whose own diagnostics raise a specific question, which is what this selects.
    """
    from replicate_tally import PAIRS, THIRD, dg
    todo, skip = [], []
    for mol, struct, _cell in PAIRS:
        if only and struct not in only:
            continue
        if dg(mol, struct + "_r2", "mmgbsa_replicates.csv") is None:
            continue                                  # pair not complete yet; nothing to compare
        legs = [struct, struct + "_r2"]
        if struct in THIRD and dg(mol, struct + "_r3", "mmgbsa_replicates.csv") is not None:
            legs.append(struct + "_r3")
        for leg in legs:
            d, prod = locate(mol, leg)
            if prod is None:
                skip.append((mol, leg, "no wrapped trajectory"))
            elif os.path.exists(f"{d}/last_{ns}ns/mmgbsa/FINAL_RESULTS_MMPBSA.dat"):
                skip.append((mol, leg, "already done"))
            else:
                todo.append((mol, leg, prod))
    return todo, skip


def dg_of(path):
    body = open(path).read().split("Differences")[1]
    f = next(x for x in body.split("\n") if x.startswith("DELTA TOTAL")).split()
    return float(f[2]), float(f[3])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ns", type=int, default=10, help="trailing window length (default 10)")
    ap.add_argument("--list", action="store_true", help="print the plan and exit")
    ap.add_argument("--only", action="append", default=[], metavar="STRUCT",
                    help="restrict to this PAIRS structure (repeatable); both its runs are done")
    args = ap.parse_args()

    todo, skip = plan(args.ns, tuple(args.only))
    for mol, leg, why in skip:
        print(f"skip  {mol:11s} {leg:34s} {why}")
    for mol, leg, prod in todo:
        print(f"TODO  {mol:11s} {leg:34s} {prod}")
    print(f"{len(todo)} windows to compute, {len(skip)} skipped", flush=True)
    if args.list or not todo:
        return 0

    for i, (mol, leg, prod) in enumerate(todo, 1):
        t0 = time.time()
        env = {**os.environ, "MODE": "trailing", "WINDOWS": str(args.ns), "SYSNAME": mol}
        r = subprocess.run(["lockf", "-k", LOCK, "zsh", "code/md_window_modal.sh",
                            leg, "20", prod], env=env, capture_output=True, text=True)
        out = f"runs/{mol}/md/{leg}/last_{args.ns}ns"
        d, _ = locate(mol, leg)
        if d and d != f"runs/{mol}/md/{leg}":
            out = f"{d}/last_{args.ns}ns"
        got = f"{out}/mmgbsa/FINAL_RESULTS_MMPBSA.dat"
        if os.path.exists(got):
            v, sd = dg_of(got)
            # the slice is ~25 MB and regenerates in seconds from the production trajectory, which
            # is kept -- leaving 24 of them behind would be 600 MB of staging area for nothing
            for f in (f"{out}/traj_wrapped.xtc",):
                if os.path.exists(f):
                    os.remove(f)
            print(f"[{i}/{len(todo)}] {leg:34s} last {args.ns} ns  {v:8.2f}  sd {sd:.2f}  "
                  f"({time.time() - t0:.0f} s)", flush=True)
        else:
            print(f"[{i}/{len(todo)}] {leg:34s} FAILED rc={r.returncode}\n"
                  f"{r.stdout[-800:]}{r.stderr[-800:]}", flush=True)
    print("TRAILING_WINDOWS_DONE", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
