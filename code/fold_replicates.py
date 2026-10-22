"""Re-fold a co-fold input several times to measure how far independent draws land apart.

**Nothing here is seeded, deliberately.** `boltz_env.run_boltz` passes no `--seed` and Boltz's own
default is "None (no seeding)", so every call is an independent draw from an unseeded stream. That is
the whole point: the user ruled on 2026-10-22 that the spread across draws is the quantity to be
measured, not something to eliminate, so **do not add a seed flag and do not re-fold with a fixed
seed "for reproducibility"** -- that would destroy what this script exists to measure.

What it does NOT change: the original co-fold of each structure stays exactly where it is, under
`runs/<mol>/boltz/boltz_results_<name>/`, and remains the published fold. Replicates are written to
`runs/<mol>/folds/foldrep/<name>/rep<k>/` so no existing path, table or leg is touched. Compare them
with `compare_folds.py --replicates`, which reports the RMSDs.

This re-uses the *same yaml* the published fold was made from (`runs/<mol>/boltz/<name>.yaml`) rather
than rebuilding an input, so the only difference between a replicate and the original is the RNG.

Why this is separate from the MM/GBSA replicate campaign (`REPLICATES.md`): those re-ran dynamics
from a FIXED starting structure, varying velocities only. This varies the structure itself, which
every MD leg, dock and table cell downstream is built from, so the replicate spread does not bound
it.

    python3 code/fold_replicates.py oxybenzone --list              # the plan, nothing run
    python3 code/fold_replicates.py oxybenzone --n 3               # every co-fold, 3 draws each
    python3 code/fold_replicates.py oxybenzone --n 5 --only ox1_esm1_control
"""

import argparse
import os
import sys
import time

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(REPO, "code"))
os.chdir(REPO)


def cif_of(out_dir, name):
    return f"{out_dir}/boltz_results_{name}/predictions/{name}/{name}_model_0.cif"


def plan(mol, n, only=()):
    """[(name, k, yaml, out_dir)] for every draw still missing, plus what is already there."""
    import glob
    todo, done = [], []
    for y in sorted(glob.glob(f"runs/{mol}/boltz/*.yaml")):
        name = os.path.basename(y)[: -len(".yaml")]
        if only and name not in only:
            continue
        for k in range(1, n + 1):
            out = f"runs/{mol}/folds/foldrep/{name}/rep{k}"
            (done if os.path.exists(cif_of(out, name)) else todo).append((name, k, y, out))
    return todo, done


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("mol", help="octinoxate or oxybenzone")
    ap.add_argument("--n", type=int, default=3, help="draws per structure (default 3)")
    ap.add_argument("--only", action="append", default=[], metavar="NAME",
                    help="restrict to this co-fold (repeatable)")
    ap.add_argument("--list", action="store_true", help="print the plan and exit")
    args = ap.parse_args()

    if not os.path.isdir(f"runs/{args.mol}/boltz"):
        sys.exit(f"no runs/{args.mol}/boltz -- is {args.mol!r} a molecule with co-folds?")

    todo, done = plan(args.mol, args.n, tuple(args.only))
    for name, k, _y, _o in done:
        print(f"have  {name:34s} rep{k}")
    for name, k, _y, out in todo:
        print(f"TODO  {name:34s} rep{k}  -> {out}")
    print(f"{len(todo)} draws to fold, {len(done)} already present", flush=True)
    if args.list or not todo:
        return 0

    from boltz_env import run_boltz
    for i, (name, k, y, out) in enumerate(todo, 1):
        os.makedirs(out, exist_ok=True)
        t0 = time.time()
        rc = run_boltz(y, out, f"{out}/{name}.log")
        got = os.path.exists(cif_of(out, name))
        print(f"[{i}/{len(todo)}] {name:34s} rep{k}  "
              f"{'ok' if got else f'FAILED rc={rc}, see {out}/{name}.log'}  "
              f"({time.time() - t0:.0f} s)", flush=True)
    print("FOLD_REPLICATES_DONE", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
