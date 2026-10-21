"""Replicate tally: signed ΔΔG per leg and the stronger / weaker / unchanged counts.

The targeted replicate campaign's deliverable (NEXT_STEPS.md, "Targeted replicate strategy"): how
much does the random seed move the dynamics? That is a *distribution*, not a verdict on any one
leg -- both runs of a pair are valid samples under different seeds, so ΔΔG carries a sign and no
run is a "failure". The tally is the input to deciding between more targeted replicates and a full
campaign.

`unchanged` is |ΔΔG| < 1.5, reusing the drift threshold the `clean` diagnostic already uses so the
two sit on one scale.

Reads run 1 from runs/<mol>/md/mmgbsa_summary.csv and run 2 from mmgbsa_replicates.csv -- a
replicate is never in the canonical summary.

    python3 code/replicate_tally.py
"""
import csv
import os

PAIRS = [("octinoxate", "orig_f12_dock1", "strong·clean"),
         ("octinoxate", "s2_esm2_control_dock6", "mid·clean"),
         ("octinoxate", "bg33_4_o3cof_dock1", "mid·clean"),
         ("octinoxate", "s3_esm2_f4_dock1", "strong·flagged"),
         ("octinoxate", "shuffle_control_o3cof_dock8", "weak·flagged"),
         ("oxybenzone", "ox2_orig_f8_dock4", "strong·clean"),
         ("oxybenzone", "bgox31_3_dock1", "strong·clean"),
         ("oxybenzone", "ox2_esm1_f8_o3cof", "mid·clean"),
         ("oxybenzone", "bgox31_4_o3cof_dock7", "weak·clean"),
         ("oxybenzone", "ox1_orig_f12_rf3cof_dock1", "weak·clean"),
         ("oxybenzone", "ox1_esm1_f8_o3cof", "mid·flagged"),
         ("oxybenzone", "ox1_esm1_control", "strong·clean (out of order)"),
         # Second batch (2026-10-21), chosen to reach 10% of the campaign's 230 original legs and to
         # correct the first batch's sampling skew: nine of those twelve were docked legs against a
         # 24% population share, and the rf3cof and o3pcof lanes had no replicate at all.
         ("octinoxate", "shuffle_control_esm0_rf3cof", "strong·clean"),
         ("octinoxate", "orig_f12", "mid·flagged"),
         ("octinoxate", "s2_esm2_control_o3pcof_dock8", "mid·flagged"),
         ("oxybenzone", "bgox31_5_rf3cof_dock1", "mid·clean"),
         ("oxybenzone", "ox1_esm1_f8_rf3cof", "mid·flagged"),
         ("oxybenzone", "ox2_shuffle_rf3cof_dock2", "weak·flagged"),
         ("oxybenzone", "ox2_orig_f8_rf3cof_dock7", "weak·clean"),
         ("oxybenzone", "ox2_orig_f8_o3pcof_dock1", "strong·flagged"),
         ("oxybenzone", "ox1_orig_f12_o3pcof_dock1", "weak·flagged"),
         ("oxybenzone", "ox2_shuffle", "weak·clean"),
         ("oxybenzone", "bgox31_4_o3cof_dock1", "mid·flagged")]
UNCHANGED = 1.5
# Third runs, for the two legs whose first two runs disagreed most -- a pair cannot say which of its
# two values is the outlier, and these two carry the campaign's largest ddG (+10.71 and +16.16).
# They extend an existing PAIRS row rather than adding one: same cell, a third sample of it.
THIRD = {"orig_f12_dock1", "ox1_esm1_control"}


def dg(mol, struct, fname):
    p = f"runs/{mol}/md/{fname}"
    if not os.path.exists(p):
        return None
    for r in csv.DictReader(open(p)):
        # skip every window, leading and trailing alike: a last_<ns>ns row is a convergence
        # diagnostic (md_window_modal.sh MODE=trailing) and would be read as the whole run here
        if r["structure"] == struct and not r["leg"].startswith(("first_", "last_")):
            return float(r["dg_bind"])
    return None


def trailing(mol, leg, ns=10):
    """The leg's trailing-window dG, or None. Written by code/trailing_windows.py, deliberately
    NOT in any csv -- see that script and collect_mmgbsa_rows.py's skip of last_<ns>ns."""
    import glob
    BG = os.path.expanduser("~/python_mac/boltzgen_local/md")
    for root in (f"runs/{mol}/md", BG):
        hits = glob.glob(f"{root}/{leg}/last_{ns}ns/mmgbsa/FINAL_RESULTS_MMPBSA.dat")
        if hits:
            body = open(hits[0]).read().split("Differences")[1]
            return float(next(x for x in body.split("\n")
                              if x.startswith("DELTA TOTAL")).split()[2])
    return None


def trailing_report(ns=10):
    """Headline ddG against back-half ddG, per pair.

    Every window in this project is a LEADING cumulative mean, so a leg that settles late has its
    weak early half averaged into the whole-run value. This asks how much of each pair's ddG
    survives when both runs are compared over their last <ns> ns instead. It is a diagnostic, not a
    replacement: the whole-run value remains the leg's result.
    """
    print(f"{'leg':30s} {'run1':>7s} {'r1 last':>8s} {'run2':>7s} {'r2 last':>8s} "
          f"{'ddG':>8s} {'ddG late':>9s} {'shift':>7s}")
    rows = []
    for mol, struct, _cell in PAIRS:
        a = dg(mol, struct, "mmgbsa_summary.csv")
        b = dg(mol, struct + "_r2", "mmgbsa_replicates.csv")
        la, lb = trailing(mol, struct, ns), trailing(mol, struct + "_r2", ns)
        if None in (a, b, la, lb):
            continue
        d, dl = b - a, lb - la
        rows.append((struct, d, dl))
        print(f"{struct:30s} {a:7.2f} {la:8.2f} {b:7.2f} {lb:8.2f} "
              f"{d:+8.2f} {dl:+9.2f} {abs(dl) - abs(d):+7.2f}")
    if rows:
        grew = sum(1 for _s, d, dl in rows if abs(dl) > abs(d))
        print(f"\n{len(rows)} pairs with both trailing windows   "
              f"|ddG| grows on the back half in {grew}, shrinks in {len(rows) - grew}")
        # The PAIRED shift is the statistic here, not two independent medians: each leg is its own
        # control, and the median of |ddG| whole against the median of |ddG| late compares different
        # legs at the median position (they came out 5.20 and 3.56 on a set that grows in 11 of 13,
        # which reads backwards). Quote the per-leg shift.
        sh = sorted(abs(dl) - abs(d) for _s, d, dl in rows)
        print(f"median shift in |ddG| {sh[len(sh) // 2]:+.2f}   range {sh[0]:+.2f} to {sh[-1]:+.2f}")


def main():
    print(f"{'cell':32s} {'leg':30s} {'run1':>8s} {'run2':>8s} {'ddG':>8s}  {'run3':>8s}  verdict")
    tally = {"stronger": 0, "weaker": 0, "unchanged": 0}
    diffs, done, triples = [], 0, []
    for mol, struct, cell in PAIRS:
        a = dg(mol, struct, "mmgbsa_summary.csv")
        b = dg(mol, struct + "_r2", "mmgbsa_replicates.csv")
        c = dg(mol, struct + "_r3", "mmgbsa_replicates.csv") if struct in THIRD else None
        third = f"{c:8.2f}" if c is not None else ("pending" if struct in THIRD else "")
        if a is None or b is None:
            shown = "--" if a is None else f"{a:8.2f}"
            print(f"{cell:32s} {struct:30s} {shown:>8} {'--':>8} {'--':>8}  {third:>8}  pending")
            continue
        d = b - a
        v = "unchanged" if abs(d) < UNCHANGED else ("stronger" if d < 0 else "weaker")
        tally[v] += 1
        diffs.append(abs(d))
        done += 1
        if c is not None:
            triples.append((struct, [a, b, c]))
        print(f"{cell:32s} {struct:30s} {a:8.2f} {b:8.2f} {d:+8.2f}  {third:>8}  {v}")
    print(f"\n{done}/{len(PAIRS)} complete   stronger {tally['stronger']}   "
          f"weaker {tally['weaker']}   unchanged {tally['unchanged']}   "
          f"(unchanged = |ddG| < {UNCHANGED})")
    if diffs:
        s = sorted(diffs)
        med = s[len(s) // 2] if len(s) % 2 else (s[len(s) // 2 - 1] + s[len(s) // 2]) / 2
        print(f"median |ddG| {med:.2f}   range {min(s):.2f}-{max(s):.2f}")
    # A pair cannot say which of its two values is the outlier; a triple can bracket the leg. Quote
    # the full spread (max - min over the three runs), not a mean, because the question is how far
    # the seed moves one leg rather than where its centre lies.
    for struct, runs in triples:
        print(f"triple  {struct:30s} " + " ".join(f"{r:8.2f}" for r in runs)
              + f"   spread {max(runs) - min(runs):.2f}")


if __name__ == "__main__":
    main()
