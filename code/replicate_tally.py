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
         ("oxybenzone", "ox1_esm1_control", "strong·clean (out of order)")]
UNCHANGED = 1.5


def dg(mol, struct, fname):
    p = f"runs/{mol}/md/{fname}"
    if not os.path.exists(p):
        return None
    for r in csv.DictReader(open(p)):
        if r["structure"] == struct and not r["leg"].startswith("first_"):
            return float(r["dg_bind"])
    return None


def main():
    print(f"{'cell':32s} {'leg':30s} {'run1':>8s} {'run2':>8s} {'ddG':>8s}  verdict")
    tally = {"stronger": 0, "weaker": 0, "unchanged": 0}
    diffs, done = [], 0
    for mol, struct, cell in PAIRS:
        a = dg(mol, struct, "mmgbsa_summary.csv")
        b = dg(mol, struct + "_r2", "mmgbsa_replicates.csv")
        if a is None or b is None:
            shown = "--" if a is None else f"{a:8.2f}"
            print(f"{cell:32s} {struct:30s} {shown:>8} {'--':>8} {'--':>8}  pending")
            continue
        d = b - a
        v = "unchanged" if abs(d) < UNCHANGED else ("stronger" if d < 0 else "weaker")
        tally[v] += 1
        diffs.append(abs(d))
        done += 1
        print(f"{cell:32s} {struct:30s} {a:8.2f} {b:8.2f} {d:+8.2f}  {v}")
    print(f"\n{done}/{len(PAIRS)} complete   stronger {tally['stronger']}   "
          f"weaker {tally['weaker']}   unchanged {tally['unchanged']}   "
          f"(unchanged = |ddG| < {UNCHANGED})")
    if diffs:
        s = sorted(diffs)
        med = s[len(s) // 2] if len(s) % 2 else (s[len(s) // 2 - 1] + s[len(s) // 2]) / 2
        print(f"median |ddG| {med:.2f}   range {min(s):.2f}-{max(s):.2f}")


if __name__ == "__main__":
    main()
