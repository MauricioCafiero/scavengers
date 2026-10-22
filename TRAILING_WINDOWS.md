# Trailing convergence windows: a diagnostic, deliberately not folded into anything

**Status: initial work, 2026-10-21/22. Nothing here is in any table.** `last_<ns>ns` rows are
excluded by four separate guards (`record_leg.py`, `watch_legs.py`, `collect_mmgbsa_rows.py`,
`replicate_tally.dg`) and the figures below appear in no other document. Whether they ever should
is the user's open call.

## Why trailing windows exist

**Every convergence window in this project is a *leading cumulative* mean.** `first_5ns`,
`first_10ns`, `first_15ns` and the whole run each average from frame zero forward, so a 20 ns value
is an average over the whole trajectory including its unsettled start. A leg whose binding tightens
late therefore has its weak early half averaged into the number that gets quoted, and a leg that
loosens late has its strong early half propping the number up.

That matters for a replicate *pair*, because a pair's ΔΔG is a difference of two such averages. If
one run settles late and the other is already converged, part of the ΔΔG is the settling rather than
the seed. The back half, measured directly, separates the two.

A trailing window cuts the **last** `<ns>` nanoseconds and runs MM/GBSA (GBn2) over those frames
alone, producing `last_<ns>ns/mmgbsa/`. It is a convergence *shape* measurement, not a result: the
whole-run value remains the leg's number.

## How to compute one

```sh
# one leg, by hand
MODE=trailing WINDOWS=10 SYSNAME=<mol> zsh code/md_window_modal.sh <struct> 20 prod_20ns

# every pair that needs one, serialised, between tails
python3 code/trailing_windows.py --list                       # the plan, nothing computed
python3 code/trailing_windows.py                              # compute, skipping what is done
python3 code/trailing_windows.py --only <struct> --only <struct>   # just these pairs
```

`md_window_modal.sh` gained `MODE` and `WINDOWS` on 2026-10-21; `trailing_windows.py` gained
`--only` on 2026-10-22. Each window is ~1.5–4 min of local CPU and the driver takes the global tail
lock **per window, not for the whole run**, so a finished leg's tail gets in at each boundary. The
~25 MB trajectory slice is deleted afterwards — it reslices from `prod_20ns` in seconds.

**Simple subtraction predicts a back half almost exactly**, because
`mean(0–20) = ½·mean(0–10) + ½·mean(10–20)`. Checked on `shuffle_control_esm0_rf3cof_r2`: predicted
−18.73 against −18.72 measured. So the leading series can be trusted to *imply* a back half when
the window is not worth computing, and a trailing window is only worth cutting when the exact
value, its standard deviation, or both runs of a pair are wanted.

## What the 15 completed pairs show

Last 10 ns of each run against the whole run, and the pair's ΔΔG recomputed back-half-to-back-half.

| structure | run 1 | r1 last 10 | run 2 | r2 last 10 | ΔΔG | ΔΔG late | shift in \|ΔΔG\| |
|---|---:|---:|---:|---:|---:|---:|---:|
| `orig_f12_dock1` | −22.11 | −20.97 | −11.40 | −7.06 | +10.71 | +13.91 | +3.20 |
| `s2_esm2_control_dock6` | −17.09 | −17.41 | −15.19 | −14.67 | +1.90 | +2.74 | +0.84 |
| `bg33_4_o3cof_dock1` | −16.04 | −16.43 | −18.93 | −20.12 | −2.89 | −3.69 | +0.80 |
| `s3_esm2_f4_dock1` | −22.42 | −22.17 | −27.62 | −27.98 | −5.20 | −5.81 | +0.61 |
| `shuffle_control_o3cof_dock8` | −7.78 | −5.97 | −9.20 | −7.59 | −1.42 | −1.62 | +0.20 |
| `ox2_orig_f8_dock4` | −24.00 | −24.67 | −18.44 | −18.17 | +5.57 | +6.50 | +0.93 |
| `bgox31_3_dock1` | −22.23 | −23.12 | −22.21 | −20.97 | +0.02 | +2.15 | +2.13 |
| `ox2_esm1_f8_o3cof` | −16.68 | −15.82 | −16.38 | −14.09 | +0.31 | +1.72 | +1.42 |
| `bgox31_4_o3cof_dock7` | −6.71 | −7.61 | −12.12 | −15.78 | −5.41 | −8.17 | +2.76 |
| `ox1_orig_f12_rf3cof_dock1` | −6.40 | −5.50 | −9.20 | −8.83 | −2.80 | −3.34 | +0.54 |
| `ox1_esm1_f8_o3cof` | −17.95 | −16.66 | −12.57 | −13.11 | +5.38 | +3.56 | −1.82 |
| `ox1_esm1_control` | −26.09 | −27.10 | −9.93 | −10.18 | +16.16 | +16.93 | +0.77 |
| `shuffle_control_esm0_rf3cof` | −21.53 | −21.58 | −15.03 | −18.72 | +6.50 | +2.86 | −3.65 |
| `ox1_orig_f12_o3pcof_dock1` | −3.19 | −1.41 | −15.75 | −16.10 | −12.56 | −14.69 | +2.13 |
| `bgox31_4_o3cof_dock1` | −10.39 | −9.43 | −4.70 | −2.23 | +5.69 | +7.20 | +1.51 |

**|ΔΔG| grows on the back half in 13 of 15 pairs and shrinks in 2. Median shift +0.84, range
−3.65 to +3.20.**

The paired shift is the statistic, not two independent medians. Comparing the median |ΔΔG| whole
against the median |ΔΔG| late compares *different legs* at the median position — they came out 5.20
and 3.56 on a set that grows in 13 of 15, which reads exactly backwards. Each leg is its own
control, so quote the per-leg shift.

### What that means, and what it does not

The hypothesis that motivated this work was that late settling *inflates* pair ΔΔGs, so that
comparing settled halves would shrink them. **Across the set it does the opposite**, and the pair
that prompted the whole investigation is the largest shrinker in it:

- `shuffle_control_esm0_rf3cof` — run 1 flat (whole −21.53, last 10 ns −21.58), run 2 still rising
  (−15.03 → −18.72). Headline ΔΔG +6.50, back-half +2.86, **shift −3.65**: on this pair 56% of the
  ΔΔG really was the settle, and run 2's back half is also *tighter* than run 1's (sd 2.73 vs 4.66),
  i.e. a slow settle rather than a different bound state.
- Generalising from it was wrong. Most shifts are small — |shift| < 1.0 in 8 of 15 and < 2.0 in 11 —
  and **no pair's verdict changes**: nothing crosses the ±1.5 unchanged threshold in either
  direction. Mean |shift| is 1.51 against a mean |ΔΔG| of 4.55.
- The direction is nonetheless real rather than noise (sign test on 13 of 15, p = 0.022) — it is
  just small. The honest summary is that **convergence is not contaminating the pair comparisons
  much**, which is the opposite of the reason the diagnostic was built.

Two pairs were cut on 2026-10-22 because their own diagnostics raised a specific question, and both
**widened**:

- `ox1_orig_f12_o3pcof_dock1` — the campaign's second-largest move (−12.56) and its flattest series.
  Late ΔΔG **−14.69**. Two seeds in genuinely different bound states 14.7 apart, with run 2's back
  half the tighter of the two (sd 1.56 vs 3.23). Not a convergence artefact.
- `bgox31_4_o3cof_dock1` — run 2 releases 38.75% of frames with 69% of it late. Its back half is
  **−2.23** against a whole-run −4.70, so the whole-run value is propped up by the half where the
  ligand was still attached. The reading is "left the site late", not "bound differently throughout".

## Pending

**18 windows across 9 pairs plus 2 third runs**, listed by `python3 code/trailing_windows.py --list`:

`bgox31_5_rf3cof_dock1`, `ox1_esm1_f8_rf3cof`, `ox2_shuffle_rf3cof_dock2`,
`ox2_orig_f8_rf3cof_dock7`, `ox2_orig_f8_o3pcof_dock1`, `ox2_shuffle`,
`s2_esm2_control_o3pcof_dock8`, `orig_f12` (both runs of each), plus the third runs
`orig_f12_dock1_r3` and `ox1_esm1_control_r3`.

At ~1.5–4 min a window that is roughly 40–70 minutes of local CPU. **Not recommended on the current
evidence** — 15 pairs already agree that the back half moves |ΔΔG| by a median 0.84 and changes no
verdict, so the remaining 18 would confirm a pattern rather than test one. The two that have the
clearest claim if anyone wants them anyway are `ox1_esm1_control_r3`, the only leg of the second
batch still tightening at 20 ns (−16.81 → −19.71, drift −2.90, zero release, so its −19.71 is a
ceiling and it is the middle value of a triple), and `orig_f12_dock1_r3`, whose −30.01 sits outside
both of its partners.

(An earlier figure of "19 pending" quoted in session was arithmetic error; the tool's own count is
18, and two of those became eligible only when their third runs recorded.)

## Relationship to the replicate campaign

See `REPLICATES.md`. Trailing windows were built to ask whether that campaign's spread is partly a
convergence artefact. The answer so far is mostly no.
