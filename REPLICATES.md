# MM/GBSA replicates: how far the random seed moves a 20 ns leg

**Status: parked 2026-10-22 at 23 structures.** The campaign will almost certainly be resumed, so
this document is the holding place for what exists — not a finished analysis, and not a policy. No
number here has been folded into `RESULTS.md`, `FOLDS.md` or
`runs/folds_md/cofold_binding_tables.md`, and none should be until the user decides how replicates
are to be reported **across every structure that has one**, rather than case by case.

## What a replicate is here

The same system, byte-identical `system.xml`, re-run with a different RNG seed. Velocities are the
only difference. Both trajectories are valid samples: neither is the true value, run 1 is not the
reference, and **no run is a failure**. A pair measures how much chance moves the dynamics, which is
why every number below carries a sign.

Rows live in `runs/<mol>/md/mmgbsa_replicates.csv`, never in the canonical
`mmgbsa_summary.csv` — `record_leg.py` routes any `_r<N>` leg there. `make_gnina_bundle.py`'s
`REFERENCE` dict is deliberately **not** extended with replicate values.

`unchanged` is |ΔΔG| < 1.5 kcal/mol, reusing the threshold the `clean` drift diagnostic already
uses so the two sit on one scale. Method is MM/GBSA with GBn2 (`igb=8`, 0.15 M salt) over all
20,000 frames; no Poisson–Boltzmann anywhere in this project.

## The 23

10% of the campaign's 228 original MD legs. Twelve were chosen 2026-10-19 as a targeted set; eleven
more on 2026-10-21 to correct that set's sampling skew (nine of the first twelve were docked legs
against a 24% population share, and the `rf3cof` and `o3pcof` lanes had no replicate at all). The
`cell` column is the parent leg's descriptor at selection time, not an identity.

| structure | cell | run 1 | run 2 | run 3 | signed ΔΔG | verdict |
|---|---|---:|---:|---:|---:|---|
| `orig_f12_dock1` | strong·clean | −22.11 | −11.40 | −30.01 | +10.71 | weaker |
| `s2_esm2_control_dock6` | mid·clean | −17.09 | −15.19 |  | +1.90 | weaker |
| `bg33_4_o3cof_dock1` | mid·clean | −16.04 | −18.93 |  | −2.89 | stronger |
| `s3_esm2_f4_dock1` | strong·flagged | −22.42 | −27.62 |  | −5.20 | stronger |
| `shuffle_control_o3cof_dock8` | weak·flagged | −7.78 | −9.20 |  | −1.42 | unchanged |
| `ox2_orig_f8_dock4` | strong·clean | −24.00 | −18.44 |  | +5.57 | weaker |
| `bgox31_3_dock1` | strong·clean | −22.23 | −22.21 |  | +0.02 | unchanged |
| `ox2_esm1_f8_o3cof` | mid·clean | −16.68 | −16.38 |  | +0.31 | unchanged |
| `bgox31_4_o3cof_dock7` | weak·clean | −6.71 | −12.12 |  | −5.41 | stronger |
| `ox1_orig_f12_rf3cof_dock1` | weak·clean | −6.40 | −9.20 |  | −2.80 | stronger |
| `ox1_esm1_f8_o3cof` | mid·flagged | −17.95 | −12.57 |  | +5.38 | weaker |
| `ox1_esm1_control` | strong·clean (out of order) | −26.09 | −9.93 | −19.71 | +16.16 | weaker |
| `shuffle_control_esm0_rf3cof` | strong·clean | −21.53 | −15.03 |  | +6.50 | weaker |
| `orig_f12` | mid·flagged | −13.71 | −21.07 |  | −7.36 | stronger |
| `s2_esm2_control_o3pcof_dock8` | mid·flagged | −12.44 | −17.24 |  | −4.80 | stronger |
| `bgox31_5_rf3cof_dock1` | mid·clean | −14.93 | −15.22 |  | −0.29 | unchanged |
| `ox1_esm1_f8_rf3cof` | mid·flagged | −12.12 | −10.39 |  | +1.73 | weaker |
| `ox2_shuffle_rf3cof_dock2` | weak·flagged | −9.29 | −9.22 |  | +0.07 | unchanged |
| `ox2_orig_f8_rf3cof_dock7` | weak·clean | −4.11 | −9.25 |  | −5.14 | stronger |
| `ox2_orig_f8_o3pcof_dock1` | strong·flagged | −20.42 | −18.50 |  | +1.92 | weaker |
| `ox1_orig_f12_o3pcof_dock1` | weak·flagged | −3.19 | −15.75 |  | −12.56 | stronger |
| `ox2_shuffle` | weak·clean | −5.30 | −4.44 |  | +0.86 | unchanged |
| `bgox31_4_o3cof_dock1` | mid·flagged | −10.39 | −4.70 |  | +5.69 | weaker |

Two structures have a third run, taken because their first two disagreed most. **Neither third run
resolved its pair**:

| structure | run 1 | run 2 | run 3 | spread |
|---|---:|---:|---:|---:|
| `orig_f12_dock1` | −22.11 | −11.40 | −30.01 | 18.61 |
| `ox1_esm1_control` | −26.09 | −9.93 | −19.71 | 16.16 |

`orig_f12_dock1`'s third run went **further out than both** of the first two; `ox1_esm1_control`'s
landed interior, 6.4 from run 1 and 9.8 from run 2, leaving the spread unchanged. So on these two
cells a third sample widened or filled the range rather than identifying an outlier — which is a
reason to be cautious about assuming more runs per cell converge on an answer.

## What the set shows so far

Read as the state of a parked, partial campaign. **The user has not drawn conclusions from it.**

- Directions are near-balanced (stronger 8, weaker 9, unchanged 6; mean signed +0.39 against a sd
  of 6.23), so **run 1 carries no systematic optimism** — the result is spread, not bias.
- Median |ΔΔG| 4.80, mean 4.55, range 0.02–16.16. Within 1.5: 6 of 23; within 3: 11; within 6: 18.
- That spread is the resolution at which the existing tables can be read. It does not invalidate
  them: designs against nulls (−18.84 vs −7 to −8) clears it comfortably; fine ranking among design
  arms does not.
- **Within-trajectory standard errors (0.01–0.10) are not uncertainty on a comparison.** They
  measure frames inside one run. Never write "this difference is N standard errors" of a pair.
- **Convergence windows do not predict reproducibility**, and this set contains the sharpest case:
  `ox1_orig_f12_o3pcof_dock1_r2` has the flattest series in the campaign (−15.32 / −15.40 / −15.59
  / −15.75, drift 0.43, se 0.01), zero release and 100% residence — and sits 12.56 kcal/mol from
  its partner. Perfect internal convergence says nothing about where an independent seed lands.
- Three predictors were tested post hoc and **none held**: run-1 depth, release fraction, and
  within-run window drift. Do not revive any of them without new data.

## Open decisions (the user's, not settled here)

1. **How a replicated cell is reported in the published tables.** The obvious instance is
   `ox1_esm1_control`, still carried as its run-1 **−26.09** in `cofold_binding_tables.md:64`,
   `RESULTS.md:422` and `RESULTS.md:509` though it is now n=3 across 16.2 kcal/mol. Deliberately
   left as is: a form chosen for one peptide has to apply to every structure with a replicate.
2. **Scope.** Stop at 10%, or run the remaining ~202 structures.
3. Whether the trailing-window diagnostic ever enters a table — see `TRAILING_WINDOWS.md`.

## Reproducing the table

```sh
python3 code/replicate_tally.py          # the live tally; PAIRS + THIRD define the set
```

`replicate_tally.py` is the source of truth for membership and arithmetic; this document is a
snapshot of it taken when the campaign was parked.
