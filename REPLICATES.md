# MM/GBSA replicates: how far the random seed moves a 20 ns leg

**Status: running 2026-10-25 — 69/101 legs complete (stronger 22, weaker 24, unchanged 23;
median |ΔΔG| 2.37, range 0.02–16.16; batch 6, 9 legs, committed 2026-10-25).** Batches 5–6 are
submitted and queued (rep5_oct_3/4/5 pending in queue, rep5_oct_6 waiting on the feeder,
rep5_oct_7 and the three rep6 legfiles after). No Monday time anchor and no end-date estimate is
quoted: the queue runs whenever an H100 frees (see the 2026-10-23 incident note below). The
campaign is mid-flight after its 2026-10-22 pause; this document is the holding place for what
exists — not a finished analysis, and not a policy. **Weekly racc allocation near its cap as of
2026-10-25 — the campaign may pause mid-batch; the legfiles and the tally above are the
resume point.** No number here has been folded into `RESULTS.md`, `FOLDS.md` or
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

## Batch 3: eleven legs queued 2026-10-22, all recorded 2026-10-23

Chosen by information content, not staging convenience — five of the eleven had to be built
(`<leg>_r2` staged from the parent's `system/`, `protein_fixed.pdb`, `ligand_prepped.sdf`,
`<leg>_ligand.sdf`). Four legfiles `runs/folds_md/rep3_{oxb_1,oxb_2,oct_1,oct_2}.legs` (BATCH=1,
one job each; `rep3_oxb_1` ran first). The rule each leg satisfies:

1. **The design-vs-null margins, whose weak side had no replicate.** The ESM2 pair now has both
   sides (`ox2_esm1_f8` −20.81 and `ox2_shuffle_esm0` −12.04 — the null whose 94.7–97.4% holding
   broke `sim`); the glycine pair its design side (`ox2_orig_f8` −24.21, the strongest oxybenzone
   cell); the octinoxate design-vs-null pair both sides (`s3_orig_f12` −24.33 and
   `shuffle_control` −15.13).
2. **Bolded column maxima with no replicate.** `orig_f12_rf3cof` (−19.51, RF3-cofold column max),
   `ox2_esm1_f8_o3pcof_dock6` (−18.74, pocket-guided column max on oxybenzone).
3. **A close comparison between lanes, 3.3 apart** — inside the 4.80 median spread, so a
   published ranking that a replicate test can actually overturn: both sides,
   `s3_esm2_f4_rf3cof_dock1` (−20.39) and `s3_esm2_f4_o3pcof_dock1` (−17.12). The latter is also
   only the pocket-guided lane's third replicate across 26 legs.
4. **A weak-end extreme.** `ox2_shuffle_esm0_o3cof` (−0.92, 1.1% holding) — the null-collapse cell;
   its design side already has a replicate (+0.31). If a second run moves it to −5 or −8, the
   "null collapses on OpenFold3" story is seed noise.
5. **Shell coverage.** `ox3_orig_f8` (−18.84) — the ox3 shell had zero replicates among the 23.

Results, each against `mmgbsa_summary.csv`'s whole-run value:

| structure | rule | run 1 | run 2 | signed ΔΔG | verdict |
|---|---|---:|---:|---:|---|
| `s3_esm2_f4_o3pcof_dock1` | 3 | −17.12 | −17.60 | −0.48 | unchanged |
| `s3_esm2_f4_rf3cof_dock1` | 3 | −20.39 | −19.09 | +1.30 | unchanged |
| `s3_orig_f12` | 1 | −24.33 | −17.54 | +6.79 | weaker (still bound: 89.6% residence, 0 releases) |
| `orig_f12_rf3cof` | 2 | −19.51 | −21.45 | −1.95 | stronger |
| `shuffle_control` | 1 | −15.13 | −9.30 | +5.84 | weaker (10.4% release in 21 episodes) |
| `ox2_shuffle_esm0_o3cof` | 4 | −0.92 | −4.03 | −3.11 | stronger (null-collapse cell; still the weakest leg) |
| `ox2_esm1_f8_o3pcof_dock6` | 2 | −18.74 | −12.34 | +6.40 | weaker |
| `ox3_orig_f8` | 5 | −18.84 | −10.19 | +8.65 | weaker |
| `ox2_esm1_f8` | 1 | −20.81 | −19.32 | +1.49 | unchanged |
| `ox2_shuffle_esm0` | 1 | −12.04 | −13.80 | −1.76 | stronger |
| `ox2_orig_f8` | 1 | −24.21 | −21.24 | +2.97 | weaker |

Batch-3 takeaways: the ESM2 design-vs-null margin **survives** (8.77 → 5.52, both sides now
n=2); the glycine margin survives (18.9 → 16.8); the octinoxate `s3_esm2_f4` RF3>pocket
direction replicates; `s3_orig_f12` deflates 6.79 but stays bound. The null-collapse cell
(`ox2_shuffle_esm0_o3cof`) moved −0.92 → −4.03 — four times stronger but still at the weak end.

**Batch 4 is done** (2026-10-23, 10/10 legs recorded): one large deflation (`ox1_orig_f12_o3pcof`
+8.72 — the bolded pocket-guided cell), one moderate (`s3_orig_f12_dock1` +4.79), and the rest
quiet (−1.80, −1.21, −0.80, +0.06), consistent with the batch's design — its pairs were the missing
sides of high-information parents.

**Batches 5 and 6 are queued but BLOCKED** (2026-10-23): batch 5, the tiered sweep of the
remainder — 41 legs in 13 batched legfiles (`rep5_{oxy_1..6,oct_1..7}.legs`, per-leg tier
rationale in each file's header) — and batch 6, the window top-up, 16 legs in 4 legfiles
(`rep6_{oct_1..3,oxy}.legs`). All 17 batches were staged and are submitted through
`racc_feed.sh` BATCH=1, driven in order by a detached `runs/folds_md/feed_rest.sh` that keeps the
queue at 3 jobs (1 running + 2 pending — standing policy). `code/racc_chain.sh` and
`code/racc_topup.sh` were **removed** (user instruction 2026-10-23) after the first chain
submission burned all 17 batches as
empty no-op jobs in 40 s (no LEGS/RUNREL injected, no env/cwd in the generated slurm, hook-masked
exit codes), left the queue empty ~2 h, and a competing 5-day-limit job took the freed GPU head
position — which can idle the whole campaign until it or another H100 finishes. `racc_feed.sh`
never had any of those bugs; filling the queue belongs to it alone.

Tally state at close of batch 4 (2026-10-23 22:00): **41/101 pairs complete — stronger 12, weaker
16, unchanged 13; median |ΔΔG| 2.89, range 0.02–16.16.** The tally's PAIRS lists all 85+16
parents; pending entries are batches 5–6. Refresh with `python3 code/replicate_tally.py`.

## What the set shows so far

Read as the state of a partial campaign (41 pairs complete of 101 registered). **The user has not
drawn conclusions from it.**

- Directions stay near-balanced (stronger 12, weaker 16, unchanged 13), so **run 1 carries no
  systematic optimism** — the result is spread, not bias.
- Median |ΔΔG| 2.89, max |ΔΔG| 16.16, range 0.02–16.16. (Distribution figures as of batch 4's
  close; refresh with the tally before quoting them.)
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
   Scale of the edit, measured 2026-10-22: 100 table rows mention one of the 23 structures, 44 carry
   a free energy, and **34 print the structure's own run-1 value** -- those 34 are the cells a rule
   would rewrite.
2. **Scope.** Stop at 10%, or run the remaining ~202 structures.
3. Whether the trailing-window diagnostic ever enters a table — see `TRAILING_WINDOWS.md`.

## Reproducing the table

```sh
python3 code/replicate_tally.py          # the live tally; PAIRS + THIRD define the set
```

`replicate_tally.py` is the source of truth for membership and arithmetic; this document is a
snapshot of it taken when the campaign was parked.
