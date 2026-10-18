# Pose reproducibility at campaign scale

Written 2026-10-18, after the full docked-campaign record. Everything here is recomputable with
`code/pose_reproducibility.py` from `runs/<mol>/md/pose_reproducibility.csv` (one row per redocked
structure) and `pose_reproducibility_docked.csv` (one row per docked MD leg). Nothing existing was
overwritten; the script reads tables and, when a value is missing anywhere, recomputes residence
with the record flow's own `md_contacts.profile`.

## The question

The redocking analysis (`NEXT_STEPS.md`, 2026-09-30) found that across 8 structures, how closely the
closest redocked pose reproduces the predicted one correlates **ρ = −0.743 (p = 0.035)** with the
structure's ligand residence in 20 ns MD — a candidate free screen: dock the predicted pose back
into its own fold, and a poor reproduction flags a leg that will release, with no dynamics needed.
The docked campaign has since run ~85 docked-leg MDs across two molecules, so the screen can be
tested at 10–20× the original sample.

## Result: the headline does not survive

| set | n | repro-RMSD vs residence | vs episodes | vs ΔG |
|---|---|---|---|---|
| octinoxate (base legs) | 30 | **−0.453 (p = .012)** | +0.322 (n.s.) | +0.320 (p = .085) |
| — S-folds only | 17 | −0.581 (p = .015) | | |
| — R-folds only | 13 | −0.154 (n.s.) | | |
| oxybenzone (base legs) | 33 | −0.136 (n.s.) | +0.108 (n.s.) | +0.098 (n.s.) |
| — docked legs, pooled both molecules | 115 | −0.205 (p = .028) | | |

Reproducing the original on its own 8 structures gives exactly ρ = −0.743, p = 0.035 — the new
tables are consistent with the old analysis, and it is the larger sample, not a measurement change,
that moved the number.

Two things kill it:

- **Chirality carries part of the original effect.** Octinoxate folds mix enantiomers (17 S / 13 of
  the R placements listed per row as `ligand_cip`); the correlation lives on the S-folds
  (−0.581) and vanishes on the R-folds, for which the "predicted pose" was never the fold's own
  answer.
- **Oxybenzone's OF3/RF3 refold folds redock beautifully and still lose the ligand.** Six examples:
  `ox1_orig_f12_o3cof` 1.16 Å → 20.5% retention, `ox2_shuffle_o3cof` 1.62 → 14.8%,
  `bgox31_3_o3cof` 2.02 → 15.8%, `ox3_orig_f8_o3cof` 1.94 → 2.4%, `ox2_orig_f8_rf3cof` 2.79 →
  19.6%, `ox2_shuffle_esm0_o3cof` 2.58 → 1.1%. Of the 26 structures that redock below 3 Å only 62%
  hold ≥50% residence — barely different from the 3–5 Å bin (82%) and the >5 Å bin (50%). The
  static pocket imprint does not decide survival; the shell does (the same lesson the MM/GBSA
  matrix keeps giving).

Release episodes never track redocking quality; the suggestive +0.32 on octinoxate is inflated by
the extreme `bg33_4` (1.30 Å, 100%) / `bg33_3` (4.47 Å, 41%) pair that supplied both ends of the
original correlation.

## What does survive as a cheap signal: the score, not the geometry

In the docked-leg analysis, a pose's own Vina score predicts that pose's MD residence — pooled
ρ = −0.306 (p < .001), oxybenzone −0.434, octinoxate −0.192 — and its post-MD ΔG (ρ ≈ +0.4 on both
molecules). Thresholds: redocked-pose score ≤ −5.5 kcal/mol → 83% of legs hold ≥50% residence;
≥ −4.5 → 48%. In the first analysis the redock's best score likewise beats repro-RMSD as a
residence screen on octinoxate (−0.618 vs −0.453).

**Practical verdict.** Re-docking the predicted pose into its own fold is only a coarse, one-sided
screen, not a residence predictor: agreement with the predicted pose is not evidence of retention,
and disagreement flags risk mostly on base-design lanes. Any cheap pass on a new fold should read
the score, not the RMSD; deciding whether the ligand is still there after the run remains 20 ns of
MD, which is why it was made the default tail step.

## Caveats worth remembering

- Octinoxate cross-structure comparisons mix enantiomers by construction of the folds; oxybenzone
  is achiral and is the cleaner molecule for any new per-structure statistic.
- Residence values come with provenance (`residence_source`): parsed from the per-leg record line,
  from FOLDS.md, or recomputed from the trajectory with `md_contacts.profile` for the 17 base legs
  predating the record flow. The recomputed values match every published figure they were
  revalidated against (`bg33_1_o3cof` 56.4%, `bg33_3` 41.0%/18 episodes).
- Apo lanes (`o3apo`/`rf3apo`/`esmapo`/`bzapo`) are out of scope by construction: their folds have
  no predicted pose to redock against.