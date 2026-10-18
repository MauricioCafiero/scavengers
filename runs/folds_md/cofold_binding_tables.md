**Table 1. Octinoxate: MM/GBSA binding free energy, ΔG (kcal/mol), and ligand retention (%) from 20 ns MD started at each folding model's co-folded pose.**

| peptide | pose | Boltz / BoltzGen co-fold | OpenFold3 co-fold | RoseTTAFold3 co-fold |
|---|---|---|---|---|
| `orig_f12` | co-fold | −13.71 (77.3%) ↑ | −5.90 (7.4%) ↑ | **−19.51** (99.9%) ↑ |
|  | Vina p1 | **−22.11** (100.0%) ↑ | −9.65 (50.2%) ↑ | **−18.71** (89.2%) ↑ |
|  | GNINA pick (p5/—/5) | −17.17 (83.3%) ↑ | — | −11.10 (71.0%) ↑ |
| `s2_esm2_control` | co-fold | −16.25 (59.6%) • | −10.21 (33.9%) ↓ | −14.98 (57.2%) ↓ |
|  | Vina p1 | −15.49 (0.0%) ↓ | −17.27 (95.2%) ↓ | −15.86 (78.9%) • |
|  | GNINA pick (p6/5/4) | −17.09 (50.3%) • | −6.84 (9.4%) ↓ | −15.41 (73.7%) ↓ |
| `s3_esm2_f4` | co-fold | **−21.08** (99.6%) • | −13.76 (84.8%) ↑ | −16.18 (88.2%) ↑ |
|  | Vina p1 | **−22.42** (99.9%) • | −14.53 (12.6%) • | **−20.39** (99.4%) • |
|  | GNINA pick (p6/8/3) | **−19.55** (59.6%) ↓ | −12.79 (6.9%) ↑ | **−20.79** (99.3%) ↑ |
| `s3_orig_f12` | co-fold | **−24.33** (100.0%) ↑ | −10.10 (44.8%) ↑ | −12.10 (68.9%) ↑ |
|  | Vina p1 | **−20.79** (98.7%) ↑ | −8.92 (25.2%) • | **−19.16** (90.0%) ↓ |
|  | GNINA pick (p9/6/5) | −13.23 (70.5%) ↑ | −11.17 (22.9%) ↑ | −8.43 (7.2%) ↑ |
| `shuffle_control` | co-fold | −15.13 (75.9%) ↑ | −9.39 (25.7%) ↓ | −11.13 (38.7%) • |
|  | Vina p1 | −16.15 (57.8%) | −10.89 (23.3%) • | −5.37 (21.8%) ↑ |
|  | GNINA pick (p7/8/7) | −8.40 (22.4%) ↑ | −7.78 (12.6%) ↑ | −13.86 (60.3%) ↓ |
| `shuffle_control_esm0` | co-fold | −14.32 (80.2%) ↓ | −15.62 (92.5%) ↓ | **−21.53** (99.3%) • |
|  | Vina p1 | −15.50 (76.1%) ↑ | −12.27 (1.9%) ↓ | −12.47 (44.1%) • |
|  | GNINA pick (p9/3/3) | −14.76 (90.5%) • | **−18.19** (98.5%) • | −12.03 (36.9%) ↑ |
| `bg33_1` | co-fold | −11.31 (14.5%) ↓ | −14.90 (56.4%) ↓ | −12.43 (68.3%) ↓ |
|  | Vina p1 | −13.47 (44.7%) ↓ | −13.16 (22.6%) • | −14.02 (60.3%) ↓ |
|  | GNINA pick (p2/6/8) | −12.55 (42.5%) ↓ | −12.93 (44.5%) ↑ | −11.78 (62.6%) • |
| `bg33_2` | co-fold | −8.02 (50.9%) ↑ | −11.58 (70.2%) • | −11.52 (58.7%) ↓ |
|  | Vina p1 | −11.82 (43.1%) ↑ | −11.55 (14.8%) • | −10.95 (81.5%) • |
|  | GNINA pick (p4/—/8) | −10.44 (57.2%) ↓ | — | −12.19 (36.7%) • |
| `bg33_3` | co-fold | −11.86 (41.0%) ↓ | −10.10 (14.9%) • | −12.44 (1.9%) • |
|  | Vina p1 | −10.84 (25.7%) ↑ | −14.89 (99.7%) • | −11.01 (0.2%) ↑ |
|  | GNINA pick (p6/6/5) | −13.09 (77.6%) ↑ | −14.92 (99.6%) • | −11.47 (0.1%) ↑ |
| `bg33_4` | co-fold | **−19.66** (100.0%) ↓ | −17.11 (93.0%) ↓ | −16.13 (81.2%) ↑ |
|  | Vina p1 | −14.81 (59.3%) ↑ | −16.04 (99.2%) • | **−18.59** (97.6%) ↓ |
|  | GNINA pick (p6/7/4) | −16.53 (97.2%) • | −16.24 (82.8%) ↑ | **−19.06** (99.6%) • |

Values are ΔG in kcal/mol with residence in parentheses (share of the 20 ns run with the ligand within 10 Å of the peptide). Each peptide carries three sub-rows: the **co-fold** pose each model produced, **Vina p1** (AutoDock rank 1 docked into that same fold), and the **GNINA pick** (max `CNNaffinity` over poses ≥1; its pose number per lane is given in the row label, and `—` means the pick was pose 1, so the Vina-p1 row is that leg). The trailing marker compares the 20 ns ΔG with the first 10 ns: • stable (within 0.5 kcal/mol), ↓ still decreasing (binding strengthening), ↑ still increasing (binding weakening). Bold marks ΔG ≤ −18. The Boltz / BoltzGen column is the original co-fold: BoltzGen for the bg* peptides, Boltz for the rest.

**Restructured 2026-10-18** from one row per peptide to three, so a row now holds one peptide × one starting pose across the three folding models, and a column is one model across every pose. ΔG and the 10 ns markers are rebuilt from `runs/<mol>/md/mmgbsa_summary.csv` (and the BoltzGen legs' own `FINAL_RESULTS_MMPBSA.dat`, which no summary CSV holds); residence is collated from the per-leg tables in FOLDS.md and RESULTS.md. Every co-fold cell was cross-checked against the values published here before the restructure, with no mismatches. **The caveat that motivated the campaign still applies to these two tables**: the Boltz column is mostly hint-forced and the OF3 and RF3 columns are not, so read Table 3 for the like-for-like comparison.

**Table 2. Oxybenzone: MM/GBSA binding free energy, ΔG (kcal/mol), and ligand retention (%) from 20 ns MD started at each folding model's co-folded pose.**

| peptide | pose | Boltz / BoltzGen co-fold | OpenFold3 co-fold | RoseTTAFold3 co-fold |
|---|---|---|---|---|
| `ox1_esm1_f8` | co-fold | −10.70 (75.0%) ↑ | −17.95 (100.0%) ↑ | −12.12 (35.2%) ↓ |
|  | Vina p1 | −16.87 (98.5%) ↑ | −10.60 (89.9%) ↑ | −7.26 (27.3%) ↑ |
|  | GNINA pick (p—/2/4) | — | −12.79 (66.9%) • | −14.92 (79.8%) ↓ |
| `ox1_orig_f12` | co-fold | −10.31 (49.1%) ↑ | −3.96 (20.5%) ↑ | −12.27 (89.1%) ↓ |
|  | Vina p1 | −8.80 (51.0%) ↑ | −4.93 (25.5%) ↑ | −6.40 (37.5%) ↑ |
|  | GNINA pick (p4/4/—) | −10.41 (77.2%) ↑ | −5.84 (33.0%) ↓ | — |
| `ox2_esm1_f8` | co-fold | **−20.81** (100.0%) ↑ | −16.68 (93.2%) ↑ | −10.90 (90.4%) ↑ |
|  | Vina p1 | −15.22 (80.5%) ↓ | −16.87 (93.9%) ↓ | −11.47 (79.1%) ↓ |
|  | GNINA pick (p3/3/5) | −7.40 (33.9%) ↑ | −10.58 (77.1%) ↑ | −9.27 (62.4%) ↓ |
| `ox2_orig_f8` | co-fold | **−24.21** (100.0%) • | −7.23 (49.5%) ↑ | −4.95 (19.6%) ↑ |
|  | Vina p1 | −17.22 (99.9%) ↓ | −9.05 (47.2%) ↑ | −4.22 (15.5%) • |
|  | GNINA pick (p4/2/7) | **−24.00** (100.0%) ↓ | −5.43 (37.1%) ↑ | −4.11 (8.2%) • |
| `ox2_shuffle` | co-fold | −5.30 (7.1%) • | −2.98 (14.8%) ↑ | −13.35 (51.0%) • |
|  | Vina p1 | −5.86 (26.2%) • | −5.91 (27.3%) ↓ | −8.60 (72.1%) ↓ |
|  | GNINA pick (p4/—/2) | −2.09 (1.2%) • | — | −9.29 (32.0%) ↓ |
| `ox2_shuffle_esm0` | co-fold | −12.04 (94.7%) • | −0.92 (1.1%) • | −11.44 (32.6%) • |
|  | Vina p1 | −12.78 (97.4%) • | −9.44 (60.4%) ↓ | −9.18 (54.8%) ↑ |
|  | GNINA pick (p9/6/—) | −11.87 (95.5%) • | −6.84 (56.3%) ↑ | — |
| `ox3_orig_f8` | co-fold | **−18.84** (100.0%) ↑ | −2.06 (2.4%) ↓ | −7.41 (45.6%) ↓ |
|  | Vina p1 | **−19.98** (94.2%) ↑ | −4.12 (11.8%) ↑ | −4.32 (19.5%) • |
|  | GNINA pick (p—/4/8) | — | −5.99 (33.2%) • | −4.90 (10.0%) ↑ |
| `bgox31_2` | co-fold | −10.89 (33.7%) ↑ | −7.45 (65.9%) ↑ | −9.06 (84.6%) ↓ |
|  | Vina p1 | −2.11 (8.8%) ↑ | −5.85 (10.0%) ↓ | −7.36 (47.9%) ↑ |
|  | GNINA pick (p6/—/3) | −8.50 (2.6%) • | — | −9.57 (52.9%) • |
| `bgox31_3` | co-fold | −15.44 (99.7%) ↑ | −5.82 (15.8%) ↑ | −14.36 (99.2%) • |
|  | Vina p1 | **−22.23** (100.0%) ↓ | −10.99 (68.5%) ↑ | −13.48 (95.7%) ↑ |
|  | GNINA pick (p3/4/—) | −17.52 (89.9%) ↑ | −12.48 (88.9%) • | — |
| `bgox31_4` | co-fold | −14.05 (100.0%) ↑ | −5.74 (17.9%) ↑ | −11.02 (72.0%) • |
|  | Vina p1 | −9.20 (36.5%) ↑ | −10.39 (49.9%) ↑ | −15.86 (99.9%) • |
|  | GNINA pick (p6/7/3) | −4.31 (4.9%) ↑ | −6.71 (11.3%) ↓ | −13.91 (100.0%) • |
| `bgox31_5` | co-fold | −16.78 (100.0%) • | −11.55 (70.7%) ↓ | −14.69 (96.4%) ↑ |
|  | Vina p1 | −12.81 (70.3%) • | −13.72 (93.4%) • | −14.93 (99.8%) • |
|  | GNINA pick (p2/7/8) | −14.28 (99.5%) • | −14.84 (99.5%) ↑ | −14.76 (94.9%) • |

Values are ΔG in kcal/mol with residence in parentheses (share of the 20 ns run with the ligand within 10 Å of the peptide). Each peptide carries three sub-rows: the **co-fold** pose each model produced, **Vina p1** (AutoDock rank 1 docked into that same fold), and the **GNINA pick** (max `CNNaffinity` over poses ≥1; its pose number per lane is given in the row label, and `—` means the pick was pose 1, so the Vina-p1 row is that leg). The trailing marker compares the 20 ns ΔG with the first 10 ns: • stable (within 0.5 kcal/mol), ↓ still decreasing (binding strengthening), ↑ still increasing (binding weakening). Bold marks ΔG ≤ −18. The Boltz / BoltzGen column is the original co-fold: BoltzGen for the bg* peptides, Boltz for the rest.

**Restructured 2026-10-18** from one row per peptide to three, so a row now holds one peptide × one starting pose across the three folding models, and a column is one model across every pose. ΔG and the 10 ns markers are rebuilt from `runs/<mol>/md/mmgbsa_summary.csv` (and the BoltzGen legs' own `FINAL_RESULTS_MMPBSA.dat`, which no summary CSV holds); residence is collated from the per-leg tables in FOLDS.md and RESULTS.md. Every co-fold cell was cross-checked against the values published here before the restructure, with no mismatches. **The caveat that motivated the campaign still applies to these two tables**: the Boltz column is mostly hint-forced and the OF3 and RF3 columns are not, so read Table 3 for the like-for-like comparison.

---

**Table 3. The strict comparison, static half: one sequence, four co-folding treatments.** Added
2026-10-18. Tables 1 and 2 set a mostly **hint-forced** Boltz column against **unforced** OpenFold3
and RoseTTAFold3 columns, so a difference down a row mixes "which model" with "forced or not". The
hint-free Boltz twin of each designed sequence (its `_control` fold) is now in dynamics; its
starting-pose geometry is already measurable, and so is the like-for-like comparison. Enclosed
fraction and ligand-centroid separation are read off each co-folded pose's own coordinates
(`runs/<mol>/dock/dock_poses.csv`, pose 0) with the same `geometry()` as `check_fold.py` — no
energy, no dynamics.

| peptide | Boltz forced | Boltz unforced | OF3 unforced | RF3 unforced |
|---|---|---|---|---|
| orig_f12 (oct) | 0.960 / 4.2 Å | 0.495 / **20.9 Å** | 0.680 / 9.4 Å | **0.800** / 5.8 Å |
| s3_esm2_f4 (oct) | 0.795 / 7.0 Å | 0.360 / 13.6 Å | **0.580** / 7.3 Å | 0.460 / 9.3 Å |
| s3_orig_f12 (oct) | 0.965 / 3.6 Å | **0.735** / 9.5 Å | 0.430 / 10.4 Å | 0.720 / 8.5 Å |
| ox1_esm1_f8 | 0.940 / 4.1 Å | **0.765** / 5.9 Å | 0.580 / 6.1 Å | 0.530 / 9.3 Å |
| ox1_orig_f12 | 0.995 / 1.2 Å | 0.520 / 7.2 Å | **0.810** / 5.7 Å | 0.550 / 6.1 Å |
| ox2_esm1_f8 | 0.665 / 6.4 Å | 0.610 / 7.6 Å | **0.745** / 5.3 Å | 0.570 / 6.4 Å |
| ox2_orig_f8 | 0.975 / 2.7 Å | 0.535 / 9.0 Å | **0.655** / 7.5 Å | 0.595 / 6.2 Å |
| ox3_orig_f8 | 0.995 / 2.1 Å | **0.615** / 5.8 Å | 0.540 / 10.3 Å | 0.465 / 11.4 Å |
| **mean** | **0.911 / 3.9 Å** | 0.579 / 9.9 Å | 0.628 / 7.8 Å | 0.586 / 7.9 Å |

Each cell is enclosed fraction / ligand-centroid separation. **Bold** marks the best enclosure among
the three *unforced* lanes in that row (Boltz forced is excluded from that contest — it is the
treatment under test, not a competitor).

**What the strict comparison says, and it is not what Tables 1 and 2 imply.** The three unforced
lanes are indistinguishable on this measure: mean enclosure 0.579 (Boltz), 0.628 (OF3), 0.586 (RF3),
with per-peptide wins splitting 3 / 4 / 1 across eight sequences and every lane's range overlapping
the others' almost completely. **Forced Boltz is in a different regime from all three** — 0.911
enclosed at 3.9 Å, against 0.58–0.63 at 7.8–9.9 Å — and no unforced lane's best cell (0.810) reaches
forced Boltz's mean. So the Boltz column's advantage in Tables 1 and 2 is attributable to the
**contact hints**, not to Boltz being the better co-folder: strip the hints and it places the ligand
no better than the other two, and on octinoxate's `orig_f12` it is the worst of the three by a wide
margin (0.495 at 20.9 Å, the furthest starting pose in the whole campaign).

**Two things this does not yet say.** It is a starting-geometry comparison, not a binding one: the
unforced Boltz legs' MM/GBSA is still running (23 legs on racc, `runs/folds_md/bzfree_*.legs`), and
the dynamics three-way — unforced Boltz vs the OF3 and RF3 ΔG already in Tables 1 and 2 — completes
when they land. And because two octinoxate unforced legs start unbound (20.9 and 13.6 Å), their ΔG
will be partly a bound-versus-unbound difference; oxybenzone, where all five unforced folds start in
contact at 5.8–9.0 Å, will be the cleaner arm. Forced OF3 (pocket-guided) and forced RF3 remain the
next campaign's business, and RF3 has no contact input at all — see NEXT_STEPS.md, "Forced vs
unforced co-folds".

*(`cofold_binding_tables.pdf` predates this table and has not been regenerated.)*
