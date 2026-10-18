**Table 1. Octinoxate: MM/GBSA binding free energy, ΔG (kcal/mol), and ligand retention (%) from 20 ns MD started at each folding model's co-folded pose.**

| Peptide | Boltz / BoltzGen co-fold | OpenFold3 co-fold | RoseTTAFold3 co-fold |
|---|---|---|---|
| orig_f12 | −13.71 (77.3%) ↑ | −5.90 (7.4%) ↑ | −19.51 (99.9%) ↑ |
| s2_esm2_control | −16.25 (59.6%) • | −10.21 (33.9%) ↓ | −14.98 (57.2%) ↓ |
| s3_esm2_f4 | −21.08 (99.6%) • | −13.76 (84.8%) ↑ | −16.18 (88.2%) ↑ |
| s3_orig_f12 | −24.33 (100.0%) ↑ | −10.10 (44.8%) ↑ | −12.10 (68.9%) ↑ |
| shuffle_control | −15.13 (75.9%) ↑ | −9.39 (25.7%) ↓ | −11.13 (38.7%) • |
| shuffle_control_esm0 | −14.32 (80.2%) ↓ | −15.62 (92.5%) ↓ | −21.53 (99.3%) • |
| bg33_1 | −11.31 (14.5%) ↓ | −14.90 (56.4%) ↓ | −12.43 (68.3%) ↓ |
| bg33_2 | −8.02 (50.9%) ↑ | −11.58 (70.2%) • | −11.52 (58.7%) ↓ |
| bg33_3 | −11.86 (41.0%) ↓ | −10.10 (14.9%) • | −12.44 (1.9%) • |
| bg33_4 | −19.66 (100.0%) ↓ | −17.11 (93.0%) ↓ | −16.13 (81.2%) ↑ |

Values are ΔG in kcal/mol with retention in parentheses. ΔG is the MM/GBSA binding free energy over the whole 20 ns production run; retention is the percentage of the run with the ligand within 10 Å of the peptide. The trailing marker compares the 20 ns ΔG with the ΔG over the first 10 ns: • stable (within 0.5 kcal/mol), ↓ still decreasing (binding strengthening), ↑ still increasing (binding weakening). The Boltz / BoltzGen column is the original co-fold: BoltzGen for the bg* peptides, Boltz for the rest.

**Table 2. Oxybenzone: MM/GBSA binding free energy, ΔG (kcal/mol), and ligand retention (%) from 20 ns MD started at each folding model's co-folded pose.**

| Peptide | Boltz / BoltzGen co-fold | OpenFold3 co-fold | RoseTTAFold3 co-fold |
|---|---|---|---|
| ox1_esm1_f8 | −10.70 (75.0%) ↑ | −17.95 (100.0%) ↑ | −12.12 (35.2%) ↓ |
| ox1_orig_f12 | −10.31 (49.1%) ↑ | −3.96 (20.5%) ↑ | −12.27 (89.1%) ↓ |
| ox2_esm1_f8 | −20.81 (100.0%) ↑ | −16.68 (93.2%) ↑ | −10.90 (90.4%) ↑ |
| ox2_orig_f8 | −24.21 (100.0%) • | −7.23 (49.5%) ↑ | −4.95 (19.6%) ↑ |
| ox2_shuffle | −5.30 (7.1%) • | −2.98 (14.8%) ↑ | −13.35 (51.0%) • |
| ox2_shuffle_esm0 | −12.04 (94.7%) • | −0.92 (1.1%) • | −11.44 (32.6%) • |
| ox3_orig_f8 | −18.84 (100.0%) ↑ | −2.06 (2.4%) ↓ | −7.41 (45.6%) ↓ |
| bgox31_2 | −10.89 (33.7%) ↑ | −7.45 (65.9%) ↑ | −9.06 (84.6%) ↓ |
| bgox31_3 | −15.44 (99.7%) ↑ | −5.82 (15.8%) ↑ | −14.36 (99.2%) • |
| bgox31_4 | −14.05 (100.0%) ↑ | −5.74 (17.9%) ↑ | −11.02 (72.0%) • |
| bgox31_5 | −16.78 (100.0%) • | −11.55 (70.7%) ↓ | −14.69 (96.4%) ↑ |

Values are ΔG in kcal/mol with retention in parentheses. ΔG is the MM/GBSA binding free energy over the whole 20 ns production run; retention is the percentage of the run with the ligand within 10 Å of the peptide. The trailing marker compares the 20 ns ΔG with the ΔG over the first 10 ns: • stable (within 0.5 kcal/mol), ↓ still decreasing (binding strengthening), ↑ still increasing (binding weakening). The Boltz / BoltzGen column is the original co-fold: BoltzGen for the bg* peptides, Boltz for the rest.

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
