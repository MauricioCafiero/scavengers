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
