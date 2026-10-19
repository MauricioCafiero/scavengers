**Table 1. Octinoxate: MM/GBSA binding free energy, ΔG (kcal/mol), and ligand residence (%) from 20 ns MD, started from each of three poses in each of three models' UNFORCED co-folds, plus the pocket-guided OpenFold3 lane.**

| peptide | pose | Boltz / BoltzGen co-fold | OpenFold3 co-fold | RoseTTAFold3 co-fold | OpenFold3 pocket-guided |
|---|---|---|---|---|---|
| `orig` | co-fold | −7.24 (29.9%) • | −5.90 (7.4%) ↑ | **−19.51** (99.9%) ↑ | −8.06 (21.1%) ↓ |
|  | Vina p1 | −2.81 (0.6%) ↑ | −9.65 (50.2%) ↑ | **−18.71** (89.2%) ↑ | −12.27 (51.0%) • |
|  | GNINA pick (p1/p1/p5/p4) | see Vina p1 | see Vina p1 | −11.10 (71.0%) ↑ | *running* |
| `s2_esm2_control` | co-fold | −16.25 (59.6%) • | −10.21 (33.9%) ↓ | −14.98 (57.2%) ↓ | −16.28 (100.0%) • |
|  | Vina p1 | −15.49 (0.0%) ↓ | −17.27 (95.2%) ↓ | −15.86 (78.9%) • | *running* |
|  | GNINA pick (p6/p5/p4/p8) | −17.09 (50.3%) • | −6.84 (9.4%) ↓ | −15.41 (73.7%) ↓ | *running* |
| `s3_esm2` | co-fold | −15.70 (29.4%) ↑ | −13.76 (84.8%) ↑ | −16.18 (88.2%) ↑ | *running* |
|  | Vina p1 | −9.65 (2.4%) ↓ | −14.53 (12.6%) • | **−20.39** (99.4%) • | *running* |
|  | GNINA pick (p7/p8/p3/p6) | −16.64 (0.1%) • | −12.79 (6.9%) ↑ | **−20.79** (99.3%) ↑ | *running* |
| `s3_orig` | co-fold | −13.37 (48.0%) ↑ | −10.10 (44.8%) ↑ | −12.10 (68.9%) ↑ | −9.78 (43.0%) ↑ |
|  | Vina p1 | −9.41 (36.8%) • | −8.92 (25.2%) • | **−19.16** (90.0%) ↓ | −11.50 (58.0%) ↓ |
|  | GNINA pick (p9/p6/p5/p9) | −11.51 (37.5%) ↑ | −11.17 (22.9%) ↑ | −8.43 (7.2%) ↑ | *running* |
| `shuffle_control` | co-fold | −15.13 (75.9%) ↑ | −9.39 (25.7%) ↓ | −11.13 (38.7%) • | — |
|  | Vina p1 | −16.15 (57.8%) | −10.89 (23.3%) • | −5.37 (21.8%) ↑ | — |
|  | GNINA pick (p7/p8/p7/—) | −8.40 (22.4%) ↑ | −7.78 (12.6%) ↑ | −13.86 (60.3%) ↓ | — |
| `shuffle_control_esm0` | co-fold | −14.32 (80.2%) ↓ | −15.62 (92.5%) ↓ | **−21.53** (99.3%) • | — |
|  | Vina p1 | −15.50 (76.1%) ↑ | −12.27 (1.9%) ↓ | −12.47 (44.1%) • | — |
|  | GNINA pick (p9/p3/p3/—) | −14.76 (90.5%) • | **−18.19** (98.5%) • | −12.03 (36.9%) ↑ | — |
| `bg33_1` | co-fold | −11.31 (14.5%) ↓ | −14.90 (56.4%) ↓ | −12.43 (68.3%) ↓ | — |
|  | Vina p1 | −13.47 (44.7%) ↓ | −13.16 (22.6%) • | −14.02 (60.3%) ↓ | — |
|  | GNINA pick (p2/p6/p8/—) | −12.55 (42.5%) ↓ | −12.93 (44.5%) ↑ | −11.78 (62.6%) • | — |
| `bg33_2` | co-fold | −8.02 (50.9%) ↑ | −11.58 (70.2%) • | −11.52 (58.7%) ↓ | — |
|  | Vina p1 | −11.82 (43.1%) ↑ | −11.55 (14.8%) • | −10.95 (81.5%) • | — |
|  | GNINA pick (p4/p1/p8/—) | −10.44 (57.2%) ↓ | see Vina p1 | −12.19 (36.7%) • | — |
| `bg33_3` | co-fold | −11.86 (41.0%) ↓ | −10.10 (14.9%) • | −12.44 (1.9%) • | — |
|  | Vina p1 | −10.84 (25.7%) ↑ | −14.89 (99.7%) • | −11.01 (0.2%) ↑ | — |
|  | GNINA pick (p6/p6/p5/—) | −13.09 (77.6%) ↑ | −14.92 (99.6%) • | −11.47 (0.1%) ↑ | — |
| `bg33_4` | co-fold | **−19.66** (100.0%) ↓ | −17.11 (93.0%) ↓ | −16.13 (81.2%) ↑ | — |
|  | Vina p1 | −14.81 (59.3%) ↑ | −16.04 (99.2%) • | **−18.59** (97.6%) ↓ | — |
|  | GNINA pick (p6/p7/p4/—) | −16.53 (97.2%) • | −16.24 (82.8%) ↑ | **−19.06** (99.6%) • | — |

Values are ΔG in kcal/mol with residence in parentheses — the share of the **whole 20 ns run** with the ligand within 10 Å of the peptide. Each peptide carries three sub-rows: the **co-fold** pose each model produced, **Vina p1** (AutoDock rank 1 docked into that same fold), and the **GNINA pick** (max `CNNaffinity` over poses ≥1). The pick's pose in each lane is in the row label, in column order; where it is **p1 the pick coincides with Vina's rank 1**, so no separate leg was ever run and the cell says *see Vina p1* rather than repeating the number. The trailing marker compares the 20 ns ΔG with the ΔG over the first 10 ns: • stable (within 0.5 kcal/mol), ↓ still decreasing (binding strengthening), ↑ still increasing (binding weakening) — a convergence flag on ΔG, not on residence. Bold marks ΔG ≤ −18. Every cell is a finished 20 ns leg; the campaign completed 2026-10-19.

**Residence is a centre-of-mass measure, so read it with the fold's shape in mind.** It is the share of frames with the ligand's centroid within 10 Å of the peptide's, which penalises an elongated fold: a ligand gripped at one end of a rod is never released and still scores low, because the far half of the peptide drags the centroid away. `s3_esm2` co-fold is the case in point — 29.4% residence with **0.15%** of frames detached, i.e. essentially permanent atomic contact. Where retention is the question, the release statistics in `md_contacts.csv` (episode count, longest episode, late share) are the shape-independent measure.

**Every column is now an UNFORCED fold. Corrected 2026-10-18, and this is the point of the table.** The version published until today put hint-forced Boltz legs (`_f4`/`_f8`/`_f12`, that many `contact ... force: true` constraints) against OpenFold3 and RoseTTAFold3 columns that were run with no constraints at all, so a difference down a row mixed *which model* with *forced or not* — it was not a model comparison. The Boltz column now holds the hint-free fold of the same sequence (`<row>_control`), run through the identical MD and MM/GBSA protocol. Peptides that were already unforced are unchanged and their published numbers stand: the shuffle nulls, `s2_esm2_control`, and every BoltzGen `bg*` design, which that generator produces from length and ligand SMILES alone with no hints.

Row labels are the **sequence**, so the `_fN` suffix is dropped where a row was remapped: a row's Boltz cell is `<row>_control` and its OF3/RF3 cells are `<row>_fN_o3cof` / `<row>_fN_rf3cof`. The `_fN` survives in those leg names only because they were named after the forced Boltz structure whose sequence they refold — no constraint was applied to them. The forced Boltz numbers this column used to carry are not lost: they remain in RESULTS.md's headline tables, and Table 3 below sets forced against unforced directly on the static read.

Built by `code/cofold_tables.py` from `mmgbsa_summary.csv` (plus the BoltzGen legs' own `FINAL_RESULTS_MMPBSA.dat`), the per-leg residence tables in FOLDS.md and RESULTS.md, and `runs/folds_md/dock_legs.txt` for which pose each GNINA pick was. It prints any cell it cannot source.

**Every column is now an UNFORCED fold. Corrected 2026-10-18, and this is the point of the table.** The version published until today put hint-forced Boltz legs (`_f4`/`_f8`/`_f12`, that many `contact ... force: true` constraints) against OpenFold3 and RoseTTAFold3 columns that were run with no constraints at all, so a difference down a row mixed *which model* with *forced or not* — it was not a model comparison. The Boltz column now holds the hint-free fold of the same sequence (`<row>_control`), run through the identical MD and MM/GBSA protocol. Peptides that were already unforced are unchanged and their published numbers stand: the shuffle nulls, `s2_esm2_control`, and every BoltzGen `bg*` design, which that generator produces from length and ligand SMILES alone with no hints.

Row labels are the **sequence**, so the `_fN` suffix is dropped where a row was remapped: a row's Boltz cell is `<row>_control` and its OF3/RF3 cells are `<row>_fN_o3cof` / `<row>_fN_rf3cof`. The `_fN` survives in those leg names only because they were named after the forced Boltz structure whose sequence they refold — no constraint was applied to them. The forced Boltz numbers this column used to carry are not lost: they remain in RESULTS.md's headline tables, and Table 3 below sets forced against unforced directly on the static read.

Built by `code/cofold_tables.py` from `mmgbsa_summary.csv` (plus the BoltzGen legs' own `FINAL_RESULTS_MMPBSA.dat`) and the per-leg residence tables in FOLDS.md and RESULTS.md; it prints any cell it cannot source.

**Fourth column added 2026-10-19: the pocket-guided OpenFold3 lane.** The first three columns are unforced — no constraint of any kind — and a difference across them is a model difference. The fourth is OpenFold3 run with a `pocket_constraint` built from the peptide's **original designed shell**, not from the Boltz hint subsets: the shell's full non-glycine slot set, with `max_distance` taken as the maximum over that structure's `_f*` contact distances. It is therefore forced, and a difference between the OpenFold3 and pocket-guided columns is the effect of the restraint within one model. OpenFold3's constraint is a single residue set scored by its centroid, so there is no dose analogue to Boltz's 4/8/12 contacts — a subset would move the pocket rather than weaken it, which is why each peptide gets one guided fold rather than three.

The lane covers only the nine rows with a designed shell to take a pocket from. The nulls (`shuffle_control`, `shuffle_control_esm0`, `ox2_shuffle`, `ox2_shuffle_esm0`) get none: a shuffled sequence has no designed pocket, and two ways of manufacturing one were considered and declined. The BoltzGen `bg*` designs get none because the generator makes them without hints by construction. Those cells read `—`, meaning not attempted by design, not missing. `*running*` marks a leg that is built and queued but has not finished.

Every GNINA pick in this column was the same pose in all three rescoring passes but **none was corroborated by `minimizedAffinity`**, whose argmin is pose 1 in seven of nine. By the standard used elsewhere in this campaign these picks are unsupported, and the consistency is not corroboration.

**Restructured 2026-10-18** from one row per peptide to three, so a row now holds one peptide × one starting pose across the three folding models, and a column is one model across every pose. ΔG and the 10 ns markers are rebuilt from `runs/<mol>/md/mmgbsa_summary.csv` (and the BoltzGen legs' own `FINAL_RESULTS_MMPBSA.dat`, which no summary CSV holds); residence is collated from the per-leg tables in FOLDS.md and RESULTS.md. Every co-fold cell was cross-checked against the values published here before the restructure, with no mismatches. **That caveat has since been removed from these two tables**: the Boltz column was replaced with the unforced `_control` legs on 2026-10-18 and the last of them landed 2026-10-19, so all three columns are now unforced and a difference down a row is a model difference. Table 3 keeps the forced-versus-unforced static read.

**Table 2. Oxybenzone: MM/GBSA binding free energy, ΔG (kcal/mol), and ligand residence (%) from 20 ns MD, started from each of three poses in each of three models' UNFORCED co-folds, plus the pocket-guided OpenFold3 lane.**

| peptide | pose | Boltz / BoltzGen co-fold | OpenFold3 co-fold | RoseTTAFold3 co-fold | OpenFold3 pocket-guided |
|---|---|---|---|---|---|
| `ox1_esm1` | co-fold | **−26.09** (100.0%) ↓ | −17.95 (100.0%) ↑ | −12.12 (35.2%) ↓ | −10.87 (95.3%) • |
|  | Vina p1 | −10.23 (68.6%) ↑ | −10.60 (89.9%) ↑ | −7.26 (27.3%) ↑ | *running* |
|  | GNINA pick (p4/p2/p4/p8) | −10.11 (56.2%) ↑ | −12.79 (66.9%) • | −14.92 (79.8%) ↓ | *running* |
| `ox1_orig` | co-fold | −8.34 (63.8%) ↑ | −3.96 (20.5%) ↑ | −12.27 (89.1%) ↓ | **−18.26** (100.0%) • |
|  | Vina p1 | −9.34 (61.6%) • | −4.93 (25.5%) ↑ | −6.40 (37.5%) ↑ | *running* |
|  | GNINA pick (p6/p4/p1/p4) | −4.79 (17.0%) ↑ | −5.84 (33.0%) ↓ | see Vina p1 | *running* |
| `ox2_esm1` | co-fold | −12.70 (71.0%) ↑ | −16.68 (93.2%) ↑ | −10.90 (90.4%) ↑ | −9.21 (31.8%) • |
|  | Vina p1 | −4.80 (21.2%) ↑ | −16.87 (93.9%) ↓ | −11.47 (79.1%) ↓ | *running* |
|  | GNINA pick (p9/p3/p5/p6) | −11.09 (86.2%) ↓ | −10.58 (77.1%) ↑ | −9.27 (62.4%) ↓ | *running* |
| `ox2_orig` | co-fold | −8.51 (43.9%) ↑ | −7.23 (49.5%) ↑ | −4.95 (19.6%) ↑ | *running* |
|  | Vina p1 | −6.71 (44.0%) • | −9.05 (47.2%) ↑ | −4.22 (15.5%) • | *running* |
|  | GNINA pick (p9/p2/p7/p5) | −3.50 (17.0%) ↑ | −5.43 (37.1%) ↑ | −4.11 (8.2%) • | *running* |
| `ox2_shuffle` | co-fold | −5.30 (7.1%) • | −2.98 (14.8%) ↑ | −13.35 (51.0%) • | — |
|  | Vina p1 | −5.86 (26.2%) • | −5.91 (27.3%) ↓ | −8.60 (72.1%) ↓ | — |
|  | GNINA pick (p4/p1/p2/—) | −2.09 (1.2%) • | see Vina p1 | −9.29 (32.0%) ↓ | — |
| `ox2_shuffle_esm0` | co-fold | −12.04 (94.7%) • | −0.92 (1.1%) • | −11.44 (32.6%) • | — |
|  | Vina p1 | −12.78 (97.4%) • | −9.44 (60.4%) ↓ | −9.18 (54.8%) ↑ | — |
|  | GNINA pick (p9/p6/p1/—) | −11.87 (95.5%) • | −6.84 (56.3%) ↑ | see Vina p1 | — |
| `ox3_orig` | co-fold | −10.13 (72.0%) ↑ | −2.06 (2.4%) ↓ | −7.41 (45.6%) ↓ | *running* |
|  | Vina p1 | −14.33 (97.9%) ↓ | −4.12 (11.8%) ↑ | −4.32 (19.5%) • | *running* |
|  | GNINA pick (p6/p4/p8/p1) | −15.44 (100.0%) ↑ | −5.99 (33.2%) • | −4.90 (10.0%) ↑ | see Vina p1 |
| `bgox31_2` | co-fold | −10.89 (33.7%) ↑ | −7.45 (65.9%) ↑ | −9.06 (84.6%) ↓ | — |
|  | Vina p1 | −2.11 (8.8%) ↑ | −5.85 (10.0%) ↓ | −7.36 (47.9%) ↑ | — |
|  | GNINA pick (p6/p1/p3/—) | −8.50 (2.6%) • | see Vina p1 | −9.57 (52.9%) • | — |
| `bgox31_3` | co-fold | −15.44 (99.7%) ↑ | −5.82 (15.8%) ↑ | −14.36 (99.2%) • | — |
|  | Vina p1 | **−22.23** (100.0%) ↓ | −10.99 (68.5%) ↑ | −13.48 (95.7%) ↑ | — |
|  | GNINA pick (p3/p4/p1/—) | −17.52 (89.9%) ↑ | −12.48 (88.9%) • | see Vina p1 | — |
| `bgox31_4` | co-fold | −14.05 (100.0%) ↑ | −5.74 (17.9%) ↑ | −11.02 (72.0%) • | — |
|  | Vina p1 | −9.20 (36.5%) ↑ | −10.39 (49.9%) ↑ | −15.86 (99.9%) • | — |
|  | GNINA pick (p6/p7/p3/—) | −4.31 (4.9%) ↑ | −6.71 (11.3%) ↓ | −13.91 (100.0%) • | — |
| `bgox31_5` | co-fold | −16.78 (100.0%) • | −11.55 (70.7%) ↓ | −14.69 (96.4%) ↑ | — |
|  | Vina p1 | −12.81 (70.3%) • | −13.72 (93.4%) • | −14.93 (99.8%) • | — |
|  | GNINA pick (p2/p7/p8/—) | −14.28 (99.5%) • | −14.84 (99.5%) ↑ | −14.76 (94.9%) • | — |

Values are ΔG in kcal/mol with residence in parentheses — the share of the **whole 20 ns run** with the ligand within 10 Å of the peptide. Each peptide carries three sub-rows: the **co-fold** pose each model produced, **Vina p1** (AutoDock rank 1 docked into that same fold), and the **GNINA pick** (max `CNNaffinity` over poses ≥1). The pick's pose in each lane is in the row label, in column order; where it is **p1 the pick coincides with Vina's rank 1**, so no separate leg was ever run and the cell says *see Vina p1* rather than repeating the number. The trailing marker compares the 20 ns ΔG with the ΔG over the first 10 ns: • stable (within 0.5 kcal/mol), ↓ still decreasing (binding strengthening), ↑ still increasing (binding weakening) — a convergence flag on ΔG, not on residence. Bold marks ΔG ≤ −18. Every cell is a finished 20 ns leg; the campaign completed 2026-10-19.

**Residence is a centre-of-mass measure, so read it with the fold's shape in mind.** It is the share of frames with the ligand's centroid within 10 Å of the peptide's, which penalises an elongated fold: a ligand gripped at one end of a rod is never released and still scores low, because the far half of the peptide drags the centroid away. `s3_esm2` co-fold is the case in point — 29.4% residence with **0.15%** of frames detached, i.e. essentially permanent atomic contact. Where retention is the question, the release statistics in `md_contacts.csv` (episode count, longest episode, late share) are the shape-independent measure.

**Every column is now an UNFORCED fold. Corrected 2026-10-18, and this is the point of the table.** The version published until today put hint-forced Boltz legs (`_f4`/`_f8`/`_f12`, that many `contact ... force: true` constraints) against OpenFold3 and RoseTTAFold3 columns that were run with no constraints at all, so a difference down a row mixed *which model* with *forced or not* — it was not a model comparison. The Boltz column now holds the hint-free fold of the same sequence (`<row>_control`), run through the identical MD and MM/GBSA protocol. Peptides that were already unforced are unchanged and their published numbers stand: the shuffle nulls, `s2_esm2_control`, and every BoltzGen `bg*` design, which that generator produces from length and ligand SMILES alone with no hints.

Row labels are the **sequence**, so the `_fN` suffix is dropped where a row was remapped: a row's Boltz cell is `<row>_control` and its OF3/RF3 cells are `<row>_fN_o3cof` / `<row>_fN_rf3cof`. The `_fN` survives in those leg names only because they were named after the forced Boltz structure whose sequence they refold — no constraint was applied to them. The forced Boltz numbers this column used to carry are not lost: they remain in RESULTS.md's headline tables, and Table 3 below sets forced against unforced directly on the static read.

Built by `code/cofold_tables.py` from `mmgbsa_summary.csv` (plus the BoltzGen legs' own `FINAL_RESULTS_MMPBSA.dat`), the per-leg residence tables in FOLDS.md and RESULTS.md, and `runs/folds_md/dock_legs.txt` for which pose each GNINA pick was. It prints any cell it cannot source.

**Every column is now an UNFORCED fold. Corrected 2026-10-18, and this is the point of the table.** The version published until today put hint-forced Boltz legs (`_f4`/`_f8`/`_f12`, that many `contact ... force: true` constraints) against OpenFold3 and RoseTTAFold3 columns that were run with no constraints at all, so a difference down a row mixed *which model* with *forced or not* — it was not a model comparison. The Boltz column now holds the hint-free fold of the same sequence (`<row>_control`), run through the identical MD and MM/GBSA protocol. Peptides that were already unforced are unchanged and their published numbers stand: the shuffle nulls, `s2_esm2_control`, and every BoltzGen `bg*` design, which that generator produces from length and ligand SMILES alone with no hints.

Row labels are the **sequence**, so the `_fN` suffix is dropped where a row was remapped: a row's Boltz cell is `<row>_control` and its OF3/RF3 cells are `<row>_fN_o3cof` / `<row>_fN_rf3cof`. The `_fN` survives in those leg names only because they were named after the forced Boltz structure whose sequence they refold — no constraint was applied to them. The forced Boltz numbers this column used to carry are not lost: they remain in RESULTS.md's headline tables, and Table 3 below sets forced against unforced directly on the static read.

Built by `code/cofold_tables.py` from `mmgbsa_summary.csv` (plus the BoltzGen legs' own `FINAL_RESULTS_MMPBSA.dat`) and the per-leg residence tables in FOLDS.md and RESULTS.md; it prints any cell it cannot source.

**Fourth column added 2026-10-19: the pocket-guided OpenFold3 lane.** The first three columns are unforced — no constraint of any kind — and a difference across them is a model difference. The fourth is OpenFold3 run with a `pocket_constraint` built from the peptide's **original designed shell**, not from the Boltz hint subsets: the shell's full non-glycine slot set, with `max_distance` taken as the maximum over that structure's `_f*` contact distances. It is therefore forced, and a difference between the OpenFold3 and pocket-guided columns is the effect of the restraint within one model. OpenFold3's constraint is a single residue set scored by its centroid, so there is no dose analogue to Boltz's 4/8/12 contacts — a subset would move the pocket rather than weaken it, which is why each peptide gets one guided fold rather than three.

The lane covers only the nine rows with a designed shell to take a pocket from. The nulls (`shuffle_control`, `shuffle_control_esm0`, `ox2_shuffle`, `ox2_shuffle_esm0`) get none: a shuffled sequence has no designed pocket, and two ways of manufacturing one were considered and declined. The BoltzGen `bg*` designs get none because the generator makes them without hints by construction. Those cells read `—`, meaning not attempted by design, not missing. `*running*` marks a leg that is built and queued but has not finished.

Every GNINA pick in this column was the same pose in all three rescoring passes but **none was corroborated by `minimizedAffinity`**, whose argmin is pose 1 in seven of nine. By the standard used elsewhere in this campaign these picks are unsupported, and the consistency is not corroboration.

**Restructured 2026-10-18** from one row per peptide to three, so a row now holds one peptide × one starting pose across the three folding models, and a column is one model across every pose. ΔG and the 10 ns markers are rebuilt from `runs/<mol>/md/mmgbsa_summary.csv` (and the BoltzGen legs' own `FINAL_RESULTS_MMPBSA.dat`, which no summary CSV holds); residence is collated from the per-leg tables in FOLDS.md and RESULTS.md. Every co-fold cell was cross-checked against the values published here before the restructure, with no mismatches. **That caveat has since been removed from these two tables**: the Boltz column was replaced with the unforced `_control` legs on 2026-10-18 and the last of them landed 2026-10-19, so all three columns are now unforced and a difference down a row is a model difference. Table 3 keeps the forced-versus-unforced static read.

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

**Two things this does not yet say.** It was a starting-geometry comparison when written; the
23 unforced legs (`runs/folds_md/bzfree_*.legs`) all completed 2026-10-19, so the dynamics three-way —
unforced Boltz vs the OF3 and RF3 ΔG in Tables 1 and 2 — is now in those tables. And because two octinoxate unforced legs start unbound (20.9 and 13.6 Å), their ΔG
will be partly a bound-versus-unbound difference; oxybenzone, where all five unforced folds start in
contact at 5.8–9.0 Å, will be the cleaner arm. Forced OF3 (pocket-guided) and forced RF3 remain the
next campaign's business, and RF3 has no contact input at all — see NEXT_STEPS.md, "Forced vs
unforced co-folds".

*(`cofold_binding_tables.pdf` is rebuilt from this file by `code/render_cofold_pdf.py`; regenerated 2026-10-19 with the completed campaign. The markdown is the source of truth — re-run that script after editing a table.)*
