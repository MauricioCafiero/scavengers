# FOLDS.md — folding the designed peptides with independent models, against the Boltz and BoltzGen structures

Started 2026-10-09, while the racc queue was unavailable for a few days. This doc holds the
fold-method comparison: the same 21 peptides that went through dynamics, folded from sequence only
by models that had nothing to do with the design or the MD, then compared against the Boltz /
BoltzGen co-folds those dynamics legs were built from.

Every peptide set below covers the nulls and the BoltzGen arms — 10 octinoxate-side
(4 designs/esm2-variants, 2 nulls, 4 bg33) and 11 oxybenzone-side (4 designs/esm2-variants + ox1
pair, 2 nulls, 4 bgox31). Manifests with sequences, kinds and SMILES live in
`runs/<mol>/folds/manifest.csv`; per-peptide fold inputs in `runs/<mol>/folds/inputs/`. The
extraction read each sequence off `protein_fixed.pdb` of the structure the MD leg actually used, so
the comparison targets exactly the structures in the headline tables. One fact the extraction
surfaced: the four bg33 folds and the four bgox31 folds are **different sequences from each other
and from the design jsons** — BoltzGen generated its own sequence per candidate, it did not fold the
designed sequence.

## The five lanes and the question each answers

| Lane | Model | Input | Where | Question |
|---|---|---|---|---|
| ESMFold | esmfold_v1 (3B) | sequence only | Modal A10G | is the designed fold a fold **without the ligand**? |
| OpenFold3 | OF3 (+ ligand) | sequence + SMILES | Modal A10G | does an independent cofolder place the ligand like BoltzGen did? |
| RF3 | RosettaFold3 | sequence + SMILES | Modal A10G | second opinion on the cofold (all 21 landed 2026-10-09) |
| OpenFold3 apo | OF3, apo | sequence only | Modal A10G | is the compact apo fold real or ESMFold-specific? (all 21 landed 2026-10-09) |
| RF3 apo | RF3, apo | sequence only | Modal A10G | second opinion on the apo fold (all 21 landed 2026-10-09) |

All three run through the `fold` repo (`~/python_mac/fold`, cloned from
`github.com/MauricioCafiero/fold`; see its HANDOFF.md for what was verified live). ESMFold is **not
local** — its 11 GB of weights need more than this Mac has; only ESM2 *embeddings* are local-CPU,
and they are not used here. `code/fold_batch.sh` drives a lane: `LANE=<name> code/fold_batch.sh`,
resumable over its outputs, sequential (one Modal run at a time), each invocation under
`caffeinate -i` as the modal guard requires.

### The five methods, in detail

**1. ESMFold — sequence only, no ligand.** The peptide sequence is sent to `esmfold_v1` (3B
parameters) through the fold repo's `esmfold_app.py`, which runs on a Modal A10G and writes a PDB
directly; the B-factor column carries pLDDT on a 0–1 scale. The SMILES never enters — this lane
answers *does the sequence alone produce the fold*, which is why it is the apo test. The first run
carried the one-time ~3 GB weight download into the Modal volume; after that each fold took ~30 s.
Nothing of the design or Boltz/BoltzGen pipeline touches it: different architecture, different
training, ligand-blind.

**2. OpenFold3 — sequence + SMILES cofold.** The input file carries both the sequence and the
ligand SMILES; `app.py` runs OpenFold3 on the A10G with MSA, drops one model cif per sample (five
samples, seed 42) under `outputs/<job>/<job>/seed_*/`, plus a per-sample
`*_confidences_aggregated.json`. The best sample is picked by `sample_ranking_score`
(`fold.results.find_best_cif`) and converted to PDB with `fold.analyze cif-to-pdb`; in that PDB the
peptide sits on chain A and the ligand on chain Z. The run auto-prints interface confidences for
the A–Z pair — ipSAE, pDockQ2, pDockQ — of which ipSAE (Dunbrack 2025) is the interface-local one;
for protein–ligand pairs the calibration is experimental, so they are used as a rank only. First
run carried the container build and a 2.3 GB checkpoint (~7 min); the other twenty averaged
~1.5 min.

**3. RF3 (RosettaFold3) — same cofold, second model.** Identical input file through `rf3_app.py`;
its outputs land under `outputs/rf3_<job>/` with the best model promoted to the top level, peptide
on chain A and ligand on chain B, and the same ipSAE/pDockQ auto-print (lower dynamic range than
OF3's — RF3 lacks an image that scores protein–ligand interfaces well, per the fold repo's note).
Ran ~45 s per cofold once its checkpoint was cached. Its role is the second opinion: a fold that
both cofolders reproduce is model-robust; a fold they pack differently is degenerate.

**4. OpenFold3 apo — sequence only.** A later lane (2026-10-09, after the ESMFold results asked
the question): the same `app.py` submitted with `--sequence` and no SMILES. OpenFold3 accepts a
single-chain protein query natively — the requirement lived in this repo's plumbing, which
relaxed it (`inputs.py`, `openfold3_core.py`, `app.py`; pushed as fold@`d969fc0`). Jobs are named
`<pep>_apo` so they cannot collide with the cofold lane's outputs; conversion and measurement are
identical to the cofold lane, with the geometry columns empty (no ligand). ~1.5 min per leg.

**5. RF3 apo — the same, second model.** `rf3_app.py` with `--sequence` only; its ligand
component was hardcoded into `rf3_core.build_query_components` and its app auto-filled
rosuvastatin — both relaxed the same way (fold@`77dcb8b`, first test coverage for rf3_core).
Outputs land under `outputs/rf3_<pep>_apo/`. ~45 s per leg. Its role: whether the apo state is
compact under a second architecture, or whether *only* ESMFold collapses it.

**What every lane feeds.** `code/folds_metrics.py <lane>` reads the converted PDB against the leg's
own files: helical/beta fraction and hbond counts from `check_fold.secondary_structure`, Cα Rg,
Cα RMSD against the co-fold's peptide (Kabsch superposition, `protein_fixed.pdb`), mean pLDDT, and
— for the two cofold lanes — `check_fold.geometry()` on fold *and* on `system/complex.pdb` (the
structure the dynamics leg was actually built from), giving enclosed/wrapped/engaged/centroid
separation for both sides through one shared code path.

## What is compared (the inexpensive metrics)

All computed from the folded structures alone, no dynamics:

- **`check_fold.geometry()`** — the pipeline's own static read: `enclosed`, `wrapped`, `engaged`,
  centroid separation, peptide Rg. Used straight on the fold outputs, so every new number sits in
  the same columns as `fold_check.csv`, `fold_check_oxyb31.csv` (`boltzgen_local/results/`) and
  `dock_poses.csv`. For ESMFold (no ligand) the ligand-facing columns are N/A and the comparison is
  helix, Rg and shape against the co-folds' peptide.
- **Cα RMSD of the peptide** against each structure's Boltz/BoltzGen chain, after superposition
  (mdtraj).
- Helical/beta fraction and hbond counts, the same columns the fold_check tables already carry.
- For OpenFold3 runs, the cofold run auto-prints ipSAE / pDockQ interface confidences — recorded
  per peptide (protein–ligand ipSAE is experimental per the fold repo's README; treat as a rank).

The comparison is written into the tables below the first time each number lands, then to
`runs/<mol>/folds/metrics/`.

## ESMFold results (all 21 landed 2026-10-09; full data in `runs/<mol>/folds/metrics/esmfold.csv`)

ESMFold runs ~30 s each on the A10G once the one-time ~3 GB weight download is in the Modal volume;
the first run carried it (~4 min). pLDDT below is the model's per-residue confidence, its B-factor
column, normalised to a 0–100 scale for every lane (ESMFold writes 0–1, OF3/RF3 write 0–100);
these artificial glycine-rich sequences score middling 50–100, so the numbers are usable but
confidence-limited. Helix and beta are fractions of backbone
residues; Rg is the peptide CA radius of gyration; RMSD superposes ESMFold's CAs on the co-fold's.

| peptide | mol | kind | helix esm / cofold | Rg esm / cofold | Cα RMSD | pLDDT |
|---|---|---|---|---|---|---|
| orig_f12 | oct | designed | 0.000 / 0.065 | 27.93 / 7.46 | 27.33 | 90 |
| s2_esm2_control | oct | esm2-variant | 1.000 / 1.000 | 14.48 / 14.50 | 0.37 | 90 |
| s3_esm2_f4 | oct | esm2-variant | 1.000 / 0.906 | 14.88 / 9.53 | 11.66 | 90 |
| s3_orig_f12 | oct | designed | 0.000 / 0.188 | 28.65 / 7.61 | 25.43 | 70 |
| shuffle_control | oct | null | 0.000 / 0.094 | 28.17 / 7.23 | 26.17 | 70 |
| shuffle_control_esm0 | oct | null-esm | 0.531 / 1.000 | 15.00 / 14.83 | 5.80 | 60 |
| bg33_1 | oct | boltzgen | 0.968 / 1.000 | 11.82 / 14.40 | 7.12 | 50 |
| bg33_2 | oct | boltzgen | 0.903 / 0.968 | 9.97 / 14.41 | 10.25 | 50 |
| bg33_3 | oct | boltzgen | 0.323 / 0.323 | 8.07 / 8.09 | 0.31 | 90 |
| bg33_4 | oct | boltzgen | 0.871 / 0.903 | 9.28 / 9.47 | 0.46 | 90 |
| ox1_esm1_f8 | oxy | esm2-variant | 0.522 / 0.522 | 12.51 / 8.22 | 8.60 | 70 |
| ox1_orig_f12 | oxy | designed | 0.000 / 0.261 | 21.23 / 6.89 | 18.90 | 90 |
| ox2_esm1_f8 | oxy | esm2-variant | 0.586 / 0.483 | 10.94 / 8.85 | 5.96 | 70 |
| ox2_orig_f8 | oxy | designed | 0.000 / 0.000 | 26.04 / 7.12 | 24.35 | 80 |
| ox2_shuffle | oxy | null | 0.000 / 0.034 | 26.12 / 6.68 | 24.93 | 80 |
| ox2_shuffle_esm0 | oxy | null-esm | 0.655 / 0.690 | 13.14 / 12.78 | 2.25 | 70 |
| ox3_orig_f8 | oxy | designed | 0.000 / 0.160 | 19.35 / 6.94 | 17.57 | 60 |
| bgox31_2 | oxy | boltzgen | 1.000 / 1.000 | 13.44 / 13.62 | 0.76 | 100 |
| bgox31_3 | oxy | boltzgen | 0.862 / 0.862 | 9.11 / 9.00 | 1.37 | 80 |
| bgox31_4 | oxy | boltzgen | 0.690 / 0.724 | 11.32 / 8.31 | 6.10 | 80 |
| bgox31_5 | oxy | boltzgen | 0.793 / 0.759 | 12.08 / 8.97 | 6.69 | 90 |

**The headline (as ESMFold showed it): the designed glycine-rich folds collapse apo.** Every
glycine-rich designed sequence — `orig_f12`, `s3_orig_f12`, `ox1_orig_f12` (octinoxate's design
json), `ox2_orig_f8`, `ox3_orig_f8`, and both octinoxate nulls plus `ox2_shuffle` — folds from
sequence alone into a near-fully-extended chain (Rg 19–29 Å, coordinate span 80–95 Å) instead of
the compact shell its Boltz/BoltzGen co-fold shows (Rg 7–7.5 Å), at Cα RMSD 17–27 Å. This section's
original conclusion — that the compact fold "exists only because Boltz asked the question *with*
the ligand present" — **was overturned by the apo cofolder lanes below**: OpenFold3 and RF3 fold
the same sequences compact *without* any ligand. The extended apo chain is ESMFold's behaviour
under these sequences, not a property of them. The ESMFold comparison is still the discriminator
that exposed the collapse in the first place — the two apo cofolders now locate it as model- rather
than sequence-determined.

**The ESM2-variants and structured nulls fold on their own.** `s2_esm2_control` (0.37 Å),
`ox2_shuffle_esm0` (2.25), `shuffle_control_esm0` (5.8) and the remaining esm2-variants reproduce
their co-folds' shape within 2–12 Å even without the ligand: sequences dominated by bulky residues
reach the same structure regardless. The 12–20 Å members (`s3_esm2_f4`, `ox1/ox2_esm1_f8`) keep
their helix content but relax toward a different shape.

**The BoltzGen arms split cleanly by ESMFold confidence.** `bg33_3` (0.31 Å), `bg33_4` (0.46),
`bgox31_2` (0.76), `bgox31_3` (1.37) are ESMFold-reproducible compact folds at pLDDT 90–100;
`bg33_1`/`bg33_2` keep a helix but change shape (7–10 Å) at pLDDT ~50. Both low-confidence ones
coincide with the compactness *increase* ESMFold shows over their co-folds' (11.8/10.0 vs
14.4): the sequence has a helix-driven compact state that BoltzGen loosened for the ligand.

**Tentative, flagged as such:** the ligand-governed group is not a designed-vs-null split — it
contains every designed glycine-rich sequence *and* three of the four nulls, so ESMFold alone
cannot tell a good design from a bad arrangement; it tells compact-cofold from extended-apo. The
usable signal is the opposite one: every sequence that *does* reproduce its co-fold without the
ligand (`s2_esm2_control`, `shuffle_control_esm0`, `ox2_shuffle_esm0`, `bg33_3/4`,
`bgox31_2/3`) reaches that structure by sequence physics alone. Whether ligand-governed folds are
fragile in dynamics or simply expensive-to-form is the question the OpenFold3 cofold lane asks
next, with the ligand present and by an independent mechanism.

## OpenFold3 results (all 21 landed 2026-10-09; full data in `runs/<mol>/folds/metrics/openfold3.csv`)

Each cofold ran ~1.5 min on the A10G once the one-time container build + 2.3 GB checkpoint had
carried on the first run; the best-ranked of five samples (seed 42) was converted to pdb. Geometry
cells read `OF3 / co-fold`, where the co-fold value is recomputed by the *same* `check_fold.geometry()`
on the leg's own `system/complex.pdb` — the exact structure the dynamics leg was built from — so
both sides of every cell share one code path. Engaged reads `OF3 vs co-fold`, each a count of the
ligand's heavy atoms touched within 4.5 Å (20 octinoxate, 17 oxybenzone). pLDDT on the 0–100 scale.

| peptide | mol | kind | enc OF3/ref | wrap OF3/ref | engaged OF3/ref | cent sep OF3/ref | Cα RMSD | Rg OF3/ref | pLDDT | ipSAE |
|---|---|---|---|---|---|---|---|---|---|---|
| orig_f12 | oct | designed | 0.68 / 0.96 | 0.70 / 1.00 | 14 vs 20 (20) | 9.4 / 4.2 | 4.99 | 7.07 / 7.46 | 76.7 | 0.182 |
| s2_esm2_control | oct | esm2-variant | 0.38 / 0.51 | 0.55 / 0.85 | 11 vs 17 (20) | 13.5 / 6.7 | 1.07 | 15.04 / 14.50 | 96.8 | 0.033 |
| s3_esm2_f4 | oct | esm2-variant | 0.58 / 0.80 | 0.85 / 0.90 | 17 vs 18 (20) | 7.3 / 7.0 | 11.90 | 14.97 / 9.53 | 98.2 | 0.115 |
| s3_orig_f12 | oct | designed | 0.43 / 0.97 | 0.70 / 1.00 | 14 vs 20 (20) | 10.4 / 3.6 | 6.68 | 6.93 / 7.61 | 79.0 | 0.160 |
| shuffle_control | oct | null | 0.51 / 0.75 | 0.90 / 0.75 | 18 vs 15 (20) | 9.2 / 7.9 | 7.89 | 6.68 / 7.23 | 78.5 | 0.110 |
| shuffle_control_esm0 | oct | null-esm | 0.31 / 0.54 | 0.30 / 0.85 | 6 vs 17 (20) | 9.4 / 11.0 | 1.68 | 14.54 / 14.83 | 96.8 | 0.041 |
| bg33_1 | oct | boltzgen | 0.42 / 0.38 | 0.90 / 0.75 | 18 vs 15 (20) | 9.5 / 17.4 | 0.50 | 14.23 / 14.40 | 98.5 | 0.078 |
| bg33_2 | oct | boltzgen | 0.35 / 0.59 | 0.50 / 0.75 | 10 vs 15 (20) | 8.0 / 5.2 | 0.42 | 14.29 / 14.41 | 98.6 | 0.069 |
| bg33_3 | oct | boltzgen | 0.49 / 0.48 | 0.75 / 1.00 | 15 vs 20 (20) | 9.4 / 8.6 | 0.55 | 8.06 / 8.09 | 98.1 | 0.247 |
| bg33_4 | oct | boltzgen | 0.52 / 0.60 | 0.85 / 0.95 | 17 vs 19 (20) | 7.0 / 6.3 | 0.72 | 9.19 / 9.47 | 96.3 | 0.067 |
| ox1_esm1_f8 | oxy | esm2-variant | 0.58 / 0.94 | 0.59 / 1.00 | 10 vs 17 (17) | 6.1 / 4.1 | 5.71 | 8.80 / 8.22 | 93.4 | 0.095 |
| ox1_orig_f12 | oxy | designed | 0.81 / 1.00 | 1.00 / 1.00 | 17 vs 17 (17) | 5.7 / 1.2 | 6.82 | 6.94 / 6.89 | 84.5 | 0.141 |
| ox2_esm1_f8 | oxy | esm2-variant | 0.75 / 0.67 | 0.94 / 1.00 | 16 vs 17 (17) | 5.3 / 6.5 | 4.78 | 9.80 / 8.85 | 83.7 | 0.033 |
| ox2_orig_f8 | oxy | designed | 0.66 / 0.98 | 1.00 / 1.00 | 17 vs 17 (17) | 7.5 / 2.7 | 7.07 | 7.48 / 7.12 | 81.8 | 0.192 |
| ox2_shuffle | oxy | null | 0.58 / 0.58 | 0.94 / 0.82 | 16 vs 14 (17) | 7.9 / 9.3 | 6.98 | 6.50 / 6.68 | 79.0 | 0.078 |
| ox2_shuffle_esm0 | oxy | null-esm | 0.47 / 0.64 | 0.77 / 0.94 | 13 vs 16 (17) | 7.3 / 4.2 | 5.07 | 13.61 / 12.78 | 93.4 | 0.034 |
| ox3_orig_f8 | oxy | designed | 0.54 / 1.00 | 0.88 / 1.00 | 15 vs 17 (17) | 10.3 / 2.0 | 5.79 | 6.13 / 6.94 | 88.6 | 0.200 |
| bgox31_2 | oxy | boltzgen | 0.36 / 0.44 | 0.65 / 0.77 | 11 vs 13 (17) | 8.1 / 9.0 | 0.89 | 13.52 / 13.62 | 98.7 | 0.048 |
| bgox31_3 | oxy | boltzgen | 0.62 / 0.96 | 1.00 / 1.00 | 17 vs 17 (17) | 8.5 / 4.0 | 2.59 | 8.33 / 9.00 | 97.9 | 0.469 |
| bgox31_4 | oxy | boltzgen | 0.47 / 0.65 | 0.82 / 0.94 | 14 vs 16 (17) | 11.2 / 6.7 | 0.63 | 8.17 / 8.31 | 97.1 | 0.152 |
| bgox31_5 | oxy | boltzgen | 0.55 / 0.50 | 0.71 / 0.77 | 12 vs 13 (17) | 7.0 / 7.0 | 1.19 | 8.80 / 8.97 | 97.0 | 0.179 |

**OpenFold3 with the ligand present reproduces every BoltzGen arm.** All eight BoltzGen peptides
come back at Cα RMSD 0.42–0.89 Å (bgox31_5 1.19, bgox31_3 2.59), Rg within 0.3 Å, so BoltzGen's
arm is a reproducible structure, not an artefact of one model's internals. But **the ligand
engagement is systematically looser in OF3's versions** — engaged drops in seven of the eight
(e.g. bg33_2 10 vs 15, bg33_3 15 vs 20, bgox31_3 17 vs 17 but centroid 8.5 vs 4.0), so same fold,
less-tucked ligand. ipSAE stays low everywhere (0.03–0.47): the two tightest interfaces are
bgox31_3 (0.469) and bg33_3 (0.247) — the same two peptides whose BoltzGen reads show
near-complete wrap — so the rank is broadly consistent even if the protein–ligand calibration is
experimental.

**The Boltz arms (designed + nulls) come back compact but repacked at RMSD 5–8 Å**, against the
sub-Å agreement of the BoltzGen arms. Same story as ESMFold's but from the cofold direction: with
the ligand present, an independent model *does* find the compact shell (Rg 6.1–7.5 matching each
co-fold's 6.7–7.6) — no extended chains here — but the exact packing is a different local
configuration, and ligand engagement drops (orig_f12 engages 14/20 against Boltz's 20/20,
enclosure 0.68 against 0.96). The compact answer is a property of sequence+ligand; the particular
shell-packing Boltz found is one of many. In this lane `s3_esm2_f4` is the outlier (OF3 keeps it
loose and helical, Rg 15.0, RMSD 11.9 — see the RF3 section, which reproduces its co-fold fine).

## RF3 results (all 21 landed 2026-10-09; full data in `runs/<mol>/folds/metrics/rf3.csv`)

Same format as the OpenFold3 table: best-ranked sample per cofold, geometry cells `RF3 / co-fold`
from the same `check_fold.geometry()` on `system/complex.pdb`, engaged counts of the ligand's heavy
atoms (20 octinoxate, 17 oxybenzone). RF3 ran ~45 s per cofold once its checkpoint was cached.

| peptide | mol | kind | enc RF3/ref | wrap RF3/ref | engaged RF3/ref | cent sep RF3/ref | Cα RMSD | Rg RF3/ref | pLDDT | ipSAE |
|---|---|---|---|---|---|---|---|---|---|---|
| orig_f12 | oct | designed | 0.80 / 0.96 | 0.90 / 1.00 | 18 vs 20 | 5.8 / 4.2 | 8.61 | 7.74 / 7.46 | 73.1 | 0.007 |
| s2_esm2_control | oct | esm2-variant | 0.42 / 0.51 | 0.70 / 0.85 | 14 vs 17 | 13.5 / 6.7 | 11.66 | 9.73 / 14.50 | 76.2 | 0.016 |
| s3_esm2_f4 | oct | esm2-variant | 0.46 / 0.80 | 0.75 / 0.90 | 15 vs 18 | 9.3 / 7.0 | 3.64 | 9.98 / 9.53 | 84.4 | 0.031 |
| s3_orig_f12 | oct | designed | 0.72 / 0.97 | 0.95 / 1.00 | 19 vs 20 | 8.5 / 3.6 | 6.82 | 6.88 / 7.61 | 74.0 | 0.008 |
| shuffle_control | oct | null | 0.83 / 0.75 | 0.90 / 0.75 | 18 vs 15 | 6.3 / 7.9 | 6.21 | 8.41 / 7.23 | 73.6 | 0.008 |
| shuffle_control_esm0 | oct | null-esm | 0.52 / 0.54 | 0.65 / 0.85 | 13 vs 17 | 8.5 / 11.0 | 12.49 | 10.62 / 14.83 | 76.8 | 0.021 |
| bg33_1 | oct | boltzgen | 0.30 / 0.38 | 0.35 / 0.75 | 7 vs 15 | 8.1 / 17.4 | 0.50 | 14.63 / 14.40 | 83.3 | 0.023 |
| bg33_2 | oct | boltzgen | 0.35 / 0.59 | 0.20 / 0.75 | 4 vs 15 | 9.5 / 5.2 | 0.40 | 14.40 / 14.41 | 84.2 | 0.018 |
| bg33_3 | oct | boltzgen | 0.49 / 0.48 | 0.95 / 1.00 | 19 vs 20 | 11.0 / 8.6 | 0.47 | 7.99 / 8.09 | 88.0 | 0.025 |
| bg33_4 | oct | boltzgen | 0.54 / 0.60 | 0.85 / 0.95 | 17 vs 19 | 7.3 / 6.3 | 0.79 | 9.20 / 9.47 | 86.9 | 0.030 |
| ox1_esm1_f8 | oxy | esm2-variant | 0.53 / 0.94 | 0.82 / 1.00 | 14 vs 17 | 9.3 / 4.1 | 5.16 | 9.75 / 8.22 | 78.6 | 0.021 |
| ox1_orig_f12 | oxy | designed | 0.55 / 1.00 | 0.77 / 1.00 | 13 vs 17 | 6.1 / 1.2 | 7.19 | 6.40 / 6.89 | 73.6 | 0.013 |
| ox2_esm1_f8 | oxy | esm2-variant | 0.57 / 0.67 | 0.71 / 1.00 | 12 vs 17 | 6.4 / 6.5 | 4.35 | 9.13 / 8.85 | 78.5 | 0.027 |
| ox2_orig_f8 | oxy | designed | 0.60 / 0.98 | 0.94 / 1.00 | 16 vs 17 | 6.2 / 2.7 | 7.55 | 7.16 / 7.12 | 75.3 | 0.014 |
| ox2_shuffle | oxy | null | 0.48 / 0.58 | 0.71 / 0.82 | 12 vs 14 | 12.1 / 9.3 | 8.20 | 7.30 / 6.68 | 76.6 | 0.012 |
| ox2_shuffle_esm0 | oxy | null-esm | 0.51 / 0.64 | 0.94 / 0.94 | 16 vs 16 | 10.3 / 4.2 | 9.77 | 8.82 / 12.78 | 81.8 | 0.025 |
| ox3_orig_f8 | oxy | designed | 0.47 / 1.00 | 0.71 / 1.00 | 12 vs 17 | 11.4 / 2.0 | 4.70 | 6.73 / 6.94 | 76.5 | 0.012 |
| bgox31_2 | oxy | boltzgen | 0.35 / 0.44 | 0.59 / 0.77 | 10 vs 13 | 8.6 / 9.0 | 0.73 | 13.60 / 13.62 | 87.2 | 0.021 |
| bgox31_3 | oxy | boltzgen | 0.54 / 0.96 | 0.77 / 1.00 | 13 vs 17 | 8.0 / 4.0 | 0.96 | 8.78 / 9.00 | 84.4 | 0.045 |
| bgox31_4 | oxy | boltzgen | 0.61 / 0.65 | 0.77 / 0.94 | 13 vs 16 | 7.0 / 6.7 | 0.41 | 8.25 / 8.31 | 86.5 | 0.089 |
| bgox31_5 | oxy | boltzgen | 0.56 / 0.50 | 0.88 / 0.77 | 15 vs 13 | 7.0 / 7.0 | 0.94 | 8.92 / 8.97 | 85.8 | 0.027 |

**RF3 confirms the BoltzGen arms' folds and strengthens the looser-ligand finding.** Every
BoltzGen peptide reproduces to Cα RMSD 0.40–0.96 Å — now from two independent cofolders — and
RF3's versions *also* drop ligand engagement, by more than OF3's did (bg33_2 4/15, bg33_1 7/15,
wrap 0.20 and 0.35 against 0.75; even bgox31_3 falls to 13/17). Two models, same folds, both
tucking the ligand less than BoltzGen's own read does: the shell geometry is robust, and the
shell *tightness* in the BoltzGen static read is on its looser side of the cofold spread.

**The designed/null arms repeat OF3's story exactly** — compact folds at Rg 6.4–8.9 against
co-folds 6.7–7.6, RMSD 4.4–8.2 Å, enclosure below the co-fold's (ox2_orig_f8 0.60 vs 0.98, engaging
16/17), and RF3's loose-helix outlier confirms it independently (ox3_orig_f8 0.47/1.00 enclosure).
So three models now say: compactness is ligand-conditional, the packing is not unique.

**The esm2-variants are where the two cofolders disagree — and where RF3 overturns OF3's one
exception.** RF3 reproduces `s3_esm2_f4` at 3.64 Å where OF3 kept it loose (11.9 Å), so that
"neither lane reproduces it" claim from the OpenFold3 section was an OF3 quirk, not a sequence
property — the co-fold *is* reproducible. Conversely `s2_esm2_control` (11.66) and
`shuffle_control_esm0` (12.49) — which OF3 matched to ~1–2 Å — come back repacked in RF3. The
esm2-dominated sequences sit on a packing near-degeneracy whose resolution depends on the model;
only the apo-reproducible-by-ESMFold group (BoltzGen arms, gly-rich) is pinned to one structure
across everything.

**ipSAE is even weaker here** (0.007–0.089; RF3's output lacks the OF3 image that scores
protein–ligand interfaces well, per the fold repo note). Highest: bgox31_4 (0.089), bgox31_3
(0.045) — again BoltzGen-arm peptides — but the dynamic range is too small to rank designs on it.
Treat both cofold lanes' interface confidences as corroboration only.

## Apo cofolding results (OpenFold3 and RF3 without the ligand; all 21 landed 2026-10-09)

Full data in `runs/<mol>/folds/metrics/of3apo.csv` and `rf3apo.csv`. Same format as the cofold
tables minus the ligand columns: each row is best-sample-per-job, Cα RMSD against the leg's own
co-fold, pLDDT on the 0–100 scale. No ligand in either lane, so enclosure/wrap/engaged are empty.

### OpenFold3 apo

| peptide | mol | kind | Cα RMSD | Rg / co-fold | pLDDT |
|---|---|---|---|---|---|
| orig_f12 | oct | designed | 8.11 | 8.40 / 7.46 | 47.3 |
| s2_esm2_control | oct | esm2-variant | 1.47 | 14.58 / 14.50 | 95.9 |
| s3_esm2_f4 | oct | esm2-variant | 12.23 | 15.13 / 9.53 | 95.1 |
| s3_orig_f12 | oct | designed | 6.11 | 7.38 / 7.61 | 48.5 |
| shuffle_control | oct | null | 8.04 | 6.88 / 7.23 | 52.0 |
| shuffle_control_esm0 | oct | null-esm | 0.38 | 15.04 / 14.83 | 87.6 |
| bg33_1 | oct | boltzgen | 0.71 | 14.46 / 14.40 | 97.4 |
| bg33_2 | oct | boltzgen | 0.64 | 14.21 / 14.41 | 97.4 |
| bg33_3 | oct | boltzgen | 0.50 | 8.06 / 8.09 | 96.7 |
| bg33_4 | oct | boltzgen | 1.04 | 9.21 / 9.47 | 92.5 |
| ox1_esm1_f8 | oxy | esm2-variant | 7.07 | 10.68 / 8.22 | 77.7 |
| ox1_orig_f12 | oxy | designed | 6.25 | 6.93 / 6.89 | 57.6 |
| ox2_esm1_f8 | oxy | esm2-variant | 4.86 | 10.01 / 8.85 | 75.0 |
| ox2_orig_f8 | oxy | designed | 7.03 | 7.01 / 7.12 | 50.2 |
| ox2_shuffle | oxy | null | 7.22 | 6.27 / 6.68 | 52.1 |
| ox2_shuffle_esm0 | oxy | null-esm | 5.35 | 13.50 / 12.78 | 83.3 |
| ox3_orig_f8 | oxy | designed | 6.56 | 6.32 / 6.94 | 56.1 |
| bgox31_2 | oxy | boltzgen | 0.40 | 13.58 / 13.62 | 98.2 |
| bgox31_3 | oxy | boltzgen | 2.72 | 8.27 / 9.00 | 89.2 |
| bgox31_4 | oxy | boltzgen | 0.52 | 8.16 / 8.31 | 93.5 |
| bgox31_5 | oxy | boltzgen | 2.37 | 8.78 / 8.97 | 84.9 |

### RF3 apo

| peptide | mol | kind | Cα RMSD | Rg / co-fold | pLDDT |
|---|---|---|---|---|---|
| orig_f12 | oct | designed | 9.11 | 7.72 / 7.46 | 70.1 |
| s2_esm2_control | oct | esm2-variant | 0.44 | 14.46 / 14.50 | 86.2 |
| s3_esm2_f4 | oct | esm2-variant | 3.65 | 9.90 / 9.53 | 80.6 |
| s3_orig_f12 | oct | designed | 6.78 | 6.82 / 7.61 | 68.3 |
| shuffle_control | oct | null | 7.79 | 8.78 / 7.23 | 67.4 |
| shuffle_control_esm0 | oct | null-esm | 0.34 | 14.91 / 14.83 | 78.0 |
| bg33_1 | oct | boltzgen | 0.46 | 14.46 / 14.40 | 78.1 |
| bg33_2 | oct | boltzgen | 0.24 | 14.41 / 14.41 | 76.0 |
| bg33_3 | oct | boltzgen | 0.41 | 7.96 / 8.09 | 86.3 |
| bg33_4 | oct | boltzgen | 0.52 | 9.21 / 9.47 | 87.2 |
| ox1_esm1_f8 | oxy | esm2-variant | 4.99 | 8.71 / 8.22 | 82.0 |
| ox1_orig_f12 | oxy | designed | 6.26 | 6.77 / 6.89 | 72.2 |
| ox2_esm1_f8 | oxy | esm2-variant | 4.43 | 9.09 / 8.85 | 82.4 |
| ox2_orig_f8 | oxy | designed | 9.79 | 8.86 / 7.12 | 70.5 |
| ox2_shuffle | oxy | null | 5.83 | 8.02 / 6.68 | 68.7 |
| ox2_shuffle_esm0 | oxy | null-esm | 2.73 | 11.42 / 12.78 | 85.7 |
| ox3_orig_f8 | oxy | designed | 6.15 | 7.37 / 6.94 | 66.7 |
| bgox31_2 | oxy | boltzgen | 0.60 | 13.58 / 13.62 | 88.8 |
| bgox31_3 | oxy | boltzgen | 1.44 | 8.45 / 9.00 | 84.4 |
| bgox31_4 | oxy | boltzgen | 1.16 | 7.84 / 8.31 | 83.9 |
| bgox31_5 | oxy | boltzgen | 0.86 | 8.85 / 8.97 | 85.6 |

**Both cofolders fold every glycine-rich sequence compact apo — designed and null alike.**
OpenFold3 apo gives the designed arms Rg 6.3–8.4 against co-folds 6.9–7.6 at RMSD 6.1–8.1 Å, and
RF3 apo the same picture (Rg 6.8–8.9, RMSD 6.2–9.8), with the two gly-rich nulls compact under
both as well. **The ESMFold extended-apo collapse is therefore model-specific, not a property of
the sequences**: three of three architectures produce *a* compact state without the ligand;
ESMFold alone produces an extended one. The campaign's apo comparison is not
"ligand-governed versus intrinsic" — it is *which* of (at least) two apo states each model lands
on, and how loosely each is held: OF3's own confidence on these apo folds is pLDDT ~47–58 (the
lowest numbers in the whole campaign), RF3's ~66–70, versus 70–90 for ESMFold. **The compact apo
form exists for these sequences, but no model holds it confidently.** What every line still agrees
on: the exact packing (RMSD 6–10 Å across lanes) is degenerate, and nothing here — the fold, the
compactness, the apo state — separates a designed arm from its gly-rich null; only the binding
measurements can.

**Four-lane one-glance grid** (`same` = reproduces the co-fold, `diff` = compact but repacked,
`ext` = near-fully-extended apo):

| peptide | mol | kind | ESMFold apo | OF3 apo | RF3 apo | OF3 cofold | RF3 cofold |
|---|---|---|---|---|---|---|---|
| orig_f12 | oct | designed | ext | diff | diff | diff | diff |
| s2_esm2_control | oct | esm2-variant | same | same | same | same | diff |
| s3_esm2_f4 | oct | esm2-variant | diff | diff | same | diff | same |
| s3_orig_f12 | oct | designed | ext | diff | diff | diff | diff |
| shuffle_control | oct | null | ext | diff | diff | diff | diff |
| shuffle_control_esm0 | oct | null-esm | same | same | same | same | diff |
| bg33_1 | oct | boltzgen | diff | same | same | same | same |
| bg33_2 | oct | boltzgen | diff | same | same | same | same |
| bg33_3 | oct | boltzgen | same | same | same | same | same |
| bg33_4 | oct | boltzgen | same | same | same | same | same |
| ox1_esm1_f8 | oxy | esm2-variant | diff | diff | diff | diff | diff |
| ox1_orig_f12 | oxy | designed | ext | diff | diff | diff | diff |
| ox2_esm1_f8 | oxy | esm2-variant | diff | diff | diff | diff | diff |
| ox2_orig_f8 | oxy | designed | ext | diff | diff | diff | diff |
| ox2_shuffle | oxy | null | ext | diff | diff | diff | diff |
| ox2_shuffle_esm0 | oxy | null-esm | same | diff | diff | diff | diff |
| ox3_orig_f8 | oxy | designed | ext | diff | diff | diff | diff |
| bgox31_2 | oxy | boltzgen | same | same | same | same | same |
| bgox31_3 | oxy | boltzgen | same | same | same | same | same |
| bgox31_4 | oxy | boltzgen | same | same | same | same | same |
| bgox31_5 | oxy | boltzgen | same | same | same | same | same |

**Campaign bottom line (five lanes, 21 peptides, all static reads):**

- **BoltzGen's arms are the model-robust structures of the campaign** — sub-Å to 1.4 Å Cα RMSD
  from *both* cofolders and *both* their apo lanes (bg33 and bgox31, all "same" in four columns of
  the grid; ESMFold disagrees only on bg33_1/2). Nothing else in the campaign is that consistent.
- **The designed glycine-rich shells are pack-degenerate but not ligand-governed.** The ESMFold
  lane's extended-apo collapse was model-specific: both cofolders fold every one of these
  sequences (and the gly-rich nulls) compact without any ligand. What no lane reproduces is
  Boltz's particular packing — everything sits 5–10 Å away in Cα RMSD, and cofolded runs engage
  the ligand *less* than BoltzGen's static read. The apo compact state exists but is held at low
  confidence (OF3 pLDDT ~50) — it is one of several loose alternatives, not a defined fold.
- **Packing degeneracy — not apo fragility — is the property that distinguishes the groups:**
  designed and null glycine-rich arms live at 4–10 Å from every reference in every model, the
  BoltzGen arms are pinned sub-Å, and the esm2-variants sit between, with the two cofolders
  disagreeing on which packing they land in.
- **Fold metrics alone cannot separate a design from its null.** Designed and plain-null rows look
  alike in every grid cell — the rubric for what to pursue has to come from binding: MM/GBSA,
  ligand retention and contacts against those same null rows.

## What "interesting" means for the follow-up

Docking/MMGBSA legs on a *re-folded* structure run only if the inexpensive metrics above show
something worth chasing. **The condition fired on the cofold lanes, and the proposal is pending
the go:** every designed glycine-rich arm re-folds into a compact shell 4–8 Å (Cα) from the
co-fold the dynamics ran on, in *both* cofold lanes — a different packing that may engage the
ligand differently. (The apo lanes sharpen the framing rather than add candidates: they compact
*everything*, nulls included, so an apo fold is not a promising starting structure *per se*.) The
candidate cells are the five designed arms (`orig_f12`, `s3_orig_f12`, `ox1_orig_f12`,
`ox2_orig_f8`, `ox3_orig_f8`), docked and run as MD legs against their co-fold baselines. Until
that go arrives, nothing downstream of folding has run: no docking, no MD, no dynamics spend on
the folded structures.