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

**Apo grid** (`same` = reproduces the co-fold, `diff` = compact but repacked, `ext` =
near-fully-extended apo). Note the apo folds are compared against the *co-fold's* packing, so
`diff` does not mean wrong — with no ligand in the query, no model can know Boltz's packing:

| peptide | mol | kind | ESMFold apo | OF3 apo | RF3 apo |
|---|---|---|---|---|---|
| orig_f12 | oct | designed | ext | diff | diff |
| s2_esm2_control | oct | esm2-variant | same | same | same |
| s3_esm2_f4 | oct | esm2-variant | diff | diff | same |
| s3_orig_f12 | oct | designed | ext | diff | diff |
| shuffle_control | oct | null | ext | diff | diff |
| shuffle_control_esm0 | oct | null-esm | same | same | same |
| bg33_1 | oct | boltzgen | diff | same | same |
| bg33_2 | oct | boltzgen | diff | same | same |
| bg33_3 | oct | boltzgen | same | same | same |
| bg33_4 | oct | boltzgen | same | same | same |
| ox1_esm1_f8 | oxy | esm2-variant | diff | diff | diff |
| ox1_orig_f12 | oxy | designed | ext | diff | diff |
| ox2_esm1_f8 | oxy | esm2-variant | diff | diff | diff |
| ox2_orig_f8 | oxy | designed | ext | diff | diff |
| ox2_shuffle | oxy | null | ext | diff | diff |
| ox2_shuffle_esm0 | oxy | null-esm | same | diff | diff |
| ox3_orig_f8 | oxy | designed | ext | diff | diff |
| bgox31_2 | oxy | boltzgen | same | same | same |
| bgox31_3 | oxy | boltzgen | same | same | same |
| bgox31_4 | oxy | boltzgen | same | same | same |
| bgox31_5 | oxy | boltzgen | same | same | same |

**Cofold grid** (same key; here the fold carries the ligand, so `same`/`diff` say whether an
independent cofolder reproduces Boltz's packing and ligand engagement):

| peptide | mol | kind | OF3 cofold | RF3 cofold |
|---|---|---|---|---|
| orig_f12 | oct | designed | diff | diff |
| s2_esm2_control | oct | esm2-variant | same | diff |
| s3_esm2_f4 | oct | esm2-variant | diff | same |
| s3_orig_f12 | oct | designed | diff | diff |
| shuffle_control | oct | null | diff | diff |
| shuffle_control_esm0 | oct | null-esm | same | diff |
| bg33_1 | oct | boltzgen | same | same |
| bg33_2 | oct | boltzgen | same | same |
| bg33_3 | oct | boltzgen | same | same |
| bg33_4 | oct | boltzgen | same | same |
| ox1_esm1_f8 | oxy | esm2-variant | diff | diff |
| ox1_orig_f12 | oxy | designed | diff | diff |
| ox2_esm1_f8 | oxy | esm2-variant | diff | diff |
| ox2_orig_f8 | oxy | designed | diff | diff |
| ox2_shuffle | oxy | null | diff | diff |
| ox2_shuffle_esm0 | oxy | null-esm | diff | diff |
| ox3_orig_f8 | oxy | designed | diff | diff |
| bgox31_2 | oxy | boltzgen | same | same |
| bgox31_3 | oxy | boltzgen | same | same |
| bgox31_4 | oxy | boltzgen | same | same |
| bgox31_5 | oxy | boltzgen | same | same |

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

## Vina docking into the folded structures (all 105, 2026-10-09)

Every fold structure was Vina-docked (`code/dock_folds.py`, which stages the fold's peptide as
receptor and calls `vina_redock.redock()` unchanged: 22 Å box, exhaustiveness 16, seed 42, 9
modes, project cutoffs). 21 peptides × 5 lanes = 105 docks in ~25 min of Mac Vina; the rows are
appended additions-only to the committed `dock_summary.csv`/`dock_poses.csv` (50 + 55 summary rows,
500 + 550 pose rows; structure names carry the lane suffix, e.g. `s3_orig_f12_o3cof`). Full scratch
tables in `runs/<mol>/folds/metrics/dock_folds_{summary,poses}.csv`.

Box definition is the interesting part. A **cofold** fold contains the ligand: its atoms are
perceived from geometry and the reference molecule rigidly aligned onto them, so the search box is
centred on exactly where the cofolder put the ligand and pose 0 is that placement. An **apo** fold
contains no ligand, so there is no local information to centre on: the OF3-cofold reference pose is
carried over onto the apo fold by Cα superposition, which makes "look here, where this peptide
packs its ligand" a site definition rather than a prediction.

Perception trap worth recording: Open Babel loses the alkene on strained cofold geometries
(`C=Cc1ccc(OC)cc1` read back as `CCc1...`), which fails the template substructure match; the fix
is RDKit proximity bonding + `AllChem.AssignBondOrdersFromTemplate`, which imposes the template's
connectivity wherever the perceived graph admits it. Obabel remains the fallback.

What the docks say:

**Best Vina score (kcal/mol) / pose-1 in-place RMSD (Å) to the structure's reference placement** —
cofold lanes' reference is the fold's own ligand (pose 0); apo lanes' reference is the
superposed-site pose. Pose-1 RMSD is how far Vina's best pose lands from where the cofolder (or
site definition) put the ligand. Full per-pose data in `runs/<mol>/folds/metrics/dock_folds_poses.csv`.

| peptide | mol | kind | o3cof | rf3cof | esmapo | o3apo | rf3apo |
|---|---|---|---|---|---|---|---|
| orig_f12 | oct | designed | -4.2 / 10.2 | -5.6 / 6.2 | -2.7 / 8.2 | -3.6 / 7.4 | -3.4 / 8.2 |
| s3_orig_f12 | oct | designed | -4.3 / 8.0 | -4.9 / 3.0 | -2.3 / 7.4 | -3.9 / 10.7 | -3.9 / 8.8 |
| ox1_orig_f12 | oxy | designed | -5.7 / 1.2 | -4.6 / 4.6 | -3.4 / 6.8 | -6.1 / 12.4 | -4.3 / 13.2 |
| ox2_orig_f8 | oxy | designed | -6.3 / 6.7 | -4.5 / 9.3 | -3.3 / 9.0 | -5.0 / 8.9 | -4.6 / 5.2 |
| ox3_orig_f8 | oxy | designed | -4.6 / 5.8 | -4.4 / 6.0 | -4.5 / 9.6 | -5.0 / 8.1 | -4.0 / 9.5 |
| s2_esm2_control | oct | esm2-variant | -4.4 / 12.6 | -4.8 / 10.6 | -5.2 / 11.0 | -4.4 / 6.4 | -4.1 / 9.2 |
| s3_esm2_f4 | oct | esm2-variant | -5.1 / 6.9 | -4.5 / 9.3 | -4.9 / 11.2 | -4.9 / 11.0 | -4.1 / 10.8 |
| ox1_esm1_f8 | oxy | esm2-variant | -5.0 / 1.5 | -4.4 / 7.0 | -4.2 / 12.8 | -4.7 / 7.7 | -4.6 / 12.0 |
| ox2_esm1_f8 | oxy | esm2-variant | -6.5 / 5.9 | -5.0 / 13.0 | -4.5 / 10.4 | -5.2 / 8.1 | -5.5 / 7.2 |
| shuffle_control | oct | null | -3.8 / 13.6 | -4.9 / 8.1 | -2.9 / 10.7 | -3.8 / 10.2 | -4.6 / 10.1 |
| ox2_shuffle | oxy | null | -5.5 / 5.9 | -3.6 / 11.2 | -3.4 / 10.7 | -4.6 / 8.9 | -4.4 / 8.6 |
| shuffle_control_esm0 | oct | null-esm | -4.5 / 10.7 | -5.3 / 4.1 | -4.9 / 9.7 | -4.9 / 7.6 | -4.9 / 8.2 |
| ox2_shuffle_esm0 | oxy | null-esm | -5.5 / 10.3 | -5.1 / 10.9 | -5.7 / 13.7 | -4.7 / 12.2 | -5.3 / 8.8 |
| bg33_1 | oct | boltzgen | -3.3 / 7.6 | -3.6 / 10.0 | -3.4 / 9.8 | -3.6 / 10.9 | -3.2 / 9.9 |
| bg33_2 | oct | boltzgen | -3.4 / 10.0 | -2.9 / 12.9 | -4.4 / 11.1 | -2.9 / 8.0 | -3.3 / 9.5 |
| bg33_3 | oct | boltzgen | -4.2 / 9.5 | -4.0 / 10.4 | -4.0 / 12.4 | -4.1 / 9.6 | -4.4 / 5.5 |
| bg33_4 | oct | boltzgen | -4.5 / 10.6 | -4.1 / 4.6 | -4.1 / 9.8 | -4.2 / 5.3 | -4.0 / 12.5 |
| bgox31_2 | oxy | boltzgen | -4.8 / 7.8 | -4.1 / 10.9 | -4.0 / 6.9 | -4.2 / 7.9 | -4.4 / 6.9 |
| bgox31_3 | oxy | boltzgen | -5.2 / 7.2 | -5.1 / 6.4 | -5.3 / 11.1 | -4.7 / 10.9 | -4.8 / 11.7 |
| bgox31_4 | oxy | boltzgen | -4.6 / 7.5 | -4.9 / 6.1 | -4.7 / 9.2 | -4.7 / 10.0 | -4.0 / 13.4 |
| bgox31_5 | oxy | boltzgen | -4.6 / 5.9 | -5.0 / 8.0 | -4.5 / 5.6 | -5.5 / 7.3 | -5.2 / 9.5 |

- **Vina does not reproduce cofolder packing in these shells.** Pose 1 sits 1.2–13.7 Å in-place
  RMSD from the fold's own ligand placement, with a median around 9 Å: in 103 of 105 structures
  the independent pose search finds different packing than the cofolder did. The two agreements —
  `ox1_orig_f12` and `ox1_esm1_f8` in the OF3 cofold lane at 1.16 and 1.51 Å — are the oxybenzone
  structures whose OF3 cofolds sit tightest to Boltz in the cofold tables. This is the
  packing-degeneracy story from the fold geometry, seen from a third direction.
- **But the docked poses stay inside the envelope.** Pose-1 wrapped fraction is 0.77–1.00
  everywhere and engaged counts run 13–20. A Vina search 9 Å from the cofold placement still
  leaves the ligand packed against the chain — these shells are pockets, not grooves, for every
  peptide including the nulls.
- **Scores say nothing separating.** Best scores span −2.3 to −6.6 kcal/mol across all 105, with
  no designed-vs-null gap in any lane. Consistent with "fold metrics cannot separate a design
  from its null"; the expectation is the same for these docks.
- **The 9 MD legs run the OF3-cofold structures' OWN ligand placement** — the cofold pose,
  pose 0 in the dock tables, not a Vina pose. The cofolder's placement is the structure being
  tested; docking it again (2026-10-09 showed Vina moves it a median ~9 Å) would have measured
  Vina's packing instead. `run_dock_pose_md.sh POSE=cofold` stages the leg-ready pose from
  `runs/<mol>/folds/dock_stage/openfold3/<pep>_o3cof_ligand_cofold.sdf` (the fold's perceived
  ligand with template bond orders, template atom order).

## What "interesting" means for the follow-up

Docking/MMGBSA legs on a *re-folded* structure run only if the inexpensive metrics above show
something worth chasing. The condition fired, and **the go arrived 2026-10-09:** nine legs are
running on Modal — the five designed arms (`orig_f12`, `s3_orig_f12`, `ox1_orig_f12`,
`ox2_orig_f8`, `ox3_orig_f8`) plus the four esm-variants (`s2_esm2_control`, `s3_esm2_f4`,
`ox1_esm1_f8`, `ox2_esm1_f8`), each as `<pep>_o3cof`: an OF3-cofold fold run with the ligand
where the cofolder put it (POSE=cofold), 20 ns and scored against the same peptides' Boltz
co-fold baselines. (One docked-pose leg, `orig_f12_o3cof_dock1`, runs first as a separate
comparison cell; the campaign itself is the cofold poses.) ~$10, one
leg at a time, tails done as each leg lands. Every one of the 105 fold structures is now
docked-ready at `runs/<mol>/dock/<name>/` should more legs be wanted — the apo-fold legs would
need that superposed-site definition, which is a site choice, not a prediction.

### The MD legs, as they land

Every leg records: whole-run + window ΔG (`mmgbsa_summary.csv`), contacts/residence
(`md_contacts.csv`), end-of-run frames, slide RMSD (`code/ligand_slide.py`), and a
`protein_stability.csv` row. Reference legs: `orig_f12` Boltz pose **−19.75**;
`orig_f12_dock1` (Boltz shell, Vina pose 1) **−31.40** — a Vina pose can read strong in a
good shell, so a weak re-fold number is not a docking artifact.

| leg | structure | ligand | ΔG 20 ns | 5/10/15 ns | slide Å | Cα drift (final) | retention |
|---|---|---|---|---|---|---|---|
| orig_f12_o3cof_dock1 | OF3 refold of orig_f12 | Vina p1 (10.2 Å from cofold placement) | −9.65 ± 0.11 | −9.13 / −10.86 / −11.28 | 2.34 | 4.6 | **released ~18 ns** (29→4.4 contacts, 15.8 Å) |
| orig_f12_o3cof | OF3 refold of orig_f12 | cofold's own placement | −5.90 ± 0.09 | −9.62 / −8.31 / −7.34 | 1.26 | 5.5 | **released** (21→1.8 contacts, 11→18.3 Å) |
| s3_orig_f12_o3cof | OF3 refold of s3_orig_f12 | cofold's own placement | −10.10 ± 0.09 | −14.04 / −10.92 / −10.51 | 1.88 | 3.9 | partial loss (33→12.5 contacts) |
| s2_esm2_control_o3cof | OF3 refold of s2_esm2_control | cofold's own placement | −10.21 ± 0.11 | −7.32 / −8.52 / −9.86 | 3.20 | 2.5 | full exit + return (18.7→0→25.6→13.0) |
| s3_esm2_f4_o3cof | OF3 refold of s3_esm2_f4 | cofold's own placement | −13.76 ± 0.08 | −13.58 / −14.61 / −14.52 | 1.30 | 1.9 | best of the refolds (11–21 contacts held, ~6–10 Å) |
| ox1_orig_f12_o3cof | OF3 refold of ox1_orig_f12 | cofold's own placement | −3.96 ± 0.10 | −7.60 / −5.17 / −4.29 | 2.94 | 2.9 | worst refold: bouncing (37→0→18.5→0.5, excursions to 22 Å) |
| ox2_orig_f8_o3cof | OF3 refold of ox2_orig_f8 | cofold's own placement | −7.23 ± 0.14 | −14.87 / −11.98 / −8.10 | 2.32 | 4.1 | walked out after 10 ns: 56/8 → 0 contacts at 12–16 ns, partial return |
| ox3_orig_f8_o3cof | OF3 refold of ox3_orig_f8 | cofold's own placement | −2.06 ± 0.06 | −0.43 / −1.44 / −1.91 | 3.12 | 5.3 | never bound: 18–19 Å separation the whole run, ~1 contact, mid-run drift to 10–13 then back out |
| ox1_esm1_f8_o3cof | OF3 refold of ox1_esm1_f8 | cofold's own placement | **−17.95 ± 0.09** | −19.54 / −19.23 / −18.92 | 3.00 | 4.2 | **the one exception**: stronger than its Boltz baseline (−10.70) by 7.2; flat profile, 28–48 contacts at ~5–7 Å the whole run |
| ox2_esm1_f8_o3cof | OF3 refold of ox2_esm1_f8 | cofold's own placement | −16.68 ± 0.11 | −19.08 / −17.54 / −15.70 | 1.20 | 4.0 | strong: weaker than its Boltz baseline (−20.81) by 4.1 but held — 23–46 contacts at ~5.6–6.7 Å, 0 released frames, best slide in the campaign (1.20 Å) |

## The cofold campaign: the five headline tables

The campaign question, answered across all ten legs: **is the OF3 re-fold worth running dynamics on,
or does the Boltz co-fold remain the start the pipeline should use?** Eight of nine re-folds read
4–18 kcal/mol *weaker* than their Boltz baseline with degraded or lost retention; the two esm1
re-folds of ox1/ox2 (`ox1_esm1_f8_o3cof`, `ox2_esm1_f8_o3cof`) are the exception — the first beats its
baseline by 7.2. Spend: ten legs ≈ $8.2 of produce (plus the $0.78 erroneous dock-pose leg), one at a
time on Modal.

### Co-folds, ΔG and retention: Boltz pose vs OF3 re-fold

| peptide | mol | Boltz co-fold ΔG (ret.) | OF3 re-fold ΔG (ret.) | Δ (refold − Boltz) |
|---|---|---|---|---|
| `orig_f12` | oct | **−13.71 (77.3%)** | −5.90 (7.4%) | +7.8 weaker, released |
| `orig_f12` (Vina p1 in the OF3 shell) | oct | −13.71 (77.3%) | −9.65 (50.2%) | +4.1 weaker, partial loss |
| `s3_orig_f12` | oct | **−24.33 (100.0%)** | −10.10 (44.8%) | +14.2 weaker, partial loss |
| `s2_esm2_control` | oct | **−16.25 (59.6%)** | −10.21 (33.9%) | +6.0 weaker, one long exit-and-return |
| `s3_esm2_f4` | oct | **−21.08 (99.6%)** | −13.76 (84.8%) | +7.3 weaker but held (best oct refold) |
| `ox1_orig_f12` | oxy | −10.31 (49.1%) | −3.96 (20.5%) | +6.4 weaker, bouncing |
| `ox2_orig_f8` | oxy | **−24.21 (100.0%)** | −7.23 (49.5%) | +17.0 weaker, walked out at 10 ns |
| `ox3_orig_f8` | oxy | **−18.84 (100.0%)** | −2.06 (2.4%) | +16.8 weaker, never bound |
| `ox1_esm1_f8` | oxy | −10.70 (75.0%) | **−17.95 (100.0%)** | **−7.2 stronger, held all run** |
| `ox2_esm1_f8` | oxy | **−20.81 (100.0%)** | −16.68 (93.2%) | +4.1 weaker but held, 0 releases |

Eight of nine re-folds are worse by 4–17 kcal/mol and most of those lose the ligand; only the two
ox esm1 re-folds retain, and `ox1_esm1_f8`'s re-fold is the single cell where the OF3 pose beats the
model that predicted the baseline. The Boltz co-fold stays the campaign's default start.

### Co-folds — window convergence

| leg | 5 ns | 10 ns | 15 ns | 20 ns |
|---|---|---|---|---|
| `orig_f12_o3cof_dock1` | −9.13 | −10.86 | −11.28 | −9.65 |
| `orig_f12_o3cof` | −9.62 | −8.31 | −7.34 | **−5.90** |
| `s3_orig_f12_o3cof` | −14.04 | −10.92 | −10.51 | −10.10 |
| `s2_esm2_control_o3cof` | −7.32 | −8.52 | −9.86 | −10.21 |
| `s3_esm2_f4_o3cof` | −13.58 | −14.61 | −14.52 | −13.76 |
| `ox1_orig_f12_o3cof` | −7.60 | −5.17 | −4.29 | −3.96 |
| `ox2_orig_f8_o3cof` | −14.87 | −11.98 | −8.10 | **−7.23** |
| `ox3_orig_f8_o3cof` | −0.43 | −1.44 | −1.91 | −2.06 |
| `ox1_esm1_f8_o3cof` | −19.54 | −19.23 | −18.92 | −17.95 |
| `ox2_esm1_f8_o3cof` | −19.08 | −17.54 | −15.70 | −16.68 |

Every release leg shows the decay signature (`orig_f12_o3cof` and `ox2_orig_f8_o3cof` give back
~7 and ~8 over their windows); the retainers (`s3_esm2_f4`, `ox1/ox2_esm1_f8`) hold a flat
−13 to −19 band the whole run, so their numbers are settled poses rather than starting-pose
artifacts.

### Co-folds — residence and release

Residence within 10 Å of the peptide centroid at 10 ps sampling; `late` is the share of released
frames in the second half of the run.

| leg | residence | released frames | episodes | longest | late |
|---|---|---|---|---|---|
| `orig_f12_o3cof_dock1` | 50.2% | 86 (4.3%) | 31 | 230 ps | 92% |
| `orig_f12_o3cof` | 7.4% | 442 (22.1%) | 53 | 2140 ps | 86% |
| `s3_orig_f12_o3cof` | 44.8% | 28 (1.4%) | 14 | 80 ps | 50% |
| `s2_esm2_control_o3cof` | 33.9% | 276 (13.8%) | 8 | **2680 ps** | 1% |
| `s3_esm2_f4_o3cof` | 84.8% | 7 (0.35%) | 5 | 30 ps | 14% |
| `ox1_orig_f12_o3cof` | 20.5% | 777 (38.9%) | 40 | 1800 ps | 66% |
| `ox2_orig_f8_o3cof` | 49.5% | 624 (31.2%) | 16 | **3010 ps** | 99% |
| `ox3_orig_f8_o3cof` | 2.4% | 1059 (53.0%) | 80 | 1470 ps | 35% |
| `ox1_esm1_f8_o3cof` | 100.0% | 0 | 0 | — | — |
| `ox2_esm1_f8_o3cof` | 93.2% | 0 | 0 | — | — |

The two esm1 re-folds and `s3_esm2_f4` are the only genuine retainers (0–7 released frames in
20 ns). `s2_esm2_control`'s single 2.68 ns excursion with 1% late is a mid-run visit that
returned — a different failure from the `ox2_orig_f8` 99%-late progressive walk-out of 3.01 ns.

### Co-folds — peptide structural stability

Cα RMSD final/max and Rg start→end, Å (`protein_stability.csv`).

| leg | Cα final | Cα max | Rg start→end |
|---|---|---|---|
| `orig_f12_o3cof_dock1` | 4.6 | 6.0 | 8.8 → 8.6 |
| `orig_f12_o3cof` | 5.5 | 6.5 | 9.0 → 8.7 |
| `s3_orig_f12_o3cof` | 3.9 | 4.6 | 8.4 → 8.7 |
| `s2_esm2_control_o3cof` | 2.5 | 3.8 | 15.3 → 14.7 |
| `s3_esm2_f4_o3cof` | 1.9 | 3.9 | 15.2 → 15.4 |
| `ox1_orig_f12_o3cof` | 2.9 | 3.8 | 8.6 → 7.8 |
| `ox2_orig_f8_o3cof` | 4.1 | 4.8 | 9.0 → 8.1 |
| `ox3_orig_f8_o3cof` | 5.3 | 6.2 | 7.8 → 8.6 |
| `ox1_esm1_f8_o3cof` | 4.2 | 4.7 | 9.5 → 9.6 |
| `ox2_esm1_f8_o3cof` | 4.0 | 5.7 | 10.0 → 10.4 |

No re-fold's Cα drift exceeds the Boltz legs' own band; failure here is ligand retention, not fold
collapse. The long-Rg legs (`s2_esm2_control`, `s3_esm2_f4` at ~15 Å) are the elongated designs and
their shells stayed the same size through dynamics.

### The pre-dynamics static read and the pairs the re-folds realise

Static half: `dock_poses.csv` on each placement's own coordinates — **enclosed** (fraction of
directions out of the ligand that meet peptide), **wrapped** (fraction of ligand atoms with peptide
within 4.5 Å), engaged atoms, centroid separation — pose 0 of the same structures the legs ran from.
Trajectory half: `pair_contacts.py` on each re-fold leg's 20 ns run — mean simultaneous designed
side chains on the ligand (`sim`) and realised pairs against the [n−1, n(n−1)/2] band.

Pair outputs are persisted as `runs/{octinoxate,oxybenzone}/md/pair_contacts_o3cof_shellN.csv`
(the `orig_f12` pair under `pair_contacts_o3cof_design.csv`); static half from
`runs/<mol>/dock/dock_poses.csv` pose 0.

| peptide | pose | enclosed | wrapped | engaged | cent sep | sim | pairs | band% |
|---|---|---|---|---|---|---|---|---|
| `orig_f12` | Boltz co-fold | 0.96 | 1.00 | 20 | 4.2 Å | | 22/66 | |
| | OF3 re-fold | 0.68 | 0.70 | 14 | 9.4 Å | 0.88 | 12/66 | 2 |
| | OF3 re-fold, Vina p1 | 0.76 | 0.95 | 19 | 6.5 Å | 1.37 | 18/66 | 13 |
| `s3_orig_f12` | Boltz co-fold | 0.965 | 1.00 | 20 | 3.6 Å | | 53/66 | |
| | OF3 re-fold | 0.43 | 0.65 | 13 | 10.4 Å | 1.39 | 18/66 | 13 |
| `s2_esm2_control` | Boltz co-fold | 0.51 | 0.85 | 17 | 6.7 Å | | 10/66 | |
| | OF3 re-fold | 0.38 | 0.50 | 10 | 13.5 Å | 1.08 | 13/66 | 4 |
| `s3_esm2_f4` | Boltz co-fold | 0.795 | 0.90 | 18 | 7.0 Å | | 32/66 | |
| | OF3 re-fold | 0.58 | 0.60 | 12 | 7.3 Å | 1.63 | 9/66 | −4 |
| `ox1_orig_f12` | Boltz co-fold | 0.995 | 1.00 | 17 | 1.2 Å | | 42/45 | |
| | OF3 re-fold | 0.81 | 0.94 | 16 | 5.7 Å | 0.76 | 15/45 | 17 |
| `ox2_orig_f8` | Boltz co-fold | 0.975 | 1.00 | 17 | 2.7 Å | | 10/55 | |
| | OF3 re-fold | 0.655 | 1.00 | 17 | 7.5 Å | 1.48 | 15/55 | 11 |
| `ox3_orig_f8` | Boltz co-fold | 0.995 | 1.00 | 17 | 2.1 Å | | 26/45 | |
| | OF3 re-fold | 0.54 | 0.77 | 13 | 10.3 Å | 0.46 | 8/45 | −3 |
| `ox1_esm1_f8` | Boltz co-fold | 0.94 | 1.00 | 17 | 4.1 Å | | 30/45 | |
| | OF3 re-fold | 0.58 | 0.47 | 8 | 6.1 Å | 2.36 | 12/45 | 8 |
| `ox2_esm1_f8` | Boltz co-fold | 0.665 | 1.00 | 17 | 6.4 Å | | 36/55 | |
| | OF3 re-fold | 0.745 | 1.00 | 17 | 5.3 Å | 2.15 | 19/55 | 20 |

The static read priced the re-folds before the GPU did: the Boltz co-folds enclose 0.51–0.995 at
1.2–6.7 Å, the OF3 placements enclose 0.38–0.81 at 5.3–13.5 Å, and every leg that read weaker in the
ΔG table sat on the low half of that drop. The one high-energy counterexample, `ox2_esm1_f8`, is
also the placement that matched its Boltz pose best statically (0.745 enclosed, 17 atoms engaged,
ligand *inside* its Boltz centroid distance at 5.3 Å) — and in dynamics it realises the campaign's
best pair band (19/55) with the longest `held` list of the ten legs. Re-fold quality is legible in
the placement's enclosure fraction at zero dynamics cost.