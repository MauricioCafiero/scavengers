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

## The six lanes and the question each answers

| Lane | Model | Input | Where | Question |
|---|---|---|---|---|
| ESMFold | esmfold_v1 (3B) | sequence only | Modal A10G | is the designed fold a fold **without the ligand**? — ⚠ **retired from downstream analysis 2026-10-11**, see below |
| OpenFold3 | OF3 (+ ligand) | sequence + SMILES | Modal A10G | does an independent cofolder place the ligand like BoltzGen did? |
| RF3 | RosettaFold3 | sequence + SMILES | Modal A10G | second opinion on the cofold (all 21 landed 2026-10-09) |
| OpenFold3 apo | OF3, apo | sequence only | Modal A10G | is the compact apo fold real or ESMFold-specific? (all 21 landed 2026-10-09) |
| RF3 apo | RF3, apo | sequence only | Modal A10G | second opinion on the apo fold (all 21 landed 2026-10-09) |
| Boltz apo | Boltz-2, apo | sequence only | local MPS | does the *originating* model keep its own co-fold's shell without the ligand? (all 21 landed 2026-10-11) |

Five of the six run through the `fold` repo (`~/python_mac/fold`, cloned from
`github.com/MauricioCafiero/fold`; see its HANDOFF.md for what was verified live). ESMFold is **not
local** — its 11 GB of weights need more than this Mac has; only ESM2 *embeddings* are local-CPU,
and they are not used here. `code/fold_batch.sh` drives a lane: `LANE=<name> code/fold_batch.sh`,
resumable over its outputs, sequential (one Modal run at a time), each invocation under
`caffeinate -i` as the modal guard requires. The sixth lane, `boltzapo`, is the exception: it is
the same local Boltz-2 that produced the co-folds in the first place, so it runs on this Mac's MPS
through `code/boltz_env.py` and costs nothing.

### The six methods, in detail

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

**6. Boltz apo — the originating model, asked without its ligand.** The lane that closes the
circle (2026-10-11). Every co-fold in this project came from Boltz-2 *with* the ligand and forced
contact constraints; this asks the identical model for the sequence alone. `LANE=boltzapo
code/fold_batch.sh` writes a minimal yaml — one protein chain, `msa: empty` exactly as the co-fold
yamls in `runs/<mol>/boltz/` carry, no ligand block, no `constraints:` section, no
`--use_potentials` — and folds it through `boltz_env.run_boltz`, the same resolver the design
pipeline uses, which lands on the local MPS wrapper (`boltz_local/code/boltz_mps.py`). Its raw
yaml, log and Boltz output tree stay under `runs/<mol>/folds/boltzapo/raw/`; the best model's cif
converts to pdb through the same `fold.analyze cif-to-pdb` the Modal lanes use. 45 s per peptide,
~15 min for all 21, no Modal spend. Its role is the control the other apo lanes cannot be: when
OF3 and RF3 disagree with ESMFold about the apo state, neither of them is the model whose answer
the co-folds — and therefore every MD leg — were built on.

**What every lane feeds, downstream of this repo.** Three consumers, all of which must be
extended when a lane is added (see CLAUDE.md's "Downstream analysis repos"): `code/dock_folds.py`
here (Vina into every fold, then `make_gnina_bundle.py` for the CNN rescoring packet),
`~/python_mac/pepgraph` (inverse-folding perplexity and z per lane — add the lane to its `LANES`
and run with `--lanes <lane>`, never a full re-run), and `~/python_mac/morph` (apo → co-fold strain
and surface patches — add the pair to its `MATCHED_PAIRS`).

**What every lane feeds here.** `code/folds_metrics.py <lane>` reads the converted PDB against the leg's
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
the ligand present" — **was overturned by the apo cofolder lanes below, then partly reinstated by
the Boltz apo lane**: OpenFold3 and RF3 fold these sequences compact without any ligand, but
Boltz-2 itself, the model that made the co-folds, extends `orig_f12`, `ox1_orig_f12`,
`ox1_esm1_f8` and `ox2_shuffle`. So the extended apo chain is not ESMFold's behaviour alone — on
those four it is also the originating model's — while for the other seventeen sequences the
compactness holds across three apo models and ESMFold is the outlier. Read the per-peptide grid
rather than any one lane's headline; the Rg magnitudes still differ (ESMFold 19–29 Å against
Boltz's 11–15 Å on the same arms), so ESMFold extends furthest wherever it extends at all.

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

### ESMFold's standing, after the inverse-folding check (2026-10-11)

**ESMFold is the least reliable lane of the six and is no longer used downstream.** Scoring every
lane's structures with a PDB-trained inverse-folding model (the `pepgraph` repo,
`scripts/compare_lanes.py`: perplexity = the effective number of residue types the geometry is
choosing between at each position, 21 = no information) puts it last by a wide margin — mean 20.4
against 10.7–15.7 for the other five lanes, i.e. at the no-information limit. It owns **7 of the
11** structures in 202 whose backbone carries no sequence-order information at all (z ≤ 0), and
**6 of those have exactly zero residue contacts**. For these ligand-dependent designs it emits
extended chains with no secondary structure, so a measurement that starts from an ESMFold
structure is measuring a model failure rather than an apo state.

Consequence: the `morph` repo's apo → co-fold strain and surface-patch analysis dropped it on
2026-10-11. That analysis originally morphed ESMFold, OF3-apo and RF3-apo all to the one Boltz
co-fold, because no Boltz apo fold existed; now that every co-folder has an apo lane it runs three
model-matched pairs instead (Boltz apo → Boltz, OF3 apo → OF3, RF3 apo → RF3), where the apo→bound
difference is the ligand's effect rather than a change of model. ESMFold has no co-fold of its own,
so it has no matched pair — and the perplexity result makes that a reason to retire it rather than
a gap to fill. The superseded set is kept runnable at `morph`'s `--set legacy`.

**What stays true from this lane.** Its folds and their docking are committed and unchanged, and
the collapse it exposed is still what prompted the apo cofolder lanes. Note also that ESMFold is
not alone on every arm: Boltz-2 apo independently extends `orig_f12`, `ox1_orig_f12`,
`ox1_esm1_f8` and `ox2_shuffle` (see the Boltz apo section), so on those the extended state is the
originating model's read too — ESMFold simply extends further, and everywhere, including the
seventeen sequences three other models fold compact.

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

## Apo results (OpenFold3, RF3 and Boltz-2 without the ligand; the first two landed 2026-10-09, Boltz 2026-10-11)

Full data in `runs/<mol>/folds/metrics/of3apo.csv`, `rf3apo.csv` and `boltzapo.csv`. Same format as the cofold
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

### Boltz apo

The same model, and the same local wrapper, that produced every co-fold in this project — asked for
the sequence alone. `apo` below classes the fold as the grid does.

| peptide | mol | kind | Cα RMSD | Rg / co-fold | pLDDT | apo |
|---|---|---|---|---|---|---|
| orig_f12 | oct | designed | 12.56 | 14.72 / 7.46 | 67.2 | ext |
| s2_esm2_control | oct | esm2-variant | 0.37 | 14.55 / 14.50 | 92.4 | same |
| s3_esm2_f4 | oct | esm2-variant | 3.73 | 9.41 / 9.53 | 80.4 | diff |
| s3_orig_f12 | oct | designed | 6.08 | 7.34 / 7.61 | 50.8 | diff |
| shuffle_control | oct | null | 7.24 | 8.46 / 7.23 | 58.1 | diff |
| shuffle_control_esm0 | oct | null-esm | 0.29 | 14.97 / 14.83 | 90.0 | same |
| bg33_1 | oct | boltzgen | 0.38 | 14.46 / 14.40 | 97.7 | same |
| bg33_2 | oct | boltzgen | 0.54 | 14.45 / 14.41 | 98.0 | same |
| bg33_3 | oct | boltzgen | 0.33 | 8.08 / 8.09 | 95.7 | same |
| bg33_4 | oct | boltzgen | 0.61 | 9.21 / 9.47 | 96.8 | same |
| ox1_esm1_f8 | oxy | esm2-variant | 9.02 | 12.50 / 8.22 | 85.5 | ext |
| ox1_orig_f12 | oxy | designed | 8.39 | 11.06 / 6.89 | 82.5 | ext |
| ox2_esm1_f8 | oxy | esm2-variant | 1.94 | 8.83 / 8.85 | 67.2 | diff |
| ox2_orig_f8 | oxy | designed | 7.65 | 9.79 / 7.12 | 60.5 | diff |
| ox2_shuffle | oxy | null | 11.96 | 13.57 / 6.68 | 68.8 | ext |
| ox2_shuffle_esm0 | oxy | null-esm | 3.03 | 11.57 / 12.78 | 84.8 | diff |
| ox3_orig_f8 | oxy | designed | 6.67 | 6.46 / 6.94 | 63.7 | diff |
| bgox31_2 | oxy | boltzgen | 0.36 | 13.62 / 13.62 | 98.2 | same |
| bgox31_3 | oxy | boltzgen | 1.38 | 8.53 / 9.00 | 92.6 | same |
| bgox31_4 | oxy | boltzgen | 1.14 | 7.93 / 8.31 | 96.6 | same |
| bgox31_5 | oxy | boltzgen | 0.93 | 8.92 / 8.97 | 95.0 | same |

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

**Boltz-2 apo splits the difference, and it is the vote that counts most.** The model that
produced every co-fold in this project does *not* uniformly keep its own shells without the
ligand. It extends octinoxate's lead design `orig_f12` outright — Rg 14.72 against the co-fold's
7.46, Cα RMSD 12.56 Å, the campaign's largest apo departure for that peptide — and extends
`ox1_orig_f12` (11.06 / 6.89), `ox1_esm1_f8` (12.50 / 8.22) and the null `ox2_shuffle` (13.57 /
6.68). But it holds `s3_orig_f12` (7.34 / 7.61), `shuffle_control` (8.46 / 7.23), `ox3_orig_f8`
(6.46 / 6.94) and `ox2_esm1_f8` (8.83 / 8.85) compact, the last at 1.94 Å. So four of twenty-one
sequences lose their shell under their own model, and the rest keep a compact state that is still
repacked 6–8 Å away from where the ligand put it.

**This re-opens what the Modal apo lanes closed, for the peptides it names.** The ESMFold collapse
was called model-specific because OF3 and RF3 both folded everything compact; Boltz-2 agreeing
with ESMFold on `orig_f12`, `ox1_orig_f12` and `ox2_shuffle` makes the tally 2–2 on those three,
with the originating model on the extended side. The honest reading: for most of these sequences a
compact apo state exists and the packing is degenerate (unchanged), but for the lead octinoxate
design the shell really does look ligand-dependent — the model that built it will not rebuild it
from sequence alone. Note also that the extended arms carry *higher* confidence than the compact
ones here (pLDDT 82.5 and 85.5 for the two ox1 arms, 67.2 for `orig_f12`, against 50.8–63.7 for the
compact designed arms): Boltz is not hedging on the extension the way OF3 hedges on compactness.

**Boltz apo confirms the BoltzGen arms a fourth time.** All eight are "same" — 0.33–1.38 Å at
pLDDT 92.6–98.2, the highest apo confidence in the campaign — so those folds now reproduce under
four independent apo reads plus two cofolders.

**Apo grid** (`same` = reproduces the co-fold, `diff` = compact but repacked, `ext` =
near-fully-extended apo). Note the apo folds are compared against the *co-fold's* packing, so
`diff` does not mean wrong — with no ligand in the query, no model can know Boltz's packing:

| peptide | mol | kind | ESMFold apo | OF3 apo | RF3 apo | Boltz apo |
|---|---|---|---|---|---|---|
| orig_f12 | oct | designed | ext | diff | diff | ext |
| s2_esm2_control | oct | esm2-variant | same | same | same | same |
| s3_esm2_f4 | oct | esm2-variant | diff | diff | same | diff |
| s3_orig_f12 | oct | designed | ext | diff | diff | diff |
| shuffle_control | oct | null | ext | diff | diff | diff |
| shuffle_control_esm0 | oct | null-esm | same | same | same | same |
| bg33_1 | oct | boltzgen | diff | same | same | same |
| bg33_2 | oct | boltzgen | diff | same | same | same |
| bg33_3 | oct | boltzgen | same | same | same | same |
| bg33_4 | oct | boltzgen | same | same | same | same |
| ox1_esm1_f8 | oxy | esm2-variant | diff | diff | diff | ext |
| ox1_orig_f12 | oxy | designed | ext | diff | diff | ext |
| ox2_esm1_f8 | oxy | esm2-variant | diff | diff | diff | diff |
| ox2_orig_f8 | oxy | designed | ext | diff | diff | diff |
| ox2_shuffle | oxy | null | ext | diff | diff | ext |
| ox2_shuffle_esm0 | oxy | null-esm | same | diff | diff | diff |
| ox3_orig_f8 | oxy | designed | ext | diff | diff | diff |
| bgox31_2 | oxy | boltzgen | same | same | same | same |
| bgox31_3 | oxy | boltzgen | same | same | same | same |
| bgox31_4 | oxy | boltzgen | same | same | same | same |
| bgox31_5 | oxy | boltzgen | same | same | same | same |

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

**Campaign bottom line (six lanes, 21 peptides, all static reads):**

- **BoltzGen's arms are the model-robust structures of the campaign** — sub-Å to 1.4 Å Cα RMSD
  from *both* cofolders and *all three* of their apo lanes, Boltz-2's own included (bg33 and
  bgox31, all "same" in five columns of the grid; ESMFold disagrees only on bg33_1/2). Nothing
  else in the campaign is that consistent.
- **The designed glycine-rich shells are pack-degenerate, and ligand-governed for some arms but
  not as a class.** The ESMFold lane's extended-apo collapse looked purely model-specific while
  only OF3 and RF3 answered — both fold these sequences, nulls included, compact without a ligand.
  Adding the originating model (2026-10-11) qualifies that: Boltz-2 apo extends four of the
  twenty-one, among them octinoxate's lead design `orig_f12` (Rg 14.7 vs 7.5 co-folded), so for
  that peptide the shell does not survive the ligand's removal under the model that produced it.
  For the remaining seventeen the compact apo state stands, and what no lane reproduces is Boltz's
  particular packing — everything sits 5–10 Å away in Cα RMSD, cofolded runs engage the ligand
  *less* than BoltzGen's static read, and the compact apo state is held at low confidence
  (OF3 pLDDT ~50) rather than as a defined fold.
- **Packing degeneracy — not apo fragility — is the property that distinguishes the groups:**
  designed and null glycine-rich arms live at 4–10 Å from every reference in every model, the
  BoltzGen arms are pinned sub-Å, and the esm2-variants sit between, with the two cofolders
  disagreeing on which packing they land in.
- **ESMFold is retired from downstream analysis** (2026-10-11): last of six lanes on
  inverse-folding perplexity (20.4, the no-information limit), owner of 7 of the 11
  zero-information structures, 6 of them with zero contacts. Its rows stay in the tables above as
  the lane that exposed the apo collapse; nothing new should be built on them.
- **Fold metrics alone cannot separate a design from its null.** Designed and plain-null rows look
  alike in every grid cell — the rubric for what to pursue has to come from binding: MM/GBSA,
  ligand retention and contacts against those same null rows.

## Vina docking into the folded structures (all 126, 2026-10-09; the Boltz apo lane 2026-10-11)

Every fold structure was Vina-docked (`code/dock_folds.py`, which stages the fold's peptide as
receptor and calls `vina_redock.redock()` unchanged: 22 Å box, exhaustiveness 16, seed 42, 9
modes, project cutoffs). 21 peptides × 6 lanes = 126 docks; the first five lanes took ~25 min of
Mac Vina and the Boltz apo lane 5 min more. The rows are appended additions-only to the committed
`dock_summary.csv`/`dock_poses.csv` (60 + 66 summary rows, 600 + 660 pose rows; structure names
carry the lane suffix, e.g. `s3_orig_f12_o3cof`, `orig_f12_bzapo`). Full scratch
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

| peptide | mol | kind | o3cof | rf3cof | esmapo | o3apo | rf3apo | bzapo |
|---|---|---|---|---|---|---|---|---|
| orig_f12 | oct | designed | -4.2 / 10.2 | -5.6 / 6.2 | -2.7 / 8.2 | -3.6 / 7.4 | -3.4 / 8.2 | -3.9 / 9.7 |
| s3_orig_f12 | oct | designed | -4.3 / 8.0 | -4.9 / 3.0 | -2.3 / 7.4 | -3.9 / 10.7 | -3.9 / 8.8 | -3.9 / 12.7 |
| ox1_orig_f12 | oxy | designed | -5.7 / 1.2 | -4.6 / 4.6 | -3.4 / 6.8 | -6.1 / 12.4 | -4.3 / 13.2 | -4.1 / 10.9 |
| ox2_orig_f8 | oxy | designed | -6.3 / 6.7 | -4.5 / 9.3 | -3.3 / 9.0 | -5.0 / 8.9 | -4.6 / 5.2 | -5.3 / 7.0 |
| ox3_orig_f8 | oxy | designed | -4.6 / 5.8 | -4.4 / 6.0 | -4.5 / 9.6 | -5.0 / 8.1 | -4.0 / 9.5 | -4.1 / 9.2 |
| s2_esm2_control | oct | esm2-variant | -4.4 / 12.6 | -4.8 / 10.6 | -5.2 / 11.0 | -4.4 / 6.4 | -4.1 / 9.2 | -4.7 / 9.9 |
| s3_esm2_f4 | oct | esm2-variant | -5.1 / 6.9 | -4.5 / 9.3 | -4.9 / 11.2 | -4.9 / 11.0 | -4.1 / 10.8 | -5.5 / 9.0 |
| ox1_esm1_f8 | oxy | esm2-variant | -5.0 / 1.5 | -4.4 / 7.0 | -4.2 / 12.8 | -4.7 / 7.7 | -4.6 / 12.0 | -4.0 / 9.6 |
| ox2_esm1_f8 | oxy | esm2-variant | -6.5 / 5.9 | -5.0 / 13.0 | -4.5 / 10.4 | -5.2 / 8.1 | -5.5 / 7.2 | -5.0 / 8.8 |
| shuffle_control | oct | null | -3.8 / 13.6 | -4.9 / 8.1 | -2.9 / 10.7 | -3.8 / 10.2 | -4.6 / 10.1 | -3.9 / 7.9 |
| ox2_shuffle | oxy | null | -5.5 / 5.9 | -3.6 / 11.2 | -3.4 / 10.7 | -4.6 / 8.9 | -4.4 / 8.6 | -4.2 / 10.8 |
| shuffle_control_esm0 | oct | null-esm | -4.5 / 10.7 | -5.3 / 4.1 | -4.9 / 9.7 | -4.9 / 7.6 | -4.9 / 8.2 | -4.4 / 11.2 |
| ox2_shuffle_esm0 | oxy | null-esm | -5.5 / 10.3 | -5.1 / 10.9 | -5.7 / 13.7 | -4.7 / 12.2 | -5.3 / 8.8 | -5.2 / 7.5 |
| bg33_1 | oct | boltzgen | -3.3 / 7.6 | -3.6 / 10.0 | -3.4 / 9.8 | -3.6 / 10.9 | -3.2 / 9.9 | -3.3 / 7.7 |
| bg33_2 | oct | boltzgen | -3.4 / 10.0 | -2.9 / 12.9 | -4.4 / 11.1 | -2.9 / 8.0 | -3.3 / 9.5 | -2.9 / 13.8 |
| bg33_3 | oct | boltzgen | -4.2 / 9.5 | -4.0 / 10.4 | -4.0 / 12.4 | -4.1 / 9.6 | -4.4 / 5.5 | -4.4 / 10.2 |
| bg33_4 | oct | boltzgen | -4.5 / 10.6 | -4.1 / 4.6 | -4.1 / 9.8 | -4.2 / 5.3 | -4.0 / 12.5 | -4.4 / 12.6 |
| bgox31_2 | oxy | boltzgen | -4.8 / 7.8 | -4.1 / 10.9 | -4.0 / 6.9 | -4.2 / 7.9 | -4.4 / 6.9 | -4.5 / 8.6 |
| bgox31_3 | oxy | boltzgen | -5.2 / 7.2 | -5.1 / 6.4 | -5.3 / 11.1 | -4.7 / 10.9 | -4.8 / 11.7 | -4.8 / 9.7 |
| bgox31_4 | oxy | boltzgen | -4.6 / 7.5 | -4.9 / 6.1 | -4.7 / 9.2 | -4.7 / 10.0 | -4.0 / 13.4 | -5.3 / 10.0 |
| bgox31_5 | oxy | boltzgen | -4.6 / 5.9 | -5.0 / 8.0 | -4.5 / 5.6 | -5.5 / 7.3 | -5.2 / 9.5 | -4.6 / 6.2 |

- **Vina does not reproduce cofolder packing in these shells.** Pose 1 sits 1.2–13.8 Å in-place
  RMSD from the fold's own ligand placement, with a median around 9 Å: in 124 of 126 structures
  the independent pose search finds different packing than the cofolder did. The two agreements —
  `ox1_orig_f12` and `ox1_esm1_f8` in the OF3 cofold lane at 1.16 and 1.51 Å — are the oxybenzone
  structures whose OF3 cofolds sit tightest to Boltz in the cofold tables. This is the
  packing-degeneracy story from the fold geometry, seen from a third direction.
- **But the docked poses stay inside the envelope.** Pose-1 wrapped fraction is 0.77–1.00
  everywhere and engaged counts run 13–20. A Vina search 9 Å from the cofold placement still
  leaves the ligand packed against the chain — these shells are pockets, not grooves, for every
  peptide including the nulls.
- **The Boltz apo lane behaves like the other apo lanes** (added 2026-10-11, same superposed-site
  box since Boltz apo folds carry no ligand): best scores −2.9 to −5.5, pose-1 RMSD 6.3–13.8 Å,
  pose-1 wrapped fraction 0.71–1.00 and engaged 12–20. Its four extended-apo arms do not stand out
  as worse-scoring — `ox1_orig_f12_bzapo` is −4.1 against −4.1 for its compact `ox3_orig_f8`
  sibling — but it gives the lane's lowest wrapped fraction (0.71) and engaged count (12) of all
  126 structures, the expected consequence of docking into a chain that has lost its shell.
- **Scores say nothing separating.** Best scores span −2.3 to −6.6 kcal/mol across all 126, with
  no designed-vs-null gap in any lane. Consistent with "fold metrics cannot separate a design
  from its null"; the expectation is the same for these docks.
- **The 9 MD legs run the OF3-cofold structures' OWN ligand placement** — the cofold pose,
  pose 0 in the dock tables, not a Vina pose. The cofolder's placement is the structure being
  tested; docking it again (2026-10-09 showed Vina moves it a median ~9 Å) would have measured
  Vina's packing instead. `run_dock_pose_md.sh POSE=cofold` stages the leg-ready pose from
  `runs/<mol>/folds/dock_stage/openfold3/<pep>_o3cof_ligand_cofold.sdf` (the fold's perceived
  ligand with template bond orders, template atom order).

## GNINA rescoring of the fold-lane docks (105 structures, all three passes, retrieved 2026-10-12)

The 105 non-ESMFold docks (bundled 2026-10-11, `runs/<mol>/gnina/README.md`) were rescored on the
GNINA side (gnina v1.3.3, CPU) and committed to `gnina_local` under `rescoring/fold_lane/<mol>/`.
Here they sit at `runs/<mol>/gnina/fold_lane/results_<mode>/` — a separate directory beside each
molecule's co-fold `results_<mode>/`, so those committed tables were not rewritten — holding every
structure's `*_scored.sdf` and `.log`, `gnina_version.txt`, and the `gnina_scores.csv` /
`collect.txt` that `collect_gnina.py` writes from them (500 octinoxate and 550 oxybenzone poses per
pass). All 105 are present in all three passes.

**Per structure, `receptor_score_only`: pose-0 CNNaffinity (its rank of 10) · the GNINA-pick docked
pose.** Pose 0 is the cofolder's own ligand placement in the cofold lanes and the superposed OF3-cofold
site pose in the apo lanes. The pick is the max-CNNaffinity pose among Vina's 1–9; `pA/B/C` means
score_only / receptor_h / minimize disagree, and **bold** marks a pick that is the same in all three
passes *and* is the argmin of `minimizedAffinity` in the minimize pass — the full CLAUDE.md criterion,
which only 15 of 105 meet (77 are consistent across the passes on CNNaffinity alone).

| peptide | mol | kind | o3cof | rf3cof | o3apo | rf3apo | bzapo |
|---|---|---|---|---|---|---|---|
| orig_f12 | oct | designed | 4.01 (3) · p1 | 4.24 (6) · p5/5/2 | 1.72 (10) · p8 | 2.44 (10) · p9 | 4.29 (1) · p7/7/8 |
| s2_esm2_control | oct | esm2-variant | 4.27 (5) · p5 | 4.00 (10) · p4/4/9 | 2.55 (10) · p8 | 4.31 (8) · **p2** | 2.71 (10) · p6 |
| s3_esm2_f4 | oct | esm2-variant | 4.67 (2) · p8 | 4.92 (6) · **p3** | 1.95 (10) · p8 | 4.59 (1) · p9 | 1.88 (10) · p4 |
| s3_orig_f12 | oct | designed | 4.74 (1) · p6 | 3.79 (4) · p5 | 1.79 (10) · p3 | 4.85 (1) · p9 | 3.34 (10) · p9/9/5 |
| shuffle_control | oct | null | 4.01 (1) · p8 | 4.30 (1) · p7/7/4 | 1.56 (10) · p1 | 3.91 (10) · p4 | 2.83 (10) · p6 |
| shuffle_control_esm0 | oct | null-esm | 4.82 (1) · **p3** | 4.30 (1) · p3/8/4 | 2.92 (10) · p9 | 2.36 (10) · p7/7/1 | 3.89 (8) · p3 |
| bg33_1 | oct | boltzgen | 4.98 (3) · p6 | 4.38 (7) · p8 | 4.60 (7) · p8 | 1.71 (10) · **p2** | 1.95 (10) · **p7** |
| bg33_2 | oct | boltzgen | 4.81 (2) · **p1** | 5.14 (1) · **p8** | 4.33 (7) · p4/4/5 | 4.33 (7) · **p3** | 1.82 (10) · **p2** |
| bg33_3 | oct | boltzgen | 3.79 (7) · p6 | 4.14 (2) · p5 | 1.28 (10) · **p2** | 1.59 (10) · p2/4/2 | 1.31 (10) · p4 |
| bg33_4 | oct | boltzgen | 5.00 (1) · p7 | 4.53 (4) · **p4** | 4.46 (7) · p2/2/8 | 4.61 (3) · p2 | 2.41 (10) · **p9** |
| ox1_esm1_f8 | oxy | esm2-variant | 3.69 (5) · p2 | 3.79 (6) · p4/5/4 | 1.06 (10) · p5 | 3.21 (10) · p3 | 1.00 (10) · p7 |
| ox1_orig_f12 | oxy | designed | 3.37 (2) · p4 | 3.75 (1) · **p1** | 1.10 (10) · p7 | 1.58 (10) · p4 | 2.77 (9) · p9 |
| ox2_esm1_f8 | oxy | esm2-variant | 3.91 (4) · p3 | 3.90 (4) · p5 | 2.17 (10) · p5/5/7 | 2.26 (10) · p6/8/6 | 2.08 (10) · p6 |
| ox2_orig_f8 | oxy | designed | 3.57 (6) · p2/2/6 | 3.64 (2) · p7 | 1.12 (10) · p9 | 3.45 (2) · p8 | 3.15 (10) · p2/9/7 |
| ox2_shuffle | oxy | null | 3.51 (4) · **p1** | 3.49 (2) · p2 | 2.22 (10) · p9 | 1.12 (10) · p1/1/4 | 1.15 (10) · p5 |
| ox2_shuffle_esm0 | oxy | null-esm | 4.05 (5) · p6/4/6 | 4.04 (9) · p1/1/2 | 2.83 (10) · p4/4/6 | 0.89 (10) · p9 | 3.14 (10) · p4 |
| ox3_orig_f8 | oxy | designed | 3.45 (2) · p4 | 3.59 (2) · p8/8/1 | 1.00 (10) · p8 | 2.12 (10) · p4/4/9 | 2.58 (10) · p7 |
| bgox31_2 | oxy | boltzgen | 4.35 (3) · p1/9/1 | 4.45 (3) · p3/4/3 | 3.15 (10) · p9 | 3.79 (7) · p3/3/8 | 2.42 (10) · p3/9/9 |
| bgox31_3 | oxy | boltzgen | 3.36 (4) · p4 | 4.44 (1) · **p1** | 1.24 (10) · p7 | 1.03 (10) · p5/5/1 | 0.91 (10) · p7 |
| bgox31_4 | oxy | boltzgen | 4.07 (1) · p7 | 3.42 (8) · p3/3/6 | 2.56 (10) · p6 | 3.13 (8) · p7 | 2.19 (10) · p9 |
| bgox31_5 | oxy | boltzgen | 4.01 (7) · p7 | 4.03 (4) · p8 | 1.90 (10) · p8 | 3.14 (1) · p4/5/8 | 2.09 (10) · p2 |

- **The CNN likes the cofolder's placement better than Vina did, but only as given.** In the two
  cofold lanes pose 0 ranks median 3rd of 10 on CNNaffinity under score_only (1st in 10 of 42,
  top three in 22), where Vina's scoring put the predicted pose 6th–9th in every original co-fold.
  After `--minimize` that advantage is gone: cofold pose 0 drops to median 6th (1st in 4 of 42),
  the same as the apo lanes.
- **Apo-lane pose 0 is a clash, not a CNN verdict.** The OF3-cofold pose carried onto an apo fold by
  Cα superposition overlaps the apo chain: its Vina-term affinity under score_only is positive in
  20/21 `o3apo`, 15/21 `rf3apo` and 21/21 `bzapo` structures (medians +92, +14, +90 kcal/mol, against
  −1.4 and −1.1 for the cofold lanes), so it ranks 10th of 10 almost everywhere (median 10, 1st in 4
  of 63). Minimizing relieves every clash but one each in `o3apo` and `rf3apo` and lifts apo pose-0 CNNaffinity from
  2.3–2.9 to 3.8, level with the cofold lanes' 3.9. Read apo pose-0 ranks as the site definition
  failing to fit, nothing more.
- **Best-pose CNNaffinity does not separate cofold from apo shells, or designs from nulls.** The mean
  of each structure's best CNNaffinity is 4.42 over the cofold lanes and 4.31 over the apo lanes;
  in the OF3 cofold lane designs + esm-variants average 4.31, nulls 4.26, BoltzGen 4.59; across the five
  lanes the design-vs-null gap is within 0.2 either way (designs ahead in two lanes, behind in three). BoltzGen's shells score highest in every lane by a few tenths,
  and octinoxate runs ~0.6 above oxybenzone throughout (a ligand-size effect on an affinity
  regressor). Same conclusion as Vina's scores and the inexpensive fold metrics: a static score
  of a static shell does not tell a design from its null.
- **No MD has been queued from these picks** (2026-10-12). Should GNINA-pick legs be wanted, the
  bold cells are the ones that meet the full pick criterion; for the 13 `_o3cof` structures with
  cofold-pose MD legs, a pick leg would be `run_dock_pose_md.sh STRUCT=<pep>_o3cof POSE=<K>`.

**Parser fix found doing this.** `collect_gnina.py` read a blank tag value (pose 0 has no
`vina_score`) by skipping to the next non-blank line, which is the next tag, so it swallowed
`minimizedAffinity` and every pose 0 came out with that column empty. Fixed 2026-10-12 and the
fold-lane tables built with the fix. The previously committed co-fold `gnina_scores.csv` tables
(octinoxate, `bg33_1_bg33_2`, oxybenzone) still carry the blank — their CNN columns are unaffected,
and they have not been regenerated.

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

**How a cofold leg stages, end to end** (the campaign legs and the four nulls both went through
this): `code/build_cofolds.sh` runs the `BUILD_ONLY=1` loop over its leg list
(`run_dock_pose_md.sh STRUCT=<leg> POSE=cofold SYSNAME=<mol>`). A leg whose protein was never
staged by a campaign build fails at `omd build` with a missing `protein_fixed.pdb`; prepping it
needs the *protonated* receptor from the docking grid —
`omd prep-protein --pdb runs/<mol>/dock/<leg>/<leg>_protein.pdb --out runs/<mol>/md/<leg>/protein_fixed.pdb`
(`--out` takes the FILE path; a directory gives `IsADirectoryError`). The dock-dir protein is
coordinate-identical to the campaign legs' `protein_fixed.pdb` (verified 0.0000 Å), so this is the
same prep, just done per leg. Then push → probe → produce on Modal (`code/modal_md.py`), and the
local tail is `run_dock_pose_md.sh` again with a `prod_L1_modal → prod_20ns` symlink, which the
build guard skips straight past to the analysis.

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
| `shuffle_control_o3cof` | OF3 refold of shuffle_control (null) | cofold's own placement | −9.38 ± 0.11 | −9.51 / −11.60 / −12.91 | 7.59 | 7.9 | released late: built to −12.9 by 15 ns then gave 3.5 back; 83% of releases in the back half, slide 7.59 Å |
| `shuffle_control_esm0_o3cof` | OF3 refold of shuffle_control_esm0 (null) | cofold's own placement | **−15.62 ± 0.09** | −12.46 / −14.55 / −15.05 | 10.62 | 4.2 | **null counterexample 2**: stronger than its Boltz null baseline (−14.32) by 1.3 and retention higher (92.5% vs 80.2%); ligand pulls IN 11.3→4.5 Å, 0% of releases late |
| `ox2_shuffle_o3cof` | OF3 refold of ox2_shuffle (null) | cofold's own placement | −2.98 ± 0.09 | −7.55 / −4.39 / −3.28 | 1.49 | 5.9 | shell opens: slide only 1.49 Å yet separation 8.9→15.4 Å, residence 14.8%, longest 2.9 ns detached |
| `ox2_shuffle_esm0_o3cof` | OF3 refold of ox2_shuffle_esm0 (null) | cofold's own placement | −0.92 ± 0.05 | −1.53 / −0.82 / −0.55 | 13.14 | 2.4 | effectively unbound: OF3 places the ligand 18.6 Å out (9.3 ns detached) while the Boltz baseline of the same null held at 94.7% |

## The cofold campaign: the headline tables, one set per molecule

Each molecule carries its **own** set of five tables below, the same five and in the same order as
`RESULTS.md` gives each ligand, rather than one pooled set with a `mol` column. The RF3 re-folds are
a different fold lane, not a second pose within this one, and are reported separately in **The RF3
re-folds** at the end of this file.

> **Status: data collection, not analysis.** Read everything below as a record of numbers as they
> arrived. Every cell is **one 20 ns trajectory with no replicate**; the ± beside each ΔG is a
> within-trajectory standard error and is **not** uncertainty on the difference between two legs,
> where run-to-run variation is larger than many of the gaps compared here. A signed difference
> between a re-fold leg and its Boltz baseline is written "+N weaker" or "−N stronger" and means only
> that two runs of two structures gave different numbers — nothing is "lost" or "gained" by the
> re-folding. 21 of the campaign's 29 legs were still unrun as of 2026-10-14, so every ordering,
> grouping and apparent predictor here is an **interpretation of the data in hand, revisable by the
> next leg** — not a rule the remaining legs are expected to satisfy. Predicting from a reading is
> how it gets tested and is encouraged; what is not allowed is letting the reading decide what a leg
> means. Where a claim written mid-collection has since failed, it is marked withdrawn in place
> rather than quietly deleted, so the record shows what was believed when.
>
> **Read the slide column before reading any cell.** Equilibration restrains protein heavy atoms
> only, so the ligand moves before production starts, and `runs/<mol>/md/ligand_slide.csv` records
> how far. Three re-fold cells moved so far that production did not begin from the predicted pose
> and they are **not tests of that pose**: `ox2_shuffle_esm0_o3cof` (13.14 Å),
> `shuffle_control_esm0_o3cof` (10.62 Å) and `shuffle_control_o3cof` (7.59 Å) — all three nulls,
> which is itself expected, since a shuffled sequence has no site to hold the ligand through
> equilibration. Seven more sit in a 2.5–3.8 Å grey band and are flagged in place. Those tables
> were written and discussed for weeks before the slide was tabulated for the octinoxate lane at
> all, which is how a 10.62 Å cell came to be cited as a counterexample; see the correction below.
>
> **Baselines slide too, and one of them badly.** The Δ column compares a re-fold leg to a *Boltz
> co-fold leg*, and that leg has its own equilibration slide. Twenty of the twenty-one baselines
> are under 3.5 Å, but **`s3_esm2_f4` slid 7.81 Å**, so its −21.08 is the ΔG of a pose nearly 8 Å
> from the Boltz prediction. Every Δ for that peptide — OF3's "+7.3 weaker, best oct refold" and
> the RF3 cell — is therefore **not** a measurement of re-folding against the Boltz pose, and that
> peptide's rows are marked ⚠. Found 2026-10-14 by checking the baselines after the re-fold cells;
> I had only checked the re-folds.

The campaign question, answered first across the ten Modal legs and now being extended on racc:
**is the OF3 re-fold worth running dynamics on, or does the Boltz co-fold remain the start the
pipeline should use?** Of the twenty-two OF3 cells that have run, eighteen read 2–18 kcal/mol
*weaker* than their Boltz baseline with degraded or lost retention. Spend: ten Modal legs at ≈ $8.2 of
produce (plus the $0.78 erroneous dock-pose leg), one at a time; the fold-lane extension runs free
on racc.

**The counterexamples, and a correction forced by the slide table (2026-10-14).** Four cells read
stronger than the model that predicted their baseline: `ox1_esm1_f8` (a design, oxybenzone, −7.2),
`shuffle_control_esm0` (a **null**, octinoxate, −1.3), `bg33_1` and `bg33_2` (BoltzGen, octinoxate,
−3.6 each). **Only three of them are readable as statements about the OF3 pose.**

`ligand_slide.py` had never been written to a table for the octinoxate lane, and back-filling it
shows `shuffle_control_esm0_o3cof`'s ligand moved **10.62 Å** relative to its peptide during
equilibration. Production therefore did not start from the OF3 pose, and per this repo's own rule —
record the slide, it is what says whether a cell is worth reading — that cell **cannot support the
claim that the OF3 pose beat the Boltz model**. It remains a real −15.62 at 92.5% residence for
whatever pose it equilibrated into; it is not evidence about OF3's prediction. The earlier gloss
here, that "one of the four is a null, so a strong OF3-cofold ΔG is legible on any sequence,
designed or shuffled", **rested on that cell and is withdrawn.**

What survives is narrower and no longer has a null in it: the three readable counterexamples are
`bg33_1` (slide 1.10 Å), `bg33_2` (1.33 Å) and `ox1_esm1_f8` (3.00 Å, itself borderline), i.e. two
BoltzGen shells on octinoxate and one design on oxybenzone, and **still none is a shell design on
octinoxate**.

**The four BoltzGen cells come out in baseline order — a reading to carry into the RF3 legs, not a
settled one.**
`bg33_1` (−11.31) and `bg33_2` (−8.02), the two weakest co-folds in the octinoxate set, each gained
3.6 kcal/mol on re-folding. `bg33_3` has a *better* baseline (−11.86) and its re-fold reads **+1.8
weaker** with retention falling 41.0% → 14.9%. That was the stated test of whether the difference tracked
the arm or the weak baselines, set before the leg ran, and it came out on the side of the baselines: the
re-fold helps where the Boltz co-fold had already failed to enclose the ligand and hurts where it
had not. So the arm-level claim is not supportable; the most that can be said from four
single trajectories is that **the weaker Boltz co-folds here happened to be the ones whose re-folds
read stronger** — which is the same pattern the oxybenzone set shows, where the charged shell `ox1` is the
one whose re-fold improves.

`bg33_4` was then run as a pre-registered test of that statement: at −19.66 it is the strongest
BoltzGen baseline, so the reading predicted its re-fold would lose the most of the four. It does —
**+2.55 weaker, the largest difference in the arm** — and the four legs rank in baseline order with no
inversions: −8.02 → −3.6, −11.31 → −3.6, −11.86 → +1.8, −19.66 → +2.55. Monotone in four points is
weak evidence on its own, but the direction was called in advance rather than fitted afterwards.

**What `bg33_4` costs the simple story is the "re-folding breaks the site" reading.** Its weaker ΔG
is *not* accompanied by a retention failure: residence is **93.0%** with **zero** released frames in 20 ns, the only
clean-zero release in the octinoxate set, and the fold barely moves (Cα 2.0 Å final, Rg 9.5 → 9.2).
It is also the most *engaged* leg in the arm by a wide margin — 1.93 side chains on the ligand on
average against 0.99–1.22 for the other three, ≥2 in 65% of frames against 24–34%, and 8 pairs held
≥5% of the run against 1–4. So the re-fold keeps the ligand, keeps its shape and makes **more**
contact than any of its peers, and still reads 2.55 kcal/mol weaker. The difference is in the quality of the
enclosure the Boltz pose had (enclosed 0.595 → 0.52, engaged 19 → 14), not in whether the complex
survives. Like `bg33_1` it is **still strengthening at 20 ns** (−14.90 → −17.11 across its windows,
+2.21 and not converged), so −17.11 is an upper bound and the true gap to −19.66 is smaller than
+2.55 — which is the one caveat that could flatten the baseline-order ranking above.

**The four cofold nulls (2026-10-10, $3.90 of produce)** asked the campaign question again on
shuffle-control sequences. The split is not "nulls release" — it is **compact vs elongated**. The
two compact nulls repeat the design-refold pattern in different directions: `shuffle_control` reads
+5.7 weaker with 83%-late release, and `ox2_shuffle`'s shell opens under it (slide 1.5 Å yet
separation 8.9 → 15.4 Å). The two elongated nulls split instead: `ox2_shuffle_esm0`'s OF3 placement
misses the site entirely (−0.92, 18.6 Å, 9.3 ns detached) against a Boltz baseline that held at
94.7%, while `shuffle_control_esm0` pulls the ligand *in* (11.3 → 4.5 Å) and becomes the null
counterexample above. The `ox2_shuffle` pairing also shows the baseline itself can bind nowhere
(7.1% residence, −5.30) on a structure whose static geometry is as good as the designs'.

### Octinoxate — co-folds, ΔG and retention: Boltz pose vs OF3 re-fold

| peptide | Boltz co-fold ΔG (ret.) | OF3 re-fold ΔG (ret.) | Δ (refold − Boltz) |
|---|---|---|---|
| `orig_f12` | **−13.71 (77.3%)** | −5.90 (7.4%) | +7.8 weaker, released |
| `orig_f12` (Vina p1 in the OF3 shell) | −13.71 (77.3%) | −9.65 (50.2%) | +4.1 weaker, partial loss |
| `s3_orig_f12` | **−24.33 (100.0%)** | −10.10 (44.8%) | +14.2 weaker, partial loss |
| `s2_esm2_control` | **−16.25 (59.6%)** | −10.21 (33.9%) | +6.0 weaker, one long exit-and-return |
| `s3_esm2_f4` | −21.08 (99.6%) ⚠ | −13.76 (84.8%) | +7.3 weaker but held (best oct refold); **Δ compromised — baseline slid 7.81 Å** |
| `shuffle_control` (null) | **−15.13 (75.9%)** | −9.39 (25.7%) | +5.7 weaker, late progressive release |
| `shuffle_control_esm0` (null) | −14.32 (80.2%) | **−15.62 (92.5%)** | **−1.3 stronger on the OF3 fold, retention higher** |
| `bg33_1` (boltzgen) | −11.31 (14.5%) | **−14.90 (56.4%)** | **−3.6 stronger, retention 4× higher** |
| `bg33_2` (boltzgen) | −8.02 (50.9%) | **−11.58 (70.2%)** | **−3.6 stronger, retention higher** |
| `bg33_3` (boltzgen) | −11.86 (41.0%) | −10.10 (14.9%) | +1.8 weaker, retention down |
| `bg33_4` (boltzgen) | **−19.66 (100.0%)** | −17.11 (93.0%) | +2.55 weaker but held, no release |

Eight of the eleven octinoxate cells are weaker, and **no designed shell re-fold on this ligand is
stronger** — the three that beat their baseline are a null (`shuffle_control_esm0`) and two of the
four BoltzGen designs, which are now complete (`bg33_1`, `bg33_2`, against `bg33_3` at +1.8 and
`bg33_4` at +2.55 weaker). `s3_esm2_f4` is the best-behaved shell re-fold
here: +7.3 weaker but it held the ligand at 84.8%, where `orig_f12`'s own re-fold released at 7.4%.
The Boltz co-fold stays the default start for octinoxate's shell designs.

**The two BoltzGen legs gain the same 3.6 kcal/mol by different routes**, which is why they are
worth reading side by side rather than pooling. `bg33_1_o3cof` (−14.90 against −11.31, retention
14.5% → 56.4%) gets there by **folding up**: the extended rod compacts over the run (Rg 14.9 → 9.6
Å) and pulls the ligand in with it (centroid separation 12.1 → 6.9 Å), the same behaviour the
`bg33_1`/`bg33_2` GNINA-docked legs showed in `RESULTS.md`, and it is **still strengthening at
20 ns** (+2.58 across its windows), so −14.90 is an upper bound rather than a converged value.
`bg33_2_o3cof` (−11.58 against −8.02, retention 50.9% → 70.2%) does the opposite: it **holds the
shape it was given** — Cα RMSD 1.4 Å final, Rg flat at 14.9 Å — and its windows are level from 5 ns
(−11.01 → −11.58, +0.57), so that number *is* converged. One leg improved by rearranging, the other
without rearranging at all.

A shared caution: both beat the weakest Boltz baselines in the set, so this is a low bar cleared.
**`bg33_3` is the third route and the negative one** — +1.8 weaker, retention 41.0% → 14.9%, the
ligand drifting out rather than in (contacts 20.4 → 12.2, centroid separation 9.5 → 13.4 Å) while
the peptide itself barely moves at all (Cα 1.1 Å final, 1.7 Å max, the most rigid fold in the
campaign). It is also the only one of the three whose input→production slide is flagged, at 2.63 Å.
Taken with the other two it says the re-fold's effect on this arm follows the baseline it starts
from, not the arm.

### Octinoxate — co-folds, window convergence

| leg | 5 ns | 10 ns | 15 ns | 20 ns |
|---|---|---|---|---|
| `orig_f12_o3cof_dock1` | −9.13 | −10.86 | −11.28 | −9.65 |
| `orig_f12_o3cof` | −9.62 | −8.31 | −7.34 | **−5.90** |
| `s3_orig_f12_o3cof` | −14.04 | −10.92 | −10.51 | −10.10 |
| `s2_esm2_control_o3cof` | −7.32 | −8.52 | −9.86 | −10.21 |
| `s3_esm2_f4_o3cof` | −13.58 | −14.61 | −14.52 | −13.76 |
| `shuffle_control_o3cof` | −9.51 | −11.60 | −12.91 | **−9.39** |
| `shuffle_control_esm0_o3cof` | −12.46 | −14.55 | −15.05 | **−15.62** |
| `bg33_1_o3cof` | −12.32 | −12.64 | −13.50 | **−14.90** |
| `bg33_2_o3cof` | −11.01 | −11.34 | −11.62 | **−11.58** |
| `bg33_3_o3cof` | −11.27 | −10.26 | −10.05 | −10.10 |
| `bg33_4_o3cof` | −14.90 | −16.42 | −16.73 | **−17.11** |

Every octinoxate release leg shows the decay signature — `orig_f12_o3cof` gives back ~4 over its
windows and `shuffle_control_o3cof` ~3.5 after peaking at 15 ns — so those endpoints are poses
coming apart, not settled ensembles. `s3_esm2_f4_o3cof` holds a flat −13.6 to −14.6 band the whole
run, which is what a settled pose looks like.

`bg33_1_o3cof` is the only leg in either molecule's set that **strengthens monotonically** (+2.58
across the four windows) and it has not levelled off by 20 ns, so its endpoint is an upper bound:
the fold is still closing on the ligand when the run stops. `bg33_2_o3cof` is its opposite and the
flattest series in the campaign (−11.01 / −11.34 / −11.62 / −11.58, a 0.6 kcal/mol spread): a 5 ns
window would have reported its 20 ns answer to within half a kcal/mol. The two legs gain the same
3.6 kcal/mol over their baselines, so **window shape, not the ΔG, is what distinguishes a leg that
has settled from one still moving** — and only one of these two has. `bg33_3_o3cof` sits between
them and decays mildly (−11.27 → −10.10, giving back 1.2 over the run, flat after 10 ns), the
signature of a starting pose that was slightly over-packed rather than one coming apart.

### Octinoxate — co-folds, residence and release

| leg | residence | released frames | episodes | longest | late |
|---|---|---|---|---|---|
| `orig_f12_o3cof_dock1` | 50.2% | 86 (4.3%) | 31 | 230 ps | 92% |
| `orig_f12_o3cof` | 7.4% | 442 (22.1%) | 53 | 2140 ps | 86% |
| `s3_orig_f12_o3cof` | 44.8% | 28 (1.4%) | 14 | 80 ps | 50% |
| `s2_esm2_control_o3cof` | 33.9% | 276 (13.8%) | 8 | **2680 ps** | 1% |
| `s3_esm2_f4_o3cof` | 84.8% | 7 (0.35%) | 5 | 30 ps | 14% |
| `shuffle_control_o3cof` | 25.7% | 84 (4.2%) | 18 | 430 ps | 83% |
| `shuffle_control_esm0_o3cof` | 92.5% | 49 (2.5%) | 11 | 160 ps | 0% |
| `bg33_1_o3cof` | 56.4% | 1 (0.05%) | 1 | 10 ps | 0% |
| `bg33_2_o3cof` | 70.2% | 7 (0.35%) | 4 | 40 ps | 0% |
| `bg33_3_o3cof` | 14.9% | 25 (1.25%) | 11 | 110 ps | 64% |
| `bg33_4_o3cof` | **93.0%** | 0 (0.0%) | 0 | — | — |

`s3_esm2_f4_o3cof` (7 released frames in 20 ns) and `shuffle_control_esm0_o3cof` (49, none late) are
the only octinoxate retainers. `s2_esm2_control_o3cof`'s single 2.68 ns excursion at 1% late is a
mid-run visit that returned — a different failure from a progressive walk-out, and the reason this
table carries the episode structure rather than residence alone.

`bg33_1_o3cof` is the case this table exists to catch, in the opposite direction from the rest: its
residence reads a middling 56.4% but it has **one** released frame in the whole run, so the ligand
never actually leaves. It sits just outside the 10 Å centroid cutoff while the rod is still extended
and comes inside as the fold compacts. Residence and the release record disagree here, and the
release record is the one describing the physics — the same pairing flagged the surface-roamer mode
for `bgox31_5_dock1` in `RESULTS.md`.

`bg33_2_o3cof` is the same story with a milder reading: 70.2% residence, 7 released frames in 4
episodes of at most 40 ps, none of it late. Both BoltzGen legs are **held ligands whose residence
percentage understates them**, because on an elongated peptide the centroid is far from the binding
site, and a 10 Å centroid cutoff is the wrong ruler for a rod. That is a measurement caveat to
state once in the paper rather than a result about these two designs.

`bg33_3_o3cof` is the one of the three where low residence **does** mean loss: 14.9% with 25 released
frames in 11 episodes and **64% of the detached time in the second half**, against 0% late for the
other two. Progressive drift, not a rod geometry artefact — and the `late` column is what separates
the two cases, which is the argument for carrying it.

### Octinoxate — co-folds, peptide structural stability

| leg | Cα final | Cα max | Rg start→end |
|---|---|---|---|
| `orig_f12_o3cof_dock1` | 4.6 | 6.0 | 8.8 → 8.6 |
| `orig_f12_o3cof` | 5.5 | 6.5 | 9.0 → 8.7 |
| `s3_orig_f12_o3cof` | 3.9 | 4.6 | 8.4 → 8.7 |
| `s2_esm2_control_o3cof` | 2.5 | 3.8 | 15.3 → 14.7 |
| `s3_esm2_f4_o3cof` | 1.9 | 3.9 | 15.2 → 15.4 |
| `shuffle_control_o3cof` | 7.9 | 9.2 | 8.8 → 9.4 |
| `shuffle_control_esm0_o3cof` | 4.2 | 5.3 | 15.3 → 14.9 |
| `bg33_1_o3cof` | 11.1 | 12.1 | 14.9 → 9.6 |
| `bg33_2_o3cof` | 1.4 | 6.6 | 14.9 → 14.9 |
| `bg33_3_o3cof` | 1.1 | 1.7 | 8.8 → 8.7 |
| `bg33_4_o3cof` | 2.0 | 3.7 | 9.5 → 9.2 |

No shell-pipeline re-fold's Cα drift exceeds the Boltz legs' own band; failure in that group is
ligand retention, not fold collapse. The long-Rg legs (`s2_esm2_control`, `s3_esm2_f4` at ~15 Å) are
the elongated designs, and their shells stayed the same size through dynamics.

`bg33_1_o3cof` breaks that pattern and is the largest fold change in the campaign: Cα RMSD 11.1 Å
final with Rg dropping 14.9 → 9.6 Å. This is **not** the loosening this table was built to detect —
it is the reverse, an extended rod folding up around the ligand, and it is why that leg's ΔG
improves rather than decays. The caveat for reading it against the shell rows: a structure that
rearranges this much is no longer the pose that was predicted, so its ΔG belongs to the fold the
dynamics found, not to the OF3 re-fold as placed.

`bg33_2_o3cof` is the control that caveat needed, and it arrived in the same batch: **the steadiest
leg in the campaign** (Cα 1.4 Å final, Rg 14.9 → 14.9 Å) and it gains the same 3.6 kcal/mol. So the
BoltzGen improvement does **not** require the rearrangement — `bg33_1`'s compaction is one route to
it, not the explanation for it, and `bg33_2`'s number is attributable to the OF3 pose exactly as
placed.

`bg33_3_o3cof` is more rigid still (Cα 1.1 Å final, 1.7 Å max, Rg 8.8 → 8.7 Å) and its ΔG got
*worse*. Across the three legs, fold stability and ΔG change are therefore **uncorrelated**: the
most mobile leg gained, the steadiest gained the same amount, and the next-steadiest lost. This
table rules out fold collapse as an explanation for any of the three, and offers nothing beyond
that.

### Octinoxate — the pre-dynamics static read and the pairs the re-folds realise

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
| `shuffle_control` | Boltz co-fold | 0.75 | 0.75 | 15 | 7.9 Å | | 18/66 | 13 |
| | OF3 re-fold (null) | 0.51 | 0.90 | 18 | 9.2 Å | 1.70 | 33/66 | 40 |
| `shuffle_control_esm0` | Boltz co-fold | 0.54 | 0.85 | 17 | 11.0 Å | | 12/66 | 2 |
| | OF3 re-fold (null) | 0.31 | 0.30 | 6 | 9.4 Å | 1.26 | 4/66 | −13 |
| `bg33_1` (boltzgen) | Boltz co-fold | 0.38 | 0.75 | 15 | 17.3 Å | | 8/55 | |
| | OF3 re-fold | 0.415 | 0.65 | 13 | 9.5 Å | 1.22 | 15/55 | 11 |
| `bg33_2` (boltzgen) | Boltz co-fold | 0.585 | 0.75 | 15 | 5.2 Å | | 8/55 | |
| | OF3 re-fold | 0.35 | 0.40 | 8 | 8.0 Å | 0.99 | 12/55 | 4 |
| `bg33_3` (boltzgen) | Boltz co-fold | 0.48 | 1.00 | 20 | 8.6 Å | | 3/55 | |
| | OF3 re-fold | 0.485 | 0.65 | 13 | 9.4 Å | 1.17 | 14/55 | 9 |
| `bg33_4` (boltzgen) | Boltz co-fold | 0.595 | 0.95 | 19 | 6.3 Å | | 14/55 | |
| | OF3 re-fold | 0.52 | 0.70 | 14 | 7.0 Å | **1.93** | 14/55 | 9 |

`bg33_1` is the one row where the OF3 re-fold reads *better* than the Boltz co-fold on the static
half — enclosure essentially unchanged (0.38 → 0.415) but the ligand centroid pulled from **17.3 Å
to 9.5 Å**, i.e. the Boltz co-fold had placed the ligand outside the fold altogether. That is the
pre-dynamics counterpart of the ΔG result, and it is consistent with the rest of the table rather
than an exception to it: the drop that predicts a weak leg is enclosure *and* separation together,
and here only one of the two is poor. Its pairs come from `design.json`'s 11 non-glycine slots
(band [10, 55], the same convention as the committed `bg33_3`/`bg33_4` rows): 15/55 realised, 2
slots never in contact, no long-range pair —
`runs/octinoxate/md/pair_contacts_o3cof_boltzgen.csv`.

**`bg33_2` is the clearest failure of the static read anywhere in the campaign, and it should be
reported as such.** Its OF3 re-fold is *worse* than its Boltz co-fold on every static column —
enclosure 0.585 → **0.35**, wrapped 0.75 → 0.40, engaged 15 → **8**, centroid separation 5.2 → 8.0 Å
— and it nonetheless binds **3.6 kcal/mol stronger** with higher retention and a converged window
series. A read that predicts weak legs from low enclosure gets this one exactly backwards. Both
BoltzGen legs also realise fewer pairs than any shell design (15/55 and 12/55, 3–4 dead slots, no
long-range pair in either), so the pair count does not rescue the prediction either. The honest
statement for the paper is that the static read prices *shell* re-folds well and has no demonstrated
skill on the BoltzGen arm.

`bg33_3` makes that concrete from the other direction: its static read is the one that barely moves
on re-folding (enclosure 0.48 → 0.485, separation 8.6 → 9.4 Å, though wrapped drops 1.00 → 0.65 and
engaged 20 → 13) and it is the leg whose ΔG got *worse*. So across the four BoltzGen legs the
static read moved sharply down and the leg improved (`bg33_2`), stayed flat and the leg worsened
(`bg33_3`), improved on separation alone and the leg improved (`bg33_1`), and moved moderately down
with the leg moderately worse (`bg33_4`, enclosure 0.595 → 0.52, engaged 19 → 14, ΔG +2.55) —
**one of four points the right way, and the three before it in three different directions**. Notably
`bg33_3` realises 14/55 pairs with the arm's only **long-range pairs** (3/10) despite being the
weakest-binding of the four, and `bg33_4` realises the same 14/55 while binding the *strongest* at
−17.11, which tells against reading the pair count as a quality measure on this arm too: the same
pair total spans 7 kcal/mol here.

What does separate `bg33_4` from the rest is the *simultaneity* column rather than any static
one — 1.93 side chains on the ligand at once against 0.99–1.22 — and that is a dynamics
measurement, available only after the leg has run. So the arm's one usable quality signal is not in
the pre-dynamics read at all, which is the finding to report rather than a repaired static rule.

### Oxybenzone — co-folds, ΔG and retention: Boltz pose vs OF3 re-fold

| peptide | Boltz co-fold ΔG (ret.) | OF3 re-fold ΔG (ret.) | Δ (refold − Boltz) |
|---|---|---|---|
| `ox1_orig_f12` | −10.31 (49.1%) | −3.96 (20.5%) | +6.4 weaker, bouncing |
| `ox2_orig_f8` | **−24.21 (100.0%)** | −7.23 (49.5%) | +17.0 weaker, walked out at 10 ns |
| `ox3_orig_f8` | **−18.84 (100.0%)** | −2.06 (2.4%) | +16.8 weaker, never bound |
| `ox1_esm1_f8` | −10.70 (75.0%) | **−17.95 (100.0%)** | **−7.2 stronger, held all run** |
| `ox2_esm1_f8` | **−20.81 (100.0%)** | −16.68 (93.2%) | +4.1 weaker but held, 0 releases |
| `ox2_shuffle` (null) | −5.30 (7.1%) | −2.98 (14.8%) | +2.3 weaker; both legs poor, the baseline binds nowhere either |
| `ox2_shuffle_esm0` (null) | **−12.04 (94.7%)** | −0.92 (1.1%) | +11.1 weaker; **slide 13.14 Å — not a test of the OF3 pose** |
| `bgox31_2` (boltzgen) | −10.89 (33.7%) | −7.45 (65.9%) | +3.45 weaker **yet retention doubled**; one 2.0 ns exit and return |
| `bgox31_3` (boltzgen) | **−15.44 (99.7%)** | −5.82 (15.8%) | **+9.62 weaker, retention collapsed** 99.7% → 15.8% |
| `bgox31_4` (boltzgen) | **−14.05 (100.0%)** | −5.74 (17.9%) | +8.31 weaker, retention 100% → 17.9%; bound first, then lost it (slid 3.83 Å) |
| `bgox31_5` (boltzgen) | **−16.78 (100.0%)** | −11.55 (70.7%) | +5.23 weaker but **held, 0 releases**; strongest baseline, second-smallest difference |

Ten of the eleven oxybenzone cells are weaker, and this is where the campaign's **largest** differences
sit: `ox2_orig_f8` +17.0 and `ox3_orig_f8` +16.8, the two best Boltz co-folds on this molecule,
both of which lose the ligand under the re-fold (49.5% and 2.4%). The one cell that goes the other
way, `ox1_esm1_f8` at −7.2 and 100% retention, is the campaign's strongest counterexample and the
only one that is a designed shell.

Read against the ΔG table in `RESULTS.md`, the pattern is that the re-fold costs most where the
Boltz co-fold was best: the charged shell `ox1`, weak and releasing at baseline, is the one whose
re-fold improves.

**`bgox31_2_o3cof` is the first oxybenzone BoltzGen cell, and it splits ΔG from retention — they
move in opposite directions.** The re-fold reads **+3.45 weaker** (−7.45 against −10.89) while
residence **doubles**, 33.7% → 65.9%. Both halves of that are real and neither is the summary. The
ΔG decays in the windows (−9.88 → −7.45, with the whole drop between 10 and 15 ns) and the contacts
trace dates the cause exactly: centroid separation runs 7.6–9.3 Å for the first 10 ns, blows out to
**28.7 Å at 14 ns**, then comes back to 6.9 Å by 20 ns. So this is one 2.04 ns excursion and return,
98% of the detached time late, inside a leg whose fold never moves (Cα 2.4 Å, Rg flat at 13.6 Å) and
whose ligand barely slid at equilibration (1.52 Å).

The reason residence still *rises* is the baseline, not the re-fold: `bgox31_2` at 33.7% was already
a poor retainer whose own detachment was 92% late. Both legs lose the ligand late; the re-fold simply
loses it for less of the run but pays more in energy while it is gone, because a 28.7 Å excursion
contributes near-zero binding to the mean where the baseline's shallower wandering did not. **This is
the cleanest case in the campaign for not collapsing a leg to one number** — the pair of (ΔG,
residence) disagrees in sign, and the `late` column plus the separation trace are what make it
legible. It also means the octinoxate rule that a weak Boltz baseline has more to gain does *not*
transfer here: `bgox31_2` is a weak baseline (−10.89, 33.7%) and its re-fold still lost 3.45.

**`bgox31_3_o3cof` goes the other way within the same arm.** `bgox31_3` is the arm's strong baseline
(−15.44 at **99.7%** residence, enclosure **0.955** with all seventeen ligand atoms engaged at 4.0 Å
centroid separation). Its re-fold reads **−5.82, i.e. 9.62 kcal/mol weaker, with residence at
15.8%** — the third-largest difference on this molecule after `ox2_orig_f8` and `ox3_orig_f8`, which
are also the two strongest Boltz co-folds here. The direction was noted before the tail finished, on
the two orderings this file had been tracking; it came out that way, which makes three cells
consistent with the reading. Stating the direction first is what made that worth anything; it is
still three single trajectories.

Two consequences for how this arm gets reported. First, **`bgox31_2`'s sign disagreement does not
generalise**: here ΔG and retention move together, both sharply down, so the opposite-sign case
remains a single well-characterised cell rather than a property of the BoltzGen arm. Second, this cell has the arm's largest static
change (enclosure 0.955 → 0.615, separation 4.0 → 8.5 Å) alongside its largest ΔG difference. I wrote
that up as evidence the static read "works on oxybenzone where it failed on octinoxate"; **`bgox31_5`
later contradicted that and the claim is withdrawn** — see the OUTCOME note below, where six
candidate orderings over four single trajectories are recorded as bookkeeping rather than signal.

**`bgox31_4_o3cof` tests a hypothesis the static read raises and refutes it.** Its OF3 pose starts at
**11.2 Å** centroid separation — the worst in the arm and already *outside* the 10 Å residence
cutoff before dynamics begin — which invites the reading that the leg never bound at all and its
17.9% residence is an artefact of where the pose was placed. The contacts trace refutes that: the
ligand pulls **in** from 11.65 Å to **8.78 Å** by 6 ns with contacts rising 14.3 → 18.9, so it did
bind, and only then opens to 23.1 Å at 14 ns before a partial return (11.7 Å, 19.4 contacts at
20 ns). With **90% late** detachment this is genuine progressive release, and the residence figure
reports a real loss rather than a placement artefact. The peptide meanwhile is the most rigid in the
whole campaign (Cα **0.7 Å** final, 1.9 Å max, Rg 8.8 → 9.1).

Its ΔG, **−5.74 against −14.05 (+8.31 weaker)**, lands between `bgox31_2` and `bgox31_3` where the
two orderings jointly suggested: strong on residence at baseline (100%, zero releases) but only mid-range on
enclosure (0.65, against `bgox31_3`'s 0.955), so it had less enclosure to lose and lost somewhat less
energy. **One caveat belongs on this cell specifically**: `ligand_slide.py` flags it at **3.83 Å**,
the largest equilibration slide in the arm, so production started from a pose appreciably displaced
from the OF3 prediction — the leg is evidence about *a* pose in that shell, and less directly about
the predicted one, than its three peers are.

The mechanism on `bgox31_3` is progressive loss with a partial return, not a clean release: contacts run 22.2 → 0.0
by 12 ns as separation opens 9.1 → **23.6 Å**, then recover to ~12.6 Å and ~13 contacts by 20 ns, for
476 released frames in 31 episodes (61% late). The fold is the most rigid in the oxybenzone set
(Cα 1.5 Å final, 2.2 Å max, Rg 8.6 → 8.5 Å) and the ligand slid 2.11 Å at equilibration, so once
again nothing structural fails — the peptide holds its shape and simply stops holding the ligand.

### Oxybenzone — co-folds, window convergence

| leg | 5 ns | 10 ns | 15 ns | 20 ns |
|---|---|---|---|---|
| `ox1_orig_f12_o3cof` | −7.60 | −5.17 | −4.29 | −3.96 |
| `ox2_orig_f8_o3cof` | −14.87 | −11.98 | −8.10 | **−7.23** |
| `ox3_orig_f8_o3cof` | −0.43 | −1.44 | −1.91 | −2.06 |
| `ox1_esm1_f8_o3cof` | −19.54 | −19.23 | −18.92 | −17.95 |
| `ox2_esm1_f8_o3cof` | −19.08 | −17.54 | −15.70 | −16.68 |
| `ox2_shuffle_o3cof` | −7.55 | −4.39 | −3.28 | −2.98 |
| `ox2_shuffle_esm0_o3cof` | −1.53 | −0.82 | −0.55 | **−0.92** |
| `bgox31_2_o3cof` | −9.88 | −9.68 | −7.50 | −7.45 |
| `bgox31_3_o3cof` | −9.15 | −6.44 | −5.20 | −5.82 |
| `bgox31_4_o3cof` | −8.80 | −7.64 | −5.70 | −5.74 |
| `bgox31_5_o3cof` | −11.35 | −9.998 | −11.19 | **−11.55** |

`ox2_orig_f8_o3cof` gives back ~7.6 across its windows and `ox1_orig_f12_o3cof` ~3.6 — the decay
signature of a starting pose coming apart. The two esm1 re-folds hold a flat −16 to −19 band the
whole run, so their numbers are settled poses rather than starting-pose artifacts.
`ox3_orig_f8_o3cof` is flat at −0.4 to −2.1, which is not convergence but a leg that never bound.

### Oxybenzone — co-folds, residence and release

| leg | residence | released frames | episodes | longest | late |
|---|---|---|---|---|---|
| `ox1_orig_f12_o3cof` | 20.5% | 777 (38.9%) | 40 | 1800 ps | 66% |
| `ox2_orig_f8_o3cof` | 49.5% | 624 (31.2%) | 16 | **3010 ps** | 99% |
| `ox3_orig_f8_o3cof` | 2.4% | 1059 (53.0%) | 80 | 1470 ps | 35% |
| `ox1_esm1_f8_o3cof` | 100.0% | 0 | 0 | — | — |
| `ox2_esm1_f8_o3cof` | 93.2% | 0 | 0 | — | — |
| `ox2_shuffle_o3cof` | 14.8% | 1092 (54.6%) | 34 | **2900 ps** | 54% |
| `ox2_shuffle_esm0_o3cof` | 1.1% | 1667 (83.4%) | 25 | **9270 ps** | 52% |
| `bgox31_2_o3cof` | 65.9% | 347 (17.35%) | 18 | 2040 ps | **98%** |
| `bgox31_3_o3cof` | 15.8% | 476 (23.8%) | 31 | **2360 ps** | 61% |
| `bgox31_4_o3cof` | 17.9% | 414 (20.7%) | 39 | 1430 ps | **90%** |
| `bgox31_5_o3cof` | 70.7% | **0 (0.0%)** | 0 | — | — |

The two esm1 re-folds are the only genuine retainers in this set, both with **zero** released
frames. At the other end, `ox2_shuffle_esm0_o3cof` spends a single **9.27 ns** episode detached and
`ox2_orig_f8_o3cof` a 3.01 ns one at 99% late — progressive walk-outs, not flicker, and their ΔG is
computed partly on an unbound ligand. `ox3_orig_f8_o3cof`'s 80 episodes at 2.4% residence is the
opposite failure: a ligand that never settled anywhere.

### Oxybenzone — co-folds, peptide structural stability

| leg | Cα final | Cα max | Rg start→end |
|---|---|---|---|
| `ox1_orig_f12_o3cof` | 2.9 | 3.8 | 8.6 → 7.8 |
| `ox2_orig_f8_o3cof` | 4.1 | 4.8 | 9.0 → 8.1 |
| `ox3_orig_f8_o3cof` | 5.3 | 6.2 | 7.8 → 8.6 |
| `ox1_esm1_f8_o3cof` | 4.2 | 4.7 | 9.5 → 9.6 |
| `ox2_esm1_f8_o3cof` | 4.0 | 5.7 | 10.0 → 10.4 |
| `ox2_shuffle_o3cof` | 5.9 | 6.6 | 8.0 → 8.6 |
| `ox2_shuffle_esm0_o3cof` | 2.4 | 5.2 | 14.1 → 14.1 |
| `bgox31_2_o3cof` | 2.4 | 5.1 | 13.6 → 13.6 |
| `bgox31_3_o3cof` | 1.5 | 2.2 | 8.6 → 8.5 |
| `bgox31_4_o3cof` | 0.7 | 1.9 | 8.8 → 9.1 |
| `bgox31_5_o3cof` | 2.8 | 3.6 | 9.3 → 9.3 |

No oxybenzone re-fold unravels either: Cα drift is 2.4–5.9 Å with Rg flat to within 1 Å in every
leg, including the two that lost the ligand outright. Fold stability and ligand retention are
decoupled on this molecule — which is the point this table makes for both sets, and the reason a
weak ΔG here cannot be explained away as a peptide falling apart.

### Oxybenzone — the pre-dynamics static read and the pairs the re-folds realise

| peptide | pose | enclosed | wrapped | engaged | cent sep | sim | pairs | band% |
|---|---|---|---|---|---|---|---|---|
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
| `ox2_shuffle` | Boltz co-fold | 0.58 | 0.82 | 14 | 9.3 Å | 1.21 | 29/55 | 42 |
| | OF3 re-fold (null) | 0.58 | 0.94 | 16 | 7.9 Å | 0.59 | 7/55 | −7 |
| `ox2_shuffle_esm0` | Boltz co-fold | 0.64 | 0.94 | 16 | 4.2 Å | 2.21 | 18/55 | 18 |
| | OF3 re-fold (null) | 0.47 | 0.77 | 13 | 7.3 Å | 0.13 | 3/55 | −16 |
| | OF3 re-fold (null) | 0.58 | 0.94 | 16 | 7.9 Å | 0.59 | 7/55 | −7 |
| `bgox31_2` (boltzgen) | Boltz co-fold | 0.44 | 0.77 | 13 | 9.1 Å | | — | — |
| | OF3 re-fold | 0.36 | 0.77 | 13 | 8.1 Å | — | — | — |
| `bgox31_3` (boltzgen) | Boltz co-fold | **0.955** | 1.00 | 17 | **4.0 Å** | | — | — |
| | OF3 re-fold | 0.615 | 0.82 | 14 | 8.5 Å | — | — | — |
| `bgox31_4` (boltzgen) | Boltz co-fold | 0.65 | 0.94 | 16 | 6.7 Å | | — | — |
| | OF3 re-fold | 0.47 | 0.65 | 11 | **11.2 Å** | — | — | — |
| `bgox31_5` (boltzgen) | Boltz co-fold | 0.50 | 0.77 | 13 | 7.0 Å | | — | — |
| | OF3 re-fold | **0.55** | 0.71 | 12 | **7.0 Å** | — | — | — |

The static read priced these re-folds before the GPU did: the Boltz co-folds enclose 0.58–0.995 at
1.2–9.3 Å, the OF3 placements enclose 0.47–0.81 at 5.3–10.3 Å, and every leg that read weaker in the
ΔG table sat on the low half of that drop. The high-energy counterexample, `ox2_esm1_f8`, is also
the placement that matched its Boltz pose best statically (0.745 enclosed, 17 atoms engaged, ligand
*inside* its Boltz centroid distance at 5.3 Å) — and in dynamics it realises the campaign's best
pair band (19/55) with the longest `held` list of the ten Modal legs. Re-fold quality is legible in
the placement's enclosure fraction at zero dynamics cost.

### How to read both sets

Columns and provenance are the same in each molecule's set. ΔG is 20 ns MM/GBSA in kcal/mol with
ligand retention (residence within 10 Å of the peptide centroid) in brackets. Window columns are
leading 5/10/15 ns slices of the same trajectory, the 20 ns column being the headline value.
Residence columns come from `md_contacts.py` at 10 ps sampling, `late` being the share of released
frames in the second half. Stability is Cα RMSD final/max and Rg start→end in Å, from
`protein_stability.csv`. The static half of the fifth table is `runs/<mol>/dock/dock_poses.csv` pose
0 — **enclosed** (fraction of directions out of the ligand that meet peptide), **wrapped** (fraction
of ligand atoms with peptide within 4.5 Å), engaged atoms, centroid separation — and its trajectory
half is `pair_contacts.py` on each re-fold leg: mean simultaneous designed side chains on the ligand
(`sim`) and realised pairs against the [n−1, n(n−1)/2] band. Pair outputs persist as
`runs/{octinoxate,oxybenzone}/md/pair_contacts_o3cof_shellN.csv`, with the `orig_f12` pair under
`pair_contacts_o3cof_design.csv`, the nulls under `pair_contacts_o3cof_nulls.csv` and the BoltzGen
legs under `pair_contacts_o3cof_boltzgen.csv`.

**Pre-registered prediction for `bgox31_5_o3cof`, written 2026-10-14 while the leg was at 7.4M of
10M steps.** This is the first cell in the arm where the file's two working rules *disagree*, so it
discriminates between them rather than adding another confirmation.

- `bgox31_5` is the arm's **strongest and only non-decaying baseline**: −16.78 at 100.0% retention,
  with windows that *strengthen* (−15.62 → −16.45 → −16.59 → −16.78) where every other bgox31
  baseline decays. The rule that a good co-fold has more to lose therefore predicts the **largest
  loss in the arm**, beyond `bgox31_3`'s +9.62.
- But its static read is the **flattest transition anywhere in the campaign**: enclosure 0.50 → 0.55
  (it rises), centroid separation 7.0 → 7.0 Å (unchanged), engaged 13 → 12. The rule that the
  enclosure-and-separation drop predicts the loss therefore predicts **little or no loss**.

Both cannot hold. If the loss is large, the baseline's own quality is doing the work and the static
read adds nothing even on oxybenzone, which would withdraw the molecule-level claim made above for
`bgox31_3`. If the loss is small, the static read is the better predictor and `bgox31_3`'s big loss
should be attributed to its enclosure collapse (0.955 → 0.615) rather than to its strong baseline.

**OUTCOME (recorded after the leg landed, 2026-10-14 14:55).** `bgox31_5_o3cof` read **−11.55**,
i.e. 5.23 kcal/mol weaker than its Boltz baseline, with residence 100% → 70.7% and **zero released
frames**, and windows flat to strengthening (−11.35 / −10.00 / −11.19 / −11.55). The expectation
stated above — that the arm's strongest baseline would show the largest difference — was wrong: this
is the second-smallest difference in the arm.

**How to hold this.** The orderings below are a *reading* of the data in hand, offered as such and
revisable by the next leg — not a rule the remaining legs are expected to obey. The distinction
matters because it decides what happens when a leg disagrees: a reading gets revised, an imposed rule
makes the disagreeing leg look like an anomaly. Each cell is **one 20 ns trajectory with no
replicate**, and the ±0.02–0.05 beside each ΔG is a within-trajectory standard error, not uncertainty
on a difference between two runs of two structures, where run-to-run spread is larger than several of
the gaps compared here. So these are worth stating and worth predicting from, at four cells per arm,
provided the prediction is a test and not a commitment:

| ordered by | cells, in that order | signed ΔΔG, in that order |
|---|---|---|
| oxybenzone, baseline strength first | 5, 3, 4, 2 | +5.23, +9.62, +8.31, +3.45 |
| oxybenzone, enclosure change first | 3, 4, 2, 5 | +9.62, +8.31, +3.45, +5.23 |
| oxybenzone, separation change first | 3, 4, 5, 2 | +9.62, +8.31, +5.23, +3.45 |
| octinoxate, baseline strength first | 4, 3, 1, 2 | +2.55, +1.80, −3.60, −3.60 |
| octinoxate, enclosure change first | 2, 4, 3, 1 | −3.60, +2.55, +1.80, −3.60 |
| octinoxate, separation change first | 2, 3, 4, 1 | −3.60, +1.80, +2.55, −3.60 |

Two of the six come out monotone — separation on oxybenzone, baseline strength on octinoxate — and
the reason to be careful is that six orderings over four points will often hand you one or two by
arithmetic alone, so the monotonicity is suggestive rather than evidential. Taken as a reading: the
two arms' numbers do not line up the same way, and separation change is the more promising thing to
watch on oxybenzone. What **is** withdrawn is the stronger form I committed earlier — that the static
read "carries signal on oxybenzone where it carried none on octinoxate" — because `bgox31_5`
contradicted it and because four single trajectories cannot carry a claim about a molecule. The
21 RF3 legs are the test; predictions made from this reading should be recorded before each leg
lands, as `bg33_4`'s and `bgox31_5`'s were, which is what made them informative when one of them
failed.

The caveat on this cell: `ligand_slide.py` flags **2.68 Å**, so like `bgox31_4` it is slid, though far
less. Its separation trace excursions to 12.7–13.4 Å at 6–8 ns without a single *released* frame also
confirm what `bg33_4` showed — the release detector counts sustained detachment, so brief excursions
beyond 10 Å depress residence without registering an episode, and residence and release have to be
read together.

The original octinoxate precedent favoured the first branch — `bg33_3` also had a near-flat static read and
its leg got *worse* — so the stated expectation is a **large loss with retention falling from
100%**, and a small loss would be the more informative outcome because it breaks the octinoxate
pattern.

**All eight new OF3 BoltzGen legs have landed** (2026-10-14): octinoxate's `bg33_1/2/3/4_o3cof` and
oxybenzone's `bgox31_2/3/4/5_o3cof`. The two OF3 table sets are therefore complete for this
campaign; what remains is the RF3 lane below. Rows are added as each leg's drain completes, not in a batch at the
end.

**The oxybenzone *BoltzGen* rows carry no pair data — a gap in the inputs, and it applies to the
`bgox31` arm only.** The shell designs (`ox1`/`ox2`/`ox3`) all have a design json and do get pair
numbers; on 2026-10-14 I wrongly extended this caveat to `ox1_esm1_f8_rf3cof` and left its pair
columns dashed, which was my omission and not a missing input. Corrected the same day. `pair_contacts.py` needs a design json naming the designed slots, and oxybenzone has
`design_shell1/2/3.json` only — there is no BoltzGen equivalent of octinoxate's `design.json` for
the `bgox31` arm. The committed baseline rows in `RESULTS.md` already show `—` in those columns for
the same reason. Deriving slots from the sequence's non-glycine positions would be inventing a
design definition, so the columns stay dashed until a real `bgox31` design json exists; the
consequence is that the arm's one useful quality signal on octinoxate — simultaneity — cannot be
read on oxybenzone at all.

## The RF3 re-folds: two further sets of five tables

RoseTTAFold3 co-folds of the same 21 sequences, kept **separate from the OF3 sets above** rather
than added as further columns: a different predictor is a different fold lane, and the question it
answers is whether the Boltz co-fold survives re-folding *in general* or only under OF3. Each
molecule gets the same five tables it gets for OF3, so this file carries **twenty** campaign tables
in all — five per fold lane per molecule. Column definitions and provenance are as in
[How to read both sets](#how-to-read-both-sets) above; pair outputs will persist as
`runs/{octinoxate,oxybenzone}/md/pair_contacts_rf3cof_*.csv`.

**Status: no RF3 leg has completed yet, so all ten tables are headers only.** Twenty-one legs are
built and queued on racc — **ten octinoxate** (`fold_oct_2`, job 512262 submitted 2026-10-14, plus
`fold_oct_3` and `fold_oct_4`) and **eleven oxybenzone** (`fold_oxb_1` through `fold_oxb_3`). Rows
land as each leg's drain completes, not in a batch at the end, and the Boltz-baseline column is
copied from the same source as the OF3 sets so the two lanes stay directly comparable.

### Octinoxate — RF3 co-folds, ΔG and retention

Ten legs, each as `<pep>_rf3cof`: `orig_f12`, `s3_orig_f12`, `s2_esm2_control`, `s3_esm2_f4`,
`shuffle_control`, `shuffle_control_esm0`, `bg33_1`–`bg33_4`.

**Prediction for `orig_f12_rf3cof`, written 2026-10-14 16:45 with the leg drained and its tail
still running, so before any ΔG was read.** Recorded as a reading under test, not a rule.

`orig_f12`'s OF3 cell is the campaign's worst retention failure: −13.71 (77.3%) → −5.90 (7.4%),
ligand effectively gone, windows decaying −9.62 → −5.90. Two readings in this file disagree about
what RF3 should do with the same sequence:

- **Sequence-level** — a peptide behaves similarly across re-folders, which `ox1_esm1_f8` supports
  by reading stronger under both OF3 (−7.2) and RF3 (−1.42). This predicts `orig_f12_rf3cof` also
  does badly.
- **Pose-level** — the static read says RF3 barely disturbed this site, where OF3 wrecked it:
  enclosure 0.96 → **0.80** against OF3's 0.68, engaged 20 → **18** against 14, centroid separation
  4.2 → **5.8 Å** against 9.4 Å. That is the gentlest re-fold transition in the octinoxate set and
  predicts it holds the ligand far better than 7.4%, with ΔG nearer −13.71 than −5.90.

**Stated expectation: the pose-level one** — it holds, well above 7.4%. The reason for preferring
it is that the static read's failures on octinoxate were all on `bg33` legs, BoltzGen shells with
weak Boltz baselines, whereas this is a designed shell with a near-perfect Boltz pose, a different
regime. That is itself a reading of four cells and is the kind I have over-trusted today, so if
the leg goes the other way the sequence-level reading gains and this note should not be used to
explain the result away.

**OUTCOME (16:50, after the tail).** The pose-level expectation was right in direction and
**understated in size**: `orig_f12_rf3cof` reads **−19.51 at 99.9% residence with zero released
frames**, i.e. 5.79 kcal/mol *stronger* than its own Boltz baseline. The prediction said "nearer
−13.71 than −5.90"; it went past −13.71 entirely, which the prediction did not anticipate.

**The same sequence and the same ligand, re-folded by two predictors, differ by 13.6 kcal/mol and by
92 points of residence** — OF3 −5.90 at 7.4%, RF3 −19.51 at 99.9%. That is the largest
predictor-to-predictor gap in the campaign and it is the clearest evidence so far that **what a
re-fold leg reports is a property of the pose it starts from, not of the sequence**: nothing about
`orig_f12` itself changed between the two legs. The pair data says the same thing — the RF3 pose
realises **35/66** pairs with **no slot ever out of contact**, against the Boltz co-fold's own 22/66
and OF3's 12/66, so RF3 did not merely preserve the site, it presented more of the designed surface
to the ligand than the pose the design was selected on.

This weakens the sequence-level reading that `ox1_esm1_f8`'s double gain had supported; one peptide
behaving alike under both predictors now looks like coincidence rather than a property. It does not
license the reverse claim that the static read is reliable — one cell, and the four `bg33` cells where
it pointed four ways are unchanged. The useful form is narrower: **a large static change between
Boltz pose and re-fold pose has, so far, gone with a large ΔG change in the same direction whenever
the Boltz pose was good** (`orig_f12` here, `bgox31_3`, `ox2_orig_f8`, `ox3_orig_f8`), and carried no
direction at all on the weak-baseline `bg33` shells. Nine cells now sit on one side of that split and
four on the other, which is a reading to take into the remaining 20 legs, not a result.

One caveat on this cell: the peptide moves more than any other RF3 leg so far (Cα **9.2 Å** final,
12.4 Å max, Rg 8.8 → 9.9), so the 99.9% residence is a ligand held by a fold that is itself
rearranging, and the windows decay mildly (−21.63 → −19.51). Slide 1.29 Å.

| peptide | Boltz co-fold ΔG (ret.) | RF3 re-fold ΔG (ret.) | Δ (refold − Boltz) |
|---|---|---|---|
| `orig_f12` | −13.71 (77.3%) | **−19.51 (99.9%)** | **−5.79 stronger, 0 releases**; OF3's same cell read −5.90 at 7.4% |
| `s2_esm2_control` | **−16.25 (59.6%)** | −14.98 (57.2%) | +1.27 weaker, residence held; OF3's cell read −10.21 at 33.9% |
| `s3_esm2_f4` | −21.08 (99.6%) ⚠ | −16.18 (88.2%) | +4.90 weaker, 0 releases; **Δ compromised — baseline slid 7.81 Å**; OF3's cell read −13.76 |
| `s3_orig_f12` | **−24.33 (100.0%)** | −12.10 (68.9%) | +12.2 weaker; OF3's cell read −10.10 at 44.8%, so RF3 better by only 2.0 |

**Prediction for `s2_esm2_control_rf3cof`, written 2026-10-14 17:58 with the leg drained and its
tail running, so before any ΔG was read.** Unlike `orig_f12`, where RF3 was far gentler than OF3,
here **both re-folders put the ligand outside the fold** and the pose-level reading predicts a
weak cell, plausibly no better than OF3's −10.21:

| `s2_esm2_control` pose | enclosed | wrapped | engaged | cent sep | peptide Rg |
|---|---|---|---|---|---|
| Boltz co-fold | 0.51 | 0.85 | 17 | 6.7 Å | — |
| OF3 re-fold | 0.38 | 0.50 | 10 | 13.5 Å | 15.3 → 14.7 (extended) |
| RF3 re-fold | 0.42 | 0.65 | 13 | **13.5 Å** | **10.2 (compact)** |

Separation rises 6.7 → 13.5 Å, the largest increase in the octinoxate set, and RF3 only marginally
beats OF3 on enclosure (0.42 vs 0.38) and engagement (13 vs 10). **The sharper point is that the
13.5 Å means something different here.** On the OF3 leg this peptide was extended (Rg 15.3 Å), so
a far centroid was partly the rod artefact this file warns about — the centroid of a rod is not
where its site is. RF3 folds it **compact** at Rg 10.2 Å, so 13.5 Å of centroid separation now
means the ligand really is outside the fold rather than off the end of one. Stated expectation:
weak, low residence, and **not** rescued by the cutoff-artefact explanation that covered
`bg33_4`, `bgox31_5` and `ox1_esm1_f8_rf3cof`. If instead it binds well from 13.5 Å, that is the
second time a pose starting outside the cutoff pulled in (after `bgox31_4`) and the static
separation column loses much of its remaining value.

**OUTCOME (18:19, after the tail): wrong on both counts, and the stated consequence applies.**
`s2_esm2_control_rf3cof` reads **−14.98 at 57.2% residence with 2 released frames** — 1.27 kcal/mol
weaker than its −16.25 baseline, where OF3's cell on the same peptide read −10.21 at 33.9%. I
predicted weak and low-residence; it is neither. Residence essentially matches the baseline's 59.6%.

The ligand **pulled in from 13.5 Å**: separation 11.08 → **5.84 Å** by 12 ns with contacts 14.2 →
**32.4**, before drifting back to 11.59 Å by 20 ns. So this is the second pose to start outside the
10 Å cutoff and bind anyway, after `bgox31_4` — and because I named that consequence in advance, it
is applied rather than argued around: **the static centroid-separation column does not support
predicting a leg's outcome.** It failed here on a compact fold, where the rod-artefact excuse that
covered the OF3 version of this peptide was explicitly unavailable. Together with its four-way split
on the `bg33` legs, what survives of the static read is description of a pose, not prediction of a
leg.

Two further things this cell settles. The **pose-level reading from `orig_f12` is weakened**: there
RF3's gentle static transition went with a strong leg, which looked like the pose explaining the
result, but here RF3's static read is as bad as OF3's and the leg is 4.8 kcal/mol *better* than OF3's,
so the static pose quality is not what separated them. And the **cutoff-artefact explanation is now
load-bearing in four cells** (`bg33_4`, `bgox31_5`, `ox1_esm1_f8_rf3cof`, and the in-and-out trace
here), which is frequent enough that residence below ~60% should be treated as uninterpretable
without its separation trace rather than as evidence of release.

Slide 2.14 Å. Pairs 15/66 with `design_shell2.json`, sim 1.58, 4 dead slots, 7/15 long-range —
against `orig_f12_rf3cof`'s 35/66 and 0 dead slots, so the stronger-binding leg of the two is also
the one engaging far more of its designed surface. Cα 1.0 Å final, Rg 10.5 → 10.8: the fold barely
moves while the ligand does all the travelling.

**Pre-registration for the six legs drained 2026-10-15 06:36–08:00, written before any ΔG was
read.** One note for all six rather than six notes, and deliberately thin: every mechanistic
reading built during yesterday's collection was falsified by a later leg, so there is no model
here worth defending. What follows is what the static read *would* say, recorded only so it can be
scored rather than quietly forgotten.

| leg | Boltz | OF3 cell | RF3 static (enc / eng / sep) | what the static read implies |
|---|---|---|---|---|
| `shuffle_control_rf3cof` | −15.13 (75.9%) | −9.39 (25.7%) | **0.83 / 18 / 6.3 Å** | best RF3 read in the set, far better than OF3's — a strong cell |
| `shuffle_control_esm0_rf3cof` | −14.32 (80.2%) | −15.62 (92.5%) | 0.52 / 13 / 8.5 Å | worse than OF3's; OF3 won this one, so RF3 should lose |
| `bg33_1_rf3cof` | −11.31 (14.5%) | −14.90 (56.4%) | **0.30 / 7 / 8.1 Å** | worst read in the campaign — a weak cell |
| `bg33_2_rf3cof` | −8.02 (50.9%) | −11.58 (70.2%) | **0.345 / 4 / 9.5 Å** | only 4 of 20 atoms engaged — weaker still |
| `ox2_shuffle_esm0_rf3cof` | −12.04 (94.7%) | −0.92 (1.1%) | 0.505 / 16 / 10.3 Å | better than OF3's, from the campaign's lowest floor — a gain |
| `ox2_shuffle_rf3cof` | −5.30 (7.1%) | −2.98 (14.8%) | 0.48 / 12 / 12.1 Å | both legs poor at baseline; no strong expectation |

**The one thing worth scoring is the `bg33` pair.** On octinoxate the static read pointed four
different ways across the four OF3 `bg33` cells, which is why it was withdrawn as predictive.
Here it is unusually emphatic — 0.30 and 0.345 enclosed, 7 and 4 atoms engaged — and octinoxate is
also the molecule where RF3 has beaten OF3 4 of 4. **Those two expectations conflict**: the static
read says these should be the weakest RF3 cells yet, the molecule split says they should still beat
OF3's −14.90 and −11.58. Both cannot hold, and the molecule split is the only reading left
standing from yesterday, so this is its first real test.

### Octinoxate — RF3 co-folds, window convergence

| leg | 5 ns | 10 ns | 15 ns | 20 ns |
|---|---|---|---|---|
| `orig_f12_rf3cof` | −21.63 | −22.03 | −19.99 | **−19.51** |
| `s2_esm2_control_rf3cof` | −12.43 | −13.84 | −15.75 | **−14.98** |
| `s3_esm2_f4_rf3cof` | −17.93 | −17.98 | −17.39 | **−16.18** |
| `s3_orig_f12_rf3cof` | −11.64 | −12.85 | −12.04 | **−12.10** |
| `shuffle_control_esm0_rf3cof` | −19.48 | −21.48 | −22.25 | **−21.53** |
| `shuffle_control_rf3cof` | −12.90 | −10.74 | −10.95 | **−11.13** |
| `bg33_1_rf3cof` | −9.41 | −11.20 | −12.05 | **−12.43** |
| `bg33_2_rf3cof` | −11.09 | −11.00 | −10.94 | **−11.52** |
| `bg33_3_rf3cof` | −11.55 | −12.18 | −11.96 | **−12.44** |
| `bg33_4_rf3cof` | −16.79 | −17.81 | −16.06 | **−16.13** |

**Prediction for `s3_esm2_f4_rf3cof`, written 2026-10-14 19:06 — leg drained, tail running, no ΔG
read — and offered with low confidence, for a reason worth stating.** Both readings I would
normally use were damaged in the preceding hour: centroid separation was withdrawn as predictive
by `s2_esm2_control_rf3cof`, and the pose-level reading lost its support from the same cell. What
remains is the enclosure-change split, and on that: `s3_esm2_f4` has a good Boltz pose (**0.795**
enclosed at 7.0 Å) and RF3 changes it substantially (**0.46** at 9.3 Å, with the fold compact at
Rg 10.5 against OF3's extended 15.2), which predicts a **substantially weaker cell, perhaps near
OF3's −13.76**.

The asymmetry matters more than the guess: a weaker reading is only weak evidence for that split
given how much else failed today, whereas a **stronger** reading would be the third strike
against the static read in one evening.

**And this cell cannot settle much either way**, because its baseline slid 7.81 Å (see the note at
the head of this section). The ΔG will be sound; the Δ against −21.08 will not be a measurement
of re-folding against the Boltz pose. Recorded so that a convenient result here is not read as
support it cannot provide.

**OUTCOME (19:29): direction right, magnitude wrong, and the pre-registered limit applies.**
`s3_esm2_f4_rf3cof` reads **−16.18 at 88.2% residence with zero released frames** — 4.90 weaker than
its (compromised) baseline, so weaker as predicted, but **2.4 kcal/mol better than OF3's −13.76**
where I had guessed it would land near it. Windows decay mildly (−17.93 → −16.18) and separation
only starts drifting in the last 2 ns (10.2 Å at 18 ns, 12.3 Å at 20 ns). Slide 1.87 Å, Cα 1.2 Å,
Rg flat at 10.5 Å. Pairs 19/66, 5 dead slots.

Per the asymmetry recorded above, a weaker cell is **weak** support for the enclosure-change split,
and the 7.81 Å baseline means even that is mostly unavailable here. So this cell is logged and
little is claimed from it.

**The comparison that does not need a baseline.** RF3 against OF3 on the *same peptide* involves two
re-fold legs and no Boltz co-fold leg at all, so the baseline-slide problem does not touch it. On
octinoxate's three shell designs RF3 wins every time, by a margin that shrinks as OF3's own cell
improves:

| peptide | OF3 re-fold | RF3 re-fold | RF3 − OF3 |
|---|---|---|---|
| `orig_f12` | −5.90 (7.4%) | −19.51 (99.9%) | **−13.6** |
| `s2_esm2_control` | −10.21 (33.9%) | −14.98 (57.2%) | **−4.8** |
| `s3_esm2_f4` | −13.76 (84.8%) | −16.18 (88.2%) | **−2.4** |
| `ox1_esm1_f8` (oxybenzone) | −17.95 (100.0%) | −12.12 (35.2%) | **+5.8** |

RF3's three octinoxate cells span 4.5 kcal/mol (−14.98 to −19.51) where OF3's span 7.9 (−5.90 to
−13.76), i.e. **RF3 is both better and less variable on this ligand** — and the one oxybenzone pair
so far goes the other way. Read as a reading, not a result: four pairs, one per peptide, no
replicates, and the ordering could be coincidence. But it is the only comparison in this file that
is immune to both the baseline slide and the Boltz-pose question, which makes it the one to watch as
the remaining RF3 legs land.

**Against `pepgraph`, this comparison points the other way — and the two have not been reconciled.**
Inverse-folding perplexity lives in the paper outline, not here: see
`$ODLOC/publications/peptidebuilder/paper_outline.md` §3.5, **Table 5** (six lanes, 202 structures)
and **Table 5b** (paired across the ligand boundary), sourced from `pepgraph/RESULTS.md` §1. Do not
restate those tables in this file — one authoritative copy, or the two drift.

What they say is that **OF3 is the better co-folder on sequence–structure self-consistency**: px mean
10.7 against RF3's 12.2, lowest perplexity on 12 of 21 sequences against the next lane's 3, and on
the paired test that holds sequences fixed, `openfold3 vs of3apo` gains +1.4 px at **+3.6σ** where
`rf3 vs rf3apo` gains +2.6 at only +1.7σ.

What the table above says is that **RF3 is the better co-folder on dynamics**, on all three
octinoxate shell designs run so far. **The lane that wins on self-consistency is not the lane that
wins on binding.** Both halves are worth stating because they are measured on the same 21 sequences
by independent means, and the pepgraph side is the stronger of the two — 21 paired sequences against
four one-per-peptide dynamics pairs with no replicates. Per-peptide the two measures agree on only
2 of those 4 (`orig_f12` and `ox1_esm1_f8` agree; `s2_esm2_control` and `s3_esm2_f4` do not), and
the sharpest single cell is `s2_esm2_control`, whose RF3 backbone scores **px 21.60 — at or past the
no-information limit — and whose leg still binds 4.8 kcal/mol better than OF3's px 13.84 backbone
for the same sequence**.

This is recorded as an unresolved tension for the paper to address, not as a result either way, and
it is the reason the RF3-vs-OF3 column above is worth watching as the remaining legs land.

### Octinoxate — RF3 co-folds, residence and release

| leg | residence | released frames | episodes | longest | late |
|---|---|---|---|---|---|
| `orig_f12_rf3cof` | **99.9%** | **0 (0.0%)** | 0 | — | — |
| `s2_esm2_control_rf3cof` | 57.2% | 2 (0.1%) | 2 | 10 ps | 50% |
| `s3_esm2_f4_rf3cof` | 88.2% | **0 (0.0%)** | 0 | — | — |
| `s3_orig_f12_rf3cof` | 68.9% | 2 (0.1%) | 2 | 10 ps | 50% |
| `shuffle_control_esm0_rf3cof` | **99.3%** | **0 (0.0%)** | 0 | — | — |
| `shuffle_control_rf3cof` | 38.7% | 13 (0.65%) | 11 | 30 ps | 8% |
| `bg33_1_rf3cof` | 68.3% | 42 (2.1%) | 12 | 290 ps | 2% |
| `bg33_2_rf3cof` | 58.7% | 5 (0.25%) | 5 | 10 ps | 0% |
| `bg33_3_rf3cof` | **1.9%** | 5 (0.25%) | 5 | 10 ps | 100% |
| `bg33_4_rf3cof` | 81.2% | **0 (0.0%)** | 0 | — | — |

**Prediction for `s3_orig_f12_rf3cof`, written 2026-10-14 19:55 — production ~6 min from done,
tail not started, no ΔG in existence.** This is the sharpest test of the pepgraph-vs-dynamics
tension in the campaign, because the two measures disagree *strongly* on this one peptide and
its baseline slide is a clean 2.34 Å, so unlike `s3_esm2_f4` its Δ is interpretable.

| pose | pepgraph px | enclosed | engaged | cent sep | ΔG |
|---|---|---|---|---|---|
| Boltz co-fold | — | 0.965 | 20 | 3.6 Å | **−24.33 (100.0%)**, campaign's strongest cell |
| OF3 re-fold | **8.19** | 0.43 | 14 | 10.4 Å | −10.10 (44.8%) |
| RF3 re-fold | 13.31 | **0.72** | **19** | **8.5 Å** | ? |

pepgraph rates OF3's backbone here among its best anywhere (px 8.19 against a 10.7 lane mean)
and RF3's as poor (13.31). The static pose read says the opposite, and by a wide margin. The
RF3-beats-OF3 dynamics pattern holds 3 of 3 on octinoxate shells so far.

**Stated expectation: RF3 reads substantially better than OF3's −10.10, somewhere around −16 to
−21, and still weaker than the −24.33 baseline.** The closest precedent is `orig_f12`, where RF3's
static read was likewise far better than OF3's and the leg came in at −19.51 against OF3's
−5.90 — there it beat its baseline, but that baseline was −13.71, not −24.33.

Confidence is moderate at best: of four predictions recorded today two were right (`bg33_4`,
`orig_f12` in direction) and two wrong (`bgox31_5`, `s2_esm2_control`). If RF3 comes in at or
below OF3's −10.10, the RF3-beats-OF3 pattern is 3 of 4 rather than 4 of 4 and pepgraph's
ordering gains its first dynamics support on this ligand.

**OUTCOME (20:11): direction right, magnitude badly wrong, and one earlier pattern breaks.**
`s3_orig_f12_rf3cof` reads **−12.10 at 68.9% residence**, 2 released frames, windows flat
(−11.64 / −12.85 / −12.04 / −12.10) so converged. I predicted −16 to −21; it is −12.10 — better
than OF3's −10.10 by only **2.0**, the smallest margin of the four pairs, and **+12.2 weaker than
its −24.33 baseline**. Slide 1.67 Å, so the Δ is clean.

The RF3-beats-OF3 count is now **4 of 4 on octinoxate shells**, but the *shrinking-margin* reading I
noted after `s3_esm2_f4` — that RF3's advantage narrows as OF3's own cell improves — **does not
survive this leg**:

| peptide | OF3 ΔG | RF3 ΔG | RF3 − OF3 |
|---|---|---|---|
| `orig_f12` | −5.90 | −19.51 | −13.6 |
| `s3_orig_f12` | −10.10 | −12.10 | **−2.0** |
| `s2_esm2_control` | −10.21 | −14.98 | −4.8 |
| `s3_esm2_f4` | −13.76 | −16.18 | −2.4 |

Two nearly identical OF3 cells (−10.10 and −10.21) give margins of 2.0 and 4.8, so ordering by OF3
quality is not monotone and that reading is withdrawn after one leg — which is what a four-point
ordering should be expected to do.

**Neither re-folder preserves the campaign's best Boltz co-fold.** −24.33 at 100% residence becomes
−10.10 under OF3 and −12.10 under RF3 — 12–14 kcal/mol weaker either way. That is the strongest case
in the file for the original campaign conclusion that the Boltz co-fold is the start to keep for a
good shell design, and it is untouched by the slide caveats (baseline 2.34 Å, OF3 1.88 Å, RF3
1.67 Å — all clean).

**The pairs column fails again, in the direction that matters.** This leg has the **highest
engagement of any RF3 leg**: sim **2.17**, ≥2 side chains in **76%** of frames, **34/66** pairs, only
1 dead slot, 12/15 long-range. `orig_f12_rf3cof` has essentially the same pair count (35/66, sim
1.96, 0 dead) and binds **7.4 kcal/mol stronger**. Engagement geometry and ΔG are not tracking each
other on this arm.

On the pepgraph tension this cell is mild rather than decisive: pepgraph strongly favoured OF3 here
(px 8.19, among its best anywhere, against RF3's 13.31) and the dynamics favoured RF3 by only 2.0,
so the two disagree in direction but barely in magnitude. It neither resolves the tension nor
sharpens it.

### Octinoxate — RF3 co-folds, peptide structural stability

| leg | Cα final | Cα max | Rg start→end |
|---|---|---|---|
| `orig_f12_rf3cof` | 9.2 | 12.4 | 8.8 → 9.9 |
| `s2_esm2_control_rf3cof` | 1.0 | 3.1 | 10.5 → 10.8 |
| `s3_esm2_f4_rf3cof` | 1.2 | 3.0 | 10.5 → 10.5 |
| `s3_orig_f12_rf3cof` | 4.4 | 5.6 | 8.9 → 9.6 |
| `shuffle_control_esm0_rf3cof` | 1.8 | 3.0 | 11.2 → 11.0 |
| `shuffle_control_rf3cof` | 4.0 | 9.3 | 9.5 → 9.5 |
| `bg33_1_rf3cof` | 10.8 | 12.3 | 14.8 → 9.9 |
| `bg33_2_rf3cof` | 5.0 | 7.3 | 14.9 → 14.0 |
| `bg33_3_rf3cof` | 2.0 | 2.6 | 8.7 → 8.9 |
| `bg33_4_rf3cof` | 2.5 | 4.3 | 9.7 → 10.0 |

**`shuffle_control_esm0_rf3cof` (2026-10-15): pre-registered call wrong.** Called to lose to OF3
(−15.62); it reads **−21.53 ± 0.03 at 99.3% residence, zero released frames** — 7.21 stronger than
its −14.32 baseline and 5.91 better than OF3, on a **null**. Windows −19.48 / −21.48 / −22.25 /
−21.53, converged. Slide 2.64 Å; fold stable (Cα 1.8 final / 3.0 max); contacts 35.3 → 20.3, separation
7.8 → 8.3 Å; pairs 20/66 (band 16%), sim 3.94, 5 dead slots. The OF3 cell it was called against slid
10.62 Å, so the call rested on a broken comparator. A null re-fold out-binding its own baseline by 7
kcal/mol is the strongest RF3 cell in the octinoxate set after `orig_f12`. RF3-beats-OF3: **8 of 11**.

**`shuffle_control_rf3cof` (2026-10-15): half right.** Called the strongest RF3 cell (best static
read: 0.83 / 18 / 6.3 Å); it reads **−11.13 ± 0.03 at 38.7% residence**, 4.00 *weaker* than its
−15.13 baseline, though 1.74 better than OF3's −9.39 (whose cell slid 7.59 Å). Windows −12.90 /
−10.74 / −10.95 / −11.13. 13 released frames in 11 brief episodes (longest 30 ps, 8% late);
separation 9.9 → 10.4 Å (max 18.6), so low residence is mostly parking at the cutoff. Slide 1.69 Å;
Cα 4.0 final / 9.3 max; pairs 24/66 (band 24%), sim 2.11. **The best static read gave a weak cell
and the null `esm0` sibling with a worse read gave the strongest** — static read again does not
predict. RF3-beats-OF3: **9 of 12**, but three of the last four wins are against slid OF3 cells.

**`bg33_1_rf3cof` (2026-10-15): the pre-registered conflict resolves against the molecule split.**
It reads **−12.43 ± 0.02 at 68.3% residence** — 1.12 stronger than its −11.31 baseline but **2.47
weaker than OF3's −14.90**, and that OF3 cell slid only 1.1 Å, so this is a **clean** comparator
(unlike the three slid-OF3 wins before it). Windows strengthen steadily (−9.41 / −11.20 / −12.05 /
−12.43, not converged, so an upper bound). 42 released frames in 12 episodes, longest 290 ps, 2%
late; contacts climb 5.5 → 16.3. The fold collapses and re-forms around the ligand (Cα 10.8 final /
12.3 max, Rg 14.8 → 9.9). Slide 2.51 Å; pairs 19/55 (band 20%), sim 0.94. The static read called a
weak cell and got a middling one (not the weakest RF3 read: `s3_orig_f12` −12.10, `shuffle_control`
−11.13 sit below it). **RF3-beats-OF3 on octinoxate is now 6/7, and its first clean-comparator loss
is here.** Tally 9 of 13.

**`bg33_2_rf3cof` (2026-10-15): a tie with OF3, and the static read misses again.** It reads
**−11.52 ± 0.02 at 58.7% residence** against OF3's −11.58 (clean, slid 1.33 Å) — a 0.06 difference,
a tie at any honest uncertainty — and **3.50 stronger than its −8.02 baseline**. Windows flat
(−11.09 / −11.00 / −10.94 / −11.52). 5 single-frame releases, none late; separation closes 9.5 →
4.4 Å and contacts rise 11.7 → 20.6, so residence is low from an early park outside the cutoff, not
loss. Slide 3.33 Å; Cα 5.0 / 7.3; pairs **27/55 (band 38%), 0 dead slots** — the best pair
realisation in the RF3 octinoxate design arm, from the static read with only 4 of 20 atoms engaged.
**Both `bg33` cells:** neither beat its clean OF3 counterpart (−2.47, −0.06), and both beat their
baselines; the static read's "weakest cells yet" was wrong on both. Neither pre-registered reading
held. Tally (tie counted as no win): 9 of 14.

**`bg33_3_rf3cof` (2026-10-15).** **−12.44 ± 0.02**, 2.34 better than OF3's −10.10 (clean, slid
2.63 Å) and 0.58 stronger than its −11.86 baseline. Windows −11.55 / −12.18 / −11.96 / −12.44.
**Residence 1.9% is a cutoff artefact, not loss:** centroid separation holds at 11.4–11.5 Å the
whole run while contacts *rise* 13.2 → 29.9 and only 5 single frames count as released — the
ligand sits bound against the side of a fold whose centroid is just past 10 Å. Read this cell by
its separation trace, as the file warns. Fold rigid (Cα 2.0 / 2.6). Slide 2.86 Å; pairs 10/55
(band 0%), 6 dead slots. Octinoxate RF3-beats-OF3 7/9; tally **15 of 20**.

**`bg33_4_rf3cof` (2026-10-15), the last RF3 cofold leg.** **−16.13 ± 0.03 at 81.2% residence,
zero releases** — **0.98 weaker than OF3's −17.11** (clean, slid 2.43 Å) and 3.52 weaker than its
−19.66 baseline. Windows −16.79 / −17.81 / −16.06 / −16.13. Fold stable (Cα 2.5 / 4.3). Slide
1.28 Å; pairs 14/55 (band 9%), 5 dead slots.

**Campaign close (all 29 cofold legs in).** RF3 beats OF3 in **15 of 21** cells: octinoxate **7/10**
(every non-win is a `bg33` BoltzGen cell — `bg33_1` −2.47, `bg33_2` a tie, `bg33_4` −0.98 — all
against clean OF3 comparators), oxybenzone **8/11**. Three of the octinoxate wins and one oxybenzone
win are against OF3 cells that slid > 7 Å, so the clean-comparator count is lower than the headline.
The one regular pattern is BoltzGen-specific: every `bgox31` RF3 cell beats OF3 and falls short of
its BoltzGen baseline; on `bg33` RF3 does not reliably beat OF3 at all.

### Octinoxate — the pre-dynamics static read and the pairs the RF3 re-folds realise

| peptide | pose | enclosed | wrapped | engaged | cent sep | sim | pairs | band% |
|---|---|---|---|---|---|---|---|---|
| `orig_f12` | Boltz co-fold | 0.96 | 1.00 | 20 | 4.2 Å | | 22/66 | |
| | RF3 re-fold | 0.80 | 0.90 | 18 | 5.8 Å | **1.96** | **35/66** | 44 |
| `s2_esm2_control` | Boltz co-fold | 0.51 | 0.85 | 17 | 6.7 Å | | 10/66 | |
| | RF3 re-fold | 0.42 | 0.65 | 13 | **13.5 Å** | 1.58 | 15/66 | 7 |
| `s3_esm2_f4` | Boltz co-fold | 0.795 | 0.90 | 18 | 7.0 Å | | 32/66 | |
| | RF3 re-fold | 0.46 | 0.60 | 12 | 9.3 Å | — | 19/66 | 15 |
| `s3_orig_f12` | Boltz co-fold | 0.965 | 1.00 | 20 | 3.6 Å | | 53/66 | |
| | RF3 re-fold | 0.72 | 0.95 | 19 | 8.5 Å | **2.17** | **34/66** | 42 |

### Oxybenzone — RF3 co-folds, ΔG and retention

Eleven legs, each as `<pep>_rf3cof`: `ox1_orig_f12`, `ox1_esm1_f8`, `ox2_orig_f8`, `ox2_esm1_f8`,
`ox3_orig_f8`, `ox2_shuffle`, `ox2_shuffle_esm0`, `bgox31_2`–`bgox31_5`. One has landed.

| peptide | Boltz co-fold ΔG (ret.) | RF3 re-fold ΔG (ret.) | Δ (refold − Boltz) |
|---|---|---|---|
| `ox1_esm1_f8` | −10.70 (75.0%) | **−12.12 (35.2%)** | **−1.42 stronger**; residence is a cutoff artefact, see below |
| `ox1_orig_f12` | −10.31 (49.1%) | **−12.27 (89.1%)** | **−1.96 stronger**, residence up; OF3's cell read −3.96 at 20.5%, so RF3 better by 8.3 |
| `ox2_esm1_f8` | **−20.81 (100.0%)** | −10.90 (90.4%) | +9.91 weaker, residence held; OF3's cell read −16.68, so **RF3 is 5.78 worse** |
| `ox2_orig_f8` | **−24.21 (100.0%)** | −4.95 (19.6%) | +19.3 weaker, **the campaign's largest drop**; OF3's cell read −7.23, so RF3 is 2.28 worse |
| `ox3_orig_f8` | **−18.84 (100.0%)** | −7.41 (45.6%) | +11.4 weaker; OF3's cell read −2.06 at 2.4%, so **RF3 better by 5.35**; still strengthening at 20 ns |

**Prediction for `ox1_orig_f12_rf3cof`, 2026-10-14 20:52, tail running, no ΔG read.** Both
pre-dynamics measures favour OF3 here — pepgraph px **7.09** against RF3's 14.26, and enclosure
0.81/17 engaged at 5.7 Å against 0.55/13 at 6.1 Å — yet OF3's leg is the weakest designed-shell
cell on this ligand (−3.96 at 20.5%, bouncing), against a −10.31 (49.1%) Boltz baseline.

**Stated expectation: RF3 beats −3.96, somewhere around −8 to −14.** The reasoning is a floor
effect rather than anything about RF3: OF3's cell is near the bottom of the campaign, and the
precedent is octinoxate's `orig_f12`, whose equally poor OF3 cell (−5.90) RF3 beat by 13.6.
Low confidence — the one oxybenzone pair so far (`ox1_esm1_f8`) went the *other* way by 5.8, so
RF3-beats-OF3 is 4/4 on octinoxate and 0/1 here.

If RF3 wins, **both static measures agreed and both were wrong on the same cell**, which is a
sharper negative than either has produced alone. If OF3 holds its lead, that is the first
oxybenzone evidence that the ligand, not the re-folder, decides which lane wins.

**OUTCOME (20:58): prediction correct and in range.** `ox1_orig_f12_rf3cof` reads **−12.27 at 89.1%
residence**, 11 released frames in 8 episodes (longest 30 ps, 27% late), inside the stated −8 to −14
window. It beats OF3's −3.96 by **8.3** and its own Boltz baseline by **1.96**. Windows strengthen
(−10.28 → −12.27), slide 1.55 Å, pairs **29/45** with **0 dead slots** against OF3's 15/45.

**Both pre-dynamics measures agreed on this cell and both were wrong.** pepgraph put OF3's backbone
among its best anywhere (px 7.09 vs RF3's 14.26) and the static pose read favoured OF3 on every
column (0.81/17 at 5.7 Å vs 0.55/13 at 6.1 Å). OF3's leg is the weakest designed-shell cell on this
ligand. That is a sharper negative than either measure produced alone, because they are independent
of each other and concordant here.

RF3-beats-OF3 now stands at **5 of 6** (4/4 octinoxate, 1/2 oxybenzone), with `ox1_esm1_f8` the lone
exception:

| peptide | OF3 | RF3 | RF3 − OF3 |
|---|---|---|---|
| `orig_f12` (oct) | −5.90 | −19.51 | −13.6 |
| `ox1_orig_f12` (oxy) | −3.96 | −12.27 | −8.3 |
| `s2_esm2_control` (oct) | −10.21 | −14.98 | −4.8 |
| `s3_esm2_f4` (oct) | −13.76 | −16.18 | −2.4 |
| `s3_orig_f12` (oct) | −10.10 | −12.10 | −2.0 |
| `ox1_esm1_f8` (oxy) | −17.95 | −12.12 | **+5.8** |

The one regularity that survives all six is **floor-dependence, not ligand or lane**: the two cells
where OF3 read worst (−3.96, −5.90) are where RF3 gains most (8.3, 13.6), and the cell where OF3
read best (−17.95) is the only one it wins. Stated as a reading on six pairs with no replicates —
and note it is close to a regression-to-the-mean artefact, which is the first thing to rule out if it
holds up: RF3's six cells span −12.10 to −19.51 while OF3's span −3.96 to −17.95, so RF3 may simply
be the more consistent lane rather than the better one.

One caveat on this cell: the peptide moves a lot (Cα **8.4 Å** final, 10.9 Å max, Rg 7.9 → 8.7), so
like `orig_f12_rf3cof` this is a ligand held by a fold that is itself rearranging.

### Oxybenzone — RF3 co-folds, window convergence

| leg | 5 ns | 10 ns | 15 ns | 20 ns |
|---|---|---|---|---|
| `ox1_esm1_f8_rf3cof` | −8.77 | −9.95 | −11.47 | **−12.12** |
| `ox1_orig_f12_rf3cof` | −10.28 | −10.84 | −12.40 | **−12.27** |
| `ox2_esm1_f8_rf3cof` | −12.61 | −11.47 | −11.15 | **−10.90** |
| `ox2_orig_f8_rf3cof` | −9.40 | −8.11 | −5.55 | **−4.95** |
| `ox3_orig_f8_rf3cof` | −2.72 | −3.52 | −5.38 | **−7.41** |
| `ox2_shuffle_esm0_rf3cof` | −14.33 | −10.99 | −10.21 | **−11.44** |
| `ox2_shuffle_rf3cof` | −13.09 | −13.78 | −13.45 | **−13.35** |
| `bgox31_2_rf3cof` | −9.10 | −8.33 | −8.80 | **−9.06** |
| `bgox31_3_rf3cof` | −14.23 | −14.45 | −14.51 | **−14.36** |
| `bgox31_4_rf3cof` | −10.65 | −11.45 | −11.81 | **−11.02** |
| `bgox31_5_rf3cof` | −15.03 | −15.20 | −15.28 | **−14.69** |

**Prediction for `ox2_esm1_f8_rf3cof`, 2026-10-14 21:55, tail running, no ΔG read.** This is the
first direct test of the **floor-dependence** reading from `ox1_orig_f12`: OF3's cell here is
**strong** (−16.68 at 93.2%, against a −20.81 baseline), so floor-dependence predicts RF3 gains
little or loses — the same regime as `ox1_esm1_f8`, the one cell OF3 has won.

The other two measures split, which makes it a cleaner test than most: the static pose read
favours **OF3** (0.745 enclosed, 16 engaged at 5.3 Å against RF3's 0.57, 12 at 6.4 Å), while
pepgraph favours **RF3** for the only time in the campaign (px 11.94 against 12.57).

**Stated expectation: RF3 lands near or below OF3's −16.68, roughly −14 to −17** — no large gain.
If RF3 instead wins by a wide margin from this high floor, floor-dependence goes the way the
shrinking-margin reading already did, and the RF3-beats-OF3 count stops being explainable as
regression to the mean.

**OUTCOME (22:07): direction right, magnitude wrong again, and floor-dependence survives its first
real test.** `ox2_esm1_f8_rf3cof` reads **−10.90 at 90.4% residence**, 2 released frames, windows
decaying gently (−12.61 → −10.90), slide 1.77 Å. That is **+9.91 weaker than its −20.81 baseline and
5.78 weaker than OF3's −16.68** — no gain, as predicted, but further below my stated −14 to −17 than
I allowed for. Four of five predictions today have been directionally right and numerically loose.

**This is now the better explanation of the whole RF3-vs-OF3 column, and it is not "RF3 is better".**
Ordering all seven pairs by how OF3 did, RF3's advantage falls away and reverses:

| peptide | OF3 | RF3 | RF3 − OF3 |
|---|---|---|---|
| `ox1_orig_f12` | −3.96 | −12.27 | −8.3 |
| `orig_f12` | −5.90 | −19.51 | −13.6 |
| `s3_orig_f12` | −10.10 | −12.10 | −2.0 |
| `s2_esm2_control` | −10.21 | −14.98 | −4.8 |
| `s3_esm2_f4` | −13.76 | −16.18 | −2.4 |
| `ox2_esm1_f8` | −16.68 | −10.90 | **+5.8** |
| `ox1_esm1_f8` | −17.95 | −12.12 | **+5.8** |

**The two cells OF3 wins are its two best cells, and the two where RF3 gains most are its two worst.**
Meanwhile RF3's seven cells span −10.90 to −19.51 (8.6 kcal/mol) against OF3's −3.96 to −17.95
(14.0) — **RF3 is the more consistent lane, not the better one.** That is the regression-to-the-mean
account flagged when floor-dependence was first written down, and with seven pairs it now explains
the column better than any claim that one predictor produces better complexes. The earlier framing
— "RF3 beats OF3 4 of 4 on octinoxate" — **should not be carried into the paper**; the count is 5 of
7 overall and the pattern is a spread effect.

It also removes the cross-molecule story: the two OF3 wins are both oxybenzone, but they are also
both high-floor cells, so molecule and floor are confounded here and the oxybenzone-vs-octinoxate
split needs a high-floor octinoxate cell to separate them. None of the remaining legs is guaranteed
to provide one.

pepgraph missed this cell too: it favoured **RF3** here (px 11.94 against 12.57), the only cell in the
campaign where it does, and RF3 lost by 5.78.

### Oxybenzone — RF3 co-folds, residence and release

| leg | residence | released frames | episodes | longest | late |
|---|---|---|---|---|---|
| `ox1_esm1_f8_rf3cof` | 35.2% | 8 (0.4%) | 8 | 10 ps | 0% |
| `ox1_orig_f12_rf3cof` | 89.1% | 11 (0.55%) | 8 | 30 ps | 27% |
| `ox2_esm1_f8_rf3cof` | 90.4% | 2 (0.1%) | 2 | 10 ps | 100% |
| `ox2_orig_f8_rf3cof` | 19.6% | 622 (31.1%) | 27 | 1290 ps | **99%** |
| `ox3_orig_f8_rf3cof` | 45.6% | 332 (16.6%) | 55 | 650 ps | **1%** |
| `ox2_shuffle_esm0_rf3cof` | 32.6% | 38 (1.9%) | 2 | 370 ps | **0%** |
| `ox2_shuffle_rf3cof` | 51.0% | 1 (0.05%) | 1 | 10 ps | 100% |
| `bgox31_2_rf3cof` | 84.6% | 70 (3.5%) | 10 | 610 ps | 3% |
| `bgox31_3_rf3cof` | **99.2%** | **0 (0.0%)** | 0 | — | — |
| `bgox31_4_rf3cof` | 72.0% | 2 (0.1%) | 2 | 10 ps | 100% |
| `bgox31_5_rf3cof` | **96.4%** | **0 (0.0%)** | 0 | — | — |

**Prediction for `ox2_orig_f8_rf3cof`, 2026-10-14 22:50, tail running, no ΔG read.** `ox2_orig_f8`
has the **campaign's strongest oxybenzone baseline** (−24.21 at 100%) and a poor OF3 cell
(−7.23 at 49.5%, +17.0 weaker — one of the two largest drops anywhere). On the regression reading
a low OF3 floor is where RF3 gains most, so: **RF3 beats −7.23, roughly −10 to −16**, and still
falls well short of −24.21.

The static read mildly favours RF3 (0.595 enclosed, 16 engaged at 6.2 Å against OF3's 0.655, 17 at
7.5 Å — better separation, slightly worse enclosure) and pepgraph favours OF3 (px 7.79 against
9.35), the usual split. Neither has earned weight.

This is a **confirmation, not a test** — the regression reading already predicts it and a hit adds
little. The informative outcome is the miss: if RF3 fails to beat a −7.23 floor, the spread
account loses its main support.

**OUTCOME (22:57): the miss. Both surviving readings are dead.** `ox2_orig_f8_rf3cof` reads
**−4.95 at 19.6% residence**, 622 released frames in 27 episodes, 99% late, windows decaying
−9.40 → −4.95. Against its −24.21 baseline that is **+19.3 weaker, the largest drop in the
campaign**, and it is **2.28 worse than OF3's −7.23** — from a low floor, which is exactly where the
regression reading said RF3 should gain most.

Two accounts die on this one cell:

- **Floor-dependence.** `ox2_orig_f8` sits at OF3 −7.23, between `orig_f12` (−5.90, RF3 gained 13.6)
  and `s3_orig_f12` (−10.10, RF3 gained 2.0). RF3 lost 2.28 there. A low floor no longer predicts a
  gain.
- **Regression to the mean**, which required RF3 to be the narrower lane. **It is now the wider one:
  RF3 spans 14.56 kcal/mol (−4.95 to −19.51) against OF3's 13.99 (−3.96 to −17.95).** The entire
  argument rested on that comparison and it has reversed.

**What is left after eight pairs is one split, and it is now deconfounded.** RF3 beats OF3 on
**4 of 4 octinoxate** cells and **1 of 4 oxybenzone** cells. When that was last noted it was
confounded with floor, because both oxybenzone losses were high-floor cells; `ox2_orig_f8` is a
**low**-floor oxybenzone cell that RF3 also loses, so molecule and floor are now separated and the
molecule split survives while floor does not.

That is the only pattern standing, it rests on eight single trajectories with no replicates, and
every other reading built during collection today — shrinking margins, static separation,
floor-dependence, regression to the mean — has been falsified by a later leg. **The honest position
is that these are not yet predictable**, and the file should say so rather than offer the
molecule split as the next candidate law. The remaining 13 legs are what decide whether even that
survives.

Slide 3.37 Å, flagged. The ligand engages well before it leaves — 22/55 pairs with **0 dead slots**
against OF3's 15/55 — so this is not a pose that never bound; it is one that bound and then lost the
ligand progressively, 99% of the detachment in the second half.

**`ox3_orig_f8_rf3cof` (2026-10-14 23:44) runs the previous leg backwards, which is worth noting.**
It reads **−7.41 at 45.6% residence**, +11.4 weaker than its −18.84 baseline but **5.35 better than
OF3's −2.06**, from the lowest OF3 floor in the campaign. Three things distinguish it:

- **Its windows strengthen monotonically** (−2.72 → −3.52 → −5.38 → −7.41, +4.69 and not converged),
  so −7.41 is an **upper bound** — the opposite of `ox2_orig_f8`'s decay to −4.95.
- **Its detachment is early, not late**: 332 released frames in 55 episodes but only **1% in the
  second half**, where `ox2_orig_f8` was 99% late. Contacts climb through the run (26 → 35) and
  separation closes (5.8 → 4.8 Å). This leg is still binding at 20 ns.
- It does this from a **worse static read than its predecessor** — enclosure 0.995 → 0.465 and
  separation 2.0 → **11.4 Å**, the largest static degradation of any oxybenzone cell.

So the two consecutive `ox*_orig_f8` legs, same ligand and near-identical design family, went opposite
directions on every axis: one decayed from a good static read and lost the ligand late, the other
improved from a poor static read and was still gaining at the end. **RF3-beats-OF3 is now 6 of 9**
(octinoxate 4/4, oxybenzone 2/5). Nothing here restores any of the falsified readings — a very low
OF3 floor (−2.06) did coincide with an RF3 gain, but `ox2_orig_f8`'s low floor (−7.23) did not, so
floor-dependence remains dead rather than rehabilitated.

**`ox2_shuffle_esm0_rf3cof` (2026-10-15): pre-registered call correct, but the comparison is weak.**
It reads **−11.44 ± 0.04 at 32.6% residence**, 0.60 weaker than its −12.04 baseline and **10.5
better than OF3's −0.92** — the gain the pre-registration called. Windows settle after an early −14.33
(−10.99 / −10.21 / −11.44). Release is two early episodes (38 frames, longest 370 ps, **0% late**);
contacts rise to 44.8 and separation sits at 10.0 Å, so the low residence is the ligand parked at
the 10 Å cutoff, not a detachment. Slide 1.71 Å; pairs 17/55 (band 16%), sim 2.18, 4 dead slots.
**Caveat:** the OF3 counterpart slid 13.14 Å in equilibration, so it never tested its predicted pose
and this "win" is against a broken cell. Counted, RF3-beats-OF3 is **7 of 10** (oxybenzone 3/6), but
this tally adds no real support to the molecule split.

**`ox2_shuffle_rf3cof` (2026-10-15): a large RF3 win on a null, against a clean comparator.** No
strong expectation was registered (both cells poor at baseline). It reads **−13.35 ± 0.02 at 51.0%
residence** — **8.05 stronger than its −5.30 baseline and 10.37 better than OF3's −2.98**, whose
cell slid only 1.49 Å. Windows flat (−13.09 / −13.78 / −13.45 / −13.35). One single-frame release;
separation closes 12.2 → 9.2 Å with contacts steady near 25, so residence reflects parking at the
cutoff from a 12.1 Å start, not loss. Slide 2.18 Å; Cα 5.1 / 7.4; pairs **29/55 (band 42%)**, sim
2.71, 1 dead slot — the best realisation in the oxybenzone shell2 RF3 arm, on a shuffled sequence.
Oxybenzone RF3-beats-OF3 is now 4/7; tally **10 of 15**.

**`bgox31_2_rf3cof` (2026-10-15).** **−9.06 ± 0.02 at 84.6% residence** — 1.61 better than OF3's
−7.45 (clean, slid 1.52 Å), 1.83 weaker than its −10.89 baseline. Windows flat (−9.10 / −8.33 /
−8.80 / −9.06). 70 released frames in 10 episodes, longest 610 ps, 3% late. Fold rigid (Cα 1.6 /
4.9). Slide 3.12 Å. No pair columns (no design json for the bgox31 arm). Oxybenzone 5/8; tally
**11 of 16**.

**`bgox31_3_rf3cof` (2026-10-15).** **−14.36 ± 0.02 at 99.2% residence, zero releases** — 8.55
better than OF3's −5.82 (clean, slid 2.11 Å), 1.08 weaker than its −15.44 baseline. Windows flat
(−14.23 / −14.45 / −14.51 / −14.36). The most rigid fold in the RF3 set (Cα 1.2 / 1.8, Rg 8.9
unchanged). Slide 2.46 Å. Oxybenzone 6/9; tally **12 of 17**.

**`bgox31_4_rf3cof` (2026-10-15).** **−11.02 ± 0.02 at 72.0% residence** — 5.28 better than OF3's
−5.74 (OF3 slid 3.83 Å, moderate), 3.04 weaker than its −14.05 baseline. Windows −10.65 / −11.45 /
−11.81 / −11.02. Two single-frame releases; separation 8.0 → 9.3 Å, so the residence gap is parking
near the cutoff. Fold rigid (Cα 1.4 / 1.9). Slide 1.50 Å. Oxybenzone 7/10; tally **13 of 18**.

**`bgox31_5_rf3cof` (2026-10-15).** **−14.69 ± 0.02 at 96.4% residence, zero releases** — 3.14
better than OF3's −11.55 (clean, slid 2.68 Å), 2.10 weaker than its −16.78 baseline. Windows flat
(−15.03 / −15.20 / −15.28 / −14.69). Fold rigid (Cα 1.4 / 2.8). Slide 0.61 Å, the smallest in the
RF3 set. **All four bgox31 RF3 cells beat OF3 and all four fall short of their BoltzGen baselines**
(by 1.1–3.0) — the most regular pattern in the campaign. Oxybenzone 8/11; tally **14 of 19**.

### Oxybenzone — RF3 co-folds, peptide structural stability

| leg | Cα final | Cα max | Rg start→end |
|---|---|---|---|
| `ox1_esm1_f8_rf3cof` | 1.5 | 4.2 | 10.7 → 10.5 |
| `ox1_orig_f12_rf3cof` | 8.4 | 10.9 | 7.9 → 8.7 |
| `ox2_esm1_f8_rf3cof` | 2.2 | 3.3 | 9.6 → 9.3 |
| `ox2_orig_f8_rf3cof` | 2.7 | 8.0 | 9.3 → 8.7 |
| `ox3_orig_f8_rf3cof` | 6.7 | 7.4 | 8.5 → 9.6 |
| `ox2_shuffle_esm0_rf3cof` | 5.9 | 7.0 | 9.8 → 10.7 |
| `ox2_shuffle_rf3cof` | 5.1 | 7.4 | 8.4 → 7.9 |
| `bgox31_2_rf3cof` | 1.6 | 4.9 | 13.4 → 13.9 |
| `bgox31_3_rf3cof` | 1.2 | 1.8 | 8.9 → 8.9 |
| `bgox31_4_rf3cof` | 1.4 | 1.9 | 9.1 → 8.9 |
| `bgox31_5_rf3cof` | 1.4 | 2.8 | 9.3 → 9.4 |

### Oxybenzone — the pre-dynamics static read and the pairs the RF3 re-folds realise

| peptide | pose | enclosed | wrapped | engaged | cent sep | sim | pairs | band% |
|---|---|---|---|---|---|---|---|---|
| `ox1_esm1_f8` | Boltz co-fold | 0.94 | 1.00 | 17 | 4.1 Å | | 30/45 | |
| | RF3 re-fold | 0.53 | 0.77 | 13 | **9.3 Å** | 1.79 | 10/45 | 3 |
| `ox1_orig_f12` | Boltz co-fold | 0.995 | 1.00 | 17 | 1.2 Å | | 42/45 | |
| | RF3 re-fold | 0.55 | 0.71 | 13 | 6.1 Å | 1.63 | **29/45** | 56 |
| `ox2_esm1_f8` | Boltz co-fold | 0.665 | 1.00 | 17 | 6.5 Å | | 36/55 | |
| | RF3 re-fold | 0.57 | 0.71 | 12 | 6.4 Å | — | — | 6 dead |
| `ox2_orig_f8` | Boltz co-fold | 0.975 | 1.00 | 17 | 2.7 Å | | 10/55 | |
| | RF3 re-fold | 0.595 | 0.94 | 16 | 6.2 Å | 0.73 | 22/55 | 27 |
| `ox3_orig_f8` | Boltz co-fold | 0.995 | 1.00 | 17 | 2.0 Å | | 26/45 | |
| | RF3 re-fold | 0.465 | 0.71 | 12 | **11.4 Å** | — | — | — |

## Fold-lane docked legs: Vina pose 1 and GNINA pick (2026-10-15 →)

76 legs on racc (`runs/folds_md/dock_legs.txt`, `dk_oct_1` + `dk5_*` batches): every OF3 and RF3
co-fold structure's Vina pose 1, plus its `receptor_score_only` max-`CNNaffinity` pose where that is
not pose 1. Same receptor as the structure's cofold-pose leg, so each row compares against that leg.
Pose `GNINA ✓` = pick corroborated (same in all three passes and `minimizedAffinity` argmin; 4 of 35) —
an uncorroborated pick rests on `CNNaffinity` alone, which RESULTS.md treats as unsupported. ΔG
marker: • stable / ↓ still strengthening / ↑ still weakening (20 ns vs first-10 ns, 0.5 kcal/mol).
Rows are appended by `code/record_leg.py` as legs land.

| leg | mol | pose | ΔG 20 ns | residence | released | slide Å | cofold-pose ΔG | Δ vs cofold |
|---|---|---|---|---|---|---|---|---|
| `bg33_1_o3cof_dock1` | oct | p1 Vina | **−13.16** • | 22.6% | 1 | 1.46 | −14.90 | +1.74 |
| `bg33_3_rf3cof_dock5` | oct | p5 GNINA | **−11.47** ↑ | 0.1% | 11 | 1.41 | −12.44 | +0.96 |
| `bg33_4_o3cof_dock1` | oct | p1 Vina | **−16.04** • | 99.2% | 0 | 1.23 | −17.11 | +1.07 |
| `bg33_4_o3cof_dock7` | oct | p7 GNINA | **−16.24** ↑ | 82.8% | 0 | 2.6 | −17.11 | +0.87 |
| `bg33_4_rf3cof_dock1` | oct | p1 Vina | **−18.59** ↓ | 97.6% | 0 | 3.01 | −16.13 | −2.46 |
| `bg33_4_rf3cof_dock4` | oct | p4 GNINA ✓ | **−19.06** • | 99.6% | 0 | 2.25 | −16.13 | −2.92 |
| `bgox31_2_o3cof_dock1` | oxy | p1 Vina | **−5.85** ↓ | 10.0% | 684 | 2.05 | −7.45 | +1.59 |
| `bgox31_2_rf3cof_dock1` | oxy | p1 Vina | **−7.36** ↑ | 47.9% | 347 | 1.17 | −9.06 | +1.70 |
