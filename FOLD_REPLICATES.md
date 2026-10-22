# Fold replicates: how far independent unseeded re-folds land from each other

**3 draws per structure, 29 structures, 87 pairs. 2026-10-22.** Each structure was folded twice more
from the **same yaml** that produced the published fold, with **no seed** — `boltz_env.run_boltz`
passes none and Boltz's default is "None (no seeding)" — so the three structures are independent
samples of the same prediction. The spread is the measurement, not a failure to reproduce.

**No draw is the reference.** The published fold is simply the first sample, so every number below is
the mean over **all three pairwise comparisons** (orig–rep1, orig–rep2, rep1–rep2), with the largest
pair as `max`. Treating the published fold as the truth would privilege one arbitrary draw — the same
mistake as treating run 1 of an MD replicate pair as correct.

Both columns are RMSDs in ångströms after a best fit on the peptide's **Cα atoms** — one per residue,
paired by position since it is the same sequence folded three times:

**RMSD (pep)** — RMSD of those Cα atoms. Backbone shape difference; 0 would be the same fold.

**RMSD (lig)** — RMSD of the ligand's heavy atoms, paired by name, with the peptide in a common
frame. Whether the ligand bound the same place on the peptide. Same construction as
`ligand_slide.py`, so it reads on that scale.

## Octinoxate

| structure | len | gly% | RMSD (pep) mean | max | RMSD (lig) mean | max |
|---|---:|---:|---:|---:|---:|---:|
| **shell 1** | | | | | | |
| `orig_control` | 33 | 63.6 | 10.28 | 12.67 | 17.40 | 25.95 |
| `orig_f12` | 33 | 63.6 | 7.98 | 10.07 | 5.99 | 7.38 |
| **shell 2** | | | | | | |
| `s2_esm2_control` | 33 | 3.0 | 0.58 | 0.65 | 9.57 | 13.12 |
| **shell 3** | | | | | | |
| `s3_orig_control` | 34 | 64.7 | 6.86 | 7.46 | 12.69 | 18.48 |
| `s3_orig_f12` | 34 | 64.7 | 7.27 | 7.90 | 12.37 | 15.28 |
| `s3_esm2_control` | 34 | 0.0 | 8.80 | 13.06 | 14.29 | 18.04 |
| `s3_esm2_f4` | 34 | 0.0 | 1.19 | 1.63 | 8.98 | 12.31 |
| **controls** | | | | | | |
| `shuffle_control` | 34 | 64.7 | 5.33 | 6.18 | 8.98 | 9.85 |
| `shuffle_control_esm0` | 34 | 0.0 | 0.58 | 0.68 | 21.38 | 27.15 |
| **BoltzGen** | | | | | | |
| `bg33_1` | 33 | 6.1 | 0.78 | 1.04 | 23.42 | 33.84 |
| `bg33_2` | 33 | 3.0 | 0.66 | 0.95 | 9.40 | 11.64 |
| `bg33_3` | 33 | 12.1 | 0.26 | 0.34 | 4.02 | 5.57 |
| `bg33_4` | 33 | 6.1 | 0.37 | 0.51 | 10.53 | 15.29 |

## Oxybenzone

| structure | len | gly% | RMSD (pep) mean | max | RMSD (lig) mean | max |
|---|---:|---:|---:|---:|---:|---:|
| **shell 1** | | | | | | |
| `ox1_orig_control` | 25 | 60.0 | 6.68 | 7.77 | 9.30 | 10.78 |
| `ox1_orig_f12` | 25 | 60.0 | 6.02 | 7.01 | 7.59 | 8.49 |
| `ox1_esm1_control` | 25 | 4.0 | 6.37 | 6.94 | 11.13 | 13.59 |
| `ox1_esm1_f8` | 25 | 4.0 | 5.03 | 7.07 | 4.99 | 6.27 |
| **shell 2** | | | | | | |
| `ox2_orig_control` | 31 | 64.5 | 6.81 | 7.70 | 10.49 | 12.85 |
| `ox2_orig_f8` | 31 | 64.5 | 8.76 | 12.25 | 5.41 | 6.77 |
| `ox2_esm1_control` | 31 | 9.7 | 2.27 | 3.14 | 5.97 | 7.08 |
| `ox2_esm1_f8` | 31 | 9.7 | 0.99 | 1.20 | 1.86 | 2.07 |
| **shell 3** | | | | | | |
| `ox3_orig_control` | 27 | 63.0 | 5.45 | 6.80 | 10.71 | 13.83 |
| `ox3_orig_f8` | 27 | 63.0 | 9.12 | 10.89 | 12.32 | 15.20 |
| **controls** | | | | | | |
| `ox2_shuffle` | 31 | 64.5 | 6.68 | 7.66 | 11.39 | 13.78 |
| `ox2_shuffle_esm0` | 31 | 6.5 | 2.20 | 2.71 | 5.08 | 6.06 |
| **BoltzGen** | | | | | | |
| `bgox31_2` | 31 | 3.2 | 0.55 | 0.64 | 15.40 | 22.53 |
| `bgox31_3` | 31 | 6.5 | 0.76 | 1.03 | 3.24 | 4.70 |
| `bgox31_4` | 31 | 6.5 | 0.49 | 0.52 | 4.80 | 6.01 |
| `bgox31_5` | 31 | 6.5 | 0.67 | 0.93 | 5.06 | 6.30 |

## What the set shows

**Glycine content predicts backbone reproducibility.** Across all 29 structures, peptides below 25%
glycine have a mean RMSD (pep) of 1.91 Å against
7.27 Å for those at 60% or more. The designed
`orig`/`shuffle` arms are 60–65% glycine by construction — that is the linker-spacing scheme — and a
chain that is two-thirds glycine has little backbone preference to reproduce. There is nothing in the
data between 12.1% and 60% glycine, so the split is clean but it is two groups, not a gradient.

**The BoltzGen arms are the reproducible ones**, 0.26–0.78 Å, because the generator writes
low-glycine sequences, not because the model is more confident about their binding.

**Backbone reproducibility and ligand placement are different things.** `shuffle_control_esm0`
(0.58 Å backbone, 21.38 Å ligand) and `bg33_1` (0.78, 23.42) fold to the same shape every time and
put the ligand somewhere unrelated every time. `ox2_esm1_f8` is the one structure tight on both
(0.99 and 1.86).

**Two zero-glycine exceptions.** `s3_esm2_control` has no glycine at all and is still
irreproducible (8.80 Å mean, 13.06 max), while its forced twin `s3_esm2_f4` sits at 1.19. So forcing
can stabilise a fold that is otherwise unstable across seeds — on this one pair.

**What this does not say.** Three draws per structure bound the spread loosely, and nothing here
retracts a forced-vs-unforced finding: a fold either had its contacts forced or it did not, and the
MD legs examined the fold that was actually produced. These numbers say at what resolution
single-draw comparisons between arms can be read.

## Reproducing

```sh
python3 code/fold_replicates.py <mol> --n 2 --only <name>       # draw (unseeded; never add --seed)
.venv/bin/python code/compare_folds.py runs/<mol> --replicates  # all-pairs table
```

Draws live in `runs/<mol>/folds/foldrep/<name>/rep<k>/`; every pair is in
`runs/<mol>/folds/fold_replicate_rms.csv`. The published co-folds are untouched.
