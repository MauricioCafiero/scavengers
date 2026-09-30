# GNINA rescoring bundle

Self-contained. Copy this folder to the Linux machine, run one script, bring back one directory.

## What this is for

`vina_redock.py` docked the ligand back into each of eight folded peptides and Vina failed at the one
thing that mattered: **in all eight structures it ranked the pose closest to the predicted one 6th to
9th of nine.** Even in `bg33_4`, where a pose lands 1.30 Å from the prediction — a genuine redocking
success — Vina put it 6th.

Recognising near-native poses is exactly what GNINA's CNN is trained to do, so this is a fair test of
it rather than a fishing trip. The question is narrow: **shown the predicted pose alongside Vina's
nine, does GNINA rank the predicted one first?**

Two secondary questions come free. Whether the CNN can tell the flip-degenerate orientations apart —
Vina scored end-for-end-flipped poses within 0.22 kcal/mol of unflipped ones — and whether
`CNNaffinity` tracks MM/GBSA ΔG better than Vina's score did (ρ = +0.61, p = 0.11, n = 8).

## Layout

```
reference.csv                 every pose joined to what is already known about it
run_gnina.sh                  the Linux-side runner
collect_gnina.py              parses GNINA's output SDFs into one table
<structure>/
  receptor.pdb                heavy atoms, the file Vina scored against
  receptor_h.pdb              PDBFixer output with hydrogens, the file MD used
  poses_with_reference.sdf    10 poses: pose 0 is the PREDICTED pose, 1-9 are Vina's in rank order
```

Eight structures: `s3_orig_f12`, `s3_esm2_f4`, `s2_esm2_control`, `orig_f12`, `shuffle_control`,
`shuffle_control_esm0` from `peptidebuilder`, and `bg33_4`, `bg33_3` from `boltzgen_local`.

Every pose in `poses_with_reference.sdf` carries `structure`, `pose`, `is_reference`,
`rmsd_to_reference_A` and `vina_score` as SDF tags, so a scored output can be joined without
reference to anything else.

## Run it

```bash
./run_gnina.sh                      # heavy-atom receptor, --score_only
python collect_gnina.py results_receptor_score_only
```

`--score_only` evaluates each pose exactly as given, with no search and no minimisation. That is what
makes this a rescoring of Vina's output rather than a fresh docking — and it matters, because the
comparison is only meaningful if GNINA sees the same coordinates Vina produced.

Worth a second pass:

```bash
./run_gnina.sh receptor minimize    # local relaxation before scoring
./run_gnina.sh receptor_h           # the protonated receptor
```

CNN scores are sensitive to small clashes and these poses came out of a rigid-receptor search that
never relaxed them, so `--minimize` may rank differently. It moves the coordinates, so its results go
in a separate directory and must not be pooled with the `score_only` numbers.

Bring back the whole `results_*` directory. `collect_gnina.py` runs on either machine and needs only
Python's standard library.

## Two things to know before reading any number

**Three of the eight folds contain the mirror-image ligand.** The source geometry,
`runs/octinoxate/ligand.xyz`, has the **S** configuration at the 2-ethylhexyl carbon. Five folds kept
it — `orig_f12`, `s2_esm2_control`, `s3_esm2_f4`, `s3_orig_f12`, `shuffle_control` — and three came out
**R**: `shuffle_control_esm0`, `bg33_3` and `bg33_4`. The ligand SDFs carry the fold geometry unchanged
(verified to 0.0005 Å against the source `.cif`), so the inversion is in what the folding models
placed, not in any downstream prep, and it is not specific to BoltzGen — one of the three is from this
repository's own Boltz runs. The MD and MM/GBSA numbers in `reference.csv` were computed on whichever
enantiomer that fold contains. Within a structure everything is consistent; comparisons *across*
structures mix the two, and a CNN trained on real complexes may be more sensitive to that than Vina's
scoring function is. Grouping by configuration is worth trying if the eight look noisy.

**The ligand is an analogue, not octinoxate.** C17H24O3, one CH₂ short of the real thing
(C18H26O3), with the ester oxygen directly on the chain CH:
`CCCC[C@H](CC)OC(=O)/C=C/c1ccc(OC)cc1`. Every number in this project is for that analogue.

## The numbers to beat

Per structure, from `reference.csv`:

| structure | MM/GBSA ΔG | residence | releases | Vina best | closest pose to prediction | ligand |
|---|---|---|---|---|---|---|
| `s3_orig_f12` | −24.33 ± 0.08 | 100.0% | 0 | −7.3 | 3.49 Å (rank 8) | S |
| `s3_esm2_f4` | −21.08 ± 0.08 | 99.6% | 0 | −6.0 | 3.21 Å (rank 9) | S |
| `bg33_4` | −19.66 ± 0.02 | 100.0% | 0 | −4.7 | **1.30 Å (rank 6)** | R |
| `s2_esm2_control` | −16.25 ± 0.08 | 59.6% | 1 | −5.0 | 5.13 Å (rank 6) | S |
| `shuffle_control` (null) | −15.13 ± 0.13 | 75.9% | 5 | −5.1 | 6.37 Å (rank 9) | S |
| `shuffle_control_esm0` (null) | −14.32 ± 0.10 | 80.2% | 5 | −4.5 | 3.90 Å (rank 9) | R |
| `orig_f12` | −13.71 ± 0.04 | 77.3% | 0 | −6.0 | 5.52 Å (rank 6) | S |
| `bg33_3` | −11.86 ± 0.02 | 41.0% | 18 | −4.2 | 4.47 Å (rank 9) | R |

The measure that did best across the eight was not a score at all: how closely the closest docked pose
reproduces the predicted one correlates **−0.743 (p = 0.035)** with residence. If GNINA's CNN ranks the
predicted pose well, that is the same signal reached more cheaply.
