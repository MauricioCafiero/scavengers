# The project, in one document

One ligand, two ways of building a peptide around it, and one physical measurement that decides
between them. This file is the definite account: what the pipeline is now, every peptide either
project has produced, what the dynamics said about the eight that were simulated, and what is still
open. It covers **both** peptide families — the shell-designed peptides built here and the
generated binders in [`../boltzgen_local`](../boltzgen_local) — because the comparison between them
is the result, and until now it could not be read from any single file.

Written 2026-10-02. Where another document holds more detail, it is named; where a number came from
a specific file, that file is named. Nothing here is drafted prose for publication.

**The one caveat that governs every number.** The ligand is **C17H24O3, one CH2 short of real
octinoxate** (2-ethylhexyl 4-methoxycinnamate, C18H26O3). Its SMILES, perceived from the stored
geometry, is `CCCC[C@H](CC)OC(=O)/C=C/c1ccc(OC)cc1`. Every result in either project is for that
analogue, and no second ligand has been run, so nothing here is known to transfer.

## Contents

- [1. What is being compared](#1-what-is-being-compared)
- [2. The pipeline as it runs now](#2-the-pipeline-as-it-runs-now)
- [3. Every peptide in the project](#3-every-peptide-in-the-project)
- [4. The result: MM/GBSA and retention](#4-the-result-mmgbsa-and-retention)
- [5. What predicts what](#5-what-predicts-what)
- [6. What is established, and how firmly](#6-what-is-established-and-how-firmly)
- [7. What is not established](#7-what-is-not-established)
- [8. What runs next](#8-what-runs-next)
- [9. Scripts, by pipeline stage](#9-scripts-by-pipeline-stage)
- [10. Where everything is written down](#10-where-everything-is-written-down)

---

## 1. What is being compared

**peptidebuilder (this repository)** designs a *shell*: it screens ten side-chain analogues against
the fixed ligand, picks an arrangement of about twelve poses that surround it, searches for an
ordering a backbone can actually follow, converts the gaps to glycine spacers, and hands the
resulting sequence to Boltz-2 to co-fold with the ligand — optionally with the designed contacts
supplied as constraints. The design is an arrangement of side chains; the sequence is a way of
delivering it.

**BoltzGen (`../boltzgen_local`)** goes the other way: given the ligand alone, it generates a binder
around it, backbone and sequence together, with no designed arrangement to reproduce. It runs on
this laptop at float32 after four patches (see that repository's `MPS_FIXES.md`), about 100 s per
design.

**The null model** is a shuffle of a design's own sequence: identical composition, length, glycine
count, spacer pattern and net charge, with none of the twelve side chains left in its designed slot
(`shuffle_control.py`). It is what says whether the *arrangement* is worth anything, as opposed to
the amino acid content.

All three go through the same back end — co-fold or generate, check the geometry, dock, rescore,
simulate, decompose — so their numbers are comparable. That back end is the subject of the next
section.

## 2. The pipeline as it runs now

Eight stages. The cost column is what it actually takes on this machine or on a rented A10G.

| # | stage | script | cost | what it decides |
|---|---|---|---|---|
| 1 | pose search | `peptide_builder.py` | 70 min measured (4,220 s, `runs/octinoxate/run.log`) | which side-chain poses contact the ligand, and which shell to build from them |
| 2 | reachability sweep | `condense.py sweep`, `condense_analyze.py` | ~1 h, pure geometry | which ordered pairs of poses a backbone can connect, at what spacer count |
| 3 | ordering | `assign.py` | minutes | the path through the poses (exact Held–Karp), and the spacer budget |
| 4 | linker fill, optional | `fill_linkers.py` | minutes | glycine spacers replaced by residues ESM2 finds probable |
| 5 | co-folding | `boltz_hints.py` (or `boltz_check.py`) | ~3 min a fold | the structure, with 0, 4, 8 or 12 designed contacts forced |
| 6 | geometry check | `check_fold.py` | seconds, no GPU | enclosure, wrapping, engagement, Rg, helicity, hints honoured |
| 7 | pose interrogation | `vina_redock.py` → `make_gnina_bundle.py` → GNINA | minutes of local CPU | whether the pose is findable, and the best static predictor of ΔG available |
| 8 | dynamics + MM/GBSA | `cif_to_md.py` → `omd` → `modal_md.py` → `md_contacts.py`, `pair_contacts.py`, `md_frames.py` | $0.45 + $0.041 per 1,000 particles | the verdict: free energy, retention, release, what the complex does in water |

Stages 1–4 are the design method and the part that is novel. Stages 6–8 are the judgement, and
their order matters: stage 6 costs seconds and should run on everything, stage 7 costs minutes of
local CPU and should run on anything that might be simulated, stage 8 is the only thing that has ever
settled a disagreement between the first two.

**Nothing here is free.** Stages 1–7 rent no GPU, which is the only sense in which they are cheap;
they still cost wall-clock time on this laptop and wear on it, and stages 1 and 2 cost about an hour
each (stage 1 measured at 4,220 s). "No rental" is the right way to describe a local stage, and a local stage that takes hours is
not cheaper than a $0.70 Modal leg in any sense that matters to a day's work.

**Stage 8 lives in another environment.** MD and MM/GBSA are the `omd` CLI in
`~/python_mac/openmm`, which needs conda — GAFF2/OpenFF parameterisation is not pip-installable and
`.venv` has no `mdtraj`. `code/cif_to_md.py` is the adapter, and it assigns ligand bond orders from
the run's own SMILES as a template; without that step the aromatic ring and the ester come back as
single bonds and the ligand is parameterised as a different molecule. **Check the printed formula.**

**Stage 7's two scorers read different axes**, which is why both are kept. Vina's own score tracks
the *escape* axis (ρ = +0.732, p = 0.039 against release episodes) and is useless for choosing a
pose — in all eight structures it ranked the pose closest to the prediction 6th to 9th of nine.
GNINA's `CNNaffinity` on the predicted pose tracks the *energy* axis at ρ = −0.857, p = 0.007
against MM/GBSA ΔG, which is the best static predictor in either project. Run GNINA with
`--score_only` on the heavy-atom receptor: `--minimize` weakens both correlations, and receptor
hydrogens change `CNNscore` by at most 0.034.

### What UMA is for, and what it is no longer for

This is the change that most of the older documents have not caught up with, so it is worth stating
flatly.

**UMA still builds the designs.** It relaxes each fragment pose with the ligand fixed, supplies the
interaction energy that accepts or rejects a pose, and measures ligand strain against the single
shared reference conformer from `ligand_reference.py`. Those uses are not in question, and
`runs/<ligand>/uma_logs/` plus `boltz/structures/` are kept because each relaxation cost 45–60
minutes and is the only record of how an energy was arrived at.

**UMA is retired as a predictor.** The static scoring pass over folded complexes — `binding_energy.py`
and `modal_score.py`, 30–70 minutes a fold — is no longer run on new structures, for three measured
reasons:

1. **It ranks structures backwards.** Against MM/GBSA ΔG on the four designs simulated, the static
   sum gives ρ = −0.80 and the UMA interaction term alone −0.80. The best static score is third of
   four in water and the worst is first.
2. **It rewards net charge.** Across 36 folds the interaction energy correlates −0.56 with charged
   residue fraction; `s2_esm2_control` at charge +5 posts −74.89 kcal/mol, of which 91% evaporates
   in water (MM/GBSA −16.25). Within the `esm2` family the energy correlates +0.93 with *centroid
   separation* and +0.18 with enclosure, the wrong sign for both.
3. **It misranks pairs where charge is not the variable.** `orig_f12` and `s3_orig_f12` are
   geometric near-twins; the static sum separates them by 52 kcal/mol the wrong way, and the
   inversion survives in the interaction term alone (UMA prefers `orig_f12` by 20.73; the MM/GBSA
   gas-phase term prefers `s3_orig_f12` by 17.33).

What replaced it: `check_fold.py` geometry at stage 6, GNINA `CNNaffinity` at stage 7, and
MM/GBSA at stage 8. The shell-*selection* criterion fell with it — total fragment interaction energy
is substantially measuring charge, and on three shells it got the fold-quality order backwards, so
it should be replaced by closure quality and pose diversity.

## 3. Every peptide in the project

### peptidebuilder: 38 folds

Three shells, each giving one glycine design plus two ESM2 linker variants, each co-folded four ways
(unconstrained, and with the 4, 8 and 12 tightest designed contacts forced) — 36 folds — plus two
nulls. All in `runs/octinoxate/boltz/fold_check*.csv`.

| | shell 1 | shell 2 | shell 3 |
|---|---|---|---|
| seed | arginine-rich, richest arrangement | second richest | **tryptophan pose 2** |
| total fragment interaction energy | −113.12 | −109.15 | **−60.40** (47% weaker) |
| charged / aromatic poses | 7 (58%) / 1 | 5 (42%) / 2 | **4 (33%) / 4 (33%)** |
| glycine design | `RGGDGGKGGGGLGGGGKGIGEGWGGDGSGGEGS` | `RGEGGEGGKGGFGGDGGIGLGGSGWGGGGLGGS` | `YGGDGFGKGGGLGGGKGGGGLGGWGGEGSGGWGS` |
| residues / spacers used (budget 16) | 33 / 21 | 33 / 21 | 34 / 22 |
| poses in the path | 12/12 | 12/12 | 12/12 |
| best enclosure | 0.960 (`orig_f12`) | 0.860 (`s2_esm2_f8`) | **0.965** (`s3_orig_f12`) |
| designed positions reproduced | not measured | 1 / 144 | **8 / 144** |
| net favourable by the static sum | 9/12 | 8/12 | 5/12 |

The nulls are `shuffle_control` (a shuffle of shell 3's glycine design, folded compact at Rg 8.2 Å)
and `shuffle_control_esm0` (the same shuffle with ESM2 linkers, folded as a 100% helical rod at Rg
15.0 Å) — deliberately one of each architecture, so release can be compared against a designed
structure of matching shape.

Two facts about this set are worth carrying forward because they shaped everything after it. **All
three shells overshoot the spacer budget by the same margin**, so the overshoot is a property of
twelve-fragment shells around this ligand rather than of any arrangement. And **the weakest shell
folded best**: shell 3 reproduces eight designed positions to shell 2's one, and its `s3_orig_f12` is
a tighter twin of shell 1's flagship on every count, from two shells sharing only part of their pose
set. That convergence is the most reproducible structural result in the project.

### boltzgen_local: 8 designs

Protocol `peptide-anything`, six at 33 residues to match this project's length and two shorter, all
at float32. Each refold was measured with *this* repository's `check_fold.py` after a column-layout
conversion, so the numbers mean the same thing. Full table in
[`../boltzgen_local/README.md`](../boltzgen_local/README.md).

| structure | len | enclosed | wrapped | engaged | Rg | sequence note |
|---|---|---|---|---|---|---|
| **`bg33_4`** | 33 | **0.59** | 0.95 | 19/20 | 9.5 | near poly-alanine, one of the two taken to MD |
| `bg33_2` | 33 | 0.58 | 0.75 | 15/20 | 14.5 | 29 alanines of 33 |
| `bg33_5` | 33 | 0.56 | 0.50 | 10/20 | 9.2 | structured |
| `bg33_3` | 33 | 0.48 | **1.00** | **20/20** | 8.5 | structured, three-stranded sheet; the other MD leg |
| `bg33_1` | 33 | 0.38 | 0.75 | 15/20 | 14.4 | 2.34 Å contact — a broken interface by this project's cutoff |
| `bgA` | 13 | 0.35 | 0.75 | 15/20 | 6.2 | poly-alanine |
| `bgB` | 20 | 0.34 | 0.60 | 12/20 | 9.2 | poly-alanine |
| `bg33_0` | 33 | 0.28 | 0.55 | 11/20 | 14.6 | poly-alanine |

**Every one of the eight is less enclosed than the shuffled null** (0.74), and the best of them is
well short of this project's best fold (0.965). Five came back near poly-alanine because
`--inverse_fold_num_sequences` defaults to 1, so there was no sequence ensemble to select among —
and the strongest binder of the eight is one of the poly-alanine designs, so that is a distribution
problem rather than a verdict. Eight designs is also the regime BoltzGen's own documentation warns
against; its filtering stage had nothing to filter. Read the comparison as *what BoltzGen gives at
laptop scale*, not as a characterisation of the method.

### The 18 dynamics legs, and where they live

| legs | structures | where |
|---|---|---|
| 8 co-folded (the predicted pose) | `s3_orig_f12`, `s3_esm2_f4`, `s2_esm2_control`, `orig_f12`, `shuffle_control`, `shuffle_control_esm0` | `runs/octinoxate/md/<structure>/` |
| | `bg33_4`, `bg33_3` | **`../boltzgen_local/md/<structure>/`** |
| 10 docked-pose legs | `s3_orig_f12_dock1`, `s3_orig_f12_dock9`, `s3_esm2_f4_dock1/6`, `s2_esm2_control_dock1/6`, `bg33_4_dock1/6`, `bg33_3_dock1/6` | `runs/octinoxate/md/<structure>_dockN/` |

`runs/octinoxate/md/mmgbsa_summary.csv` holds 16 of the 18: the two BoltzGen co-folded legs were run
in the other repository and were never added. **So the headline table below cannot be rebuilt from
this repository's CSV alone**, which is one concrete form of the fragmentation this document exists to
fix. `make_gnina_bundle.py` carries the same ΔG values again in a hardcoded `REFERENCE` dict, so a
new leg has to be written into both places.

## 4. The result: MM/GBSA and retention

Lead with this. Everything else in either project is either upstream of it or methodology about how
to spend GPU time.

| peptide | source | ΔG (kcal/mol) | retention | releases | Vina | GNINA `CNNaffinity` | closest redock |
|---|---|---|---|---|---|---|---|
| `s3_orig_f12` | peptidebuilder | **−24.33** | 100.0% | 0 | −7.3 | 4.85 | 3.49 Å |
| `s3_esm2_f4` | peptidebuilder | −21.08 | 99.6% | 0 | −6.0 | 5.04 | 3.21 Å |
| `bg33_4` | boltzgen | −19.66 | 100.0% | 0 | −4.7 | **5.05** | **1.30 Å** |
| `s2_esm2_control` | peptidebuilder | −16.25 | 59.6% | 1 | −5.0 | 4.57 | 5.13 Å |
| `shuffle_control` | peptidebuilder (null) | −15.13 | 75.9% | 5 | −5.1 | 4.27 | 6.37 Å |
| `shuffle_control_esm0` | peptidebuilder (null) | −14.32 | 80.2% | 5 | −4.5 | 4.50 | 3.90 Å |
| `orig_f12` | peptidebuilder | −13.71 | 77.3% | 0 | −6.0 | 4.04 | 5.52 Å |
| `bg33_3` | boltzgen | **−11.86** | **41.0%** | **18** | −4.2 | 4.20 | 4.47 Å |

| group | n | mean ΔG | best | worst | mean retention | releases |
|---|---|---|---|---|---|---|
| **peptidebuilder (designs)** | 4 | **−18.84** | **−24.33** | −13.71 | **84.1%** | **1** |
| boltzgen | 2 | −15.76 | −19.66 | −11.86 | 70.5% | 18 |
| peptidebuilder (nulls) | 2 | −14.73 | −15.13 | −14.32 | 78.1% | 10 |

**The designs from this pipeline come out ahead on every measure available**: best binder by 4.67
kcal/mol, better mean free energy, higher mean retention, one release episode against eighteen. The
worst structure in either project on energy *and* retention simultaneously is a BoltzGen design.

Three qualifications, none of which changes the ordering. n is 4 against 2, so nothing here
characterises BoltzGen as a method. Both BoltzGen structures carry the **R** ligand where all four
peptidebuilder designs carry **S**, so pipeline and configuration are confounded —
`shuffle_control_esm0` is the only peptidebuilder structure that came out R and is the one available
lever on that. And `shuffle_control` at −15.13 beats the design `orig_f12` at −13.71, so a shuffle
of a sequence still beats a real design outright in one of four cases: "designed" is not sufficient,
only the good designs are distinguishable from their own shuffle.

**What the designed arrangement is worth, with one variable.** `s3_orig_f12` against its own shuffle
— identical composition, length, glycine count, spacer pattern and net charge, 0 of 12 side chains
in place, folded to the same compactness — is **9.2 kcal/mol**. That is the cleanest single-variable
result in the project.

**Release is what separates design from null, and averages hide it.** On mean behaviour
`shuffle_control` and `orig_f12` are the same trajectory (35.3 against 34.3 contacts, 8.64 against
8.45 Å separation, identical Rg drift and closest approach). The difference is that the null lets go
— five episodes up to 70 ps, **entirely in the second half of the run**, progressive loss rather
than thermal noise — and the signature replicated in the second null, which Boltz folds as a rod
rather than a globule: five episodes, up to 40 ps, again all in the second half. Both designed
structures in the matching architectures give zero or one isolated blip. Release tracks how the side
chains are *ordered*, not whether the peptide is a globule or a rod.

Measure it at matched frame spacing or not at all: sampled every 1 ps `orig_f12` shows 8 released
frames, at the 10 ps spacing of the others it shows zero, because all eight were isolated single
frames. `md_contacts.py` resamples to a common spacing for exactly this reason.

## 5. What predicts what

Ranked by what it costs to compute — GPU rental where there is any, otherwise local wall-clock
time — with the sign convention that + means agreement.

| measure | cost | what it predicts | strength |
|---|---|---|---|
| mean side chains simultaneously on the ligand | a trajectory | ΔG | **+0.83** (n=6), +0.90 (n=5) |
| % of frames with ≥2 side chains engaged | a trajectory | ΔG | **+0.83**, +0.90 |
| MM/GBSA `VDWAALS` | a trajectory | ΔG | +0.71 |
| GNINA `CNNaffinity`, predicted pose | **minutes** | ΔG | **−0.857, p = 0.007** |
| closest-pose redock RMSD | **minutes** | retention | **−0.743, p = 0.035** |
| Vina score | **minutes** | release episodes | **+0.732, p = 0.039** |
| pose-1 `enclosed_fraction` | seconds | release episodes | −0.713, p = 0.047 |
| `enclosed_fraction`, `contacts_under_cutoff` | seconds | ΔG | +0.30 |
| UMA interaction energy, ± ligand strain | 45–60 min a fold | ΔG | **−0.80 — backwards** |

Four things follow, and they are the practical content of this section.

**The combinatorial count is not the predictor; simultaneous engagement is.** The pipeline delivers
*n* positions chosen to contact the ligand, which is n(n−1)/2 pairwise options in the band
[n−1, n(n−1)/2] — [11, 66] for twelve slots. Pairs realised correlates at +0.49; mean simultaneous
engagement at +0.83. `s3_orig_f12` realises 53/66 against its shuffle's 18/66 and sustains 24 pairs
where the null sustains 8, *below the free floor of 11*. And engagement explains release without
further assumptions: at 3.82 side chains engaged, losing one does not detach the ligand; at 1.75 the
ligand is often held by a single contact.

**The pair framing fails on rods.** Designed and shuffled rods are indistinguishable on it — 10/66
against 12/66, 1.45 against 1.73 engaged, six dead slots each, 0/15 long-range pairs both — yet they
score −16.25 and −14.32. It diagnoses how a compact fold binds and says nothing useful about a
groove on a helix.

**Wrapping is a contact census, not a verdict.** The highest wrapped value in shell 2, 1.00, belongs
to a groove binder; `bg33_3` wraps 1.00 and engages 20 of 20 atoms and is the worst binder of the
eight. Read `peptide_rg` first to know whether a fold encapsulates (9.3–9.8 Å) or grooves (11.8–14.9
Å) — enclosure separates those two groups by 0.015, which is far too narrow to act on.

**Hint satisfaction is a third metric that disagrees with both others.** `orig_f12` honoured none of
its twelve requested contacts and is the best structure in shell 1; `s3_esm2_f4` satisfies 0/4 and
binds at −21.08. Forcing contacts works, but not by being satisfied — it works by making the fold
surround the ligand. What forcing does reliably is improve the cavity; what it does unreliably is
improve the score.

## 6. What is established, and how firmly

In descending order of how much weight each can carry.

1. **20 ns is the minimum and short windows are confidently wrong.** `orig_f12`'s leading windows
   read −25.37, −21.52, −16.86, −14.61, −13.71 at 40 ps through 20 ns: successive changes halving,
   every window tight, every short one too negative, because they are all still measuring the
   predicted pose rather than the ensemble. Replicated on `bg33_3` (2.7 kcal/mol of drift) and
   `s3_orig_f12_dock1` (+6.22, still moving at 20 ns).
2. **The designed arrangement is worth 9.2 kcal/mol against its own shuffle**, single-variable.
3. **A co-folded pose converges where a docked pose drifts.** `bg33_4` predicted: −19.18, −18.89,
   −19.43, −19.66 across the series, 0.77 kcal/mol total. `s3_orig_f12_dock1` docked: −27.00 →
   −20.79 monotonically and still moving. A docked pose starts over-contacted and decays; a
   co-folded one starts at equilibrium.
4. **The static gas-phase advantage of a charged peptide is an artifact, quantitatively.** 29.0
   kcal/mol becomes 2.5 in water; 91% evaporates, and the +5 peptide's electrostatic attraction to
   the neutral ligand is *weaker*, not stronger. What UMA was rewarding is long-range polarisation,
   which GB screens almost completely.
5. **The predicted pose beats Vina's best-scoring pose on the same fold**, 3.55 kcal/mol under an
   identical protocol, and it is the pose that tightens over 20 ns rather than loosening.
6. **The peptide grips the ligand but does not orient it.** 26 of 72 docked poses are turned end for
   end while keeping the same contact residues (Jaccard mean 0.60, 60 of 72 sharing at least half),
   and flipped poses score −4.85 against unflipped −5.07. Both scorers are blind to it: GNINA
   separates flipped from unflipped at p = 0.283.
7. **One structure out of eight actually redocks.** `bg33_4`'s pose 6 at 1.30 Å is the only hit by
   the 2 Å standard; this repository's six run 3.2–9.9 Å at best. Both folders placed the ligand
   inside the peptide, so the pocket is the ligand's own imprint and this was the easiest possible
   target — which makes agreement weak evidence and disagreement strong.
8. **Dock only when the co-folded pose is poorly retained.** Four supporting cells, no
   counterexample among legs that held their starting pose: at 100% retention co-folding wins by
   3.5–11.1 kcal/mol, at 59.6% and 41% a GNINA-chosen pose wins by 0.84 and 1.23.

## 7. What is not established

- **There is no replicate anywhere.** Every number in section 4 is one leg. The only two
  near-identical starting points ever run (`bg33_4` co-folded and `bg33_4_dock6`, 0.52 Å apart at
  placement) diverged by **3.13 kcal/mol** — larger than the 3.08 pipeline gap and comparable to the
  4.11 design-versus-null gap. That divergence bounds *pose* discrimination, not protocol noise, and
  protocol noise has never been measured.
- **MM/GBSA standard errors are not uncertainty on a comparison.** They measure frames within one
  trajectory (0.04–0.10 kcal/mol). Never write "this difference is N standard errors".
- **Every correlation in section 5 is n = 8 or smaller**, release episodes are (0, 0, 0, 0, 1, 5, 5,
  18) so that variable is nearly a design/null split, and one BoltzGen pair supplies both ends of the
  retention correlation.
- **Three of eight folds contain the mirror-image ligand.** The source is **S**;
  `shuffle_control_esm0`, `bg33_3` and `bg33_4` came out **R**, placed that way by the folding
  models (SDFs verified against the CIFs to 0.0005 Å). Within-structure numbers are unaffected;
  cross-structure score comparisons mix enantiomers. Checking the other 30 folds is an RDKit pass
  over the CIFs — perceive the CIP label, not the canonical SMILES, which carries double-bond stereo
  the CIF-derived SDF does not and flags all eight as different.
- **The protocol does not hold the pose it is given.** Equilibration restrains protein heavy atoms
  only (`openmm_md/dynamics.py:144`), so the ligand slides 1.15–15.11 Å before production. Do not
  "fix" this — every run in both projects equilibrated it freely and a restrained leg would be
  incomparable. Record the number; it is what says whether a cell is worth reading.
- **Designed side-chain positions are essentially never reproduced**: 8/225 within 3 Å in the
  shell-1 analysis, 1/144 in shell 2, 8/144 in shell 3, at ligand superposition RMSD 0.06–0.86 Å so
  the comparison is sound. The best binder in the project reproduces none of its twelve. The
  distances say something the count hides — forcing contacts moves side chains from 12–14.5 Å to
  7.6–9.8 Å from their designed positions, and the nearest side chain *of any type* reaches 3.0–4.9
  Å — so the failure is in which residue lands where, not in reachability, which the sweep shows is
  satisfiable for every pair.
- **One ligand, and an analogue of it.**

## 8. What runs next

**Before any of the pose-matrix cells is quoted** — [NEXT_STEPS.md](NEXT_STEPS.md) item 7 — window
the five matrix legs and the two Modal baselines, which were run without a convergence series
(`md_window_live.py`, `run_windows_live.sh`). MM/GBSA only, no GPU, a few minutes of local CPU per
window against trajectories already on disk. A series still walking away from its start at 20 ns was
not at equilibrium whatever its endpoint reads, and that is knowable without new dynamics.

**Two legs, ≈$1.55 of GPU together.** Each is a fold whose comparison partner is already on disk, so
one new leg buys a complete sequence-matched pair — the only form of this test that holds
composition, net charge and ligand configuration fixed and varies nothing but the fold.

1. **`s3_orig_f8` co-folded, 20 ns** — the *same sequence* as the best binder in the project,
   folded with 8 forced contacts instead of 12. Enclosure 0.88 against `s3_orig_f12`'s 0.965,
   wrapping 0.85 against 1.00, 59 contacts against 67, centroid separation 7.9 against 3.6 Å. The
   partner reads −24.33, so this says directly whether the cavity quality the design objective
   measures is what produced that number. Compact at Rg 7.8 Å, so about 5–6,000 particles, ≈$0.70.
2. **`s2_esm2_f8` co-folded, 20 ns** — the same sequence as `s2_esm2_control`, which is already
   simulated at −16.25, and the one fold in shell 2 whose helix *broke* to encapsulate: enclosure
   0.86 against the control's 0.51, wrapping 0.95 against 0.85, 31 contacts against 18, Rg 9.8 Å
   against 15.0. It is the best-enclosed structure in shell 2 and it scored poorly on the static
   energies, which is the clearest standing case of enclosure and energy disagreeing. Being compact
   it also solvates to far less than its partner's 23,121 particles, so ≈$0.85 rather than $1.40.

Together they take the designs from n = 4 to n = 6 in [section 4](#4-the-result-mmgbsa-and-retention)
and give the design objective its first two controlled tests. Both pairs are within one sequence, so
neither inherits the enantiomer confound or the charge confound that every cross-structure comparison
in this project carries. This is [NEXT_STEPS.md](NEXT_STEPS.md) item 1 — more dynamics on the
remaining shell rungs — made specific. The runner-up on that item is **a second shuffle seed**, which
would say whether the 9.2 kcal/mol design-against-null margin replicates; it is the next leg after
these two, and it is a weaker use of the first two because the enclosure pairs come with their
partners already measured.

**Neither of these is a redock.** The eight structures of interest were all docked with Vina and
rescored with GNINA on 2026-10-01 to 10-02, and none of that is repeated. A fold that has never been
prepared for MD has no receptor or reference ligand on disk, so stage 7 runs once on each of the two
new structures as part of their prep, which is the normal order of the pipeline and not a second pass
over anything. It happens to land each one in the retention correlation, whose 1.5–3.0 Å band is
currently empty with only `bg33_4` anchoring the low end.

Build with `--box-shape dodecahedron` explicitly — `omd build` defaults to a cube, and six of the
seven co-folded baselines are dodecahedral. Use `run_dock_pose_md.sh` with `N` parameterised rather
than a new driver. Record input → production ligand RMSD. Write each result into **both**
`mmgbsa_summary.csv` and the `REFERENCE` dict in `make_gnina_bundle.py`.

**Queued, deliberately not now: the first replicate.** Run-to-run spread has never been measured and
the only two near-identical starting points ever run diverged by 3.13 kcal/mol, so every comparison
in section 4 rests on single legs. That is worth fixing, but not before the single-leg set is larger
— spending the next legs on new structures buys more than spending them on a second copy of one.
When the time comes, the first replicate should be `s3_orig_f12`, which is the structure every
headline number is anchored to, and it is the cheapest replicate available: `modal run
code/modal_md.py::produce --structure s3_orig_f12 --leg L2` reuses the system already on the Modal
volume, so there is no rebuild, the box and its 5,615 particles are identical, and the integrator
takes no explicit seed so the trajectory is independent. ≈$0.68. Recorded as
[NEXT_STEPS.md](NEXT_STEPS.md) item 8, with the reasoning, so it is not lost.

**The design-method work queued beyond that** is in [NEXT_STEPS.md](NEXT_STEPS.md): fragment
condensation (A), reachability-constrained assignment (B), and the attachment-geometry piece both
need, which is the information `build_sequence` currently discards when it chains pose *centroids*
by nearest neighbour. Option C is explicitly not being pursued. Two further shells —
`tryptophan/2` and `isoleucine/2` — are specified in [HANDOFF.md](HANDOFF.md) §8. A second ligand is
the one experiment that would say whether any of this transfers.

## 9. Scripts, by pipeline stage

Nothing in `code/` needs to be rewritten to vary a parameter. **A hardcoded value is a parameter to
add, not a reason to start over** — on 2026-10-01 four purpose-written MD drivers were written
instead of parameterising `run_dock_pose_md.sh` and silently dropped three of its steps across nine
legs. README.md has the full per-file table; this is the same inventory arranged by stage, including
the seven scripts the README's tables do not mention.

| stage | scripts |
|---|---|
| design | `peptide_builder.py`, `attach_geom.py`, `condense.py`, `condense_analyze.py`, `assign.py`, `fill_linkers.py`, `check_design.py` |
| alternative route | `condense_chain.py`, `condense_strain.py`, `design_test.py` |
| folding | `boltz_hints.py`, `boltz_check.py`, `boltz_env.py`, `boltz_offline.py` |
| geometry | `check_fold.py`, `compare_folds.py`, `overlay.py` |
| static energetics (retired as a predictor) | `binding_energy.py`, `uma_binding.py`, `ligand_reference.py`, `strain_global.py`, `solvate.py`, `modal_score.py` |
| nulls | `shuffle_control.py`, `random_control.py` |
| docking and rescoring | `vina_redock.py`, `dock_compare.py`, `dock_vs_md.py`, **`dock_pose_to_sdf.py`**, **`make_gnina_bundle.py`** |
| MD and MM/GBSA | `cif_to_md.py`, `modal_md.py`, **`run_dock_pose_md.sh`**, **`md_window_live.py`**, **`run_windows_live.sh`** |
| trajectory analysis | `md_contacts.py`, `pair_contacts.py`, `md_frames.py` |
| correlation and figures | `correlate.py`, `make_figures.py` |
| shell runners | `run_shell.sh`, **`run_shell2.sh`**, **`run_shell2_forced.sh`** |

The five in bold under docking and MD are the current pipeline's own stages, and they are the ones
missing from the README's file tables — which is the documentation lag this file is here to end.
What each does, in one line:

- **`dock_pose_to_sdf.py`** — writes one docked pose as an SDF `omd prep-ligand` will accept. The
  pose is not read as a molecule: the reference SDF supplies the topology and only coordinates are
  replaced, reindexed by `vina_redock.ref_order_poses`, because Open Babel's pdbqt reader loses
  valence.
- **`make_gnina_bundle.py`** — assembles a self-contained folder for rescoring elsewhere: ten poses
  per structure in a fixed order (pose 0 the predicted one, then Vina's nine in rank order), two
  receptors, and a hardcoded `REFERENCE` dict of the MM/GBSA values.
- **`md_window_live.py`** — slices the leading *n* ns out of a trajectory that is still being
  written, by `iterload` with an explicit frame budget read from `energy.csv`, because
  `mdtraj.load` on a growing 20 ns file is ~700 MB resident on a machine already swapping.
- **`run_windows_live.sh`** — drives those windows as the run passes each mark. A window is 4–16
  minutes of one core against a trajectory holding another, so the series adds little to a run that
  is already occupying the machine — which is the argument for computing it during the run rather
  than after it.
- **`run_dock_pose_md.sh`** — the one MD driver. Stage guards so a re-run resumes, MM/GBSA invoked
  from inside the leg directory, the leading windows, `caffeinate`, and `md_contacts.py` plus
  `md_frames.py` at the end.
- **`run_shell2.sh` / `run_shell2_forced.sh`** — shell 2 end to end and its nine forced folds. The
  second exists because `--boltz-args '--use_potentials'` reaches argparse as two argv entries and
  every forced fold dies at parse time; the working form is `--boltz-args=--use_potentials`. Both
  `force: true` *and* `--use_potentials` are required — `force: false` is discarded by the
  featurizer rather than softened, and `force: true` alone is only conditioning.

**Housekeeping that is not optional.** Run MM/GBSA from inside the leg directory: MMPBSA.py scatters
`reference.frc` plus a `_MMPBSA_*` set into the working directory, ~300 MB a leg, and `reference.frc`
reached 1.9 GB once when run from the repository root. `md_contacts.py` and `md_frames.py` are
pipeline steps, not extras — a leg without them is incomplete.

## 10. Where everything is written down

**Live documents.**

| file | what it is for |
|---|---|
| **PROJECT.md** (this file) | the definite account: pipeline, both peptide families, the result, what is open |
| [README.md](README.md) | how to install and run it, what every file does, what the metrics mean, the traps, the limitations |
| [RESULTS.md](RESULTS.md) | the long-form evidence, ~1,580 lines, ordered chronologically by shell — the three worked examples, dynamics, redocking, GNINA, the pose matrix |
| [NEXT_STEPS.md](NEXT_STEPS.md) | the design-method plan: methods A and B, the attachment geometry both need, and what is rejected |
| [CLAUDE.md](CLAUDE.md) | working rules for this repository, written after a day in which each one was learned by breaking it |
| METHODS_OUTLINE.md | untracked: a bullet outline of the manuscript methods section, with `[decide]` and `[gap]` markers |
| [`../boltzgen_local/README.md`](../boltzgen_local/README.md) | the BoltzGen side: install and the four Mac patches, the eight designs, their dynamics, and the two format details that stand between the projects |

**Superseded, kept for provenance.** [PROGRESS.md](PROGRESS.md) is a session log that stops at
2026-09-25 and predates dynamics, docking and GNINA entirely. [HANDOFF.md](HANDOFF.md) is mostly
discharged — its code fixes are done and shell 2 is run — but §8 (the next two shells) and §10 (the
ligand strain definition) are still the reference for those. README_previous.md is a dead copy of an
earlier README.

**The data.** `runs/octinoxate/` holds everything for this ligand: `poses/` and `energies.csv` from
the search, `design*.json` and `condense_*.csv` from the ordering, `boltz/` with the folds and
`fold_check*.csv`, `figures*/` with readable copies and PyMOL scripts, `md/` with one directory per
leg and `mmgbsa_summary.csv`, `dock/` with the redocking tables, `gnina/` with the rescoring bundle
and `results_*/gnina_scores.csv`. `uma_logs/` and `logs/` are split by whether UMA ran, because a
UMA log is the only record of how an energy was arrived at and regenerating one means re-running the
hours that produced it.

**Full trajectories are off this drive, on iCloud.** The large full-system trajectory files — every
`traj.dcd`, `traj_ext.dcd`, the MMPBSA `traj.nc` intermediates and the `checkpoint*.chk` (13 GB, each
hours of GPU) — were moved to **iCloud Drive at `Dynamics_trajectories/peptidebuilder_POC/`**, in
subfolders mirroring the run names here (e.g. `peptidebuilder_POC/orig_f12_dock5/prod_20ns/traj.dcd`),
then evicted locally to free space. They are never in git (`.gitignore` excludes them) and were never
deleted. What stays in the repo is enough to redo any analysis without them: the wrapped solute
trajectory (`traj_wrapped.xtc`), the prmtops, `energy.csv`, `FINAL_RESULTS_MMPBSA.dat` and the frame
PDBs. To re-run MM/GBSA or re-window a leg from the full trajectory, pull its folder back from iCloud
first (in Finder, or `brctl download`).

### Why this file did not exist until now

Worth recording, because the cause will recur. Every document above was written for the *next
session* at the moment it was written, and the organising axis in all of them is chronology —
shell 1, then shell 2, then shell 3; "written 2026-09-24", "added 2026-09-30", "computed
2026-10-02". Nothing was reorganised by subject when a later finding superseded an earlier one, and
nothing was retired when it was discharged. The pipeline then changed underneath the documents:
UMA's static scoring was the whole judgement stage when README.md and PROGRESS.md were written, and
it is now retired as a predictor, while the Vina–GNINA–MD stages that replaced it are the ones
missing from the README's file tables. The two peptide families ended up in two repositories with
one README each, so the comparison that is the actual result — four designs and two nulls against
two generated binders — existed only as a table two thirds of the way down a 1,580-line file, and
the ΔG values it rests on were split across two repositories and a hardcoded dict.

Keep this file first and current: add a leg to section 3 and section 4 as it is computed, move a
claim from section 7 to section 6 when a replicate earns it, and let RESULTS.md stay what it is good
at, which is the long-form argument with the figures.
