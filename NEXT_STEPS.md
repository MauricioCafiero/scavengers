# Where this stands and what to try next

Written 2026-09-24. Everything below describes the built-in `octinoxate` ligand, which is
**C17H24O3, one CH2 short of real octinoxate** (2-ethylhexyl 4-methoxycinnamate, C18H26O3). Its
SMILES, perceived from the stored geometry, is `CCCC[C@H](CC)OC(=O)/C=C/c1ccc(OC)cc1` — the ester
oxygen sits directly on the chain CH. Every result here is for that analogue.

## The result that motivates the next step

The pipeline builds a *shell*: fragments placed all around the ligand, every side chain in contact.
A folded peptide cannot reproduce it.

| Measure | Value |
|---|---|
| Designed side-chain positions reproduced within 3 A | **8 / 225 (3.6%)** |
| Ligand superposition quality (so the comparison is sound) | 0.07-0.86 A, mean 0.54 A |
| Residues of a folded peptide touching the ligand | 3-6 of 12-21 |
| Fragments touching the ligand by construction, in a design | all of them |

This holds for every variant tried: linear, `--linker-slack` 0-3, cyclic, and ESM2 linker
substitution.

Substitution deserves its own note, because Boltz and UMA disagree about it. Boltz's dG improved in
4 of 4 pairs (mean 0.22 kcal/mol), but the UMA interaction energy on the same structures went one
better, one flat, one *positive*, and one much worse:

| design | parent | substituted |
|---|---|---|
| `DKSIFEDIKWGEGLGYGLS` | -7.21 | **-11.67** |
| `RYGLSGIKWDKSFEGDGGE` | -10.05 | -10.24 |
| `SLWYGKGDGFLGYSGEGGR` | -12.47 | **+0.78** (broken interface) |
| `cyclo-LGSGFGIGGYGWGEGGKGG` | -16.20 | **-6.66** |

Residue overlap did not improve either: 2 hits within 3 A across the four parents, 0 across the four
substituted versions. So substitution improves Boltz's score while leaving or degrading the actual
complex — better binders by one metric, not better realisations of the design. It also shifts net
charge a lot (one peptide went 0 to +4 as ESM2 inserted lysines and arginines, another -1 to -3).

Against a null model there is no evidence of added value: at length 19 the designs give UMA
interaction energies of -16.2, -12.5, -10.1, -7.2 kcal/mol (mean -11.5) and five random sequences of
the same length give -40.8, -9.5, -7.7, -5.1, +12.1 (mean -10.2). n is far too small for
significance, but the best random peptide is 2.5x better than the best design.

## What redocking already settled, and what it did not

Added 2026-09-30. The six structures with dynamics have now been redocked with AutoDock Vina
(RESULTS, [redocking](RESULTS.md#an-independent-check-on-the-pose-redocking)). Two findings bear on the
proposals below:

* **The pocket does not specify the ligand's orientation.** 18 of 54 docked poses are turned end for
  end while keeping the same contact residues, and flipped and unflipped poses score −5.18 against
  −5.30 kcal/mol. Wrapping and enclosure saturate either way round, so the design objective is blind to
  it. If either proposal below is meant to place the chromophore head specifically, the objective needs
  a term that distinguishes the two ends.
* **A static score predicts retention, though not energy.** Across all eight redocked structures —
  the six here plus BoltzGen's `bg33_4` and `bg33_3` — Vina's score correlates +0.732 (p = 0.039) with
  release episodes and +0.61 (p = 0.11) with MM/GBSA ΔG, and the docked pose's enclosure correlates
  −0.713 (p = 0.047) with release. Both release correlations fell when the BoltzGen pair was added
  (+0.845 and −0.833 at n = 6), which is the expected cost of breaking a nearly binary variable, and
  both survived.
* **GNINA's `CNNaffinity` is the first static measure here that predicts the binding energy.** Rescoring
  the eight structures' poses with GNINA v1.3.3 (`--score_only`, so the CNN sees exactly Vina's
  coordinates) gives ρ = **−0.857, p = 0.007** against MM/GBSA ΔG on the predicted pose, strengthening to
  −0.900 on the five S-configuration structures alone. Vina's score manages +0.61 (p = 0.11), nesso
  +0.14, and the static UMA energy ranks backwards. The pose-ranking improvement is smaller than it
  looks: GNINA puts the predicted pose 1st to 4th where Vina put it 6th to 9th, but in `bg33_4` it ranks
  a pose 9.3 Å away above the one 1.3 Å away, so part of the gain is that a co-folded geometry looks like
  a crystal structure. Controlling for that by using docked poses only, a real but modest position signal
  survives — `CNNaffinity` ρ = −0.400 (p = 0.0005) against RMSD, where Vina had none — and it is the CNN,
  not the empirical term, that carries it. The flip degeneracy survives both scorers: GNINA separates
  flipped from unflipped at p = 0.283.
* **The measure worth pursuing is pose reproducibility, not the score.** How closely the closest docked
  pose reproduces the predicted one correlates **−0.743 (p = 0.035) with residence**, the strongest
  single relationship found. `bg33_4` redocks to 1.30 Å — the only structure of eight that redocks at
  all by the 2 Å standard — and holds 100% residence with zero releases; `bg33_3` redocks to 4.47 Å and
  has the worst residence measured in either project, 41.0% with 18 episodes. Neither extreme existed in
  this repository's own six, so this could not have been seen before the BoltzGen pair was added. If it
  holds it is a free screen: dock the predicted pose back into its own fold and ask how well the best
  pose reproduces it, with no dynamics needed. It also separates two things `bg33_4` shows are not the
  same — Vina *scores* it 7th of 8 while MM/GBSA puts it 3rd, so the number is wrong while the pose
  agreement is right. Next test: three or four more structures with intermediate release, since three of
  the eight are still tied at zero episodes and one pair supplies both ends of the correlation.
* **Three folds of eight contain the mirror-image ligand, and it is not a BoltzGen quirk.** The source
  `ligand.xyz` is **S** at the 2-ethylhexyl carbon; `shuffle_control_esm0` (this repository's own Boltz
  run), `bg33_3` and `bg33_4` came out **R**, the other five kept S. The SDFs carry the fold geometry
  unchanged to 0.0005 Å, so the folding models placed it that way. Within-structure results are
  unaffected — every RMSD and wrapping number compares a pose to its own reference — but cross-structure
  score and ΔG comparisons mix enantiomers, and anything chirality-aware applied later will be more
  sensitive to it than Vina was. Worth checking the other 30 folds for how often this happens, which is
  free: perceive the CIP label from each `.cif` ligand and compare against the source. Check the label,
  not the canonical SMILES — the source carries `/C=C/` double-bond stereo that a cif-derived SDF does
  not, so a string comparison flags all eight as different and hides the three that are.

## The information the code currently throws away

Each fragment is a side-chain analogue capped where the backbone would attach — the capping
methyl/H marks roughly where CB, and hence CA, would sit. **Nothing in the pipeline uses that.**
`build_sequence` chains fragment *centroids* by nearest neighbour and converts gaps to glycine with
`round(d / 3.8) - 1`, so a design can require a backbone path that no peptide can follow. That is
the most likely explanation for the 3.6%, and it is what both proposals below exploit.

## A. Fragment condensation

Keep the fragments where the search put them and join them chemically instead of writing a sequence
and asking a folder to rediscover the geometry.

1. For each selected pose, identify the attachment atom and its direction (the capping methyl or H).
2. Between consecutive fragments in the chosen order, add the backbone atoms (N, CA, C, O) needed to
   connect their attachment points, inserting glycine-like spacer residues where the gap needs them.
3. Build the resulting covalent peptide in place, ligand fixed.
4. Relax with UMA and read the strain.

The strain energy after relaxation is the answer: low strain means the shell was realisable and the
structure is now in hand; high strain means it never was, and *where* the strain localises names the
fragment pair that cannot be connected. There is no sequence-to-structure step, which is the stage
we have measured failing, and UMA needs no parameterisation for whatever geometry results.

Untested. Note that linking fragments into a lead is routine in small-molecule design, so any
novelty claim rests on doing it with peptide backbone chemistry judged by an ML potential, from
ML-scored side-chain poses. Worth a literature check first.

## B. Reachability-constrained assignment

Feed A. Rather than nearest-neighbour chaining of centroids, choose which subset of the 170 poses to
use and in what order, maximising total interaction energy subject to backbone reachability between
consecutive attachment points — both the distance and the direction the attachment vectors point.
The output is a shell that is known to be traversable before any compute is spent folding it. The
current chaining has no notion of traversability at all.

## The plan

Decided 2026-09-24. The order is **not** fixed as B then A, because B needs a reachability criterion
and guessing one would bake in an unvalidated assumption. A can measure it instead:

1. **Shared first piece — attachment geometry.** Every fragment is a side-chain analogue capped with
   a methyl where the backbone would attach, so that methyl carbon marks roughly where CB sits and
   the bond direction into it points where CA would go. This is the information the current pipeline
   discards. Found by graph search over each fragment: a carbon with three hydrogens and exactly one
   heavy-atom neighbour. Verified for all ten fragments:

   | fragment | atoms | capping CH3 index |
   |---|---|---|
   | arginine | 19 | 15 |
   | lysine | 17 | 12 |
   | aspartic | 7 | 2 |
   | glutamic | 7 | 3 |
   | isoleucine | 14 | 2, 9 |
   | leucine | 14 | 0, 5, 6 |
   | serine | 6 | 2 |
   | tryptophan | 19 | 15 |
   | tyrosine | 17 | 10 |
   | phenylalanine | 15 | 11 |

   Isoleucine, leucine and valine-like fragments have several terminal methyls, so the attachment one
   is ambiguous from the graph alone and needs picking by chemistry (the one that corresponds to CB),
   not just by the CH3 test. Leucine's three candidates and isoleucine's two must be resolved before
   either method can use them.

2. **Calibrate A on pairs.** Condense two placed poses with k spacer residues between them, relax
   with UMA, and record strain against attachment-point distance, the angle between the two
   attachment vectors, and k. Pairs are small systems, so this sweep is cheap, and it yields an
   empirical reachability criterion rather than an assumed one. It is also the direct diagnostic for
   why designs do not fold as drawn.

3. **B with the measured constraint.** Solve the pose subset and ordering against the criterion that
   step 2 produced.

4. **A on the full selection.** Condense and relax the whole peptide.

A later iteration could feed A's strain results back into B's constraints, dropping fragment pairs
that prove infeasible and re-solving.

C below is **not** being pursued.

## Superseded ordering note: B feeding A

Decided 2026-09-24: build **B then A** as one pipeline. B chooses a pose subset and ordering that a
backbone can actually follow; A condenses the peptide onto those poses and relaxes it with UMA. B
must run first because A needs a connectable ordering to build along. A later iteration could feed
A's strain results back into B's constraints, excluding fragment pairs that prove infeasible and
re-solving the assignment.

C below is **not** being pursued.

## C. Grow the peptide, not the fragments (not pursued)

A third option, in the spirit of the repo's own `grow_fragments`: start from the best pose and add
one residue at a time, enumerating backbone dihedrals, scoring each candidate side-chain placement
against the fragment-derived energy map, keeping a beam. The chain is connected at every step, so
connectivity is never a post-hoc problem.

Rejected as too close to existing tools: backbone-first sampling with an external builder,
restrained MD toward the designed positions, and pharmacophore-style screening of folded candidates.

## Ligand strain: read HANDOFF.md section 10 before quoting any strain number

Redefined 2026-09-26. Strain is now `E(bound) - E_ref` against **one shared reference per ligand** --
the ligand's own lowest conformer, from `code/ligand_reference.py` (20 ETKDG embeddings, MMFF-ranked,
five lowest relaxed with UMA), stored in `runs/<ligand>/ligand_reference.json` and reused by every
shell. `code/strain_global.py` applies it to each fold as a single point on the bound ligand as the complex
relaxation left it -- that state must not be re-relaxed, not even its hydrogens; see the wrong turn
recorded in HANDOFF.md section 10. `binding_energy.py` no longer needs to compute `strain_ligand` at all.

The old definition relaxed each complex's own bound pose, so every structure was measured against a
different local minimum. It was sensitive to the force cutoff (+1.5 to +7.7 kcal/mol between fmax 0.1
and 0.01), to the step cap (a generous one lets the ligand change conformer -- up to 1.0 A of
heavy-atom drift), and to `pdbfixer`'s non-deterministic hydrogens (3.6 kcal/mol on one structure).
Interaction energy was never affected, because it compares the complex with its own parts at one fixed
geometry and the hydrogen error cancels.

All 24 existing folds were recomputed; `runs/octinoxate/boltz/strain_global.csv` holds the result.

## Queued on Modal credit: renewed October 2026

The September 2026 allocation was spent (about $4 across scoring and six 20 ns runs). Credit renewed
**2026-10-01**; item 5 is the first thing drawing on it. In rough priority order:

1. **More dynamics.** Six trajectories is what every conclusion in the README's dynamics section rests
   on, and it is too few: the best predictor found, mean simultaneous side-chain engagement, moved from
   +0.90 to +0.83 when the sixth was added and got that structure's pair backwards. The cheapest useful
   additions are the remaining shell-3 rungs and a second shuffle seed, since a single design against a
   single shuffle cannot put a number on what the ordering search is worth. Budget roughly $1 per 20 ns
   run for a compact fold, $1.10 for an extended one.
2. **A second ligand, whole pipeline, probably without scoring.** Every number in this repository comes
   from octinoxate, so nothing here is known to transfer. The case for skipping UMA scoring is now
   empirical rather than a matter of taste: across the four designs with MM/GBSA, interaction energy and
   interaction-plus-strain both rank structures at Spearman **−0.80** against the free energy, while
   costing 45–60 minutes per fold. Folding plus `check_fold.py` is free and ranks at +0.30 -- no better
   than chance, but not actively inverted, and it does reliably identify the folds where the ligand is
   not bound at all. So: design, fold, filter on geometry, and spend the compute on dynamics for the two
   or three survivors rather than on scoring all of them. Keep `binding_energy.py` for the folds that go
   to MD, where the decomposition is worth having beside the MM/GBSA one.

3. **Dynamics on a docked pose — done for `s3_orig_f12`, 2026-09-30; superseded by item 5.** The
   +3.55 below is measured between two runs whose ligands moved 2.51 Å and 2.78 Å before production
   began (item 6), so it compares two relaxed positions rather than a docked pose against a predicted
   one. Comparable drift makes this the best case in the set, not a clean one. The window series
   survives better than the endpoint. Vina's best-scoring pose was run
   for 20 ns under an identical protocol and binds **3.55 kcal/mol worse** than the predicted pose
   (−20.79 ± 0.04 against −24.33 ± 0.08), holding about eighteen fewer contacts and loosening over the
   run where the predicted pose tightens. Both keep the ligand: no releases either way. Full account in
   RESULTS, [redocking](RESULTS.md#dynamics-on-a-docked-pose-the-predicted-pose-binds-better).

   Three things that run left open. The estimate was **still moving at 20 ns** (+0.97 kcal/mol over the
   last 5 ns, 24× its standard error), so −20.79 is an upper bound and a longer run would tighten the
   comparison. The predicted pose's own value came from a Modal production leg with no windows, so its
   convergence is unknown — a 30 or 40 ns leg on both poses would fix both problems at once. And **where
   the ligand ends up cannot be measured on this system**: the peptide superposition RMSD is 5.89 Å,
   larger than the ligand distances being compared, and the control is decisive — the predicted-pose run
   reads 6.57 Å from its own starting pose. Do not spend more effort on that measure here.

4. **Run the pose GNINA picks, not the pose Vina picks.** Superseded by item 5, which widens this to
   three structures and both scorers; the reasoning below still stands for the `s3_orig_f12` cell.
   GNINA and Vina disagree about the best docked pose in **8 of 8 structures**: Vina's pick is
   always its own pose 1, GNINA's are poses 9, 6, 7, 6, 7, 9, 5 and 6. The pose already simulated,
   `s3_orig_f12` pose 1, is GNINA's **7th of 10** — so the experiment so far tested the pose the better
   scorer thinks is poor.

   For `s3_orig_f12`, GNINA's pick is **pose 9**: `CNNscore` 0.733 and `CNNaffinity` 5.166, both the best
   of the nine docked poses, against pose 1's 0.401 and 5.010. It sits 4.08 Å from the prediction, almost
   exactly as far as pose 1's 4.09 Å, but differs in how it got there — 110° of rigid-body rotation
   against pose 1's 167°, and 48° of head/tail change against 64°. So the comparison is close to
   controlled: near-identical displacement, different orientation, and the two scorers disagree about
   which is better. If pose 9 beats pose 1's −20.79 the CNN's pose preference is worth something beyond
   its ΔG correlation; if it does not, `CNNaffinity` predicts affinity without preferring poses, which is
   still useful but a narrower claim.

   Everything needed is in place: `code/dock_pose_to_sdf.py s3_orig_f12 --pose 9` writes the input, the
   receptor is the same `protein_fixed.pdb`, and `code/run_dock_pose_md.sh` needs only its `N` changed.
   Roughly 4.5 h locally, or an hour of Modal credit. Worth pairing with **pose 8** if there is budget:
   at 3.49 Å and only 27° of rotation it is the mild-perturbation control that would separate how much of
   the 3.55 kcal/mol penalty is reorientation and how much is displacement.

5. **Starting-pose selection: the two-scorer matrix — RUN 2026-10-01, results in
   [RESULTS.md](RESULTS.md#the-two-scorer-pose-matrix).** Item 4 is folded into this. Five 20 ns legs on
   Modal with MM/GBSA local. Extended the same day to `bg33_3` and `s2_esm2_control`: nine legs in all,
   **$10.10** of GPU, 23 rows in `mmgbsa_summary.csv`. Of that, ~$0.89 went to a preempted leg that
   restarted from zero and ~$1.30 to cubic boxes (see traps).

   **This is methodology, not a result about the designs.** It asks which starting pose is worth 20 ns
   for a peptide you already have. Which peptides are worth having is the MM/GBSA-and-retention
   comparison in [RESULTS](RESULTS.md#what-the-pipelines-deliver-mmgbsa-and-retention), where this
   pipeline's four designs average −18.84 kcal/mol at 84.1% retention with one release, against
   BoltzGen's two at −15.76, 70.5% and eighteen. Nothing in this item bears on that.

   Each cell measures where a given starting pose leads under 20 ns, not whether one pose geometry is
   better than another — the geometric labels below (flip, roll) describe the inputs and do not survive
   equilibration (item 6). **Five rows are now run, and ordered by the co-folded pose's ligand retention
   they give a usable rule.** Where the co-folded pose holds the ligand completely, co-folding is the
   better start; where it holds it poorly, a docked start wins, by a margin that grows as retention falls:

   | peptide | retention | best route | margin over co-folding |
   |---|---|---|---|
   | `s3_orig_f12` | 100% | co-folded | — |
   | `bg33_4` | 100% | co-folded | — |
   | `s2_esm2_control` | 59.6% | GNINA p6 | **0.84** |
   | `bg33_3` | 41% | GNINA p6 | **1.23** |

   `s3_esm2_f4` is excluded: all three of its legs moved 7.7–15.1 Å before production. So is
   `s2_esm2_control`'s Vina cell, at 12.6 Å. `s2_esm2_control`'s GNINA cell is the strongest in the set
   on its own terms — 1.58 Å of drift, the best-held docked start anywhere here — so its −17.09 against
   the co-folded −16.25 is attributable to the pose rather than to where it wandered.

   **Working rule: dock only when the co-folded pose is poorly retained.** Four supporting cells, no
   counterexample among legs that held their poses. Caveat: two structures each side, and the two
   co-folded winners are also the two best binders, so retention and binding strength are not separated.

   What the item does establish is about the scorers: GNINA's `CNNaffinity` ranks structures against
   MM/GBSA at rho = −0.857, and across five rows its pick beats Vina's in four of them. It is still not
   dependable: its one failure is on `s3_orig_f12`, the project's best binder, where its pick lands 7.56
   kcal/mol below the pose it ranked 7th of 10 — the largest gap in the table.

   The original plan, kept because the reasoning still applies to the cells:

   The picks are unambiguous. Vina chooses its own pose 1 in all three structures, which is what
   ranking by the score the poses were sorted on amounts to. GNINA's picks by `CNNaffinity` are poses
   9, 6 and 6, and they are stable: the `score_only` and `minimize` passes agree on all three, and
   `CNNscore` agrees except for `bg33_4` under `score_only`, where it prefers pose 7 by 0.012.

   `bg33_4`'s GNINA cell is a replicate, and is run anyway — deliberately, as the matrix's one
   reproducibility control. Measured against the predicted pose it sits 0.52 Å away by centroid, 8.0° in
   head-to-tail direction, and 14.2° of rigid-body rotation — the same binding mode by every measure, on
   a ligand whose long axis is 12.9 Å, so the ends move about 1.6 Å. Equilibration erases that in
   picoseconds. The 1.301 Å symmetry-minimised RMSD that made it look like a distinct pose is barely
   below the plain atom-indexed 1.43 Å, so no automorphism was hiding a difference either. What it
   therefore measures is not a pose comparison but agreement between two independent paths to the same
   geometry: `bg33_4`'s predicted pose gives −19.66 ± 0.02 kcal/mol from BoltzGen's co-folded structure
   on Modal, and this run reaches the same pose by docking, prepares it in this repository, and runs it
   under the matrix protocol. How far apart those two numbers land is the error bar that belongs on every
   other cell, and nothing else in the project supplies it. It goes last, because it is the one cell
   whose result can be anticipated.

   What the six cells are, measured against each structure's predicted pose:

   | run | centroid | head-to-tail | rigid rotation | what it tests |
   |---|---|---|---|---|
   | `s3_esm2_f4` pose 1 (Vina) | 2.48 Å | **152.2°** | 167.0° | the flip, under dynamics |
   | `bg33_4` pose 1 (Vina) | 2.27 Å | 8.6° | 104.5° | orientation alone: a pure roll |
   | `s3_esm2_f4` pose 6 (GNINA) | 2.20 Å | 34.2° | 149.6° | the flip's unflipped partner |
   | `s3_orig_f12` pose 9 (GNINA) | 3.04 Å | 26.5° | 118.4° | item 4's controlled comparison |
   | `bg33_4` pose 6 (GNINA) | 0.52 Å | 8.0° | 14.2° | replicate: protocol reproducibility |
   | `s3_orig_f12` pose 1 (Vina) | 2.57 Å | 17.8° | 164.9° | **done**: −20.79 against −24.33 |

   Run the five in that order, which is by information gained rather than by structure. The first two are
   the ones that answer something no static measure can. `s3_esm2_f4`'s two picks are an end-for-end
   flipped pose and an unflipped one at nearly the same displacement, which is a direct test of the
   degeneracy GNINA could not separate statically (p = 0.283) and that redocking showed the design
   objective is blind to. `bg33_4` pose 1 holds position and head-to-tail direction nearly fixed and
   rolls the molecule ~105° about its own long axis, which is the cleanest available test of whether
   orientation by itself costs binding energy.

   Head-to-tail is the angle between the vectors joining the ligand's two most distant atoms (indices
   6 and 17 of the 20 heavy atoms); it detects an end-for-end flip, which RMSD on a molecule this
   elongated partly absorbs. Rigid rotation is the Kabsch angle after centering, which detects a roll
   about the long axis that head-to-tail cannot see. Neither is in the bundle CSVs; both are computed
   from `runs/octinoxate/gnina/<structure>/poses_with_reference.sdf`, whose poses share atom indexing
   with the reference, so the comparison needs no matching.

   These do not match item 4's angles for the same two poses, and the discrepancy is unresolved. Item 4
   reports 110° and 167° of rigid-body rotation for `s3_orig_f12` poses 9 and 1 where the Kabsch
   measure above gives 118.4° and 164.9°, and 48° and 64° of head/tail change where the two-most-distant-
   atoms measure gives 26.5° and 17.8°. The rotation figures are close enough to be an alignment or
   hydrogen-inclusion difference; the head/tail figures are far enough apart to be a different
   definition, and item 4's was not recorded. Before either set is quoted outside this file, find item
   4's definition and keep one. The displacements agree (4.08 and 4.09 Å, both the symmetry-minimised
   RMSD in `reference.csv`), so nothing in the pose selection depends on this.

   **Split execution: dynamics on Modal, MM/GBSA local.** The trajectories go to Modal — October credit
   has renewed, and at roughly $1 per 20 ns run the four cost about what six did in September. MM/GBSA
   runs here afterwards on the downloaded trajectories, as `s3_orig_f12_dock1` did. Two things that
   follow from splitting it. The production leg has to come back whole, not just its summary: anything
   expensive to recompute gets downloaded, and a ΔG with no trajectory behind it cannot be rewindowed.
   And MMPBSA.py writes `reference.frc` plus a `_MMPBSA_*` set into the working directory, about 6 GB
   per structure, so run each in its own scratch directory and check staged files by size rather than
   by extension before cleaning up.

   One prep trap: `bg33_4`'s inputs are split across two repositories. Its docked poses are here in
   `runs/octinoxate/dock/bg33_4/`, but its predicted-pose MD lives in
   `~/python_mac/boltzgen_local/md/bg33_4/` under that repository's system preparation. Holding the
   protocol identical across the matrix means preparing it here from `dock_pose_to_sdf.py` like the
   others, not reusing that run's system.

   **Recommended first, and not yet decided: backfill windows on the two Modal baselines.** Two of the
   three predicted-pose ΔG values have no drift estimate, because Modal runs skip the time checks:
   `s3_orig_f12` (−24.33) and `s3_esm2_f4` (−21.08). `bg33_4` is not among them — it ran locally in
   `boltzgen_local` and has all four, and they matter, because they show what a converged run on this
   system looks like: **−19.18, −18.89, −19.43, −19.66**, a range of 0.77 kcal/mol over 20 ns and only
   −0.23 over the last 5 ns.

   Set that against `s3_orig_f12_dock1`: **−27.00, −24.04, −21.76, −20.79**, monotone, 6.2 kcal/mol of
   drift over 20 ns and +0.97 in the last 5. So 20 ns is not too short for this system in general, and
   the drift is a property of the docked starting pose rather than of the protocol — a docked pose starts
   over-contacted and decays toward equilibrium, a co-folded pose starts at it. That makes the shape of
   the window series a result in its own right, not just a convergence check: if the five new runs drift
   the way `dock1` did and the one replicate does not, the drift is diagnostic of a pose that was never
   at equilibrium.

   Backfilling the two Modal baselines is feasible: their wrapped trajectories hold the solute MM/GBSA
   needs (421 and 670 atoms) across 2000 frames over the full 20 ns, so it is MM/GBSA on frame subsets
   of files already on disk, no new dynamics. The frame density is lower than a local run's 15,000, which
   costs roughly 0.11 kcal/mol of standard error against `dock1`'s 0.039 — negligible beside a 3.55
   effect. `run_windows_live.sh` does it.

   **Where the rows go: six rows, one per cell.** A Modal run gets MM/GBSA on the production leg only,
   not the four time checks, so each cell contributes a single row to
   `runs/octinoxate/md/mmgbsa_summary.csv` and the matrix is six rows. The windowed four-row shape
   belongs to the local runs (`orig_f12`, `s3_orig_f12_dock1`, and `bg33_4` in `boltzgen_local`), which
   is why the existing Modal rows all read `prod_L1_modal` and stand alone. New rows take
   `prod_L1_modal` too, since the leg names where the trajectory came from.

   Structure names follow the existing `{structure}_dock{pose}` convention, giving `s3_esm2_f4_dock1`,
   `s3_esm2_f4_dock6`, `bg33_4_dock1`, `bg33_4_dock6` and `s3_orig_f12_dock9`. The pose number is in the
   name; which scorer picked it is not, and the table above is the only place that mapping is written
   down. `s3_orig_f12_dock1`'s four rows were missing from the CSV and were added 2026-10-01 from its own
   `FINAL_RESULTS_MMPBSA.dat` files; its `prod_20ns` row is the −20.7851 quoted in item 3, and it is the
   row to compare against, not its windows.

   Append each row as its run finishes rather than batching all five at the end. Nothing prevents a
   batch — the file is flat and no code parses it — but a ΔG that exists only in a MMPBSA output
   directory is one cleanup away from being gone, and the runs are hours apart.

   Because MM/GBSA runs locally here, windows are available for any of these cells after the fact, at
   the cost of three more MM/GBSA passes over a trajectory already on disk. They are not part of the
   plan, but the option survives the run in a way it does not when MM/GBSA happens on Modal.

   **The duplication trap.** Nothing reads `mmgbsa_summary.csv`. `make_gnina_bundle.py` names it in a
   comment but carries the numbers in a hardcoded `REFERENCE` dict, so the CSV and that dict are two
   hand-maintained copies of the same values and can drift silently. A new row has to be added in both
   places, or the next GNINA bundle will be built against stale ΔG values. Worth collapsing: have
   `make_gnina_bundle.py` read the CSV, and the problem disappears.

   Also from item 3, now qualified: 20 ns converges a co-folded pose on this system but did not converge
   the one docked pose tried, so read a drifting cell as a ranking rather than a free energy, and check
   the window shape before quoting any of them as converged.

6. **The ligand is unrestrained during equilibration, and what that means for reading item 5.** Found
   2026-10-01 while checking the replicate. Not a bug: `openmm_md/dynamics.py:144` builds the positional
   restraint from `protein_heavy_indices(topology)`, so the receptor is pinned and the ligand and solvent
   relax. That is a deliberate and conventional choice, and for most of this project's dynamics it is the
   right one — a docked pose is scored by a docking function, not by this force field, and freezing it
   would be worse.

   What it costs is the ability to attribute a result to the geometry of the input pose. Measured from
   each system's input to its production frame 0, superposed on the peptide:

   | run | peptide RMSD | ligand RMSD, input → production start |
   |---|---|---|
   | `bg33_4` co-folded | 0.88 Å | **1.15 Å** |
   | `bg33_2` co-folded | 0.82 Å | 1.68 Å |
   | `bg33_4_dock1` | 0.85 Å | 2.30 Å |
   | `s3_orig_f12_dock1` | 1.15 Å | 2.51 Å |
   | `s3_orig_f12` co-folded | 1.51 Å | 2.78 Å |
   | `s3_orig_f12_dock9` | 2.12 Å | 4.24 Å |
   | `bg33_4_dock6` | 1.26 Å | **4.82 Å** |
   | `s3_esm2_f4` co-folded | 1.17 Å | **7.66 Å** |
   | `s3_esm2_f4_dock1` | 1.73 Å | **11.05 Å** |
   | `s3_esm2_f4_dock6` | 2.24 Å | **15.11 Å** |

   The peptide barely moves, so this is the ligand sliding, not the receptor refolding; minimum
   ligand–peptide distance stays near 2 Å, so nothing unbinds and no periodic image jumps.

   So the matrix compares **starting points and where they lead**, which is exactly what a dock-then-MD
   workflow does and is the practically useful question. Every cell remains a valid measurement of that.
   What cannot be read off it is any claim of the form "the flipped pose binds better than the unflipped
   one", because the flip does not survive to production — the geometric labels in item 5's table
   describe the inputs, not the trajectories.

   **Do not change the equilibration to fix this.** Fifteen runs on disk all ran with the ligand free,
   and their value is being comparable to each other; a restrained leg would not belong in the same
   table. What is worth doing is recording input → production ligand RMSD on each leg as a diagnostic.
   It changes nothing about how a leg runs and tells you afterwards whether a cell is worth reading —
   which is how `s2_esm2_control_dock1` at 12.57 Å came to be discounted rather than taken at face
   value.

7. **Window the trajectories already on disk. Nice to have, after other exploratory work.** One item,
   not two: MM/GBSA over the leading 5, 10 and 15 ns of a trajectory, which is three more passes over
   frames already here and no new dynamics. It applies to the five matrix legs and to the two Modal
   baselines (`s3_orig_f12` at −24.33, `s3_esm2_f4` at −21.08) alike; `bg33_4`'s co-folded run already
   has windows, and `bg33_4_dock6` was windowed on 2026-10-01 (−16.72, −16.90, −16.36, −16.53 — flat,
   which is how we know it is a converged measurement of a *different* pose rather than one caught
   mid-drift). Those three window rows reached `mmgbsa_summary.csv` on 2026-10-02, and their directories
   were renamed `win_Nns` → `first_Nns` to match the convention `orig_f12` and `s3_orig_f12_dock1` use.

   The value is that the shape survives what the endpoint does not: a series still walking away from its
   start at 20 ns was not at equilibrium, whatever its absolute number turns out to be. No GPU and no new
   dynamics -- a few minutes of local CPU per window on the strided trajectories -- so it can wait, but
   it should happen before any matrix cell is quoted.

8. **The first MM/GBSA replicate. Not now -- after more data points.** Added 2026-10-02, and recorded
   here so the reason survives: nothing in this project has ever been run twice from the same input, so
   there is no run-to-run error bar on any MM/GBSA number. Every comparison rests on single legs and the
   differences being read are 3 to 4 kcal/mol -- four designs averaging −18.84 against BoltzGen's two at
   −15.76 and the two nulls at −14.73. The only two near-identical starting points ever run, `bg33_4`
   co-folded and `bg33_4_dock6` at 0.52 Å apart, ended **3.13 kcal/mol** apart, which bounds how finely
   two *poses* can be told apart by one leg each but says nothing about protocol noise, because those two
   genuinely relaxed into different basins (1.15 Å against 4.82 Å of pre-production drift). MM/GBSA's own
   standard errors, 0.04 to 0.10, measure frames inside one trajectory and are not uncertainty on a
   comparison.

   The legs available now are better spent on new structures, which is why this waits. When it runs, run
   it on **`s3_orig_f12`** -- the structure every headline number is anchored to, and the cheapest
   replicate available:

   ```sh
   modal run code/modal_md.py::ls      --structure s3_orig_f12   # confirm the Volume still has it
   modal run code/modal_md.py::produce --structure s3_orig_f12 --leg L2
   ```

   `_produce` reuses the `system.xml` and `complex.pdb` already on the Volume, so there is no `omd
   build`: box, water count and the 5,615 particles are identical to the original leg, which removes the
   one thing that would otherwise differ between them. `dynamics.py` passes no seed to
   `LangevinMiddleIntegrator` or to `setVelocitiesToTemperature`, so velocities and the Langevin stream
   are independent and the result is a genuine repeat rather than a rerun. About $0.68 and an hour, then
   MM/GBSA locally as with any Modal leg; the row goes in `mmgbsa_summary.csv` as `prod_L2_modal`.

   What it decides: whether a 3 kcal/mol pipeline difference is a result or a single-leg artifact. Within
   about 1 kcal/mol of −24.33 and the comparisons already written stand; 3 or more and they have to be
   restated as what one leg each showed.

9. **Test `CNNscore` as the pose picker.** Added 2026-10-04. The switch to `CNNaffinity` was made
   because it tracked the co-folded MM/GBSA ranking, but the newest legs undercut it: `shuffle_control`'s
   top-affinity pose 9 (`CNNaffinity` 3.97, best of poses 1–9) bound **worst** of that peptide's three
   starts (−11.15 at 46.4% residence, against −16.15 for Vina p1 and −15.13 co-folded, 2026-10-04), and
   `orig_f12`'s GNINA p5 leg releases more than the nulls do (RESULTS, residence table). Two cheap ways
   to test the other CNN head on data partly in hand:
   - Tabulate `CNNscore` against the 20 ns ΔG of every docked leg already simulated — **done
     2026-10-04, first half of this item** (`results_receptor_{minimize,h_score_only,score_only}/
     gnina_scores.csv` × `mmgbsa_summary.csv`). Finding: on **7 of 8** structures the `CNNscore`-top
     pose of 1–9 *is* the `CNNaffinity`-top pose in all three sets — every docked leg run so far is
     already a `CNNscore`-top pose, so the two heads are indistinguishable on existing data. Across
     peptides, `CNNscore` vs ΔG is r = −0.50 (n = 8), but that weights design-vs-null geometry, not pose
     choice; it does not license or condemn either head as a picker.
   - Run the one robust disagreement: `shuffle_control`, whose `CNNscore` best is **pose 7** in all
     three sets (0.460 / 0.461 / 0.485) while its affinity ranks near-last (3.57) — **done
     2026-10-04, second half of this item** — and its affinity-top p9 was already measured at
     −11.15 / 46.4%. The p7 leg ran on the Reading ARC (first cluster leg, job 490801): ΔG −8.40
     (windows −12.49 / −9.56 / −8.99, the wrong-pose shape), residence **22.4%**, 24 episodes,
     2.42 ns detached — worst of the structure's four starts on both metrics. **Item closed:
     `CNNaffinity` keeps the picker job, `CNNscore` is a cross-check that on this fold actively
     picks worse than affinity (not interchangeable — worse).** The failures also add one more
     "proxy to not read alone": p7 makes 31/66 pair contacts (the structure's best) while the
     ligand drifts across the peptide surface. (The minimise set's affinity-top p8, 4.085, remains
     an optional parked leg — nothing left unsettled by it that this test would change.)

10. **A new molecule, through both pipelines.** Added 2026-10-04. Everything in this repo is built on
    the octinoxate analogue, so the next question is transfer, and it needs both design paths on the same
    new ligand: the peptidebuilder shell design (peptide_builder) and the BoltzGen co-folding design
    (`~/python_mac/boltzgen_local`), on the same molecule, same target metric set as octinoxate got —
    co-fold, dock (Vina + GNINA), dynamics on the top poses by the criterion item 9 settles. This
    concretises the earlier "second ligand, probably without scoring" plan above: skip the UMA scoring
    stage (it ranks backwards), spend the budget on dynamics for the two or three folds whose geometry
    survives `check_fold.py`. The ligand choice is the one real decision: small, flexible, a couple of
    acceptors — different enough from the ester chain that transfer is being tested, similar enough that
    the docking box setup transfers. Everything else (dodecahedron, windows, MM/GBSA locally) follows the
    now-fixed recipe.

    **Ligand chosen 2026-10-04: oxybenzone** (benzophenone-3,
    `COc1cc(C(=O)c2ccccc2)ccc1O`). The first repo hard-coded four built-in ligands — octinoxate,
    octocrylene, oxybenzone, avobenzone — and the choice between them was made on Morgan fingerprints
    (radius 2, 2048 bits) against the octinoxate analogue the project actually uses
    (`CCCC[C@H](CC)OC(=O)/C=C/c1ccc(OC)cc1`): oxybenzone **0.17**, octocrylene **0.19**, avobenzone
    **0.28**. Oxybenzone wins on the transfer logic too — no ester, no long chain, rigid and planar
    (3 rotatable bonds) against octinoxate's flexible 2-ethylhexyl ester, one phenol donor alongside
    the three acceptors, MW 228 against 276 — while staying small enough that the docking box and the
    MD recipe transfer unchanged. Rigid rather than flexible, which trades half the stated criterion;
    tested transfer is the point, so the rigidity is the feature.


## Environment map

Nothing extra is installed in this repo; each external tool is called from its own environment.

| Need | Where |
|---|---|
| This repo's code | `.venv` (Python 3.12, fairchem 2.22, `uma-s-1p2p1`) |
| Boltz-2 co-folding | `~/python_mac/boltz_local`, via its `code/boltz_mps.py` (two passes, so the structure and affinity models do not share 8 GB) |
| pdbfixer (adds H to the peptide) | `~/python_mac/pocket_assist/venv`, `--fixer-venv` |
| ESM2 / transformers | `~/python_mac/GenMaskFill/.venv`, `--genmask-venv` |
| SMILES to 3D recipe this follows | `~/python_mac/mace/code/mace_calc.py` (`smiles_to_atoms`) |
| AutoDock Vina + Open Babel (redocking only) | `~/python_mac/dock_assist/dock-env`, run as `dock-env/bin/python code/vina_redock.py`. Vina is the x86_64 build vendored in `dockstring`, which needs Rosetta on Apple silicon; `obabel` comes from `openbabel-wheel` in that venv, so nothing is installed system-wide |
| **MD build / run / analyze / MM/GBSA** | `~/miniforge3/envs/openmm-md`, as `~/miniforge3/envs/openmm-md/bin/omd`. Code lives in `~/python_mac/openmm/src/openmm_md`. This **must be conda**, not a venv: `omd build` needs the OpenFF/GAFF2 stack to parameterise the ligand and `omd mmgbsa` needs AmberTools, and neither is pip-installable. Also provides the `mdtraj` used by every drift/QC calculation — `.venv` does not have it |
| **Modal GPU dynamics** | `code/modal_md.py`, entrypoints `push` / `probe` / `produce` / `fetch`. Needs a Modal token (`~/.modal.toml`) and creates the Volume `peptidebuilder-md-state` on first use. It mounts `~/python_mac/openmm/src` from the local disk, so that repo must be present even though the GPU work is remote |
| **BoltzGen designs (`bg33_*`)** | `~/python_mac/boltzgen_local`. Their receptors (`md/bg33_*/protein_fixed.pdb`) and co-folded baseline trajectories live there, not here — `prep_matrix_poses.sh` and `run_row4_md.sh` read them by absolute path |

### Resuming on a different machine

Checked 2026-10-02, when this was about to be cloned onto a laptop with a better GPU.

* **Paths are now derived, not hardcoded** (fixed 2026-10-02). Every `code/run_*.sh` and
  `prep_matrix_poses.sh` takes the repo root from its own location (`REPO=${0:A:h:h}`) and reads the
  externals from overridable variables: `OMD_ENV` (default `~/miniforge3/envs/openmm-md`),
  `BOLTZGEN_ROOT` (default `~/python_mac/boltzgen_local`), and `DOCK_ASSIST` (default
  `~/python_mac/dock_assist`, used by `dock_pose_to_sdf.py` and `make_gnina_bundle.py`). The scripts
  check `OMD_ENV` exists and exit with a readable message rather than failing midway. So a clone works
  from wherever it lands, and only the three env vars need setting if the companion repos are not in
  `~/python_mac/`.
* **Full trajectories do not travel.** `.gitignore` excludes `traj.dcd` / `traj.nc`, by design. What is
  committed is the wrapped solute trajectory (2,000 frames), the prmtops, `energy.csv` and
  `FINAL_RESULTS_MMPBSA.dat` — enough to redo any MM/GBSA or window analysis, not enough for a warm
  restart. The full `.dcd` files also still sit on the Modal Volume, so they can be re-fetched there.
* **A better GPU changes the economics.** Everything so far split dynamics onto Modal because this
  laptop has no CUDA. Measured rates for comparison: A10G ran 0.264 ms/step at 14,105 particles and
  0.362 at 23,142; this laptop's OpenCL managed 1.44 ms/step at 5,615. If the new machine beats the
  A10G, `omd run --platform CUDA` locally removes the $1-per-leg cost, the preemption risk, and the
  2,000-frame stride that the Modal return path imposes.
* **What is ready to run immediately:** `shuffle_control`, `shuffle_control_esm0` and `orig_f12` all have
  tracked docking inputs and receptors, so their Vina p1 and GNINA top cells need no new preparation.
  GNINA's picks for them are p8, p9 and p5 by `CNNaffinity`.

## Traps

**`omd build` defaults to a cube; every baseline here is dodecahedral.** Pass
`--box-shape dodecahedron` explicitly. The seven legs run 2026-10-01 inherited the cube default while
six of the seven co-folded runs they are compared against were built dodecahedral — same physics and no
effect on MM/GBSA, which strips to the solute, but ~30% more waters and about $1.30 of avoidable GPU
across them. The waste scales with how elongated the solute is and was invisible on the compact systems:
$0.09 a leg on `bg33_3` at 7,384 atoms, $0.41 on `s2_esm2_control`, whose solute spans 47.6 Å in one
dimension and built to 33,057 atoms cubic against 23,142 dodecahedral. Unlike the equilibration
restraints, switching to dodecahedral *restores* comparability with the baselines rather than breaking it,
so it is safe to fix. Done in `prep_matrix_poses.sh`, `run_row4_md.sh` and `run_row5_md.sh`.

**Modal leg progress is unobservable mid-run.** The volume's `energy.csv` is a stale early flush and the
container prints nothing between start and finish, so there is no way to read a step count while a leg
runs. Any "percent complete" is a guess. Cost is observable — `modal billing report --for today
--show-resources` gives accrued spend per app — and the leg prints its own `minutes` and `usd` when it
finishes.

- **fairchem 2.22+ refuses the UMA 1.0 checkpoint** (`uma-s-1`), telling you to install
  `fairchem-core<=2.21.0`. `uma-s-1.pt` is cached in `~/.cache/fairchem` if you ever need it; the
  current code uses `uma-s-1p2p1`.
- **fairchem does not support MPS** — `device must be either 'cpu' or 'cuda'`, so everything is CPU.
- **Boltz writes heavy atoms only**, so hydrogens must be added before any energy calculation.
- **Boltz's affinity head barely discriminates and is not geometry-consistent.** Its dG spans
  1.65 kcal/mol across 37 sequences, and it gave `RYGLGEFSDWKI` pIC50 6.96 / binder probability 0.75
  with the ligand 6.6 A away — not bound. UMA on the same coordinates gives -1.35. Always read a
  contact distance next to an affinity.
- **A positive UMA interaction energy means the predicted structure is broken**, not weakly bound.
  5 of 33 interfaces came out positive and 4 of those have a heavy-atom contact under 2.6 A (shortest
  1.63 A). Use it as a filter.
- **UMA interaction energy as computed grows with size** (Spearman -0.52 against peptide length): it
  is a rigid gas-phase energy with no desolvation or entropy. Normalise (per contact atom or buried
  area) before ranking anything by it. Boltz's dG has the opposite length bias (+0.68), which is why
  the two anti-correlate (-0.47) — mostly an artifact, not a chemical disagreement.
- **Refuted explanations**, so nobody re-tests them: the best random binder is *more* charged than
  the designs (37% vs 23%), and interaction energy does not track contact area (+0.17 vs contact
  pairs, +0.42 vs buried ligand atoms — wrong sign; two structures of near-identical burial differ by
  35 kcal/mol). At fixed length the energy is set by contact quality, and what drives that is
  unexplained.
- **ESM2 linker filling must sample, not take the argmax.** With this many masks the 35M
  checkpoint's per-position distribution is nearly flat (confidence 0.06-0.10) and the argmax
  collapses every linker to leucine. It can also change net charge (one fill took a peptide from -1
  to -3 by inserting aspartates), so consider constraining the residue set.
- **A substituted sequence has no shell of its own.** `overlay.py` needs a design structure — the
  ligand-plus-fragment-poses `.xyz` in `runs/octinoxate/sequences/` — so copy the parent's, since the
  poses are identical and only linkers changed.

## Keeping the machine awake

**An assertion is armed right now: `caffeinate -is -t 28800`, PID 4534, started 2026-09-24 23:47
local, 8 hours.** It is deliberately left running. Disarm it when the working session ends:

```sh
pgrep -f 'caffeinate -is'        # find it (ignore the Bash tool's own short -t 300 wrapper)
kill <pid>
pmset -g assertions | grep -E 'PreventSystemSleep +[01]'   # 0 once released
```

This machine sleeps mid-conversation, not only mid-job, which stalls detached work and interrupts
the session, so the assertion has to cover the whole working session rather than an individual job.
Bind it to a timeout, never to a job's PID: an assertion bound with `-w <pid>` dies when that job
ends and leaves everything afterwards unprotected, which is the mistake that let the machine sleep
earlier in this session.

`-i` prevents idle sleep, `-s` prevents system sleep and needs AC power. The Bash tool wraps its own
commands in a short `caffeinate -t 300` which expires by itself — that is not the one to kill.

## Repo

- Work is pushed to **github.com/MauricioCafiero/scavengers**, which is now the only remote. `main`
  tracks `scavengers/main`, so a bare `git push` goes to the right place. The `lucia-71/peptidebuilder`
  remote has been removed, and the worktree dance previously needed to keep the inherited history out
  is no longer necessary.
- **No inherited material is in this history.** The 12 commits authored by the original author were
  dropped when `main` was rebased onto the scavengers root (2026-09-25); the history now begins at
  your own `Initial commit` and every commit is yours. Verified by blob hash across both trees: no
  file here is byte-identical to any file of hers.
- `data/` and `notebooks/` stay out. Both are inherited: `data/` is 197 MB of the original project's
  MD trajectories, and `notebooks/frag_grow.ipynb` is her file, which one of our commits merely
  relocated. **Checking which commit added a file to our history does not establish its origin** —
  check whether the content appears in her tree. Restoring the notebook on the strength of the wrong
  test is a mistake already made once.
- `runs/octinoxate` is committed and other `runs/*` are ignored. Logs are split by engine:
  `runs/octinoxate/uma_logs/` for the 21 runs where UMA ran (the only record of how long each
  relaxation took and whether it converged), `runs/octinoxate/logs/` for the 10 with no UMA in them
  (geometry sweeps, Boltz folds, analysis). All 31 sat loose in `runs/` until 2026-09-26, kept by a
  `!runs/*.log` negation that is now only a safety net for a stray. No code writes into either
  folder — every log path but Boltz's comes from a shell redirect, so choose the right one.
- Credit to the original author stays in the README and in any publication. Keeping her commits out
  of this repository's git history is separate from, and does not affect, that attribution.
