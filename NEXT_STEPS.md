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
(README, [redocking](README.md#an-independent-check-on-the-pose-redocking)). Two findings bear on the
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

## Blocked on Modal credit: October 2026

The Modal allocation is spent for September 2026 (about $4 across scoring and six 20 ns runs). When it
renews in **October 2026**, in rough priority order:

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

3. **Dynamics on the docked poses.** This is the one experiment that would settle which pose is right,
   and it is cheap. Redocking (README, [redocking](README.md#an-independent-check-on-the-pose-redocking))
   puts the ligand 3–10 Å from where Boltz put it, mostly by turning it end for end in the same pocket,
   and the docked poses are as wrapped and as enclosed as the predicted ones. That comparison is
   circular — Vina selected them for close packing against a rigid receptor — so it cannot be resolved
   statically. Run each structure's best-scoring docked pose through the same MD and MM/GBSA path the
   Boltz poses took and compare residence, release episodes and ΔG against the numbers already in
   `md/mmgbsa_summary.csv`. Three outcomes, all informative: the docked pose holds as well or better, in
   which case the orientation genuinely does not matter and the pocket is degenerate as the static
   analysis suggests; it is released quickly, in which case Boltz's orientation is right and Vina's
   score is not reading anything real about this system; or it drifts into the Boltz pose, which would
   be the strongest possible confirmation of the prediction.

   Start with **`bg33_4`**, which is now the obvious first case: it is the only structure of eight whose
   predicted pose Vina reproduces (1.30 Å, its pose 6), so it is the one place where "run the docked
   pose" and "run the predicted pose" are nearly the same experiment and any difference in retention is
   attributable to the 0.52 Å and 7° between them rather than to a different binding mode. Then
   `s3_orig_f12` and `shuffle_control` — the design/null pair already rendered side by side, the two
   extremes of docked-pose promiscuity (1/9 flipped at Jaccard 0.70 against 6/9 at 0.42), and the pair
   whose residence numbers differ most (100% with zero releases against 75.9% with five episodes). `cif_to_md.py` will not do the input prep here, because the input is a pose in a
   `.pdbqt` rather than a Boltz `.cif`: write the docked pose out as an SDF in the reference's atom
   order with `vina_redock.ref_order_poses`, and reuse the existing `protein_fixed.pdb`, which is the
   same receptor. Two 20 ns runs, roughly $2 of Modal credit, or locally as `orig_f12`'s was.

Shell 4 (`isoleucine/2`) was the plan before this and is still unstarted; it needs no Modal credit for
design and folding, only for dynamics. Shell 1's per-fold designed-position reproduction is also still
unmeasured, and needs per-fold shell `.xyz` copies.

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

## Traps

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
