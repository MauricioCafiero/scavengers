# Where this stands and what to try next

**Active right now (2026-10-27, 19:50):** MM/GBSA replicates — tally **88/107** via
`code/replicate_tally.py` (stronger 30 / weaker 28 / unchanged 30, |ΔΔG| < 1.5; median |ΔΔG|
2.50). Batch 9 (six legs) recorded and committed: s3_orig_f12_rf3cof_dock1_r2 +7.03,
s3_esm2_f4_rf3cof_dock3_r2 −0.34, s2_esm2_control_rf3cof_dock1_r2 +2.06, s3_esm2_control_r2 −1.73,
s3_orig_f12_dock9_r2 +7.65 (vs its Modal run 1), shuffle_control_dock7_r3 (third run: −8.40 racc /
−24.75 / −24.91 Young — enrolled in `THIRD`).

- **Full sweep of the remaining unreplicated canonical legs started** (user 2026-10-27: "over 100
  runs we still need to replicate"): **129** at the start (59 octinoxate / 70 oxybenzone, read off
  both `mmgbsa_summary.csv` files minus the tally's PAIRS). First six on Young — jobs 116195
  `sweepY1` (bg33_1_o3cof_dock1/dock6), 116196 `sweepY2` (bg33_1_rf3cof, _rf3cof_dock1), 116197
  `sweepY3` (bg33_1_rf3cof_dock8, bg33_2_o3cof), `_r2` each, PENDING on Resources; each has its own
  drain on `runs/folds_md/young_<job>.legs`. **123 remain** (53 oct, 70 oxy); oxybenzone needs its
  `runs/oxybenzone/md` tree staged on Young first. Packing: first-fit in the inventory's order, ≤6 h
  of Young model estimate per job (the model over-predicted 115841 by ~1.3×, so it is a safe
  `--time` basis).
- **15:00 2026-10-28: the last 53 oxybenzone legs went to Young** (user: "load them on young" --
  the balance NOT already in racc's rep7_oxy_1-3). `~/openmm/runs/oxybenzone/md` created on Young,
  318 files md5-matched; sweepY19-Y32 = 117342-117355 (13 x 4 legs + 1, `RUNREL=runs/oxybenzone/md`,
  `DCD_DIR` scratch), three `afterany` chains behind Y16 (116337), Y17 (116338) and a fresh one
  (Y21, PENDING on Priority -- other users now ahead); one drain each with `SYSNAME=oxybenzone`,
  `DR_LEGDIR='runs/oxybenzone/md/${leg}/prod_20ns'`. **Every canonical leg is now enrolled: tally
  148/230, 0 unreplicated.** racc (stalled since its 10:40/10:48 preemptions behind another user)
  still holds rep7_oct_1 (3 legs left) + rep7_oxy_1/2 queued, rep7_oxy_3 held by feed_rest.sh.
  Measured Young per-leg wall today: 67-91 min (sacct). Read job state from sacct, not from when a
  leg's MM/GBSA lands -- that lags the job end by 30-60 min of fetch + tail.
- **18:06 2026-10-29: Young SCRATCH quota (250 GB) filled with full dcds** -> Y25, Y28, Y29, Y30,
  Y32 (17 legs) FAILED on `Errno 122`; Y26's last leg finished just before and is complete. User
  rule: never keep dcds once the wrapped solute is verified local -> 91 verified dcds deleted
  (scratch 250 -> 4 GB). **Resubmitted:** Y25=119187, Y29=119188 (heads), Y28=119189 afterany
  Y25, Y32=119190 afterany Y29, Y30=119191 afterany Y31 (117354, untouched, pending). Legfiles
  keep the old 1173xx IDs; their drains are still running. Delete each Young dcd once recorded.
- **03:15 2026-10-28: Young home quota (100 GB) filled with full dcds** -> Y7 leg 3, Y8 legs 3-4
  aborted and Y10-Y17 (116252-116259) FAILED in 1-2 s; queue empty 02:34-03:15. 25 dcds moved to
  `/scratch/mmm1166/dcd` (symlinked; nothing deleted), all submit scripts now `DCD_DIR=` scratch.
  **Resubmitted:** Y10=116330 Y11=116331 Y12=116332 (start now), **Y18=116333** (the 3 aborted
  legs: orig_f12_o3cof_dock1_r2, orig_f12_o3pcof_dock4_r2, orig_f12_rf3cof_dock5_r2), Y13=116334
  Y14=116335 Y15=116336 Y16=116337 Y17=116338 chained afterany as before. The legfiles keep the OLD
  job ids in their names (`young_116252..9.legs`, and Y18's legs sit in `young_116215/6.legs`);
  their drains key on xtc paths, not job ids, so they need no restart. See YOUNG.md.
- **Overnight cover (21:20, user: "make sure the queues don't empty out overnight"):** 31 more
  sweep legs staged and queued, **85 unreplicated remain**.
  racc: `runs/folds_md/feed_rest.sh runs/folds_md/rep7.chain` (feed_rest now takes chain files as
  args) feeds `rep7_oxy_1/2/3` (sweep R2-R4, oxybenzone bgox31_2..5, 5.4-5.8 h est each) behind
  rep7_oct_1, one per freed slot; each has its own drain (SYSNAME=oxybenzone).
  Young: sweepY4/5/6 = 116212/116213/116214 (submitted now) and sweepY7/8/9 = 116215/116216/116217
  (`--dependency=afterany` on Y4/Y5/Y6), octinoxate; `~/openmm/submit_sweepY<n>.sh`; all 84 staged
  files md5-matched; one drain per `young_<job>.legs`. Session caffeinate re-armed `-t 43200` 21:20.
  **22:55:** Young ran ~2x faster than the model on big legs (YOUNG.md: overhead-bound, ~85 min a
  leg flat on a good node, ~150 on a slow one), so Y7-Y9 would have drained by ~04:00. Queued the
  **whole rest of octinoxate** (32 legs) as sweepY10-Y17 = 116252-116259, 4 legs each,
  `--time=7:10:00` (not padded for slow nodes -- user: small requests are for backfill), in three
  `afterany` chains behind Y7/Y8/Y9; one drain each. **Octinoxate is fully enrolled; 53
  unreplicated remain, all oxybenzone** (17 in rep7_oxy_1-3). Young drains before 22:25 polled racc
  by mistake (`DR_SOCK` `:-` bug, fixed + restarted). Next racc batch at `bgox31_5_rf3cof` (oxy);
  Young's next work needs the oxybenzone tree staged there.
- **racc (21:01):** rep6_oct_3 548596 RUNNING, rep6_oxy 548701 + **rep7_oct_1 548907** PENDING.
  rep7_oct_1 = sweep R1, the first racc sweep batch (user 2026-10-27 "go ahead and submit"): the 7
  small bg33_3/bg33_4 co-fold legs, 5.6 h est, `--time=7:07:00`; its drain (new script) runs on
  `runs/folds_md/rep7_oct_1.legs`. `feed_rest.sh` is finished -- the next racc batch is packed by
  hand off `replicate_tally.py --unreplicated` (racc rate model), named `rep7_oct_<n>`/sweep R<n>.
  rep6_oct_2 all four legs MM/GBSA'd 20:04-20:55 (+2.01, -0.09, +0.97, +3.22); records queued.
- **racc (earlier):** feeder on rep6_oct_2 (RUNNING) / rep6_oxy / rep6_oct_3 (PENDING). Socket back 19:46. rep6_oct_2 legs 1–2
  (`bg33_4_rf3cof_dock1_r2`, `_dock4_r2`) done on scratch, re-fetching (dock1's first fetch died at
  16:48 leaving a 0-byte xtc, removed so the drain re-fetches).
- **Young at 19:48:** 116195 RUNNING (started 19:46:28, a minute after submission); 116196/116197
  PENDING on `Resources` — all 40 A100s busy, and they are 2 of the only 3 jobs on the partition
  with that reason (the other ~150 pending are `Priority`/`Dependency`, behind us), so they start on
  the next free GPU.
- **Local queue:** empty — every tail and record through s3_esm2_control_r2 is done.
- **Pipeline changes 2026-10-27 (committed in batch 9; described in CLAUDE.md):**
  `record_leg.py` upserts only its own leg's rows (five records ran in 30 s, not ~3 min each);
  records take `/tmp/peptidebuilder_record.lock`; `racc_drain.sh` stages fetches and moves the
  xtc in last. The racc drains started before 19:10 keep the old copy; the Young drains have both.
- **Waiters:** `watch_legs --on mmgbsa` on rep6_oct_2/3/oxy and on the three young_*.legs. Report
  off the MM/GBSA, not the record.
- **racc's rep5_oct_6 guard-missed bg33_3_r2** (scratch cleaned, so it recomputed); the leg was
  already recorded, the second trajectory stays unfetched, no double rows.
- **Next commit (batch 10):** the next ~5 recorded legs (rep6_oct_2's four + the first sweep
  legs). The Young-vs-racc ΔΔG read-out waits for ~10 Young legs (user ruling).
  The `runs/folds_md/replicate_ddG_*.png` stay uncommitted until ~10 more legs, then regenerate
  (via the `replicate_tally` import) and commit.

**Resume checklist (written 2026-10-27 19:55, before a context reset).** The drains are detached
and survive a reset; the waiters do not — they were tasks of the old session.
1. `date`; `pgrep -fl 'caffeinate -i -s'` — the session caffeinate (pid 33291, `-t 21600`) expires
   ~00:56 on 2026-10-28. Re-arm detached if it is gone:
   `~/miniforge3/bin/python -c "import subprocess; subprocess.Popen(['caffeinate','-i','-s','-t','21600'], start_new_session=True)"`.
2. Drains that should be alive (`ps ax | grep racc_drain`; short-lived duplicates are forked poll
   children): **racc** `rep6_oct_2`, `rep6_oct_3`, `rep6_oxy` (started 2026-10-23 on the OLD
   script: direct fetch, records on the tail lock — compare each fetched `traj_wrapped.xtc` size
   with the cluster's, CLAUDE.md "wedged fetch"); **Young** `young_116195/6/7` (new script). A drain
   exits by itself once all its legs are drained or flagged.
3. Re-arm the waiters (each exits after one report or 120 min; re-arm after every report):
   `~/miniforge3/envs/openmm-md/bin/python code/watch_legs.py --on mmgbsa --legfile runs/folds_md/rep6_oct_2.legs --legfile runs/folds_md/rep6_oct_3.legs --legfile runs/folds_md/rep6_oxy.legs --legfile runs/folds_md/young_116195.legs --legfile runs/folds_md/young_116196.legs --legfile runs/folds_md/young_116197.legs`
   (`--selftest` first), backgrounded. Report each leg off its whole-run MM/GBSA as it lands:
   dG, windows, signed ΔΔG vs run 1, release/residence from `tail.log`; slide and tally after the record.
4. Queue state: racc `squeue -u $USER` through the socket (`ssh -o ControlPath=~/.ssh/cm-racc
   racc.rdg.ac.uk`), Young `ssh young-ng 'squeue -u $USER'`. A Young `Resources` reason means
   "next for a free A100", not a long queue (YOUNG.md).
5. **Next Young batch** (the full sweep): `python3 code/replicate_tally.py --unreplicated` lists
   the remaining canonical legs (123 at 19:55: 53 oct, 70 oxy) with true particle counts. Take the
   next ones in order; pack first-fit to ≤6 h of `8 + 10M × 0.466 ms × (N/7327)^e` min (e 0.4 ≤13k,
   0.64 above), `--time` = estimate + 1 h 30; stage `<leg>_r2` from the parent's `system/`,
   `protein_fixed.pdb`, `ligand_prepped.sdf`, `*_ligand.sdf`; tar each leg to
   `young-ng:~/openmm/runs/<mol>/md/` and **md5-compare every file** (a dropped stream left a
   truncated `system.xml` on 2026-10-27); sbatch from `~/openmm/submit_sweepY1.sh` (copy, change
   `-J`, `--time`, `LEGS`); write `runs/folds_md/young_<jobid>.legs`, start its drain with the
   YOUNG.md `DR_*` recipe at submission; add the structures to `PAIRS` in `replicate_tally.py`
   as `sweep Y<n>`. Oxybenzone legs also need `RUNREL=runs/oxybenzone/md`, `SYSNAME=oxybenzone`,
   `DR_LEGDIR='runs/oxybenzone/md/${leg}/prod_20ns'`, and the `~/openmm/runs/oxybenzone/md` tree
   created on Young (absent so far).
6. Commit every ~5 recorded legs (tables + each leg's `md_contacts.csv` + legfiles), additions-only
   `git diff --numstat`, and **push each commit immediately** (`git push scavengers
   "${c}:refs/heads/main"`, one commit at a time). Do not let commits accumulate unpushed — seven
   did on 2026-10-27 and the user asked why on earth they were saved up.

**Third run of `shuffle_control_dock7`: produced and tailed.** r1 −8.40 (racc), r2 −24.75 (Young),
**r3 −24.91 (Young)** — the two Young draws agree within 0.16 and racc's run 1 is the outlier.
When its record lands, enroll the structure in `THIRD` (`code/replicate_tally.py:133`) so the
triple column lights up. The racc duplicates of the two Young batches (547254/547288) were
canceled the same day.

**The accounting workbook (2026-10-25).** `PB_runs_accounting.xlsx` at the repo root tracks every
run by peptide/fold-type. Its fill/audit machinery now lives in `code/workbook/`
(`j_pairs.py` → `premd_pairs.csv`; `fill_accounting.py` fills empty cells only; `audit_gihjlm.py`
and `audit_wb.py` are independent audits — run all three with `.venv/bin/python`). Protocol: fill
EMPTY cells only, overwrite only with an audit-proven mismatch; the workbook is untracked in git
(the user's Excel copy is authoritative); writes only happen while Excel has the file closed.
Column sources and semantics are in the scripts' docstrings — key facts: pre-MD pairs
denominator comes from `pair_contacts_*.csv`'s `ceiling` column (66 = C(12,2) oct, 45 = all oxy
guided rows), never from adjacent/nonadjacent cell text; GNINA pose keys are lane-aware
(`o3gninaN`/`rf3gninaN`/bare `gninaN`); generator arms (`bg33_*`/`bgox31_*`) get NO pair cells
(2026-10-22 ruling). G for the BoltzGen cofold rows comes from the **new `boltzcf` lane**
(pepgraph `LANES` now includes it; scored 2026-10-25: perp 8.0–13.3, z +2.0 to +5.9 for all 8;
pepgraph's `compare_lanes.py` rewrites `lane_comparison.json` wholesale per run, so score a new
lane to a separate `--out` and merge rows back, or full-run).

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
  **Done at campaign scale 2026-10-18 (n = 63 structures + 115 docked legs): the headline does not
  survive — pooled ρ = −0.27, and oxybenzone alone −0.14 (n.s.), because OF3/RF3 refold folds
  routinely redock closely and still lose the ligand; the score, not the geometry, is the surviving
  cheap signal.** Full account in [POSE_REPRODUCIBILITY.md](POSE_REPRODUCIBILITY.md); tables at
  `runs/<mol>/md/pose_reproducibility*.csv`, script `code/pose_reproducibility.py`.
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

## ARCHIVED 2026-10-22: the fragment-condensation direction (A, B, C) is not being pursued

**Do not re-propose any of this.** The user's ruling, 2026-10-22: *"we're way past using methods A
and B, so you can archive them in that document but don't bring them up again because we're not
doing it."* The three sections that follow -- the discarded attachment geometry, A (fragment
condensation), B (reachability-constrained assignment), the four-step plan built on them, and C
(already rejected in 2026-09) -- are kept as the record of a direction the project chose not to
take. The campaigns that followed went a different way: design, fold, filter on geometry, and
spend compute on dynamics, which is what produced the 228-leg MM/GBSA record. Nothing below is an
open task, a next step, or a suggestion to revive.

## The information the code currently throws away (archived)

Each fragment is a side-chain analogue capped where the backbone would attach — the capping
methyl/H marks roughly where CB, and hence CA, would sit. **Nothing in the pipeline uses that.**
`build_sequence` chains fragment *centroids* by nearest neighbour and converts gaps to glycine with
`round(d / 3.8) - 1`, so a design can require a backbone path that no peptide can follow. That is
the most likely explanation for the 3.6%, and it is what both proposals below exploit.

## A. Fragment condensation (archived)

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

## B. Reachability-constrained assignment (archived)

Feed A. Rather than nearest-neighbour chaining of centroids, choose which subset of the 170 poses to
use and in what order, maximising total interaction energy subject to backbone reachability between
consecutive attachment points — both the distance and the direction the attachment vectors point.
The output is a shell that is known to be traversable before any compute is spent folding it. The
current chaining has no notion of traversability at all.

## The plan (archived)

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

## Superseded ordering note: B feeding A (archived)

Decided 2026-09-24: build **B then A** as one pipeline. B chooses a pose subset and ordering that a
backbone can actually follow; A condenses the peptide onto those poses and relaxes it with UMA. B
must run first because A needs a connectable ordering to build along. A later iteration could feed
A's strain results back into B's constraints, excluding fragment pairs that prove infeasible and
re-solving the assignment.

C below is **not** being pursued.

## C. Grow the peptide, not the fragments (not pursued; archived)

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

## Queued on Modal credit: renewed October 2026 — items 1 and 2 SUPERSEDED, 3-5 done

**Status 2026-10-22.** Items 1 and 2 were overtaken by the campaigns that followed and are not
open work: item 1 ("six trajectories is too few") is answered by the **228 legs** now carrying a
canonical whole-run MM/GBSA row, and item 2 ("a second ligand, whole pipeline") is **oxybenzone**,
which has had its own full set of tables since 2026-10-14. Items 3-5 carry their own DONE markers
below. Modal is also no longer the main compute lane -- racc's `gpuscavenger` is, and it is free --
so a Modal leg is now the exception rather than the queue; per-leg spend is recorded in
`runs/<mol>/md/modal_receipts.csv`. The item text is kept because the *reasoning* in 2 (skip UMA
scoring, fold and filter on geometry, spend compute on dynamics for the survivors) is still the
operating rule.

**What 3-5 leave unanswered is NOT a next step** (user, 2026-10-22: *"something that we will
resolve while looking at the data, not necessarily a next step but just things to consider down the
line"*). Treat them as caveats to apply when reading the pose-selection results, not as work to
queue or propose: (a) the `s3_orig_f12` docked comparison was **still moving at 20 ns** (+0.97
kcal/mol over the last 5 ns, 24x its standard error), so its -20.79 is an upper bound and the 3.55
margin against the predicted pose is not a settled number; (b) the predicted pose's own value came
from a Modal leg with no convergence windows, so its own settling is unmeasured; (c) pose 8, the
mild-perturbation control that would separate how much of that penalty is reorientation and how
much is displacement, was never run. Quote the margin with (a) and (b) attached rather than
proposing legs to close them.

The original list, for the record:

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
   equilibration (item 6). **The rows now run, ordered by the co-folded pose's ligand retention, give a
   usable rule.** Where the co-folded pose holds the ligand completely, co-folding is the better start;
   where it holds it poorly, a docked start wins. But the winning docked scorer is not always GNINA: on
   the two BoltzGen extended designs the Vina rank-1 pose beats both co-folding and GNINA's pick, so the
   margin does not fall cleanly with retention the way the GNINA-only rows had suggested:

   | peptide | retention | best route | margin over co-folding |
   |---|---|---|---|
   | `s3_orig_f12` | 100% | co-folded | — |
   | `bg33_4` | 100% | co-folded | — |
   | `s2_esm2_control` | 59.6% | GNINA p6 | **0.84** |
   | `bg33_2` (racc) | 50.9% | Vina p1 | **3.80** |
   | `bg33_3` | 41% | GNINA p6 | **1.23** |
   | `bg33_1` (racc) | 14.5% | Vina p1 | **2.16** |

   `bg33_1` and `bg33_2` ran on the Reading ARC (free) rather than Modal, 2026-10-04/05. For `bg33_1`
   the co-folded start reads −11.31 at the project's worst residence (14.5%); GNINA's pick (pose 2,
   `CNNaffinity` 5.40, `minimizedAffinity` −3.65) reads −12.55 at 42.5%; and **Vina's rank-1 pose reads
   −13.47 at 44.7% residence** — the best of the three, drifting only 1.52 Å before production (table in
   item 6), so the cell is attributable to the pose. `bg33_2` runs the same order: co-folded −8.02,
   GNINA p4 −10.44, **Vina p1 −11.82 at 43.1%** (1.50 Å drift). So on both BoltzGen extended designs a
   docked start wins and Vina's rank-1 is the best docked pose — each also the steadiest fold of its
   design's three legs (Cα RMSD 1.6–2.1 Å). The `bg33_1` Vina window series still climbs over the run
   (−11.09 → −11.63 → −13.30 → −13.47), so its margin could move with a longer leg.

   `s3_esm2_f4` is excluded: all three of its legs moved 7.7–15.1 Å before production. So is
   `s2_esm2_control`'s Vina cell, at 12.6 Å. `s2_esm2_control`'s GNINA cell is the strongest in the set
   on its own terms — 1.58 Å of drift, the best-held docked start anywhere here — so its −17.09 against
   the co-folded −16.25 is attributable to the pose rather than to where it wandered.

   **Working rule: dock only when the co-folded pose is poorly retained.** Six supporting cells (bg33_2
   joined 2026-10-04: co-fold −8.02 at 50.9% vs GNINA p4 **−10.44 ± 0.05 at 57.2%**, margin 1.30), no
   counterexample among legs that held their poses. Caveat in force: the cells ran on two platforms
   (racc legs are local-unstrided; Modal legs strided), run-to-run variation exceeds MM/GBSA standard
   errors, and the three co-folded winners are also three of its best binders, so retention and binding
   strength are not separated.

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
   | `bg33_1_dock2` (GNINA p2, racc) | 0.08 Å | **0.38 Å** |
   | `bg33_2_dock4` (GNINA p4, racc) | 0.54 Å | **2.30 Å** |
   | `bg33_4` co-folded | 0.88 Å | **1.15 Å** |
   | `bg33_2` co-folded | 0.82 Å | 1.68 Å |
   | `bg33_4_dock1` | 0.85 Å | 2.30 Å |
   | `s3_orig_f12_dock1` | 1.15 Å | 2.51 Å |
   | `s3_orig_f12` co-folded | 1.51 Å | 2.78 Å |
   | `bg33_1` co-folded | 0.84 Å | 3.46 Å |
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

    **Design generation to complete (added 2026-10-05).** The peptidebuilder side is done — three shells
    (charged `ox1`, aromatic `ox2`, hydrophilic `ox3`), folded, with Vina + GNINA bundles and the
    co-fold / Vina-p1 / GNINA-pick MD legs running on RACC. Still to generate, for the full two-pipeline
    comparison on oxybenzone:
    - **Four BoltzGen candidates** — co-folded designs via `~/python_mac/boltzgen_local` (MPS,
      `precision=32` mandatory, ~100 s/design; four local patches), then taken through the same
      co-fold -> dock (Vina + GNINA) -> dynamics matrix as the shells, so BoltzGen and peptidebuilder are
      compared on identical metrics (as octinoxate's `bg33_*` were).
    - ~~**Two null candidates** — shuffle controls~~ **DONE 2026-10-06.** `ox2_shuffle` and
      `ox2_shuffle_esm0`, from `design_shell2.json` (seed 4, 0 of 11 side chains left in their designed
      slot), through the full matrix on racc: co-fold, Vina p1 and GNINA pick, six 20 ns legs, all five
      oxybenzone tables filled. **The design wins both pairs — 18.9 kcal/mol on the glycine pair and 8.8
      on the ESM2 pair, and the null loses in all six cells.** The two pairs fail differently and only the
      second isolates the arrangement: `ox2_shuffle`'s fold never enclosed the ligand to begin with
      (enclosure 0.58), while `ox2_shuffle_esm0` holds it 94.7–97.4% of the run, posts the highest
      simultaneous engagement in the whole set, and still loses by 8 — so it also **breaks the `sim`
      predictor**, which ranks that pair backwards. Use `dead_slots` instead (3–5 of 11 against 0).


## Replicate every MM/GBSA leg: run-to-run variation — PARKED 2026-10-22 at 23 structures

**The result and the current state live in `REPLICATES.md`, not here.** A targeted set of 23
structures (10% of the campaign's 228 original MD legs) ran on racc 2026-10-19 → 2026-10-22 and
the campaign was then parked by the user, with the expectation of returning to it. Headline:
stronger 8, weaker 9, unchanged 6; median |ddG| 4.80, range 0.02–16.16; two third runs, neither of
which resolved its pair. The full ~202-leg sweep described below was NOT run.

**There is exactly ONE open question here: what to do about replicates.** How many structures get
a second run and how a replicated cell is written in the published tables are **not two decisions**
-- they are one, and the user is still thinking about it (2026-10-22). Do not split them, do not
ask about them separately, and do not edit a single table cell in the meantime: a form chosen for
one peptide would set a precedent for the other 22 by accident, which is why every affected row
across RESULTS.md, FOLDS.md and cofold_binding_tables.md was left untouched. The size of that edit,
measured 2026-10-22: **100** rows mention one of the 23 structures, **44** of those carry a free
energy, and **34** actually print the structure's own run-1 value -- 34 is the number a reporting
rule would rewrite. (`orig_f12` and `ox2_shuffle` account for 64 of the 100 mentions but only 13 of
the 34, because they appear across many tables that hold no dG.) The sections below
are the original plan, kept for the recipe.

**Status 2026-10-10: paused until the racc queue reopens Monday 2026-10-12 — nothing submitted.**
One stale figure when it resumes: the
inventory was read 2026-10-08, and the fourteen legs since then (ten `*_o3cof` cofold legs + four
`_o3cof` nulls, all with MM/GBSA rows) are NOT staged as `_r2` and not in the four batch lists —
extend the inventory and stage them before the smoke leg if they are to be replicated in the same
pass.

CLAUDE.md records that "MM/GBSA standard errors are not uncertainty on a comparison" but nothing has
ever put a number on run-to-run variation. This plan runs a second, independent 20 ns trajectory of
every leg that produced an MM/GBSA row — both molecules, 62 legs — so the headline comparisons can
carry a real ΔΔG reproducibility figure. **Plan only; no leg is launched until you approve it.**

**Why a re-run is a replicate without touching any code.** `openmm_md/dynamics.py` never sets a
seed: the `LangevinMiddleIntegrator` (line 77) and `setVelocitiesToTemperature` at equilibration
(line 183) both draw fresh random numbers every instantiation. So running `omd run` again on the
*same* `system.xml` is an independent trajectory, and the replicate definition is exact: identical
starting coordinates and box, new initial velocities and new Langevin noise. No seed plumbing, no omd
change.

**The inventory (read off `system/system.xml` particle counts, not estimated): 62 legs.** Octinoxate
31 (27 in `runs/octinoxate/md/` + the four `bg33_*` co-folds in `~/python_mac/boltzgen_local/md/`),
oxybenzone 31 (27 + the four `bgox31_*` co-folds). Exclusions: `zzz_test_dock1` (leftover test build,
no system) and `orig_f12_dock5/prod_25ns` (the matrix leg is its 20 ns run — the 25 ns extension is a
different experiment, not a duplicate). Particle range 3,382–33,057; racc wall-time estimate
(using wall ≈ 8 min + 10M × 0.236 ms × (N/5838)^0.4): **~31.8 h octinoxate, ~26.0 h oxybenzone,
~58 h total ≈ 2.4 days continuous racc**. Modal by the cost model would be ≈ $55
($0.45/leg + $0.041/1,000 particles on 665,689 total) — an estimate, not a quote, and racc is free
and was where every recent leg ran anyway.

**Replicate construction: each leg gets its own directory `<leg>_r2` holding a byte-identical copy
of the original's inputs and a fresh dynamics run.** Copy from the original leg dir:
`system/`, `protein_fixed.pdb`, `ligand_prepped.sdf`, and `<N>_ligand.sdf` — nothing else (no
trajectory). Because `run_dock_pose_md.sh`'s build stage is guarded on `$M/system/system.xml`, a
staged `_r2` dir skips the build and goes straight to `omd run` into a new `prod_20ns`. Naming it as
its own leg dir (rather than `prod_20ns_r2` inside the original) is what makes every downstream glob
work unchanged — `md_stability.py`'s `runs/*/md/*/prod_20ns`, `pair_contacts.py`, `ligand_slide.py`
all pick the replicates up automatically.

**Two parameterisations are needed, both one-line (per "extend, don't reimplement"):**
1. `code/run_dock_pose_md.sh` gains `N=${LEGNAME:-${STRUCT}_dock${POSE}}` — then a replicate co-fold
   is `LEGNAME=<leg>_r2 STRUCT=<orig base> SRC=<orig dir>` and a replicate dock is
   `LEGNAME=<leg>_dock1_r2 STRUCT=<base> POSE=1 SRC=<co-fold dir>`.
2. `code/racc_drain.sh`'s dispatch parses `<leg>_dock<K>` to set STRUCT/POSE; parsing a
   `..._dock1_r2` name needs the `_r2` suffix stripped first, then passed back through. Small edit,
   flag it when the first batch is staged.

**Queue plan — racc, four batch jobs (each well inside the 24 h partition cap, one running at a
time):** octinoxate-A (~16 h) and octinoxate-B (~16 h) split the 31 legs so each batch carries a mix
of 5,000–8,000-particle legs (~46–54 min) and the few 21–33k ones (~74–87 min); same for oxybenzone-A
and -B (~13 h each). Stage with `racc_feed.sh <legfile>` in `BATCH=1` mode
(`<local dir>:<leg>` lines); fetch with `racc_drain.sh`, which also strips `traj.nc`. Preemption
restarts a leg from zero, but `racc_run.py`'s per-leg `traj_wrapped.xtc` guard means a resubmitted
batch skips completed legs and re-runs only what wasn't done. Total wall ~2.4 days if nothing
preempts; realistic figure 3–4 days.

**Order of operations.** (a) Make the two one-line driver edits. (b) Stage all 62 `_r2` dirs locally
(a find-loop of `cp` — the only new file this task needs). (c) **Smoke leg first**: run one
small leg end-to-end through run → fetch → tail (window MM/GBSA, contacts, frames, pairs) and verify
every stage before submitting any batch — the one-leg-checking rule, here applied to a whole new leg
kind. Suggested smoke: `ox2_orig_f8_r2` (5,231 particles, ~46 min racc, and its original is the best
characterised leg on the oxybenzone side, so agreement vs disagreement is interpretable at a glance).
(d) Submit batches in order, keeping the queue topped up. (e) Tail each leg as it lands —
never let the drain backlog (see [[drain-remote-completions-proactively]]).

**Per-leg recording set for a replicate** (same-day, one commit per leg — per the outstanding-work
list): whole-run + three window MM/GBSA rows; `md_contacts.py` and `md_frames.py` (pipeline steps);
`pair_contacts.py`; input→production ligand RMSD (`ligand_slide.py`); `md_stability.py` **full glob**
(the globs above pick `_r2` legs up with no change; include them in the same argv pass, never a
single-leg run). One change from a normal leg: **replicate rows go in a separate file
`runs/<mol>/md/mmgbsa_replicates.csv`, not into `mmgbsa_summary.csv`**, and `make_gnina_bundle.py`'s
`REFERENCE` dict is NOT extended — REFERENCE carries canonical-leg numbers for the GNINA bundle, and
a replicate is not a canonical leg.

**Deliverable.** Per molecule, `RESULTS.md` gains a replicate table after the MM/GBSA matrix:
`leg | ΔG run 1 | ΔG run 2 | |ΔΔG| | retention r2 (r1) |` — 62 rows, plus a median/mean |ΔΔG| per
molecule quoted in prose, which finally replaces "the standard errors are not uncertainty" with an
actual reproducibility number. If |ΔΔG| turns out small (a few tenths), the headline tables stand as
they are; if it is comparable to the design-vs-null gaps, the matrix gains a run-to-run spread column
across all rows at once. Ligand retention may differ more than ΔG does (residence is a tail event) —
record it, but treat a released-frames difference as a flagged cell, not a conclusion.

Nothing here runs until you say go; the first action on a go is the two one-line edits and the smoke
leg.

### Targeted replicate strategy (2026-10-19): answer the question with 11 legs, not 202 — DONE, see `REPLICATES.md`

**Do not run the full sweep first.** The plan above was written for a 62-leg inventory; there are
now **202** legs with a canonical whole-run MM/GBSA row and only **63** staged `_r2` directories,
because the unforced campaign, the 66-leg docked campaign, both cofold lanes and the guided lane all
landed after the 2026-10-08 staging. A full sweep is therefore ~170 h of racc across ~30 batches
(and the four staged lists `{oct,oxb}_{A,B}`, at 12.7–16.3 h each, are now illegal under the ≤6 h
rule and must be re-split whatever the scope).

**The question is not "what is the spread of every cell".** It is: **is the 16 kcal/mol spread on
`ox1_esm1_control` typical, or peculiar to the campaign's largest outlier?** Those two worlds imply
opposite actions — one invalidates every single-run cell in the headline tables, the other says the
tables stand and the outliers need an error bar — and they are distinguishable with ~11 legs.

**Stratify on two axes, because the second is the one that bites.** ΔG band (strong ≤ −18, mid −18
to −8, weak > −8) and *diagnostic quality*: `clean` = released <1% AND |20 ns − 10 ns window drift|
<1.5 AND slide <2 Å; `flagged` = any of released ≥10%, drift ≥3, slide ≥3 Å. The population is
heavily skewed toward the mid band, and the interesting cells are small:

| | clean | mixed | flagged |
|---|---|---|---|
| strong (≤ −18) | **6** | 8 | 6 |
| mid | 43 | 42 | 41 |
| weak (> −8) | 6 | 14 | 22 |

**`strong`+`clean` is the decisive cell, and it holds only six legs in the whole project** — one of
which is `ox1_esm1_control` itself. That is what makes this cheap: four of six is a near-census, not
a sample. And it is the cell that matters because run 1 *was* clean by every marker — flat monotone
windows, zero releases, 1.92 Å slide — and still failed to reproduce. So:

* **if the other strong+clean legs also fail to reproduce, no cell in the headline tables is
  trustworthy on one run**, and the full 202-leg sweep is justified;
* **if they reproduce, `ox1_esm1_control` is a genuine outlier**, the tables stand, and the spread
  column is needed only for the extremes — which is a presentation fix, not a week of queue.

Either way the answer arrives in ~10 h instead of ~170, and the full sweep stays available.

**The eleven legs** (`ox1_esm1_control_r2` is already run and is the twelfth; picked widest-|ΔG|
first within each cell, no two legs sharing a parent structure so a single bad build cannot carry a
whole stratum):

| cell | leg | ΔG | drift | released | slide | `_r2` staged |
|---|---|---|---|---|---|---|
| strong·clean | `ox2_orig_f8_dock4` (oxy) | −24.00 | −0.66 | 0.0% | 1.7 | yes |
| strong·clean | `bgox31_3_dock1` (oxy) | −22.23 | −0.89 | 0.0% | 1.6 | yes |
| strong·clean | `orig_f12_dock1` (oct) | −22.11 | +1.14 | 0.0% | 0.9 | yes |
| mid·clean | `s2_esm2_control_dock6` (oct) | −17.09 | −0.31 | 0.0% | 1.1 | yes |
| mid·clean | `ox2_esm1_f8_o3cof` (oxy) | −16.68 | +0.86 | 0.0% | 1.2 | no |
| mid·clean | `bg33_4_o3cof_dock1` (oct) | −16.04 | −0.39 | 0.0% | 1.2 | no |
| weak·clean | `bgox31_4_o3cof_dock7` (oxy) | −6.71 | −0.91 | 0.0% | 0.8 | no |
| weak·clean | `ox1_orig_f12_rf3cof_dock1` (oxy) | −6.40 | +0.90 | 0.0% | 1.4 | no |
| strong·flagged | `s3_esm2_f4_dock1` (oct) | −22.42 | +0.25 | 0.0% | **10.9** | yes |
| mid·flagged | `ox1_esm1_f8_o3cof` (oxy) | −17.95 | +1.28 | 0.0% | **3.0** | no |
| weak·flagged | `shuffle_control_o3cof_dock8` (oct) | −7.78 | +1.81 | 0.0% | **3.2** | no |

Particle counts 4,069–23,142; estimate **607 min = 10.1 h**, so **two batches** under the 6 h cap
(the 23,142-particle `s2_esm2_control_dock6` at 103 min and `s3_esm2_f4_dock1` at 77 are the only
large ones). Five are already staged from 2026-10-08; the other six need their `_r2` directory
copied per the construction recipe above, which is `system/`, `protein_fixed.pdb`,
`ligand_prepped.sdf` and the leg's `*_ligand.sdf` — and nothing else.

**What to report.** Per leg, ΔΔG **with its sign** and the change in released share, since retention
moved as violently as ΔG did on the first replicate (0% → 9.3%). The deliverable is the **tally of
stronger / weaker / unchanged** across the set, plus median |ΔΔG| per cell — a measure of how much
the random seed moves the dynamics, which is then the basis for deciding between more targeted
replicates and a full campaign. Take **unchanged = |ΔΔG| < 1.5**, reusing the drift threshold the
`clean` definition already uses, so the two are on one scale.

**Two framing errors to avoid, both made and corrected on 2026-10-20.**

*A replicate that lands elsewhere is not a "failure to reproduce".* Both trajectories are valid
samples of the same system under different seeds; neither is the true value and the first is not the
reference. The pair measures seed sensitivity, and every run is informative — the language should be
"moved by N kcal/mol", not "failed".

*The convergence windows are not a reproducibility predictor and must not be scored as one.* They
measure whether **one** trajectory's running mean has settled over its own 20 ns. A leg can be
perfectly converged internally and still sit far from an independent run of the same system; those
are different quantities, and faulting the windows for not anticipating the second is a category
error. (An earlier version of this paragraph asked for exactly that test and called the n=1 result a
negative; that was wrong.)

**Replicate rows go to `mmgbsa_replicates.csv`, never the canonical summary, and the REFERENCE dict
is not extended** — a replicate is not a canonical leg.

**One leg of this campaign has now run, out of order, and it landed hard (2026-10-19).**
`ox1_esm1_control_r2` was run on its own as the highest-value single leg rather than as part of a
batch, because the leg it replicates carried a headline claim. **|ΔΔG| = 16.16 kcal/mol**
(−26.09 → −9.93), with retention collapsing too: 0% released → 9.3%, 54.3 mean contacts → 23.4.
So the answer to this plan's own question — "if |ΔΔG| turns out small (a few tenths), the headline
tables stand as they are; if it is comparable to the design-vs-null gaps, the matrix gains a
run-to-run spread column" — is, on n=1 of 62, emphatically the second. One replicate does not
establish the distribution, and this leg may be the worst case rather than typical: it is the
campaign's largest outlier, and an outlier is exactly where regression to the mean is strongest.
But it rules out the benign reading, and a 16 kcal/mol spread on any cell is larger than most
differences the headline tables are asked to carry. Recorded in `mmgbsa_replicates.csv`; the other
61 staged `_r2` legs are untouched.

**Prepared 2026-10-08 (uncommitted; nothing launched, nothing staged to racc).** Both one-line edits
are in (`run_dock_pose_md.sh` now honours `LEGNAME`; `racc_drain.sh` parses `..._dock1_r2` by
stripping the suffix and passing `LEGNAME=$leg`). All **62 `_r2` dirs are staged local** under
`runs/<mol>/md/` — including the BoltzGen co-fold replicates, which live in `runs/` next to their
molecule's docks rather than in `boltzgen_local/md/` — each holding a verified byte-identical
`system/` (`system.xml` particle count equal to the original's, checked all 62) plus
`protein_fixed.pdb`, `ligand_prepped.sdf`, the leg's `*_ligand.sdf`/`*_protein.pdb`, and no
trajectory. The four wall-balanced leg lists are written: `runs/replicates/{oct_A,oct_B,oxb_A,oxb_B}.legs`
(`oct_A` 15 legs ~15.6 h, `oct_B` 16 legs ~16.3 h, `oxb_A` 16 legs ~13.3 h, `oxb_B` 15 legs ~12.7 h,
in dir:leg format for `racc_feed.sh BATCH=1`). The remaining actions on a go are: one smoke leg
(`ox2_orig_f8_r2`), then the four `racc_feed.sh BATCH=1` submissions with `BATCH_NAME` per batch and
`SYSNAME` set to the batch's molecule, a `racc_drain.sh` watcher per molecule, and the same-day tails.


## Co-fold determinism: are the co-folds reproducible at all? — DONE 2026-10-22, see `FOLD_REPLICATES.md`

**Answered.** 29 structures re-folded twice each, unseeded, 87 pairwise comparisons; commit
e8743bfb. The result and the tables live in `FOLD_REPLICATES.md`. Headline: glycine content
predicts backbone reproducibility (below 25% glycine ~0.9 A mean RMSD over peptide CAs, at 60%+
~6.6 A, and the designed orig/shuffle arms are 60-65% glycine by construction), and backbone
reproducibility is independent of where the ligand lands. Nothing from it is folded into any
other table. What remains open below is only the seeding policy, which the user settled: do not
seed.

**STILL TO DO, after 2026-11-01: the same replicates for OpenFold3 and RosettaFold3** (user,
2026-10-22 — gated on the Modal credit renewal, since unlike Boltz-2 these two do not run on
this Mac and bill per fold). What was measured above is Boltz-2's seed-to-seed spread only; the
cofold tables compare three co-folders, so the other two need the same treatment before any
cross-model difference can be read against each model's own noise. Do not launch before the
credits renew, and size it first: `code/fold_batch.sh LANE=openfold3` and `LANE=rf3` are the
lanes, the same structures as the Boltz set, and the per-fold cost is the thing to check against
the renewed balance rather than assumed. Whether OF3 and RF3 are even unseeded by default is
unverified — check each lane's invocation the way `boltz_env.run_boltz` was checked, because a
model that seeds by default would give a spread of zero and that would be a finding, not a bug.

**Established fact, checked 2026-10-22: no co-fold in this project has ever been seeded.**
`boltz_env.run_boltz()` builds `predict <yaml> --num_workers 0 --out_dir <dir>` plus any
`extra_args`, and nothing in the repo passes `--seed`. Boltz's own help states the default
explicitly: *"Seed to use for random number generator. Default is None (no seeding)."* Sampling
defaults are `recycling_steps 3`, `sampling_steps 200`, `diffusion_samples 1` — one diffusion sample
per fold, from an unseeded stream.

So **a second co-fold of the same peptide and ligand is an independent draw, not a reproduction**,
and we do not know how far it lands from the first, because no structure has ever been folded twice.
119 co-fold yamls exist under `runs/{octinoxate,oxybenzone}/boltz/` and each was drawn once.

**Why this matters more than it looks.** The 23-structure replicate campaign (`REPLICATES.md`)
re-ran *dynamics from a fixed starting structure*: same `system.xml`, different velocity seed. Fold
variation sits one level upstream of that — if a second co-fold puts the ligand somewhere else, the
starting structure is itself a draw and every MD leg inherits it. The two are different sources of
variation and the replicate spread does not bound this one.

**SET UP 2026-10-22, NOT RUN** (user: "set up but do not yet run folding replicates"). Two pieces,
both committed and dry-run verified; no fold has been drawn.

```sh
python3 code/fold_replicates.py <mol> --n 3 --list                  # the plan, nothing folded
python3 code/fold_replicates.py <mol> --n 3 --only <name>           # draw them
.venv/bin/python code/compare_folds.py runs/<mol> --replicates      # the RMS comparison
```

`code/fold_replicates.py` re-folds from the **same yaml the published fold was made from**
(`runs/<mol>/boltz/<name>.yaml`), so the only difference between a draw and the original is the
RNG. Output goes to `runs/<mol>/folds/foldrep/<name>/rep<k>/`, which touches no existing path,
table or leg — the published co-fold stays where it is and remains the published fold. Resumable:
a draw whose cif exists is skipped. It passes no `--seed` and must never be given one.

`compare_folds.py --replicates` compares each draw against its own original and writes
`runs/<mol>/folds/fold_replicate_rms.csv`. **Three numbers per draw, because one "RMSD" would hide
the interesting case:**

| column | superposed on | question |
|---|---|---|
| `ca_on_peptide` | the two peptides' CAs | is it even the same fold? |
| `ca_on_ligand` | the ligand | how differently is the peptide arranged around it? |
| `lig_on_peptide` | the peptide's CAs | did the ligand land in the same pocket of its own peptide? |

A small `ca_on_peptide` with a large `ca_on_ligand` is the same fold presenting a different face to
the ligand, which matters more here than a different shape would. `lig_on_peptide` is built exactly
like `ligand_slide.py`'s measure, so it reads on the same scale as a leg's slide.

**Scale, for scoping.** 81 co-fold yamls for octinoxate and 38 for oxybenzone, so n=3 over
everything is 243 + 114 = 357 draws. Boltz-2 runs on this Mac's MPS with no Modal spend; the apo
lane measured ~45 s/peptide and BoltzGen ~100 s/design, so a co-fold draw is probably 1–2 min —
**not measured, and co-folds carry a ligand, so measure one before trusting any total.** The
sensible shape is several draws of a few chosen structures (`--only`), not n=3 across 119: the
question is within-arm spread, which one peptide answers, and the per-structure choice is the
user's. `--n`, `--only` and the ordering are all still open; **do not launch it unasked**.

**DECIDED by the user 2026-10-22: do not seed.** *"I don't want to seed because I'd like to see
what the spread is from doing different replicates."* So `run_boltz` keeps passing no `--seed`,
every fold stays an independent draw, and the spread is something to *measure* rather than
something to eliminate. Do not add a seed flag, do not propose one, and do not re-fold anything
with a fixed seed "for reproducibility" — that would destroy the quantity being measured.

**Two distinct questions, and the second does not undermine the first** (user, 2026-10-22):

1. **Forced vs unforced folding** — the existing comparison, and it remains valid as it stands.
   A fold either had its contacts forced or it did not, and what we examine is the fold we
   actually got. That is a real object with real dynamics, whatever other draws might have
   existed; the cofold tables compare them on identical measures.
2. **How different the same fold would be under a different seed** — the new, separate question.
   It asks for the *within-arm* spread: re-fold the same peptide + ligand and see how far the
   second draw lands from the first.

Keeping these apart matters for how results get written. Answer 2 does not retract answer 1; it
tells you the resolution at which single-draw comparisons between arms can be read — the same
role the replicate spread plays for the MM/GBSA tables, one level upstream. Report it that way
rather than as a correction to the forced/unforced finding.

## Formal data tables for the paper (added 2026-10-22, awaiting the user's specification)

The user will specify what to build. **Do not start, and do not guess at a table list** — the
existing tables were restructured on instruction more than once (FOLDS.md's twenty campaign tables,
five per molecule per fold lane; never pool two molecules in one table; never add a fold lane as
extra columns), so inventing a paper layout would mean redoing it. What exists to draw on:
`RESULTS.md`'s ten headline tables, `FOLDS.md`'s twenty campaign tables,
`runs/folds_md/cofold_binding_tables.md` (Tables 1 and 2, spliced by `cofold_tables.py --write`,
rendered to PDF by `code/render_cofold_pdf.py`), `REPLICATES.md` and `TRAILING_WINDOWS.md`.
Note that the replicate reporting question is unresolved, and it lands squarely in this work.

## Oxybenzone ABFE set: complementary campaign, proposed 2026-10-14, not launched

**Plan only; no leg is submitted until you approve it.** The octinoxate ABFE cross-check campaign
(fep repo `RESULTS.md`) finished 9/9 points on 2026-10-14 — OpenFE 1.12 `AbsoluteBindingProtocol`,
30×10 ns complex + 14×10 ns solvent, ~10–14 h/leg on the racc H100s (dock7 took 11.6 h) — and its
working rule ("MM/GBSA dynamics screens, ABFE prices the survivors") is written. Oxybenzone already
has one ABFE point: `ox2_orig_f8` co-folded → −2.60 ± 0.61. The set below adds the cells that make
oxybenzone ligand-matched with it. Oxybenzone is achiral, so no S/R prep confound — solvent-leg
estimates are directly comparable across legs.

Core tier (one-to-one parallels of the cells that discriminated on octinoxate):

| # | leg (20 ns peptidebuilder record) | cell type | MM/GBSA / residence | what it tests |
|---|---|---|---|---|
| 1 | `ox2` GNINA p4 (aromatic, docked) | docked & **retained** −24.00/100% | −24.00 | The one ligand-matched cell where docking did not fail. On octinoxate every docked leg priced unbound (−0.60…+0.91); here the docked pose *retains*. Does it price like its co-fold sibling (−2.60), or is octi's 2–6 kcal prep-route gap about the pose, not the route? |
| 2 | `ox2_esm1_f8` co-folded | aromatic·esm2 linker −20.81/100% | −20.81 | Linker invariance, the parallel of `s3_esm2_f4` (ABFE −6.53 vs fold-parent −5.47). |
| 3 | `ox2_shuffle_esm0` co-folded | null rod that **retains** −12.04/94.7% | −12.04 | The sharpest single test: the octi null rod priced −9.68 (exterior groove). Does a *retaining* null price like the designed set or like the octi null? "ABFE ranks the site, not the design" hangs on this. |
| 4 | `ox2_shuffle` co-folded | **released** null −5.30/**7.1%** | −5.30 | The octi adsorption-control surrogate (`s2_esm2_control` → −4.52): released by MD, so ABFE prices the surface it lands on. |

Second tier (fills the remaining chemistry axes):

| # | leg | cell type | what it tests |
|---|---|---|---|
| 5 | `ox1_orig_f12` co-folded | charged shell, weak −10.31/**49.1%** | dose-response: octi showed −8.40 MM/GBSA → ABFE ≈ unbound (−1.25). Does a weak-designed cell fall the same way? |
| 6 | `bgox31_3_dock1` | BoltzGen best, docked & retained −22.23/100% | MM/GBSA says the shell beats BoltzGen by only 2.0 kcal on oxybenzone (octi BoltzGen never enclosed at all). ABFE tests whether that 2-kcal ranking is real. |
| 7 | `ox3_orig_f8` co-folded | hydrophilic shell −18.84/100% | third shell chemistry; completes the per-shell ABFE-offset spread at fixed ligand. |

Practical notes (all learned the hard way):

* Prep: all-atom ligand is only in the MD `complex.pdb` (peptidebuilder `runs/oxybenzone/md/<leg>/`);
  rebuild the SDF by the graft-Hs-to-prepped-SDF path (fep `scratch/make_allatom_lig_b2.py` lineage),
  heavy coords verified 0.000 Å.
* Boresch: apply wide settings on principle to any null/shuffle rod
  (`--rmsf-cutoff 0.3 --anchor-strategy multi-residue --host-max-distance 2.0`) and validate the
  geometry stage *before* submitting (host-pool + all-anchor sweep, fep `scratch/dock7_diagnose2.py`
  lineage) — dock7 burned a job on the defaults even with banked equils.
* Check the release column of any BoltzGen pick before submitting (5 of the 12 `bgox31` legs are
  ligand releases; their ΔG is computed partly on an unbound ligand).
* Core tier ≈ 2 days of scavenger GPU at the measured pace; submit top-up ≤3 via the
  drain/feeder pair (fep `code/racc_abfe_feed.sh` + `code/abfe_feed_legs.txt` pattern), and never
  without a validated wide-Boresch tf if the peptide is a rod.


## Forced vs unforced co-folds: making the cofold tables like-for-like — Part 1 DONE 2026-10-19, Part 2 CLOSED 2026-10-21

**Part 2 is closed, not open.** The pocket-guided OF3 lane ran to 26 legs and every cell of cofold
Tables 1 and 2 is populated (`cofold_tables.py` reports `pending: 0`); the finding was that the
pocket constraint damages ligand PLACEMENT rather than the backbone. RF3 was closed on a source
read — its ground-truth conformer is re-centred and randomly rigid-augmented, and inter-molecule
template distances are masked and never exposed — so it has no usable guided mode and no
comparable column was possible. The three "open questions" below are answered; they are kept
because the reasoning is the record of why RF3 was abandoned.

> **2026-10-19, and it changes the campaign's headline: `ox1_esm1_control` does NOT reproduce.**
> The replicate (`ox1_esm1_control_r2`, racc 529604, recorded and audited) came back **−9.93
> against the original's −26.09 — 16.16 kcal/mol apart**, from a byte-identical system differing
> only in the initial velocities `openmm_md` never seeds. Released 9.3% of frames in 11 episodes
> where run 1 released none; mean contacts 23.4 against 54.3; separation 9.3 Å mean and 19.0 Å max
> against 4.3/5.4. Windows −13.48 / −9.69 / −8.53 / −9.93 against −24.28 / −25.08 / −25.83 / −26.09.
> Component-wise the engagement simply failed to form: EEL −41.14 → −9.52 (partly offset by EGB
> 42.61 → 16.19), vdW −23.79 → −14.46.
>
> **Three things follow.** (1) Point 3 below — "the campaign's strongest leg is an unforced
> co-fold, which cuts against the premise the campaign was built on" — **is now unsupported**, and
> the sentence that called a 9.5 kcal/mol outlier on n=1 too big a claim was right. (2) Run-to-run
> variation on this system is **~16 kcal/mol** against within-trajectory standard errors of
> 0.02–0.04: CLAUDE.md's standing warning that those errors are not uncertainty on a comparison now
> has a number, and it is larger than most design-vs-null gaps in the project. (3) **Convergence
> did not catch it** — run 1's windows were flat and monotone, textbook converged. Convergence
> within a trajectory says nothing about reproducibility across them. Replicate rows live in
> `runs/<mol>/md/mmgbsa_replicates.csv` (never the canonical summary); `record_leg.py` routes any
> `_r2` leg there. This is the strongest argument yet for the 62-leg replicate campaign already
> staged below, and it is not started.
>
> Part 2 prep is done and **submitted**: nine pocket-guided OF3 folds from the designed shells,
> docked, GNINA bundles pushed, and 18 MD legs (9 cofold + 9 Vina p1) on racc as four batches
> `runs/folds_md/o3pc_{oct_1,oct_2,oxb_1,oxb_2}.legs` via `code/feed_chain.sh`. RF3 is closed.
> See "Status 2026-10-19" and "Answers, 2026-10-19" under Part 2.
>
> **Start here after a context reset.** Part 1 is finished, audited, committed and pushed (through
> `c826ff9`); nothing from it is pending. Part 2 — forced OF3 and RF3 folds — has not started, and
> its plumbing is already written and pushed in the sibling `~/python_mac/fold` repo (`fe3abed`:
> `POCKET:`/`POCKET_MAX_DISTANCE:` input fields, `build_pocket_constraint`, a `--pocket` flag) plus
> `code/pocket_inputs.py` here, which writes pocket-guided inputs from a forced fold's own yaml.
> Nothing has been *run* with any of it, so every claim about OF3's pocket behaviour below is
> code-reading, not measurement.
>
> **Three things to carry forward, each of which cost something to learn:**
>
> 1. **Read release, not residence.** Residence is a centre-of-mass measure and understates
>    retention badly on these elongated folds — `s3_esm2_control_dock7` reads 0.1% residence having
>    released *zero* frames of 20 ns. Several legs were compared on residence before this surfaced.
> 2. **Don't generalise from one triple.** GNINA's pick beat Vina p1 in 4 of 7 and lost 2. It was
>    called a clean win for rescoring after a single triple, and the next triple contradicted it.
> 3. **The campaign's strongest leg is an unforced co-fold** (`ox1_esm1_control`, −26.09, zero
>    releases). That cuts against the premise the campaign was built on, so **a replicate of that
>    one leg is the highest-value next run** — a 9.5 kcal/mol outlier carrying a headline claim
>    should not rest on n=1. It is also the leg that NaN'd once and ran clean on a plain
>    resubmission, so treat a single failure as a draw, not a verdict.

**Part 1 ran and is complete: 23 legs, all recorded, Tables 1 and 2 rebuilt.** The plan as written
below is kept for the reasoning; where it was superseded in execution, the differences are recorded
in "What actually ran" at the end of Part 1 rather than edited into the steps, except for three
places where leaving the original text would mislead (the rate exponent, the batch sizing and the
step-6 table plan, all corrected in place). `runs/folds_md/cofold_binding_tables.md`
sets the Boltz co-fold of each peptide against its OpenFold3 and RF3 co-folds, but most Boltz cells are
**hint-forced** folds (`_f4`/`_f8`/`_f12`: that many `contact` constraints in the yaml) while the OF3
and RF3 lanes were run with **no constraints at all** (`fold_batch.sh` passes sequence + SMILES only).
A column difference therefore mixes "which model" with "forced or not". Two fixes:

### Part 1 (now): MD/MM-GBSA on the Boltz *unforced* folds that already exist

Every forced fold in the tables has an unconstrained Boltz twin of the identical sequence already on
disk (`runs/<mol>/boltz/boltz_results_<x>_control`, 0 `contact` lines, checked 2026-10-18). Rows that
are already unforced need nothing new.

| table row | Boltz cell now | unforced twin (Boltz dir) | new leg |
|---|---|---|---|
| `orig_f12` (oct) | forced, 12 | `RGGDGGKGGGGLGGGGKGIGEGWGGDGSGGEGS_control` | `orig_control` |
| `s3_esm2_f4` (oct) | forced, 4 | `s3_esm2_control` | `s3_esm2_control` |
| `s3_orig_f12` (oct) | forced, 12 | `s3_orig_control` | `s3_orig_control` |
| `ox1_esm1_f8` | forced, 8 | `ox1_esm1_control` | `ox1_esm1_control` |
| `ox1_orig_f12` | forced, 9 | `ox1_orig_control` | `ox1_orig_control` |
| `ox2_esm1_f8` | forced, 8 | `ox2_esm1_control` | `ox2_esm1_control` |
| `ox2_orig_f8` | forced, 8 | `ox2_orig_control` | `ox2_orig_control` |
| `ox3_orig_f8` | forced, 8 | `ox3_orig_control` | `ox3_orig_control` |
| `s2_esm2_control`, `shuffle_control`, `shuffle_control_esm0`, `ox2_shuffle`, `ox2_shuffle_esm0` | already unforced | — | none (reuse) |
| `bg33_*`, `bgox31_*` | BoltzGen design, no hints | — | none (reuse) |

**8 new co-fold legs** (3 octinoxate, 5 oxybenzone), all on racc. Steps, each an existing script:

1. **Stage** each leg: `code/cif_to_md.py runs/<mol> --structure <boltz name> --out-dir runs/<mol>/md/<leg>`
   (model_0, the same model the forced legs used), check the printed formula (C17H24O3 / C14H12O3),
   then `omd prep-protein --pdb <leg>_protein.pdb --out protein_fixed.pdb`.
   One small parameter to add: `cif_to_md.py --name` (default = `--structure`) so the shell-1
   octinoxate fold's files come out as `orig_control_*` instead of the 33-letter sequence name.
2. **Build** with the one driver: `STRUCT=<leg> POSE=cofold SYSNAME=<mol> SRC=runs/<mol>/md/<leg>
   code/run_dock_pose_md.sh BUILD_ONLY=1`. Its ligand guard finds `<leg>_ligand.sdf` from step 1, so it
   never reads the OF3 dock stage; dodecahedral box by default. Read the particle counts. Unforced
   folds tend to be rods, so some may be large (`s2_esm2_control` was 23,142); quote wall time from
   `8 min + 39.3 min × (N/5838)^e` with **e = 0.4 only up to ~13k particles and 0.64 above** —
   measured 2026-10-18 on `orig_control` at 19,617 particles, where the 0.4 model predicted 71 min
   against ~94 actual. These rods run 18k–23k, which is exactly where the old exponent under-buys
   the `--time` request.
3. **Submit**: legfiles under `runs/folds_md/bzfree_*.legs`, `racc_feed.sh BATCH=1` each.
   **A batch is capped at ~6 h of estimated wall, so this is five legfiles, not two** — a preempted
   leg on `gpuscavenger` restarts from zero, so a smaller batch bounds the loss. Every submission
   must pass `SBATCH_OPTS="--time=<estimate + 1h30> --cpus-per-task=8"`; `racc_feed.sh` exits 2
   without both. Free (gpuscavenger), one job runs at a time under `MaxJobsPU=1`.
4. **Drain**: `racc_drain.sh` per legfile. One change is needed: a bare `<leg>` name currently goes to
   `boltzgen_local/md/run_md20.sh`, which `cd`s into boltzgen_local and would not find a
   peptidebuilder leg. Add a branch: a bare leg whose local dir is under `$REPO` tails through
   `run_dock_pose_md.sh POSE=cofold`, same as `_o3cof`/`_rf3cof` legs.
5. **Record**: `record_leg.py <mol> <leg>` per leg. No pair table holds these structures yet, so seed
   one per shell with the first leg via `pair_contacts.py` (`--run-dir runs/<mol>`). Use
   `pair_contacts_bzfree_design.csv` (`design.json`) and `_shell3` (`design_shell3.json`) for octinoxate,
   and `_shell1/_shell2/_shell3` for oxybenzone. After that `record_leg.py` finds the table itself.
   Audit line `-> OK` per leg; strip check; commit additions-only.
6. **Tables**: ~~`cofold_binding_tables.md` gains a "Boltz unforced" column beside the existing Boltz
   one~~ — **superseded, and it was wrong.** Adding a column beside the forced one leaves the
   confound in place; Tables 1 and 2 exist to vary *one* thing, so the Boltz column is **replaced**
   by the unforced legs and the forced values come out (they survive in RESULTS.md and in Table 3's
   static forced-vs-unforced read). See the table convention immediately below, which is the
   authority here. The forced-minus-unforced difference is still the useful quantity; it lives in
   Table 3 rather than side-by-side in Tables 1 and 2.

Out of scope unless you ask: pepgraph/morph (these are Boltz folds under `runs/<mol>/boltz/`, not a
new `folds/` lane).

**Table convention, fixed 2026-10-18 before the data landed. Two different rules, and which applies
depends on what the table is for.**

* **Everywhere else — RESULTS.md's ten headline tables, FOLDS.md — the new dynamics are *added*.**
  Those tables are inventories of every leg taken through dynamics, so an unforced leg is another
  row beside the forced ones, labelled as unforced. Nothing existing is rewritten or removed.
* **`runs/folds_md/cofold_binding_tables.md` Tables 1 and 2 are *replaced*, because they are a
  different comparison.** Their whole job is one sequence across three folding models, so every
  column must be the same treatment. Adding unforced rows beside forced ones there would leave the
  very confound the campaign exists to remove: the published version set hint-forced Boltz against
  unconstrained OF3/RF3, which is not a model comparison. So the Boltz column *becomes* the
  `<row>_control` legs, forced values come out (they survive in RESULTS.md and in Table 3's
  forced-vs-unforced static read), and row labels become the sequence with the `_fN` suffix
  dropped. **Getting this wrong once already cost a cycle**: on 2026-10-18 the tables were
  restructured into pose sub-rows while the forced column was left in place, which changed the
  shape and fixed nothing.

The additive half, in detail. The eight unforced peptides extend
the **existing** per-molecule RESULTS.md tables rather than getting their own set: the forced and
unforced co-folds of one sequence belong in the same column so the pair reads down it. The
distinction goes in the label column — `pb (unforced)` for octinoxate, `·unforced` appended to the
shell label for oxybenzone — and the five co-folds that were **already** unforced but unlabelled
(`s2_esm2_control`, `shuffle_control`, `shuffle_control_esm0`, `ox2_shuffle`, `ox2_shuffle_esm0`)
have been relabelled to match, with a note under the headline table saying a forced and an unforced
row are not two runs of one protocol. The pre-dynamics table's four static columns were filled for
all 23 poses the same day, from `dock_poses.csv`, which needs no trajectory; ΔG, slots, sim and
pairs cells stay blank until `record_leg.py` has the leg. Per-leg CSVs
(`mmgbsa_summary.csv`, `protein_stability.csv`, `ligand_slide.csv`, `pair_contacts_*.csv`) are
written by `record_leg.py` as each leg lands and need no manual step; RESULTS.md's ten headline
tables are the curated layer and are updated per molecule once its legs are in, not per leg.

**What the static read already says, before any dynamics** (detail under both pre-dynamics tables):
octinoxate's unforced folds do not put the ligand on the peptide — `orig_control` at **20.9 Å**
centroid separation is the furthest cell in that table — so two of its three legs start unbound and
their ΔG will be partly a bound-versus-unbound difference rather than a pose-quality one. All five
oxybenzone unforced folds do start in contact (5.8–9.0 Å) but every one is looser than its forced
twin (e.g. `ox1_orig_control` 0.52 enclosed/7.2 Å against `ox1_orig_f12`'s 0.995/1.2 Å), which
makes oxybenzone the cleaner of the two forced-vs-unforced comparisons. The unforced co-folds also
land at or below their own shuffle nulls on enclosure, so enclosure alone does not separate design
from null without the hints.

### What actually ran (Part 1, completed 2026-10-19)

**Scope grew from 8 legs to 23** at the user's instruction on 2026-10-18: each of the 8 unforced
co-folds also got its **Vina p1** and **GNINA pick** pose as MD legs, so Tables 1 and 2 compare
three models × three starting poses rather than three models × one. Vina docking ran locally into a
scratch `--out` with additions-only appends; the GNINA bundle went out to `gnina_local` and its
results came back the same day. Submitted as five `bzfree_*.legs` batches (≤6 h each) plus a
one-leg retry; every leg recorded, audited `-> OK`, and committed (`419e9e6`, `1d08d7f`).

The eight unforced co-folds, ΔG (kcal/mol) with release share — not residence, see below:

| leg | ΔG | released | note |
|---|---|---|---|
| `ox1_esm1_control` | **−26.09** | **0.0%** | strongest leg in the campaign; 46.1 mean contacts, vdW-dominated |
| `s3_esm2_control` | −15.70 | 0.15% | flat windows, converged |
| `s3_orig_control` | −13.37 | 2.5% | |
| `ox2_esm1_control` | −12.70 | 9.6% | windows decay −22.17 → −12.70 |
| `ox3_orig_control` | −10.13 | 13.5% | windows decay −15.28 → −10.13 |
| `ox2_orig_control` | −8.51 | 5.2% | |
| `ox1_orig_control` | −8.34 | 13.4% | |
| `orig_control` | −7.24 | — | starts unbound at 20.9 Å; ΔG is partly bound-vs-unbound |

**Four findings, with their confidence.**

1. **The campaign's strongest binder is an unforced co-fold** (`ox1_esm1_control`, −26.09, never
   releases, 100% residence). This is the result that bears on the original question: the published
   Boltz advantage in these tables was not simply "forcing helps". One leg and a 9.5 kcal/mol
   outlier, so it wants a replicate before it carries weight.
2. **GNINA's pick beats Vina p1 in 4 of 7 triples, loses 2, ties 1.** A lean, not a result — and
   recorded here because a single triple was read as a clean win for rescoring partway through,
   which the next triple contradicted. Both losses are `*_orig` peptides and most wins are
   esm-substituted, but `s3_orig` breaks that split.
3. **No consistent winner between co-fold and docked starts.** `ox1_orig` and `ox3_orig` are better
   from the docked poses than from their own unforced co-fold; `ox2_esm1` and `s3_esm2` the reverse.
   So the unforced co-fold is not reliably the best starting geometry.
4. **Release is late-weighted on oxybenzone and near-absent on octinoxate shell 3.** Eight of eight
   oxybenzone shell-1/2/3 legs release mostly in the second half; the shell-3 octinoxate legs barely
   release at all. The exception is `ox1_esm1_control`, which never releases — on that molecule a
   deeply enclosed pose behaves categorically differently from a surface-bound one.

**Read release, not residence.** `md_contacts.py`'s `residence_pct` is a *centre-of-mass* measure
(ligand centroid within 10 Å of the peptide centroid, line 101) while `released_pct` counts frames
with zero atom–atom contacts under 4 Å. These unforced folds are rods, so a ligand gripped at one
end scores low on residence while never letting go: `s3_esm2_control_dock7` is 0.1% residence with
**0.0%** release — perfect retention — and `s3_esm2_control` 29.4% with 0.15%. Several legs were
contrasted on residence alone before this was noticed; `CLAUDE.md`'s reporting section now carries it.

**One table bug, found while filling the dock rows.** `cofold_tables.py` sourced residence by
scanning `FOLDS.md`/`RESULTS.md` for the first `N%` cell in any leg-keyed row — and RESULTS.md's
*static fold-geometry* tables are leg-keyed too, so `s3_esm2_control` was rendered at **100%** from
a fold-geometry row (`RESULTS.md:1159`) against its real 29.4%. Fixed by reading
`runs/folds_md/record_<leg>.log` first, which is what `record_leg.py` writes from the profile
itself; two cells changed, both in the Boltz column.

**Both of Part 1's loose ends are now closed (2026-10-19) — do not redo them.** The four missing
pair tables were seeded from the kept trajectories (`pair_contacts_bzfree_shell3.csv` for octinoxate
shell 3, `_shell1/_shell2/_shell3` for oxybenzone), so all 23 legs have pairs data; and RESULTS.md's
ten headline tables were filled — 77 rows *added* across the eight dynamics tables and the 23 blank
pairs cells filled in the two static ones. Audited afterwards: every leg has three windows plus the
whole run, contacts, both frame files, a figure, zero MM/GBSA scatter, four rows in
`mmgbsa_summary.csv` and one each in `protein_stability.csv`, `ligand_slide.csv` and a pairs table.
**Nothing is missing.** If a later note claims otherwise, check the tables before believing it —
a stale "still open" line once cost this project two weeks of false "missing" work.

**The table chain, which is three steps and not one.** `cofold_tables.py` does *not* write the
markdown: it prints diagnostics and writes table *blocks* to `$OUT` (default `/tmp/tables.md`),
which are then spliced into `runs/folds_md/cofold_binding_tables.md`, after which
`code/render_cofold_pdf.py` rebuilds the PDF from that markdown. Miss the third step and the PDF
goes stale silently, which is exactly what happened for four days in October. Audit those three
tables' matching whenever a leg is added:

```sh
OUT=/tmp/t.md $OMD/bin/python code/cofold_tables.py   # blocks + "missing/pending" report
# splice the blocks into runs/folds_md/cofold_binding_tables.md
.venv/bin/python code/render_cofold_pdf.py            # -> cofold_binding_tables.pdf
```

### Part 2 (ran 2026-10-19 → 2026-10-21, CLOSED): forced OF3 and RF3 folds

The research below is what the campaign was sized from; it was written before anything ran and
its code-reading conclusions held up. OF3's pocket-guided lane ran to 26 legs; RF3 was never
launched because the read showed it has no usable guided mode.

The fold repo (`~/python_mac/fold`) has no constraint plumbing today: no `constraint`/`contact`/`pocket`
handling in its OF3 or RF3 apps (checked 2026-10-18). **Upstream support was researched 2026-10-18 by
reading the two projects' released source and issues — nothing was installed or run, so every
behavioural claim below is code-reading, not measurement.** The answer differs per model, and neither
reproduces a forced Boltz fold.

**OF3 has a pocket hint, and the pinned image already carries it.** `pocket_constraint` arrived in
OpenFold3 0.4.5 and `fold`'s image is 0.5.0, so no upgrade is needed:

```json
"pocket_constraint": {"ligand_chain_id": "Z", "pocket_residues": [["A",4],["A",24]], "max_distance": 4.0}
```

One per query (the field is singular; the docs state no maximum), residue-level, `max_distance`
default 4.0 and documented only as "the distance threshold used when scoring whether ligand atoms
contact the specified pocket". Verified against upstream 2026-10-18 — `docs/source/input_format_reference.md`,
`core/config/pocket_sampling_config.py` and `core/model/structure/pocket_constraints.py` at v0.5.0.
It is **weaker than a Boltz `contact` in two ways that matter here**:

1. **It names no ligand atom**, so each hint's `token2` is discarded; `[chain_id, residue_id]` pairs are
   the only granularity the schema has.
2. **Placement is a ranking preference, not a constraint.** Enabling the field triggers "pocket
   proposal and partial-diffusion refinement": 1,024 rigid ligand proposals (32 RDKit conformers,
   random rotation plus a translation onto the pocket centroid or a random pocket atom) are scored and
   sorted — pocket-centroid distance, in-pocket atom count, a contact penalty
   `relu(lig_to_pocket - contact_distance)²`, vdW overlap — then greedily diversity-filtered at 0.5 Å
   RMSD and handed to refinement as seeds for the last 25% of the diffusion schedule
   (`noise_frac=0.75`, 16 parents). **Nothing rejects a seed that fails to contact the pocket**: if the
   diversity filter leaves too few poses the remaining slots are filled from the top of the sorted list
   regardless of contact quality. So there is no equivalent of `force: true`'s potential.

A third concern — that the peptide fold comes from an *unconstrained* parent and so never folds
around the ligand — is **plausible but unverified, and an earlier draft of this section overstated
it**. `pocket_constraints.py` only builds seeds: it copies every non-ligand atom from the parent
unchanged and overwrites ligand coordinates. Where the parents come from (`xl_base`) and whether the
peptide moves during the refinement that follows are both outside that file. Worth settling before
the campaign, because if it holds, a pocket-guided leg cannot show a shell closing around the ligand.

**RF3 has nothing to use.** Confirmed against the `production` README and issue thread 2026-10-18.
The README states the limitation itself — "The model was trained with random pairwise contacts; our
inference API could be further extended to support specification" — and its interface-templating
section ("We have yet to extend our inference API to support this use case") is commented out, so it
is not live documentation. A fold job's documented fields are `name`, `components`, `bonds` and
`template_selection`; there is no contact, restraint or hotspot input, and covalent `bonds` is the
only restraint-like one.

Issue #160 (open, created 2026-01-08) is the request, and its thread is the useful part. A maintainer
called it "just an inference API change", then reported that on their de novo binder set **a single
pairwise interface contact "significantly hindered the model's ability to act as a filter"** —
because monomer prediction is strong, one contact collapses the degrees of freedom — and that they
were exploring softer constraints the model can ignore. Inter-chain templating was promised "within
the next 2 weeks" on 2026-02-14 and the thread never confirms a release; two separate requests for
*ligand*–protein contacts (2026-02-06, 2026-04-28) went unanswered. So an upgrade would not help:
there is nothing to upgrade to.

The `template_selection` route is not a workaround — it would need a 3D pose templated from the
forced Boltz fold, which is circular, and it fixes intra-set distances too. One thing to check rather
than assume: the README advertises `ground_truth_conformer_selection` as "useful for keeping known
ligand conformations while allowing the model to predict protein structure around them", which is
closer to what we want than the agent-report characterisation ("fixes the ligand's own conformer, not
where it sits") allowed. Whether it also fixes *placement* is unverified and worth one look before
writing RF3 off entirely.

So the plan is **OF3 only, labelled honestly**:
* add an optional `pocket_constraint` to the query built in `openfold3_core.run_predict`, taking the
  residue list from each Boltz yaml's `contact` `token1` entries and `max_distance` from the largest of
  that yaml's values (~4.2 Å) — a new input-file field, or a parser over the yaml we already have;
* new `LANE=openfold3p` in `fold_batch.sh` and leg suffix `_o3pcof`, with its own branch in
  `run_dock_pose_md.sh`'s `LANE_DIR` and in `racc_drain.sh`;
* MD reuses Part 1's flow exactly; the table gains **one** new column, named "OF3 pocket-guided" —
  **not "forced"**, because it is not the same intervention.
* pin `rc-foundry` in `rf3_app.py` while in there regardless: it installs unpinned, so the version in
  Modal's cached image is unknowable.

**Status 2026-10-18: the OF3 side is written, not run.** `fold` commit `fe3abed` adds
`pocket_constraint` support — a `POCKET:`/`POCKET_MAX_DISTANCE:` pair in the plain-text input
format, a `--pocket "4,24,29,34"` flag on both CLIs, and `build_pocket_constraint` in
`openfold3_core`, with six tests. Inputs are generated, not typed:
`code/pocket_inputs.py <mol> <forced fold names...>` reads each `_fN` yaml's `contact` block and
writes `runs/<mol>/folds/inputs_pocket/<struct>.txt` with the `token1` positions and the largest of
that yaml's `max_distance` values (so no contact is tightened relative to Boltz). It refuses an
unforced fold, which has no hint set to carry. Verified on `s3_esm2_f4` (4 residues, 4.2 Å),
`orig_f12` and `s3_orig_f12` (12 residues, 5.2 Å). What is left is a `LANE=openfold3p` in
`fold_batch.sh` and the `_o3pcof` branches, then one smoke fold.

**Status 2026-10-19: plumbing done, smoke fold run, questions 1–3 answered (below). RF3 is out;
OF3 pocket-guided is a go.**
* **The guided set is 9, one per non-BoltzGen, non-null row of Tables 1 and 2: 4 octinoxate + 5
  oxybenzone.** Counted wrong twice on 2026-10-19 (reported as 8) by deriving the list from
  `pocket_inputs.py`'s refusals instead of from the tables' own rows. The script refuses a fold
  with no `contact` block, which is correct, but a refusal means "*this fold* is unforced", not
  "this peptide has no hints" — `s2_esm2_control` is the Table 1 row whose Boltz cell has always
  been the unforced control, while `s2_esm2_f4/_f8/_f12` sit on disk with full contact sets and no
  MD leg. **Count the rows of the table the column is going into; do not infer a row list from a
  tool's error messages.** That row has no canonical dose for the same reason (it was never
  forced), so it takes **f4**, matching `s3_esm2_f4` so both esm2 variants carry the same dose —
  a choice, recorded here, not an inherited value. Its input is named `s2_esm2_control.txt` after
  the manifest row rather than after the yaml it came from, so the lane lines up with
  `folds/openfold3/s2_esm2_control.pdb`; sequence and SMILES verified identical to that fold.

  **Anchor each pocket on the variant that has MD — that is the only thing the guided cell can be
  compared against** (user, 2026-10-19). Checked against `mmgbsa_summary.csv`, and 8 of the 9 already
  satisfy it: `orig_f12`, `s3_esm2_f4`, `s3_orig_f12`, `ox1_esm1_f8`, `ox1_orig_f12`, `ox2_esm1_f8`,
  `ox2_orig_f8`, `ox3_orig_f8` each have a forced MD leg *and* a `_control` leg, so the pocket comes
  from the dose that was actually simulated. **`s2_esm2` is the lone exception and must be read
  differently: it has `s2_esm2_control` MD only — no forced MD at any dose** — so its f4 pocket is
  anchored by nothing, and its guided cell has no forced cell for the same peptide to sit beside.
  Its reference columns and its guided-vs-unguided comparator do come from MD
  (`s2_esm2_control` and `s2_esm2_control_o3cof`), so the row is sound as a guided-vs-unguided
  read and unsound as a forced-vs-guided one. Say so in the caption rather than letting the cell
  imply a dose comparison it cannot support.
* **The nulls get no guided fold (user ruling, 2026-10-19).** All four (`shuffle_control`,
  `shuffle_control_esm0`, `ox2_shuffle`, `ox2_shuffle_esm0`) have zero `contact:` lines, so no
  pocket can be read off them. Two ways to manufacture one — borrowing a designed peptide's
  positions onto a null sequence, or shuffling the positions on a designed sequence — were both
  **declined**. So the guided OF3 column ships without a null, and nothing in it may be phrased as
  separating design from noise; it compares guided against unguided OF3 on the same peptide only.
* Pocket inputs written for all **9** rows of the Part 1 table
  (`runs/{octinoxate,oxybenzone}/folds/inputs_pocket/`): `orig_f12` 12 res/5.2 Å, `s3_esm2_f4`
  4/4.2, `s3_orig_f12` 12/5.2, `ox1_esm1_f8` 8/4.8, `ox1_orig_f12` 9/4.9, `ox2_esm1_f8` 8/4.7,
  `ox2_orig_f8` 8/4.7, `ox3_orig_f8` 8/4.6. SEQUENCE and SMILES lines are byte-identical to the
  unguided lane's `inputs/`, so the pocket is the only difference.
* `fold_batch.sh` has `LANE=openfold3p` (reads `inputs_pocket/`, job `<pep>_pocket`, output
  `runs/<mol>/folds/openfold3p/<pep>.pdb`), a working `MOL=` filter (documented before, never
  implemented) and `ONLY="<names>"`, and every Modal launch now holds `caffeinate -i -s`.
  `_o3pcof` is wired into `run_dock_pose_md.sh` (stage lane `openfold3p`), `dock_folds.py`'s
  `LANES` and `record_leg.py`'s FOLDS-row regex. **`racc_drain.sh` still needs `*_o3pcof` added to
  its `*_o3cof || *_rf3cof` branch** — not done because a drain was running it (never edit a
  running zsh script).
* Remaining campaign steps: fold the other 7 (`LANE=openfold3p code/fold_batch.sh`, ~90 s and
  cents each on the A10G), `folds_metrics.py openfold3p`, `dock_folds.py` (stages
  `<pep>_o3pcof_ligand_cofold.sdf` and Vina-docks), build + racc exactly as Part 1.

### The three open questions, and how to settle them cheaply — ALL ANSWERED 2026-10-19/21

**1 — answered: the peptide is not frozen.** The `s3_orig_f12` smoke fold through `LANE=openfold3p`
refolded the peptide and raised ligand engagement, so the hint can affect more than the ligand and
the column was worth running. It ran; see the paragraph above. **2 and 3 — moot: RF3 is closed**
on the source read described above, so neither its ground-truth conformer nor
`p_provide_inter_molecule_distances` was ever going to give a comparable guided column.

Each is a yes/no answerable in one or two runs, before any campaign is sized. **The first decides
whether a pocket-guided OF3 column is worth running at all.**

1. **Does the peptide refold during OF3's refinement, or is it the unconstrained parent's fold?**
   Run `s3_orig_f12` (12 hints, the tightest designed shell) through `LANE=openfold3p`, and run the
   same sequence unguided — we already have the unguided fold. Compare the two peptides' Cα RMSD
   and radius of gyration, and the ligand's wrapped/enclosed fractions, with the tools that already
   do this (`code/compare_folds.py`, `code/folds_metrics.py`). **If the peptide Cα RMSD is ~0 and
   only the ligand has moved, the hint cannot produce a shell** — report that as the finding and do
   not run 8 legs of MD on it. If the fold does change, proceed.
2. **Does RF3's `ground_truth_conformer_selection` fix the ligand's placement or only its internal
   conformer?** Cheapest decisive test: take one co-folded structure we already have, feed its
   ligand as the ground-truth conformer, and run RF3 twice with different seeds. Then measure the
   ligand centroid's position relative to the peptide (superposing on peptide Cα first, exactly as
   `code/ligand_slide.py` does). Placement fixed → centroid lands in the same place both times and
   matches the supplied pose; conformer only → internal geometry matches but the centroid wanders.
   If placement is fixed, RF3 has a usable guided mode after all and question 3 is moot.
3. **Does RF3's `p_provide_inter_molecule_distances` override reach inference?** Only worth asking
   if 2 says no. Set it to 1.0 through a Hydra override alongside a `template_selection` naming the
   hint residues and the ligand, and check whether inter-molecule distances appear in the templated
   features at all (the quick read is whether the predicted ligand moves toward the templated
   residues versus the same run without the override). Even a yes is **not** a route to a
   comparable column: it needs a pose templated from the forced Boltz fold, which is circular, and
   it fixes intra-set distances too. Treat it as a question about the tool, not a campaign arm.

**Answers, 2026-10-19.**

1. **The peptide is not a frozen parent copy, so proceed.** One smoke fold (`s3_orig_f12`, seed 42
   like the unguided run; `pocket_sampling.enabled: true` confirmed in the run's
   `experiment_config.json`). Peptide Cα RMSD between the two best-ranked models is **7.3 Å**, not
   ~0. The ligand reading improves on every measure (`folds_metrics.py`, unguided → guided, with the
   forced Boltz fold in brackets): enclosed 0.43 → **0.785** (0.965), wrapped 0.7 → 0.9 (1.0),
   engaged 14 → 18 (20), centroid separation 10.4 → 6.4 Å (3.6). Hint residues within 5.2 Å of the
   ligand go from 3/12 to 6/12. **What this does not show:** that the guidance *shaped* the
   peptide. OF3's own sample-to-sample spread on this G-rich 34-mer is as large as the
   guided-vs-unguided gap (5 samples per run, pairwise Cα RMSD: unguided-vs-unguided median 6.50 Å,
   guided-vs-guided 5.76, across runs 6.19), and the guided peptide is no closer to the forced Boltz
   fold (Cα RMSD 7.46 vs 6.68). So the peptide changes from run to run, and the guided ligand sits
   better on it, but there is no sign of a Boltz-like shell closing. One structure: read it as a go
   <!-- superseded by the all-8 result below; kept because it is the reasoning the campaign was
        sized on, and because it is the exact shape of error the "don't generalise" rule names -->
   for the campaign, not as a result.

   **All 8 folded, 2026-10-19 14:42 — the smoke fold was the good case, and the full set is a
   coin-flip.** Every guided fold against its own unguided twin, forced Boltz in brackets:

   | peptide | enclosed ung → guid (Bz) | engaged (Bz) | centroid sep Å (Bz) |
   |---|---|---|---|
   | `orig_f12` | 0.68 → **0.94** (0.96) | 14 → 16 (20) | 9.4 → **4.7** (4.2) |
   | `s3_esm2_f4` | 0.58 → **0.69** (0.80) | 17 → 18 (18) | 7.3 → **4.9** (7.0) |
   | `s3_orig_f12` | 0.43 → **0.79** (0.96) | 14 → 18 (20) | 10.4 → **6.4** (3.6) |
   | `ox1_esm1_f8` | 0.58 → 0.48 (0.94) | 10 → 11 (17) | 6.1 → 7.4 (4.1) |
   | `ox1_orig_f12` | 0.81 → 0.52 (0.99) | 17 → 15 (17) | 5.7 → 6.0 (1.2) |
   | `ox2_esm1_f8` | 0.74 → **0.79** (0.67) | 16 → 17 (17) | 5.3 → **4.4** (6.5) |
   | `ox2_orig_f8` | 0.66 → **0.93** (0.97) | 17 → 17 (17) | 7.5 → **4.6** (2.7) |
   | `ox3_orig_f8` | 0.54 → 0.45 (0.99) | 15 → 14 (17) | 10.3 → **8.0** (2.0) |

   Guided beats unguided on enclosure **5/8**, engagement 5/8, centroid separation 6/8 — and
   the three losses on enclosure are large (−0.10, −0.29, −0.09), so this is not a uniform
   improvement with noise on top. It splits by molecule: **octinoxate 3/3 improved, oxybenzone
   2/5**, and oxybenzone's guided folds stay far from their forced Boltz twins (0.48/0.52/0.45
   against 0.94/0.99/0.99) where octinoxate's close most of the gap. Naming the designed residues
   helps OF3 when its unguided placement was poor and hurts when it was already decent
   (`ox1_orig_f12`, best unguided cell at 0.81, drops to 0.52) — the same "docking rescues a bad
   pose and degrades a good one" shape the docked campaign found, which is worth stating as a
   prediction the MD can test rather than a conclusion.

   **This does not change the go decision** — a column that improves 5/8 and is measurably
   different in all 8 is worth 20 ns — but it does change what the column is *for*: it answers
   "does naming the residues help OF3 place the ligand", and the answer already looks like
   "sometimes, and it depends on the molecule". Size the campaign on that, not on the smoke fold.
   Metrics: `runs/{octinoxate,oxybenzone}/folds/metrics/openfold3p.csv`.

<!-- the smoke-fold-only reading, left in place above, concluded: a go
   for the campaign, not as a result.
2. **Conformer only. Answered from source, so no run was needed** (rc-foundry 0.2.0, the latest on
   PyPI, plus atomworks 3.0.3). The selection sets `GroundTruthConformerPolicy.REPLACE`, which
   swaps the RDKit reference conformer for the supplied coordinates. But
   `af3_reference_molecule.py` then calls `masked_center()` ("to avoid leaking 1D information") and
   `random_rigid_augmentation()`, and RF3's inference pipeline builds
   `GetAF3ReferenceMoleculeFeatures` with that augmentation left on (the default). So the input
   position is discarded before the model sees it: internal geometry survives, placement does not.
3. **Not reachable without patching, and pointless if patched.** `p_provide_inter_molecule_distances`
   defaults to 0.0, and `RF3InferenceEngine.__init__` builds `transform_overrides` itself with no
   entry for it. It also passes `transform_overrides=` explicitly, so supplying the key through
   `**kwargs` or Hydra would collide rather than merge. With p = 0, `ground_truth_template.py`
   masks every inter-molecule distogram pair. A subclass could set it to 1.0, but the template
   would still need a peptide–ligand pose, which only the forced Boltz fold can supply (the
   circularity noted above).
   **Upstream as of `production` 2026-10-16** (Foundry 1.0 refactor #404, AtomWorks 3 #401/#405):
   no contact or interface-template input; the README still says "could be further extended".
   Issue #160 is open, last updated 2026-04-28. **RF3 has no forced or guided mode to run. Close
   it unless upstream ships one**; the comment on #160 below is the only open action.

A fourth, non-blocking one: **ask upstream.** Foundry issue #160 has two unanswered requests for
ligand–protein contacts; a comment naming this use case (short designed peptide, 4–12 residue–ligand
contacts) costs nothing and might date the feature.

**Caveats to carry into the caption, not discover later.** f4/f8/f12 collapse to pocket sets of
different sizes, which is not the same dose of restraint; OF3's guided ipTM/ipSAE are not comparable
to its unguided ones (same issue #160 warning applies to guided runs generally); and a pocket-guided
OF3 column answers "does naming the designed residues help OF3 place the ligand", not "what do the
hints buy", which only Part 1's Boltz forced-minus-unforced difference measures.


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

### Pocket-guided lane: downstream status (2026-10-19)

Six of twenty-six legs recorded. Tables 1 and 2 of `cofold_binding_tables.md` carry a fourth
`OpenFold3 pocket-guided` column, spliced by `cofold_tables.py --write` and refreshed by re-running
it as each leg lands; the PDF is rebuilt from the markdown with `.venv/bin/python
code/render_cofold_pdf.py`. RESULTS.md needs nothing — it holds no co-fold lane tables at all (zero
mentions of `rf3cof`), so that campaign lives only in the cofold markdown.

**pepgraph: done.** `openfold3p` added to `LANES` and scored — median perplexity 7.2 / mean 9.7
across the nine, against the unconstrained `openfold3` lane's 10.2 / 10.7 on the same nine
sequences (diff +1.0 ± 0.5, 7 of 9 favouring guided), no zero-information structures. The guided
backbones are markedly more collapsed: 71 contacts against 36, `frac_aR` 0.19 against 0.79, so they
are not helices. `lane_comparison.json` was **merged, not replaced**: z comes off one RNG stream
consumed in structure order, so `openfold3`'s z moved on all 21 rows between runs while perplexity
moved on none. Every existing lane keeps its committed z. Read perplexity, not z, for this lane —
its own z was drawn in a two-lane run.

**morph: deliberately not extended, and this is the reason, not an omission.** A `--set matched`
pair is apo → holo within one model so that a ΔE is the ligand's effect. OpenFold3's
`pocket_constraint` takes a `ligand_chain_id`, so a guided *apo* fold is undefined and cannot be
built; pairing `of3apo → openfold3p` would measure the ligand's effect and the restraint's at once,
the exact confound `matched` exists to remove. The comparison the lane wants is
`openfold3 → openfold3p`, both holo, differing only by the restraint — a strain between two bound
states, which needs its own set rather than a row in `MATCHED_PAIRS`. Recorded in morph's
`run_batch.py` beside the pair list.

Still open and unactioned, both the user's call: the form of the `ox1_esm1_control` table cell after
the replicate (recommended `−26.09 / −9.93 (n=2, spread 16.2)` across `cofold_binding_tables.md`,
`RESULTS.md:422` and `RESULTS.md:509`), and targeted-11 versus full-202 replicate scope.
