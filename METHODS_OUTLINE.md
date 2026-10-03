# Methods — bullet outline for the manuscript

Untracked working document. Every entry is a fact pulled from the code, the README or the run
records, with the source named where it is not obvious. Nothing here is drafted prose.

Numbers in `code` font are defaults as coded; where a run used something other than the default it
is said so. Items marked **[decide]** are choices only you can make; **[gap]** marks something the
record does not currently contain.

---

## 1. Ligand and system definition

- Ligand: built-in `octinoxate` entry, **C17H24O3** — one CH2 short of real octinoxate
  (2-ethylhexyl 4-methoxycinnamate, C18H26O3). SMILES perceived from the stored geometry:
  `CCCC[C@H](CC)OC(=O)/C=C/c1ccc(OC)cc1`; the ester oxygen sits on the chain CH.
  (`NEXT_STEPS.md`, written 2026-09-24.)
- **[decide]** whether the paper calls this an octinoxate analogue throughout, or re-runs on the
  true C18 molecule. Every number in the project is for the C17 analogue.
- Ligand net charge 0, singlet; 20 heavy atoms (the denominator of every "engaged x/20" count).
- Binding-site box: ligand bounding box padded by `--box-pad` (default `0.0` Å); poses are sampled
  uniformly inside it.

## 2. Fragment library

- Ten side-chain **analogues**, each capped where the backbone would attach — not residues.
  Isoleucine entry is n-butane, leucine entry is isobutane; neither carries its namesake branch.
  Tyrosine is p-cresol, glutamic is propionate; all ten are chemically distinct.
- One line needed for reproducibility: two entries differ from the upstream library they came from
  (tyrosine and glutamic were wrong there), so results are not comparable with runs made against it.
  `attach_geom.py` checks electron parity, so a radical declared closed-shell cannot pass.
- `attach_geom.pose_attachment` locates each fragment's capping carbon (the CB marker) and the
  positions a CA can be grafted onto it.

## 3. Pose placement and scoring

- Placement: random position inside the padded box, random rotation about x, y, z
  (`get_frag_coordinates`), `--tries` `100` attempts per fragment.
- Acceptance: unrelaxed interaction energy below `--energy-cutoff` `500` kcal/mol, then BFGS
  relaxation of the fragment with all ligand atoms fixed (`FixAtoms`), `fmax 0.05` eV/Å, `steps 200`.
- Interaction energy `IE = E(complex) − E(ligand) − E(fragment)`, all at the relaxed pose geometry,
  converted with `23.06035` kcal/mol per eV.
- Potential: **FAIRChem UMA**, checkpoint `uma-s-1p2p1`, `FAIRChemCalculator(task_name="omol")` via
  ASE (`peptide_builder.py:1079`). CPU only on Apple Silicon (FAIRChem accepts `cpu`/`cuda`, no MPS).
  Charge and spin are passed per subsystem in `atoms.info`; totals combine as
  `charge_A + charge_B` and `spin_A + spin_B − 1`.
- Yield for octinoxate: **165 accepted poses across ten side-chain types**.
- Pose-class energetics, needed later for the charge-bias argument: charged analogues average
  **−10.79**, aromatic **−4.71**, neutral aliphatic **−2.94** kcal/mol.

## 4. Shell selection

- Greedy and order-dependent: walk the fragments in a given order, each taking its lowest-IE pose
  that does not come within `--clash` `1.4` Å of any pose already selected; repeated `--copies`
  times so a fragment may be placed more than once.
- Rotating which fragment picks first (`--all-starts`) × {1, 2} copies → **20 shells**, 8–13
  fragments, drawing on 48 distinct poses.
- Independence caveat that belongs in Methods, not just Limitations: only **18 of the 20 are
  distinct**, and **9 of those are strict subsets** of another (each single-copy shell is the first
  pass of its two-copy counterpart; the lysine and aspartic walks converge). → nine independent
  arrangements plus nested subsets.
- Three shells swept and folded in full:

  | | shell 1 | shell 2 | shell 3 |
  |---|---|---|---|
  | seed / selection | richest | second richest | tryptophan pose 2, chosen deliberately weak |
  | total fragment IE (kcal/mol) | −113.12 | −109.15 | −60.40 |
  | charged poses | 7 (58%) | 5 (42%) | 4 (33%) |
  | aromatic poses | 1 (8%) | 2 (17%) | 4 (33%) |
  | shared poses with shell 1 | — | 5 / 12 | fraction |

- **[decide]** whether shell 3's selection (an a-priori test of the selection criterion) is
  presented in Methods as a designed experiment or in Results as a follow-up.

## 5. Reachability measurement (the step that distinguishes the method)

- What is bridged is **attachment point to attachment point**, not centroid to centroid — the two
  are several Å apart for a large residue like tryptophan, and only the attachment points carry a
  direction.
- Construction: hold two poses fixed, graft a CA onto each capping carbon, grow an **ideal trans
  backbone** from one CA to the other through *k* glycine spacers by NeRF.
- Ideal geometry held fixed (Å, degrees): `N–CA 1.458`, `CA–C 1.525`, `C–N 1.329`, `C=O 1.231`,
  `CA–CB 1.530`, `N–H 1.010`, `CA–HA 1.090`; `N–CA–C 111.2`, `CA–C–N 116.2`, `C–N–CA 121.7`,
  `CA–C–O 120.8`, `CB–CA–C 110.5`, `CB–CA–N 110.5`, `C–N–H 119.0`; **omega fixed at 180°**.
- Free parameters: `t1` (rotation of residue A's tripod about CB–CA), `t2` (psi of A, shifted by
  the CB reference), and φ/ψ per spacer → **2 + 2k**.
- Closure condition: residue B's amide N must land on the circle compatible with B's fixed CA and
  CB — a 111° cone of radius 1.458 Å — i.e. two conditions; residual reported in Å.
- Optimisation: `scipy.optimize.minimize`, L-BFGS-B, `--restarts 24` random starts, `--seed 0`.
- Optional steric term (`--steric`) penalising backbone–ligand and backbone–fragment overlap.
- Sweep: every **ordered** pair at `k = 0…4` (`--kmin 0 --kmax 4`); 12 poses → 132 ordered pairs.
- Tolerance used downstream: closure < `--tol` `0.2` Å counts as connectable (README's rationale:
  absorbable by real backbone flexibility).
- Measured results to report as method validation:
  - **k = 0 closed in 0 of 810 attempts**; 2% even at a deliberately generous 1.0 Å tolerance.
  - At k = 2, a partner within 30° of the residue's CA→CB direction reaches **7.5 Å**; a partner
    approaching from the side reaches **11.1 Å**.

## 6. The spacer budget

- Counting argument: for a path through *m* steps the free dihedrals number `1 + m + 2·Σk` against
  `4m` closure conditions → generically closable when **Σk ≥ (3m − 1)/2**.
- Stated as a **global** budget, not a per-step test: surplus in one segment pays a shortfall in
  the next.
- Criterion is **necessary, not sufficient** — degrees of freedom plus a steric term; one
  three-fragment path that met its budget closed while another did not.

## 7. Path assignment

- Nodes = placed poses (worth their IE); edge a→b exists only if the sweep found a connection under
  tolerance; edge cost = *k* spacers + residue b.
- **Exact** search: Held–Karp dynamic program over (visited set, current node) — 2^n·n, no greedy
  walk. `--lam` optional charge penalty per residue (`0.0` in the runs).
- Output `design.json`: ordering, per-step spacer counts, sequence, budget check.
- Reproducible structural counts across the three shells (all 132 ordered pairs connectable, path
  through 12/12 poses in each):

  | | shell 1 | shell 2 | shell 3 |
  |---|---|---|---|
  | sequence | `RGGDGGKGGGGLGGGGKGIGEGWGGDGSGGEGS` | `RGEGGEGGKGGFGGDGGIGLGGSGWGGGGLGGS` | `YGGDGFGKGGGLGGGKGGGGLGGWGGEGSGGWGS` |
  | residues | 33 | 33 | 34 |
  | spacers used / budget | 21 / 16 | 21 / 16 | 22 / 16 |
  | worst closure | 0.196 Å | 0.120 Å | 0.188 Å |

## 8. Linker substitution

- `facebook/esm2_t12_35M_UR50D`; glycine spacers masked and refilled, designed positions untouched.
- **Sampled, not argmax**: candidates with probability ≥ `--prob-cutoff` `0.05`, `--seed` `1`,
  `--variants 2` per design → `esm1`, `esm2` per shell.
- Net charge is not controlled and moves a lot — record it per variant, since it drives the
  artifact in §13:
  - shell 1: design −1; shell 2 esm1 0, esm2 **+5**; shell 3 esm1 0, esm2 **+4**.

## 9. Null control

- `shuffle_control.py`: randomise which side chain occupies which slot, **keeping the glycine
  linker pattern**, `--seed 2`, plus one ESM2-filled variant (`shuffle_control_esm0`).
- Matched to the design on composition, length, 22 glycines, spacer pattern and net charge (0);
  **0 of 12 side chains left in their designed slot**.
- Rationale to state: shuffling the whole string instead produces 4–5 adjacent side-chain pairs
  with no glycine between them, randomising spacing as well as ordering — so the control tests
  **ordering only**.
- No fragment IE is fabricated for the control (it has no placed poses); `sequences.csv` untouched.

## 10. Co-folding

- **Boltz-2 2.2.1**, one peptide chain + ligand, `msa: empty`.
- Contact constraints written from the designed shell: ligand atom ↔ residue index, `max_distance`
  = contact cutoff `--cutoff 4.5` Å + `--margin 1.0` Å; the `--top N` tightest designed contacts.
- **Both flags are required and this must be in Methods**:
  - `force: false` constraints are **discarded by the featurizer**, not softened;
  - `force: true` alone is conditioning — the sampling-time guidance is gated behind
    `--use_potentials`, off by default. Runs used `--force` **and** `--boltz-args=--use_potentials`.
- Constraint ladder per sequence: **control (none), f4, f8, f12** → 12 structures per shell,
  3 shells = **36 structures**.
- Affinity head omitted: Boltz 2.2.1 calls `AllChem.Descriptors.MolWt`, which RDKit 2026.03.6 no
  longer exposes. Independent reason to omit: it spans only 1.65 kcal/mol across 37 sequences and
  scored a complex with the ligand 6.6 Å away confidently.
- Exact command (`boltz_env.run_boltz`): `boltz predict <yaml> --num_workers 0 --out_dir <dir>`,
  plus `--use_potentials` on the forced runs. Nothing else is overridden, so the sampling parameters
  are the 2.2.1 defaults, enumerated from the pinned install
  (`boltz_local/.venv/.../boltz/main.py:1050–1056, 1231`):

  | parameter | value | note |
  |---|---|---|
  | model | `boltz2` | |
  | recycling steps | 3 | |
  | sampling steps | 200 | |
  | diffusion samples | 1 | one structure per input, hence `model_0` |
  | step scale | 1.5 | Boltz-2 default (1.638 is Boltz-1) |
  | max parallel samples | None | |
  | output format | mmcif | |
  | num workers | 0 | the one non-default, passed explicitly |
  | **seed** | **None — no seeding** | |
  | affinity head | not requested | see the RDKit incompatibility above |

- Hardware and precision, from the fold logs: **MPS**, bfloat16 automatic mixed precision, on all
  **162** logged folds (`grep "GPU available" runs/octinoxate/boltz/*.log` → 162 × `(mps), used: True`).
  No fold in the project ran on CUDA or CPU.
- Structure taken: `*_model_0.cif` (the only sample, given `diffusion_samples 1`).
- **Consequence to state plainly**: with no seed and `diffusion_samples 1`, **the folds are not
  bit-reproducible** — re-running a sequence gives a different sample from the same distribution.
  Every structure in the paper is one draw. **[decide]** whether to say this in Methods (it is
  honest and cheap) or to re-fold a subset with `--seed` to quantify the spread. The latter is the
  referee-proof option and costs 2–3 min per fold.

## 11. Structural metrics (`check_fold.py`)

- `enclosed` — fraction of **200** sampled directions out of the ligand centroid that run into
  peptide; each ray given a 3 Å tolerance out to 12 Å. *Is the ligand in a shell or on a face.*
- `wrapped` — fraction of ligand heavy atoms with a peptide heavy atom within `--wrap` `4.5` Å.
  *Contact census.*
- `engaged` — count of the 20 ligand heavy atoms in contact (`--contact` `4.0` Å).
- `centroid sep` — ligand-to-peptide centroid distance, to be read against `peptide_rg`.
- `peptide_rg` — radius of gyration. Empirically the cleanest discriminator of the design
  objective: encapsulating folds 9.3–9.8 Å vs groove binders 11.8–14.9 Å, no overlap, where
  `enclosed` separates the same two groups by 0.015 and `wrapped` inverts them.
- `helical_fraction`, `beta_fraction` — φ/ψ basin counts with deliberately generous windows
  (helix φ −63 ± 35, ψ −43 ± 35) because Boltz geometry is unrefined; termini excluded.
- `helical_hbonds` / `nonlocal_hbonds` — backbone N···O < 3.3 Å split at sequence separation 3–5
  (i,i+4) vs >5 (strand pairing). **The split is the method point**: 81% β with 6 non-local bonds
  is an extended unpaired chain, not a sheet; a combined count ranks a pure helix highest.
- `hints` — requested contacts satisfied / requested. Report it; do not treat it as success
  (uncorrelated with everything else measured).
- Overlay to the designed shell (`overlay.py`): superposition on the ligand, RMSD 0.06–0.85 Å, so
  the designed-position comparison is sound; designed position counted reproduced within 3 Å.

## 12. Static energetics

- Hydrogens added to Boltz's heavy-atom CIFs with **pdbfixer**, then **relaxed with heavy atoms
  fixed** (`relax_hydrogens`, `fmax 0.10` eV/Å, `steps 75`).
  - Justification to state: pdbfixer's placement is **not deterministic** — same structure twice
    gave −14.47 and −13.31 kcal/mol (1.2 kcal/mol from hydrogens alone); three protonations of one
    CIF differed by up to 2.3 Å, RMSD 0.78 Å, heavy atoms untouched. Relaxing hydrogens brings
    reproducibility to **0.08 kcal/mol**. Stable on a Linux container; unstable on macOS.
- `interaction` = `E(complex) − E(peptide) − E(ligand)`, all three **at the complex geometry** →
  rigid interaction energy, contains no strain.
- `ligand strain` = `E(bound ligand, exactly as the complex relaxation left it) − E_ref`, a
  **single point**, no re-relaxation of any kind.
- `E_ref` = one shared reference per ligand (`ligand_reference.py`): 20 ETKDGv3 embeddings from the
  run's own SMILES (`--seed 0xC0FFEE`), MMFF94 optimised (`maxIters 2000`) and ranked, the **5
  lowest** relaxed with UMA at `fmax 0.01` eV/Å; lowest kept. Value used throughout:
  **−557162.371 kcal/mol-equiv**. `strain_global.py` refuses to run if the stored reference's
  SMILES or model does not match the run.
- Why one shared reference and not a per-structure relaxation — **one sentence in Methods, numbers
  in SI**: relaxing each complex's own bound pose measures every structure against a different
  local minimum, and the term then moves with the relaxation settings (+1.5 to +7.7 kcal/mol
  between `fmax` 0.1 and 0.01) and with pdbfixer's non-deterministic hydrogens (3.6 kcal/mol on
  bit-identical heavy atoms) rather than with the bound geometry.
- **The bound state is not re-relaxed, in whole or in part** — including no re-placement of
  hydrogens on the isolated ligand, which would insert an intermediate state and report only the
  second leg.
- Negative strain = failed reference (a bound pose below the supposed global minimum), not a
  finding; re-run the conformer search with more embeddings.
- Terms deliberately **not** reported, each disabled by default in the code, with the reason:
  - peptide strain (`--terms strain_peptide`) — a free peptide relaxed in vacuum collapses into a
    hydrogen-bonded ball: 237.8 kcal/mol on a 33-mer, not converged. Needs implicit solvent, which
    puts it in MM/GBSA territory.
  - explicit-water cavity desolvation — ~420 atoms, still descending 1.6 kcal/mol per step at step
    102 of 200, so the number would be set by the step budget rather than the physics.
- Known biases of the static energy, to state in Methods rather than defend in Results:
  - grows with system size (Spearman **−0.52** against peptide length) → compare matched lengths;
  - rewards net charge: across 36 folds, interaction energy correlates **−0.56** with charged-
    residue fraction;
  - a **positive** interaction energy means a broken structure — check heavy-atom contacts under
    2.6 Å. Of three found here, one was an n→π* approach 70° out of plane, one a CH···O bond 17°
    off the C–H vector, one genuinely unphysical (carbonyl O 2.55 Å from an ether O, no donor).

## 13. Molecular dynamics

- Structure split for MD (`cif_to_md.py`): Boltz writes peptide + ligand in one heavy-atom CIF;
  the ligand SDF is built by **assigning bond orders from the run's SMILES as a template** before
  adding hydrogens — without it the aromatic ring and the ester come out as single bonds.
- Force fields: **ff14SB** (`amber14/protein.ff14SB.xml`), **TIP3P** (`amber14/tip3p.xml`, ions
  bundled), ligand **GAFF-2.11** via `openmmforcefields` GAFFTemplateGenerator.
- Solvation: padding **1.5 nm**, **rhombic dodecahedron** box (`--box-shape dodecahedron`), 0.15 M
  ionic strength, neutralised, pH 7.0 protonation by pdbfixer (`prep-protein`).
  - Box shape is a methods point worth one line: same minimum-image distance in ~71% of the volume
    — 33,045 → **23,121** particles for the rod, 30% saving, no change in physics.
- System: nonbonded cutoff **1.0 nm**, PME, **HBonds** constraints, rigid water.
- Protocol: L-BFGS minimisation to a **10 kJ/mol/nm** force tolerance (no iteration cap);
  **5,000 steps** restrained equilibration with **100 kJ/mol/nm²** on protein heavy atoms;
  production **10,000,000 steps × 2 fs = 20 ns**, Langevin middle integrator, **300 K**,
  friction **1 ps⁻¹**, MonteCarlo barostat **1 atm** every **25** steps; frames every 500 steps.
- Platforms: local Apple Silicon **OpenCL** (single precision; the restrained equilibration phase
  has to run on CPU there) vs rented **A10G CUDA** via `modal_md.py`, which imports
  `openmm_md.dynamics.run` so remote and local trajectories come from identical code.
  - Rates: 1.72 ms/step at 8,647 particles (local) vs 0.376 ms/step at 23,121 (A10G); 20 ns of the
    23k system ≈14 h local vs **63 min**.
  - Trap worth a footnote: the `auto` platform probe does not list CUDA and silently selects CPU on
    a GPU host; CUDA passed explicitly. Warm-up required before timing (lazy compile/autotune reads
    2–5× too slow).
- Six trajectories, 20 ns each: `orig_f12`, `s3_orig_f12`, `s3_esm2_f4`, `s2_esm2_control`
  (designs) and `shuffle_control`, `shuffle_control_esm0` (nulls).
- **[decide]** whether to report the frame spacing difference explicitly: `orig_f12` was saved at
  1 ps, the rest at 10 ps, and release statistics are resampled to a common spacing (see §15).

## 14. MM/GBSA

- **MMPBSA.py** (AmberTools 24.8), **single-trajectory** binding mode with pre-stripped
  receptor/ligand prmtops (`-cp`/`-rp`/`-lp`).
- **GBn2, igb = 8**, saltcon = 0.15 M, mbondi2 radii set explicitly (ParmEd leaves RADII at 0 when
  converting from OpenMM, which would give a NaN GB term).
- Every frame used: `startframe=1, endframe=N, interval=1` over the wrapped solute-only trajectory
  (`analyze` writes `traj_wrapped.xtc`/`.pdb`; it is a required step, not optional).
- Parameters reused from the production system, so MM/GBSA and MD are the same force field.
- **No entropy term**; single-trajectory means BOND/ANGLE/DIHED are exactly zero by construction —
  so ΔG contains **no ligand strain**, which is why interaction terms and not the static sum are
  the fair comparison against UMA.
- Convergence test to report as methods, on `orig_f12` (leading *n* ns of one 20 ns run):
  −25.37 ± 0.50 (40 ps pilot), −21.52 (5 ns), −16.86 (10 ns), −14.61 (15 ns), **−13.71 ± 0.04**
  (20 ns); successive changes 3.85, 4.66, 2.25, 0.90 — halving. Point to make: a short window gives
  a **confidently wrong** answer, not a noisy one (the pilot's ±0.50 is sampling error over 40
  correlated frames of the starting pose and overstates binding by 11.7 kcal/mol).

## 15. Trajectory analysis

- `md_contacts.py`: heavy-atom contacts at `--cutoff 4.0` Å (matched to `fold_check.csv` so counts
  are comparable), centroid separation, peptide Rg, residence within 10 Å, release episodes, in
  `--deciles 10` windows.
  - **Release statistics resampled to a common `--resample-ps 10.0`** spacing — otherwise a finely
    saved run looks worse. Demonstration: `orig_f12` shows 8 released frames at 1 ps and **zero**
    at 10 ps, all eight being isolated single frames.
- `pair_contacts.py`: how many designed side-chain **pairs** contact the ligand simultaneously.
  Slots taken from the parent design's non-glycine positions, so shells, ESM variants and shuffles
  are handled identically. `--cutoff 4.0` Å, `--persist 0.05` (pair held ≥5% of the run).
  - Framing to state in Methods, since it defines what the construction delivers: *n* slots give
    **n(n−1)/2** pairwise options, chain connectivity hands over the *n−1* adjacent pairs free, so
    the informative band is **[n−1, n(n−1)/2]** = [11, 66] at twelve slots; the informative content
    is the non-adjacent pairs.
  - Derived measures: mean side chains engaged, % frames with ≥2 engaged, pairs realised, band
    position, pairs held ≥5%, long-range pairs (gap ≥ 7), dead slots.
- `md_frames.py`: representative frames for figures — ten-model ensemble superposed on the peptide
  and the **medoid** (frame nearest the window mean) of the last `--last-fraction 0.10` of the run,
  or an explicit `--window-ns` interval.
  - Disclosure that belongs in a caption, not hidden: `orig_f12` is shown at 11.0 ns rather than at
    the end, because its last 2 ns are the only stretch where the ligand is back in the cavity.

## 16. Statistics and reporting conventions

- Correlations by `correlate.py`, joining the geometry, energy and overlay CSVs on structure name;
  Pearson and Spearman reported with *n* beside each.
- **`--group` is mandatory practice, not a convenience**: these correlations invert between binding
  mechanisms. Documented instance — interaction energy vs centroid separation runs +0.23 / +0.70 /
  +0.93 at net charge −1 / 0 / +5, and vs enclosure −0.56 / −0.76 / +0.18; pooling all twelve
  reports −0.14 and +0.42, hiding both.
- n is small everywhere and should be said so once, in Methods: four structures per sequence,
  twelve per shell, 36 static structures, **six** trajectories.
- Spearman ρ at n = 5–6 quoted as indicative; the best predictor found (mean simultaneous
  engagement) fell **+0.90 → +0.83** when the sixth structure was added and got that structure's
  pair backwards. **[decide]** whether to quote both columns as the README does.

## 17. Software and versions

| component | version | role |
|---|---|---|
| Python | 3.12 | design pipeline |
| fairchem-core (UMA `uma-s-1p2p1`) | 2.22.0 | every static energy |
| ASE | 3.29.0 | BFGS, structure handling |
| RDKit | 2026.03.6 | SMILES, ETKDGv3, MMFF94, bond perception |
| PyTorch | 2.13.0 | UMA, ESM2 |
| transformers (ESM2 `esm2_t12_35M_UR50D`) | 5.17.0 | linker substitution |
| pdbfixer / OpenMM (design env) | 1.12.0 / 8.6.1 | hydrogen addition |
| SciPy | 1.18.1 | L-BFGS-B closure search, Held–Karp support |
| Boltz | 2.2.1 | co-folding with contact constraints |
| OpenMM (MD env) | 8.4.0 | dynamics |
| openmmforcefields / openff-toolkit | 0.15.1 / 0.18.0 | GAFF2 ligand parameterisation |
| MDTraj / ParmEd | 1.11.1 / 4.3.1 | trajectory and topology handling |
| AmberTools (MMPBSA.py) | 24.8 | MM/GBSA |

- Note for the record: **FAIRChem ≥ 2.22 refuses the UMA 1.0 checkpoint**, directing to
  `fairchem-core<=2.21.0`; `uma-s-1p2p1` is used throughout.
- Two environments, called at their boundary: the design pipeline is pip-installable; GAFF2/OpenFF
  parameterisation is not, so MD/MM-GBSA lives in a conda environment (`openmm-md`).

## 18. Reproducibility and data

- Seeds fixed and recorded: placement `--seed`, closure sweep `0`, ETKDG `0xC0FFEE`, ESM2 fill `1`,
  shuffle `2`.
- Every relaxed geometry is written, not derived and discarded (`boltz/structures/`), because each
  costs 45–60 min; per-term results are checkpointed to `partial_*.json` as computed.
- UMA runs logged separately (`uma_logs/`) from non-UMA runs (`logs/`): step counts, timings and
  convergence are the only record of how an energy was arrived at.
- **[decide]** data availability statement: repository URL, and whether trajectories
  (`md/<structure>/prod*/`) or only the stripped solute + `mmgbsa_summary.csv` are deposited.

---

## Ordering suggestion for the section

1. Ligand and fragment library (§1–2) — short.
2. Placement and scoring (§3), shell selection (§4).
3. Reachability and the budget (§5–6), assignment (§7) — **this is the novel method; it deserves
   the most space and the geometry figure.**
4. Sequence construction: linkers (§8), controls (§9).
5. Co-folding and constraints (§10).
6. Metrics (§11) and static energetics (§12) — the strain-reference argument (§12) is long; it may
   belong in SI with a two-sentence pointer.
7. MD and MM/GBSA (§13–14), trajectory analysis (§15).
8. Statistics (§16), software (§17), availability (§18).

## Candidates for SI rather than Methods

- The fmax sweep and the pdbfixer determinism measurements behind §12's two justification bullets.
- The MM/GBSA convergence table (§14) — arguably Results, since "run 20 ns, do not trust 5" is a
  finding about the protocol.
- Platform/timing table (§13) — SI, unless the paper makes a cost argument.

## Scope note

This outline documents the pipeline as it currently stands. The superseded centroid-linker
construction (`build_sequence`'s `round(d / 3.8) − 1`), its twenty-odd designs, the
`--linker-slack` sweep, the cyclic variants, the random-sequence control and the abandoned
per-structure strain reference are deliberately absent: they belong in the Introduction as
motivation or nowhere, not in Methods. The measured reachability statistics that *justify* current
choices (§5's 0/810 at k = 0, §6's budget) are kept, since they are properties of the method being
described.

## Gaps to close before the section is complete

- **[gap]** Shell 3's pose overlap with shells 1 and 2 is described qualitatively but never counted;
  shell 2 vs shell 1 is 5/12.
- **[gap]** Shell 1's twelve folds were never scored for reproduced designed positions
  (`boltz/overlay.csv` holds sequence-level comparisons only), so that column is blank, not zero —
  say so or fill it.
