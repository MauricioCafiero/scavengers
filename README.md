# PeptideBuilder

Design a peptide that binds a small molecule by deciding where the side chains should go first, and
working out the backbone afterwards.

The tool places amino-acid side-chain analogues around a ligand, scores each placement with a machine
learning interatomic potential, and then — this is the part that distinguishes it — measures whether a
real peptide backbone can actually reach from one placed side chain to the next before committing to a
sequence. Designs that ask for geometry no backbone can follow are rejected at that stage instead of
being discovered to be wrong after an expensive folding calculation.

The result is a sequence with the right number of spacer residues in the right places, which can then
be co-folded with the ligand, and the fold checked against what was asked for.

**Start with [PROJECT.md](PROJECT.md)** if you want the state of the project rather than the tool: it
is the one document that holds the pipeline as it runs now, every peptide in both this repository and
[`../boltzgen_local`](../boltzgen_local), the MM/GBSA and retention comparison between them, and what
is still open. This README is the tool — install it, run it, what each file does, what the metrics
mean, and where it goes wrong.

---

## Contents

- [The problem this solves](#the-problem-this-solves)
- [End-to-end workflow](#end-to-end-workflow)
- [Installation](#installation)
- [Repository layout](#repository-layout)
- [What each file does](#what-each-file-does)
- [Where output goes](#where-output-goes)
- [Results](#results) — summary; the full record is in [RESULTS.md](RESULTS.md)
- [What the metrics mean](#what-the-metrics-mean)
- [Traps](#traps)
- [Limitations](#limitations)
- [Credits and provenance](#credits-and-provenance)
- [Citation](#citation)
- [Licence](#licence)

---

## The problem this solves

Fragment-based placement is good at finding where a side chain wants to sit. It puts an arginine
against a carbonyl, a tryptophan across an aromatic face, and each placement is individually
favourable. The difficulty is that the placements are found independently, so the resulting set is a
*shell* of side chains surrounding the ligand with no regard for whether one polypeptide chain can
visit all of them.

Earlier versions of this pipeline chained the placed fragments by the distance between their centroids
and inserted `round(d / 3.8) - 1` glycines between consecutive pairs. That estimate is wrong in two
ways. It measures the wrong points — what has to be bridged is the attachment point where the backbone
joins each side chain, not the fragment's centre of mass, and for a large residue like tryptophan those
are several ångströms apart. And it has no notion of direction: two side chains the same distance apart
can be easy or impossible to connect depending on which way their attachment points face.

The consequence was measurable. Of 225 designed side-chain positions across twenty designs, folding
reproduced 8 — about 3.6%. The shells were not being lost by the folding model; they were being
transcribed into sequences that could not hold them.

This version measures reachability directly. For every pair of placed side chains it grows an ideal
trans peptide backbone out of one attachment point and tries to land it on the other, with bond lengths,
bond angles and planar peptide bonds all held at their proper values, varying only the rotatable
dihedrals. What comes back is a number in ångströms: how far the connection misses by. A miss under a
couple of tenths of an ångström is absorbable by the small flexibility real backbones have. A miss of
several ångströms is not.

Three results follow from that measurement, and they shape the whole workflow.

**Two independently placed side chains essentially cannot be adjacent residues.** With zero spacers
between them the connection succeeded in 0 of 810 attempts, and still only 2% at a deliberately
generous 1.0 Å tolerance. Pinning a residue's CA and CB leaves just two rotatable dihedrals against
four conditions needed to place the next residue's CA where its pose demands, so the system is
over-determined. Two fragments placed independently around a ligand will not satisfy it by luck.

**Reachability depends on approach angle as much as on distance.** At two spacers, a partner sitting
behind its own fragment — within 30° of that residue's CA→CB direction — can be reached out to 7.5 Å,
while one approaching from the side reaches 11.1 Å. Half again the distance for the same linker cost.
The backbone leaves CA pointing away from CB, so a partner behind it requires the chain to double back.

**The spacer requirement is a budget over the whole path, not a test on each step.** For a path through
*m* steps the free dihedrals number `1 + m + 2·Σk` against `4m` closure conditions, so a path is
generically closable when

```
Σk ≥ (3m − 1) / 2
```

Surplus freedom in one segment pays for a shortfall in the next, which is why the condition is global.
This reproduces every closure rate measured, and it explains the older behaviour: a ten-fragment design
needs about thirteen spacers, the centroid formula supplied two, and adding one glycine per gap
(`--linker-slack 1`) took the total to eleven — which is why that setting improved the fold and why
adding more made it worse.

---

## End-to-end workflow

```
peptide_builder.py run     place side-chain analogues around the ligand, score each with UMA
        │
        ├─ poses/*.xyz, energies.csv
        │
condense.py sweep          grow a real backbone between every pair of placed poses, at 1–4 spacers,
        │                  and record how far each connection misses
        ├─ condense_<shell>.csv
        │
assign.py --save           choose the order through the side chains and the spacers each step needs,
        │                  subject to the measured reachability and the global budget
        ├─ design.json  (path, spacer counts, sequence, budget check)
        │
fill_linkers.py            optional: replace the glycine spacers with residues ESM2 finds probable
        │
boltz_hints.py --run       co-fold each sequence with Boltz-2, supplying the designed shell as
        │                  contact constraints, plus an unconstrained control
        ├─ boltz/boltz_results_*/…/*_model_0.cif
        │
check_fold.py              is the ligand bound, how much of it is enclosed and wrapped, and were
        │                  the requested contacts honoured
        ├─ boltz/fold_check.csv
        │
binding_energy.py          interaction energy for each folded complex, and the bound-state geometry
        │                  the strain step needs
        ├─ boltz/binding_*.csv, boltz/structures/*_ligand_bound.xyz
        │
ligand_reference.py        the ligand's lowest conformer, as the one reference every strain is
        │                  measured against. Computed once per ligand and reused by every shell
        ├─ ligand_reference.json
        │
strain_global.py           ligand strain per fold: E(bound) − E(reference), a single point on the
        │                  bound ligand exactly as the complex relaxation left it
        ├─ boltz/strain_*.csv
        │
make_figures.py            flat filenames, a manifest joining structures to their numbers, and a
                           PyMOL script that loads everything superposed on the ligand
```

Two optional branches sit alongside this. `condense_chain.py` builds the entire path as one covalent
chain rather than folding a sequence, which verifies that the chosen spacer counts really are sufficient
for the whole chain instead of only pairwise; it is expensive and off by default. `overlay.py`,
`compare_folds.py` and `random_control.py` provide comparisons against the designed shell, against each
other, and against random sequences of matched length.

A single command runs steps 3 through 8:

```bash
python code/design_test.py runs/octinoxate \
    --poses "arginine:3,lysine:5,aspartic:21,..." \
    --pairs-csv condense_shell12.csv --variants 2
```

---

## Installation

One environment runs everything except the folding, and Boltz-2 is optional and discoverable rather
than assumed at a fixed path.

### This repository

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

That gives `fairchem-core`, `ase`, `numpy`, `torch`, `rdkit`, `py3Dmol`, `pdbfixer` and
`transformers`. `pdbfixer` pulls `openmm`, which it is built on. Python 3.12 is what this has been
run on.

If the environment was made with `uv`, note that it has no `pip` in it, so install with
`uv pip install --python .venv/bin/python -r requirements.txt` — `.venv/bin/python -m pip` fails with
"No module named pip".

`transformers` currently resolves `huggingface-hub` below 2.0. That is safe here: `fairchem-core`
declares `huggingface-hub>=0.27.1` with no upper bound and uses exactly one function from it,
`hf_hub_download`. A UMA single point gives the same energy either side of the change.

The UMA potential is downloaded on first use from Hugging Face and needs a token with access to the
Meta FAIR-Chem repository:

```bash
huggingface-cli login
```

The default checkpoint is `uma-s-1p2p1`. Pass `--model` to change it.

### Boltz-2, the main external dependency

Boltz-2 does the co-folding. It is heavy and GPU-bound, so it may live in its own environment rather
than this one — but nothing here assumes where. `code/boltz_env.py` resolves how to run it, first
match winning:

1. `--boltz-cmd`, or `$PEPTIDEBUILDER_BOLTZ_CMD` — a complete command, used as given
2. `--boltz-venv`, or `$PEPTIDEBUILDER_BOLTZ_VENV` — a virtualenv holding Boltz
3. `boltz` importable in this environment — the plain `pip install boltz` case
4. `boltz` on `$PATH`
5. an MPS wrapper beside the venv, if that layout happens to exist

If none match, the error lists these options rather than failing obscurely.

The only other external tool is optional and used by the redocking scripts alone: `vina_redock.py` and
`dock_compare.py` call AutoDock Vina and Open Babel out of a sibling repository,
`~/python_mac/dock_assist`, where they already exist — Vina as the copy vendored in `dockstring`, Open
Babel from the `openbabel-wheel` package. Nothing is installed here for them, and they are run with
that repository's interpreter rather than this one's:

```bash
~/python_mac/dock_assist/dock-env/bin/python code/vina_redock.py --all
~/python_mac/dock_assist/dock-env/bin/python code/dock_compare.py --all
```

Everything else in this repository runs without them; only the [redocking
section](RESULTS.md#an-independent-check-on-the-pose-redocking) needs them.

**No Boltz at all?** `code/boltz_offline.py export` writes every input into one folder with the
command to run, you fold them anywhere — a GPU box, Colab, a cluster queue — and `import` puts the
structures back where the scoring scripts look. That route is also the practical way to use a GPU for
scoring, which is where the time actually goes.

### Hardware

Everything here runs on CPU. FAIRChem does not support Apple's MPS backend — it accepts only `cpu` or
`cuda` — so on an Apple Silicon machine the potential runs on the CPU cores while Boltz-2 uses the GPU
through its own wrapper.

That asymmetry dominates the runtime. On a six-core Apple Silicon laptop, co-folding a 33-residue
peptide with its ligand takes two to three minutes, while scoring that same complex takes 30 to 70
minutes depending on size. If a CUDA device is available, `binding_energy.py` picks
it up automatically and the scoring bottleneck largely disappears.

Memory matters more than core count. A scoring process holds 0.5 to 2.5 GB depending on what it is
doing, and running two concurrently on a 16 GB machine will push it into swap. Run them serially.

---

## Repository layout

```
peptidebuilder/
├── code/                        all scripts; see the table below
├── runs/                        pipeline output, one directory per ligand
│   └── octinoxate/              the worked example below
├── requirements.txt
├── LICENSE
└── README.md
```

## What each file does

### The design pipeline

| file | purpose |
|---|---|
| `peptide_builder.py` | ligand and fragment definitions, pose placement, UMA scoring, pose selection, sequence assembly, and the command line (`run`, `combine`, `list`) |
| `attach_geom.py` | finds each fragment's capping carbon — the atom marking where CB sits — and the positions a CA can be grafted onto it. Also audits the fragment library |
| `condense.py` | grows an ideal trans backbone between two placed poses and measures the closure error. `sweep` does this for every ordered pair at a range of spacer counts |
| `condense_analyze.py` | turns a sweep into lookup tables: what closes, at what distance, at what approach angle |
| `assign.py` | exact Held–Karp search for the best order through the placed poses, subject to reachability, reporting the global spacer budget |
| `fill_linkers.py` | replaces glycine spacers with residues ESM2 finds probable, sampling rather than taking the argmax |

### Folding and scoring

| file | purpose |
|---|---|
| `boltz_hints.py` | writes a Boltz-2 input with the designed shell as contact constraints, and co-folds it |
| `boltz_check.py` | the simpler path: co-folds sequences straight from `sequences.csv` |
| `check_fold.py` | per fold: is the ligand bound, how enclosed and how wrapped is it, were the contacts honoured |
| `binding_energy.py` | interaction energy, ligand and peptide strain, cavity desolvation. Writes every relaxed structure it produces |
| `uma_binding.py` | reads Boltz CIFs and adds hydrogens — `parse_cif`, `protonate_peptide`, `protonate_ligand`, `peptide_charge` — used by everything above. Its own command line predates `binding_energy.py` and computes the interaction term alone |
| `solvate.py` | wraps a molecule in explicit waters, for desolvation energies |

### Analysis and comparison

| file | purpose |
|---|---|
| `check_design.py` | walks an existing design's own ordering and reports which steps are physically connectable |
| `compare_folds.py` | pairwise structural comparison, superposed on the ligand and on the peptide |
| `overlay.py` | superposes a fold on the designed shell and counts reproduced side-chain positions |
| `random_control.py` | null model: random sequences of matched length |
| `ligand_reference.py` | finds the ligand's lowest conformer once — 20 ETKDG embeddings, MMFF-ranked, the five best relaxed with UMA — as the single reference every strain is measured against |
| `strain_global.py` | ligand strain for every fold against that reference. A single point on the bound ligand as the complex relaxation left it — `structures/<name>_ligand_bound.xyz`, or the step-0 energy of that fold's free-ligand block in the scoring log, which is the same state and lets folds scored before geometry saving be corrected without re-running |
| `correlate.py` | joins every results CSV on structure name and correlates the properties against each other. `--group` reports each subset separately, because these correlations invert between binding mechanisms |
| `make_figures.py` | flat filenames, a manifest, and a PyMOL loading script |
| `cif_to_md.py` | splits a Boltz complex into the protein PDB and ligand SDF the OpenMM MD pipeline wants, assigning ligand bond orders from the run's SMILES. See [dynamics and MM/GBSA](RESULTS.md#what-happens-after-this-pipeline-dynamics-and-mmgbsa) |
| `md_contacts.py` | contacts, centroid separation, peptide Rg, residence and ligand-release episodes in windows across a trajectory: whether the designed pose survived, which one ΔG cannot tell you. Heavy atoms only, so counts stay comparable with `fold_check*.csv`, and release statistics are resampled to a common frame spacing because otherwise a finely saved run looks worse than a coarsely saved one |
| `md_frames.py` | representative end-of-run frames as PDBs for figures — a ten-model ensemble superposed on the peptide, and the single frame nearest the window mean. `--window-ns` picks a different interval where the end of a run is unrepresentative |
| `pair_contacts.py` | how many designed side-chain **pairs** reach the ligand over a trajectory, against the [n−1, n(n−1)/2] band. Mean simultaneous engagement is the best predictor of MM/GBSA ΔG measured here. Slots come from the parent design, so shells, ESM variants and shuffles all work |
| `shuffle_control.py` | the null model: the same residues in a random arrangement, keeping the linker pattern, with one ESM2 variant. No scoring — it has no poses, and `sequences.csv` is left alone |
| `vina_redock.py` | redocks the ligand into its own folded peptide with AutoDock Vina, from the receptor and reference ligand the MD prep already wrote. `--md-root` adds structures prepared outside this repository — BoltzGen's two go in with `--md-root ~/python_mac/boltzgen_local/md` — so they land in the same tables and compare directly; a `source` column records where each came from. Vina and Open Babel come from `~/python_mac/dock_assist`; nothing is installed here. Reports Vina's score, symmetry-corrected RMSD to the Boltz pose, and `check_fold.py`'s wrapping and enclosure on every pose, so a docked pose and a fold are measured by identical cutoffs. See [redocking](RESULTS.md#an-independent-check-on-the-pose-redocking) |
| `dock_compare.py` | splits each pose difference into its parts — rigid-body translation and rotation from a Kabsch fit, internal torsion change as the RMSD left after that fit, and the angle between the two poses' head→tail vectors — then asks separately, by contact-residue Jaccard, whether it is even the same pocket. Writes a `compare.pml` per structure |
| `dock_vs_md.py` | superposes each MD medoid's peptide onto the docking receptor and asks whether the ligand ended nearer the Boltz pose or a docked one. Reports the peptide's own superposition RMSD, because when that is as large as the pose differences the comparison cannot settle anything |

### Running a whole shell, and renting a GPU

| file | purpose |
|---|---|
| `run_shell.sh` | the parameterised runner for one shell end to end: design, variants, folding, checking, scoring, overlay, strain, figures. Stage 0 preflights the forced-contact flags before anything expensive starts, `SKIP_SCORING=1` stops after folding so scoring can go elsewhere, and strain runs at stage 6b *after* scoring because scoring is what writes `_ligand_bound.xyz` |
| `modal_score.py` | runs the scoring step on a rented GPU, shipping the CIFs out and the logs and structures back through `runs/<ligand>/modal/<tag>/`. It calls `binding_energy.py` unmodified, so a remote score and a local one are the same computation |
| `modal_md.py` | runs the MD production leg on a rented GPU by importing `openmm_md.dynamics.run`, so a remote trajectory and a local one come from the same protocol. Keeps the full solvated trajectory on a Volume and returns only the stripped solute |

Both Modal scripts leave the local code untouched: nothing in the list above knows or cares whether it
was run on this laptop or rented hardware.

### The alternative route

| file | purpose |
|---|---|
| `condense_chain.py` | builds a whole path as one covalent peptide onto the placed poses, rather than folding a sequence. Verifies that spacer counts suffice for the entire chain |
| `condense_strain.py` | full-atom build and UMA relaxation of a condensed pair, giving the strain of holding two side chains at their designed positions |
| `design_test.py` | driver chaining ordering, substitution, folding, checking and scoring |

---

## Where output goes

Everything for one ligand lives under `runs/<ligand>/`.

```
runs/octinoxate/
├── poses/                       one .xyz per accepted pose: ligand + one fragment
├── energies.csv                 every pose's interaction energy, and which were selected
├── ligand.xyz                   the ligand alone
├── combined.xyz                 ligand + all selected fragments: the designed shell
├── sequences/                   one .xyz per design, ligand + its placed fragments
├── sequences.csv                every design's sequence, start residue, fragment count, total IE
├── ligand_reference.json        the ligand's lowest conformer and the energy every strain uses
├── condense_<shell>.csv         reachability sweep: closure per pair, spacer count, CA choice
├── design.json                  chosen path, spacer counts, sequence, budget check
├── condense/                    condensed structures from the chain route
├── uma_logs/                    logs of every run where UMA ran: BFGS step counts, timings,
│                                convergence -- the only record of how an energy was arrived at
├── logs/                        logs of runs with no UMA in them: geometry sweeps, Boltz driving
├── boltz/
│   ├── *.yaml, *.log            Boltz inputs and logs
│   ├── boltz_results_*/         folded complexes, as .cif
│   ├── structures/              every geometry a relaxation produced
│   ├── fold_check.csv           enclosure, wrapping, contacts, hints honoured
│   ├── binding_*.csv            interaction energy per complex
│   ├── strain_*.csv             ligand strain per complex, against the shared reference
│   └── partial_*.json           per-term results, saved as computed
├── figures/                     shell 1; `figures_shell2/` and `figures_shell3/` alongside
│   ├── *.cif                    every fold under a readable name
│   ├── designed_shell.xyz       the arrangement that was asked for
│   ├── manifest.csv             every structure joined to its numbers
│   ├── load_folds.pml           PyMOL script: loads all, superposed on the ligand
│   ├── style.pml                the representation, including the residue labels
│   └── *.png                    the renders reproduced in this README
├── modal/<tag>/                 staging for a rented-GPU run, before its logs and structures
│                                are copied back into the paths above
├── md/figures/                  the MD renders reproduced in this README, one per run
└── md/<structure>/              dynamics and MM/GBSA, one directory per structure
    ├── protein_fixed.pdb        and ligand_prepped.sdf: the prepared inputs
    ├── system/                  the solvated, parameterised system
    ├── prod*/                   trajectory, energy log, wrapped solute
    │   ├── frames_last.pdb      ten end-of-run models, and frame_medoid.pdb, for figures
    │   └── mmgbsa/              FINAL_RESULTS_MMPBSA.dat and its prmtops
    └── mmgbsa_summary.csv       every leg's decomposition in one table (at md/ root)
└── dock/                        redocking, one directory per structure, plus three tables
    ├── dock_summary.csv         per structure: source, score, RMSD to the predicted pose, and
    │                            wrapping and enclosure for the reference, pose 1 and closest pose
    ├── dock_poses.csv           every pose's geometry; pose 0 is the Boltz reference
    ├── pose_differences.csv     translation, rotation, internal RMSD, head/tail angle, Jaccard
    ├── md_vs_dock.csv           which pose the dynamics ended nearer
    └── <structure>/             receptor and ligand PDBQT, poses.pdbqt, poses.sdf, vina.log,
                                 and compare.pml to look at the poses against the reference
```

`load_folds.pml` loads and superposes; it does not style. Run `@style.pml` beside it for the
representation used in this README's figures, which is what draws the residue labels. `style.pml` is
written only when absent, so a hand-tuned copy in one figure directory is never overwritten — and, for
the same reason, is not inherited by the next shell's directory either.

Logs are split by which engine ran, and both are kept rather than ignored. A run that invoked UMA
goes in `runs/<ligand>/uma_logs/`: it holds each relaxation's step count, timing and convergence,
which is the only record of how an energy was arrived at, and regenerating one means re-running the
hours that produced it. A run with no UMA in it — a geometry sweep, a Boltz fold, an analysis pass —
goes in `runs/<ligand>/logs/`. Redirect a run's output to whichever applies, since the scripts print
to stdout and do not choose a path themselves:

```bash
python code/binding_energy.py runs/octinoxate … > runs/octinoxate/uma_logs/score_x.log 2>&1
python code/condense.py sweep runs/octinoxate … > runs/octinoxate/logs/condense_x.log 2>&1
```

Boltz is the one exception: it writes its own per-fold log beside its structures, at
`boltz/<name>.log`.

Three conveniences worth knowing. `partial_*.json` holds each energy term the moment it is computed,
so interrupting a long scoring run does not discard what it has already paid for. `structures/` holds
the hydrogen-relaxed complex, the bound ligand and the relaxed free ligand for every complex scored —
these cost 45 to 60 minutes each to produce, so they are written rather than derived once and thrown
away. And `figures/` exists because the folded structures are otherwise buried at
`boltz/boltz_results_<name>/predictions/<name>/<name>_model_0.cif`, which is tedious to load twelve of.

The large full-system trajectories are not here. The big `.dcd` files `.gitignore` keeps off git —
`traj.dcd` and the `traj_ext.dcd` — are also kept off this drive: they live on **iCloud Drive at
`Dynamics_trajectories/peptidebuilder_POC/<leg>/<run>/`**, mirroring the names under `md/`, evicted
locally to save space and never deleted. Only the big `.dcd` moved; the small `traj.nc` and
`checkpoint*.chk` stay here. The committed `traj_wrapped.xtc`, prmtops, `energy.csv` and
`FINAL_RESULTS_MMPBSA.dat` are enough to redo the analysis; pull a leg's folder back from iCloud only
to re-window or re-score from the full trajectory.

---

## Results

The measured results have their own file: **[RESULTS.md](RESULTS.md)**.

It holds the three worked examples, the dynamics and MM/GBSA work, the redocking and GNINA rescoring,
and the two-scorer pose matrix -- about 1,300 lines that had grown to be two thirds of this README and
were burying the part that explains how to run anything. What stays here is the pipeline: what it does,
how to install it, what each file is for, what the metrics mean, and where it goes wrong.

Shortest useful summary, with the detail and the caveats in RESULTS.md:

* **The designs from this pipeline are the best binders measured here.** Across nine structures with
  20 ns of dynamics and MM/GBSA, the four peptidebuilder designs average **-18.84 kcal/mol** with 84.1%
  mean ligand retention and one release episode between them; the three BoltzGen structures average
  **-12.71** with 51.6% retention and forty-three releases. The best design, `s3_orig_f12` at **-24.33**,
  leads the best BoltzGen structure by 4.67 kcal/mol; the worst retention is BoltzGen's `bg33_1`
  (14.5% within 10 Å — but it kept ligand contact the whole run, a surface roamer) and the worst
  energy is BoltzGen's `bg33_2` (-8.02, bound well for 15 ns then fully released). n is 4 against 4
  and the two pipelines' ligands differ in configuration, so this characterises these structures
  rather than the methods.
* **A design does not reliably beat a shuffle of itself — on octinoxate.** Its nulls average -14.73, but
  `shuffle_control` at -15.13 outscores the design `orig_f12` at -13.71. That was the sharpest standing
  criticism of the design method. **Oxybenzone answers it, twice:** two shuffles of `design_shell2.json`
  lose to their designs by **18.9 and 8.8 kcal/mol** and lose in all six cells of the pose matrix, so
  docking does not rescue them. Two matched pairs on one molecule and one shell, which does not make it
  general — but the criticism no longer stands unopposed. See
  [the oxybenzone tables](RESULTS.md#oxybenzone--mmgbsa-Δg-by-starting-pose).
* **The shell cannot be reproduced by a folded peptide.** 8 of 225 designed side-chain positions land
  within 3 A, across every variant tried.
* **Static scores do not rank these structures.** UMA interaction energy ranks against MM/GBSA at
  Spearman -0.80 -- inverted, not merely weak. GNINA's `CNNaffinity` does correlate across structures
  (rho = -0.857) and picks the better docked starting pose in four of five rows, but fails badly on the
  best binder, by 7.56 kcal/mol.
* **Dock only when the co-folded pose is poorly retained.** Across five peptides run from three
  starting poses each, co-folding wins where its pose holds the ligand 100% of the time, and a docked
  start wins where it does not -- by 0.84 kcal/mol at 59.6% retention and 1.23 at 41%. Equilibration
  leaves the ligand free, so it slides 1.2 to 15.1 A before production starts depending on the run; that
  figure is recorded per leg and is what says whether a given cell is worth reading.

---

## What the metrics mean

Two of these are easy to confuse, and the fact that they disagree is the point.

**`enclosed`** — the fraction of directions out of the ligand's centre that run into peptide. Two
hundred directions are sampled, each treated as a ray with a 3 Å tolerance out to 12 Å. This asks
whether the ligand sits in a shell or against a face.

**`wrapped`** — the fraction of the ligand's heavy atoms with a peptide heavy atom within 4.5 Å. This
asks how much of the ligand is in contact.

Both are computed by `check_fold.geometry`, which takes two coordinate arrays rather than a file, so a
Boltz fold read from a `.cif` and a docked pose read from a `.pdbqt` are measured by one definition and
one set of cutoffs. Anything that holds a peptide and a ligand in a shared frame can be measured the
same way, and the numbers stay comparable across the two.

They disagree usefully. `orig_control` contacts three quarters of the ligand's atoms while enclosing
less than half of the directions out of it, because the ligand is stuck to one face. Reporting wrapping
alone over-rates that structure considerably.

**`centroid sep`** — distance between the ligand's centroid and the peptide's. Compare it to the
peptide's radius of gyration: `orig_control`'s 20.9 Å against an Rg of 13.7 Å says the ligand is not
inside the peptide at all.

**`peptide_rg`** — the peptide's radius of gyration. Reported because it turned out to be the best
single discriminator of the design objective: across the second worked example's twelve folds, every
structure that visibly encapsulates the ligand has Rg 9.3–9.8 Å and every one that holds it in a
surface groove has 11.8–14.9 Å, with nothing in between, while `enclosed` separates the same two groups
by 0.015 and `wrapped` gets them backwards.

**`helical_fraction`** and **`beta_fraction`** — the fraction of residues whose backbone φ/ψ fall in
the α-helical and extended-β basins, using generous windows (helix φ −63±35, ψ −43±35) because Boltz's
geometry is unrefined and a tight window reports zero for folds that are plainly helical. Together
they explain binding mode: a substituted sequence often comes back as a rigid helix that can only grip
the ligand along its surface, and encapsulating then requires the helix to break.

**`helical_hbonds`** and **`nonlocal_hbonds`** — backbone N···O contacts under 3.3 Å, split by
sequence separation: 3 to 5 residues apart, which is the i,i+4 bond a helix makes with itself, against
more than 5 apart, which is where strand pairing would appear. **The split is the point.** A high
`beta_fraction` means an extended backbone, not a sheet, and only the non-local count distinguishes
them — the glycine design's unconstrained fold reaches 81% β with just 6 non-local bonds, so it is
extended and unpaired rather than pleated into a sheet. Counting all separations together hides this
and ranks a pure helix highest, since 40-odd i,i+4 bonds swamp everything while pairing nothing.

**`interaction`** — `E(complex) − E(peptide) − E(ligand)`, all three at the complex geometry, so it is a
rigid interaction energy containing no strain. Hydrogens are relaxed first with heavy atoms fixed,
which leaves Boltz's predicted geometry untouched while removing the arbitrariness of where pdbfixer
placed them.

**`ligand strain`** — `E(ligand at its geometry in the complex) − E_ref`, the conformational price the
ligand pays, where **`E_ref` is one shared reference: the ligand's own lowest conformer, computed once
per ligand** by `ligand_reference.py` and reused for every structure and every shell. It is a property
of the ligand and the potential, not of any fold, so it is computed on first use and then simply read;
`strain_global.py` refuses to run if the stored reference's SMILES or model does not match the run's,
since a reference from a different molecule would silently corrupt every strain rather than fail.

That reference matters more than it sounds. Earlier this term relaxed each complex's own bound pose to
get its free-ligand energy, which measured every structure against a *different* local minimum —
whichever basin that relaxation happened to fall into. The consequences were large:

* At the default `fmax 0.1` the free-ligand relaxation stops well short of a minimum. Tightening to
  0.01 moved individual strains by **+1.5 to +7.7 kcal/mol**, and the sweep only converged by 0.02
  (6.96 → 8.08 → 8.56 → 8.61 for one structure at 0.10, 0.05, 0.02, 0.01).
* With a loose step cap the ligand cannot travel; with a generous one it changes conformer. The
  largest shifts came with up to **1.0 Å of heavy-atom drift** — a rotor flipping, not a relaxation.
* The bound-state hydrogens were relaxed inside the complex, starting from a `pdbfixer` protonation
  that is **not deterministic**. Two scorings of one structure differed by 3.6 kcal/mol from that
  alone, while the heavy atoms stayed bit-identical.

A single reference removes all three: differences between structures now reflect only their bound
geometries. `strain_global.py` computes the term for every fold as one single point on the bound ligand
**exactly as the complex relaxation left it** — `structures/<name>_ligand_bound.xyz`, or equivalently
the step-0 energy of that fold's free-ligand block in the scoring log, which is the same state and so
lets folds scored before geometry saving existed be corrected without re-running anything.

**The bound state must not be re-relaxed, in whole or in part.** Strain is the energy released going
from the bound state to the free one, and the bound state is the ligand as it exists in the complex,
hydrogens included. Re-placing those hydrogens on the isolated ligand inserts an intermediate state —
bound, then bound-heavy-atoms-with-free-molecule-hydrogens, then the global minimum — and reports only
the second leg, discarding the first. It was tried, and it lowered all 24 strains by 0.2 to 6.1
kcal/mol. The hydrogen reorganisation is part of what the relaxation releases, and it is not
double-counted against the interaction energy, which is evaluated at one fixed geometry and says nothing
about the path to the free state.

**A negative strain is a failed reference, not a finding.** It means a bound pose sits below the
supposed global minimum, so the conformer search missed it; re-run with more embeddings.

Interaction energy is unaffected by any of this, which is why it was reliable throughout: it compares
the complex against its own parts at one fixed geometry, so hydrogen-placement error largely cancels.
Strain compares a bound state against a separately relaxed one, and nothing cancels.

**`hints`** — requested contacts satisfied, of those requested. Worth reporting and worth not trusting
as a measure of success: it is uncorrelated with everything else in the table.

`correlate.py` computes the relationships between all of these, joining the geometry, energy and
overlay CSVs on structure name and reporting Pearson and Spearman with n beside each coefficient:

```bash
python code/correlate.py runs/octinoxate --group '(esm1|esm2|orig)'
```

Use `--group`. These correlations invert between subsets, so a figure pooled over folds that bind by
different mechanisms is close to meaningless — across one set of four folds, interaction energy
correlated +0.93 with the ligand-to-peptide centroid separation and +0.18 with enclosure, while a
second set of four from the same run gave −0.80 for enclosure. Pooling the eight reports −0.21 and
hides both.

---

## Traps

These cost time to find. They are recorded here so they cost no one else any.

**The folding models sometimes place the mirror-image ligand, and nothing raises.** The source
`runs/octinoxate/ligand.xyz` has the **S** configuration at the 2-ethylhexyl carbon. Three of the eight
folds taken to dynamics came out **R**: `shuffle_control_esm0` from this repository's own Boltz runs, and
both of BoltzGen's, `bg33_3` and `bg33_4`. The other five kept S. The ligand SDFs carry the fold geometry
unchanged — verified to 0.0005 Å against the `.cif` — so this is what the folder placed, not something a
prep step did, and it is not specific to either tool.

Same formula, same constitution, same 20 heavy atoms, same bond graph, so every check that compares
composition passes and the inversion survives into the MD and the MM/GBSA numbers. Two things follow.
Check the **CIP label**, not the canonical SMILES and not the formula: the source SMILES carries `/C=C/`
double-bond stereo that a cif-derived SDF does not, so a whole-string comparison reports every structure
as different and hides the one difference that is real. And when a chirality-sensitive method is applied
across structures — a CNN scorer, a force field with improper torsions — group by configuration before
reading the spread.

**Open Babel's `pdbqt` → `sdf` loses valence, so RDKit cannot match the docked pose to its own
reference.** Carbons come back bracketed, `[C]`, with no implicit hydrogens, and `CalcRMS`,
`GetBestRMS` and even `AssignBondOrdersFromTemplate` all fail with "No sub-structure match found
between the reference and probe mol" — so every pose RMSD silently comes back empty. Do not try to
repair the perception. Vina preserves the atom order of the input ligand PDBQT in its output poses, and
that PDBQT was built from the reference SDF without touching coordinates, so the mapping is recoverable
exactly by coordinate identity with no chemical perception at all. `vina_redock.ref_order_poses` does
this and asserts the mapping is bijective.

**RMSD between a docked pose and its reference must not be superimposed.** They already share a frame,
so `GetBestRMS` would align away precisely the displacement being measured and report a ligand that
moved 8 Å as barely changed. Use `CalcRMS`, or an explicit automorphism minimum as here, and keep the
superposed value as a separate column measuring internal conformation.

**A docked pose packing better than the predicted one is not evidence it is better.** Vina selects
poses for close packing against a rigid receptor, so they score well on wrapping and enclosure by
construction. The comparison is circular, and only dynamics or an independent energy can break it.

**A Boltz contact constraint with `force: false` is discarded, not softened.** The featurizer skips it
outright, so a "hinted" fold with unforced constraints is an unconstrained fold. Set `force: true`.

**`force: true` alone is conditioning, not enforcement.** The constraint enters the pair
representation and the model may decline it. The guidance that acts on the structure during sampling is
gated behind `--use_potentials`, which is off by default. Both are needed.

**Boltz 2.2.1's affinity head fails against RDKit 2026.03.6.** It calls
`AllChem.Descriptors.MolWt`, and `AllChem` no longer exposes `Descriptors`, so requesting affinity
raises at parse time. `boltz_hints.py` omits the affinity block by default. Little is lost — that head
spans only 1.65 kcal/mol across 37 sequences and gave a confident score to a complex with the ligand
6.6 Å away.

**Boltz writes heavy atoms only.** Hydrogens must be added before any energy calculation, which is what
pdbfixer is for.

**pdbfixer's hydrogen placement is not deterministic.** The same structure scored twice gave −14.47 and
−13.31 kcal/mol, a 1.2 kcal/mol spread from hydrogen positions alone. Relaxing hydrogens with heavy
atoms fixed brings reproducibility to 0.08 kcal/mol. This is the default and should stay on.

**Peptide strain is meaningless in vacuum.** A relaxed peptide collapses into a compact
hydrogen-bonded ball, so `E(bound) − E(relaxed)` measures collapse rather than strain: 237.8 kcal/mol
on a 33-residue peptide, and not converged. It is available behind `--terms strain_peptide` and
disabled by default.

**Explicit-water cavity desolvation is not affordable on CPU.** A peptide plus a cavity water shell is
around 420 atoms and its relaxation was still descending 1.6 kcal/mol per step at step 102 of 200.
Truncating gives a number set by the step budget rather than by the physics. For same-length designs it
should largely cancel anyway, varying only with the distribution of polar residues.

**A positive interaction energy means the structure is broken, not weakly bound.** Check for heavy-atom
contacts under 2.6 Å. Not every short contact is a defect, though — inspect the geometry before
concluding. Of three sub-2.6 Å contacts found here, one was an n→π\* interaction approaching a carbonyl
70° out of plane, one a CH···O hydrogen bond 17° off the C–H vector, and only the third — a backbone
carbonyl oxygen 2.55 Å from an ether oxygen, with no donor available to either — was genuinely
unphysical.

**Ligand strain needs a shared reference, and the defaults do not give one.** Relaxing each complex's
own bound pose measures every structure against a different local minimum. At `fmax 0.1` the relaxation
stops 1.5–7.7 kcal/mol short; raise the step cap and the ligand changes conformer instead (up to 1.0 Å
of heavy-atom drift). Compute one reference conformer per run with `ligand_reference.py` and use
`strain_global.py`. Interaction energy is immune, because it compares the complex with its own parts at
one geometry and the hydrogen-placement error cancels.

**`pdbfixer` is non-deterministic on some platforms.** Three protonations of one CIF on the same
machine gave hydrogen positions differing by up to 2.3 Å, RMSD 0.78 Å, with heavy atoms untouched. On a
Linux container the same call was stable to 0.01 kcal/mol across runs. This is why strain must not
depend on hydrogens optimised inside the complex.

**Interaction energy grows with system size** (Spearman −0.52 against peptide length) and does not track
contact area. Compare designs of matched length, or normalise.

**FAIRChem refuses the UMA 1.0 checkpoint** on version 2.22 and above, directing you to
`fairchem-core<=2.21.0`. Use `uma-s-1p2p1`.

**The fragment library holds analogues, not residues.** Each entry is a side-chain fragment capped
where the backbone would attach, so the isoleucine entry is n-butane and the leucine entry is
isobutane. Neither carries its namesake's branch. Two entries were found to be outright wrong and have
been repaired: tyrosine was a C₇H₇ radical — 59 electrons at neutral charge, declared as a singlet —
carrying a detached water and no phenol oxygen, and glutamic was acetate, the same molecule as
aspartic. They are now p-cresol and propionate, and all ten entries are chemically distinct.
`attach_geom.py` checks electron parity so a radical declared as a closed shell is caught rather than
assumed.

---

## Limitations

**The ligand's orientation in the pocket is unconstrained, and nothing in the design objective asks
about it.** Redocking turns the ligand end for end in 18 of 54 poses while keeping the same contact
residues, and the flipped and unflipped poses score the same to within 0.12 kcal/mol. Wrapping and
enclosure saturate either way round, so both the design criterion and the fold check are blind to
which way the chromophore faces. If head placement matters for a photostabiliser, it is not currently
being designed for.

**The one redocking success is a single pose in a single structure.** `bg33_4` is the only one of eight
whose predicted pose Vina reproduces (1.30 Å), and the retention correlation that follows from it is
anchored at both ends by the same BoltzGen pair. It is a lead, not a result, and the promising screen it
suggests — redock and measure agreement, no dynamics needed — has been tested on eight structures with
release episodes at (0, 0, 0, 0, 1, 5, 5, 18), three of them tied at zero.

**Three of the eight structures carry the mirror-image ligand** — `shuffle_control_esm0`, `bg33_3` and
`bg33_4` are R at the 2-ethylhexyl carbon where the source geometry and the other five are S — so any
comparison across the eight mixes enantiomers. Within-structure numbers are unaffected, and that
includes every RMSD, wrapping and flip measurement here, since each pose is compared only against its
own reference. What it touches is the cross-structure score and ΔG columns.

**The redocking rests on a rigid receptor and one scoring function.** AutoDock Vina 1.1.2 was
parameterised on globular protein–ligand complexes, not on 12–33-mer peptides with this much exposed
surface, and no side chains were left flexible. That the docked poses have never been relaxed or run is
the reason the pose question is settled by the dynamics rather than by the docking.

**One ligand.** Every number here comes from octinoxate. Whether the spacer thresholds, the enclosure
relationship or the constraint behaviour transfer to other ligands is untested.

**Three shells.** Twenty shells were enumerated and three have been swept and folded in full, giving
thirty-six structures. That is enough to show that several single-shell conclusions did not
generalise — the effect of forcing contacts reverses between them — and enough to establish two that
did: the 21–22 spacer overshoot against a budget of 16, and the fully forced glycine design, which two
independent shells converged on at 0.96 enclosure. The remaining seventeen are untested.

**The twenty shells are less independent than they sound.** Only 18 of the 20 are distinct, and nine
of those are strict subsets of another: each single-copy shell is the first pass of its two-copy
counterpart, and the `lysine` and `aspartic` walks converge on identical selections. The family offers
nine independent arrangements plus nested subsets.

**Shell selection is biased toward charge by the same artifact that distorts the energies.** Greedy
selection ranks shells by total fragment interaction energy, and that energy strongly favours charged
fragments: across all 165 poses, charged side-chain analogues average −10.79 kcal/mol against −4.71
for aromatics and −2.94 for neutral aliphatics, on a gas-phase potential with no desolvation term.
Across the twenty shells, interaction energy per fragment correlates −0.56 with charged fraction and
+0.63 with aromatic fraction. The first two shells tested are therefore among the most charge-dominated
in the family — the first is 58% charged with a single aromatic among twelve fragments — which is
backwards for a ligand that is a methoxyphenyl chromophore on a branched alkyl chain.

**That bias has now been tested, and the criterion fails.** Shell 3 was chosen deliberately weak and
deliberately aromatic — 47% weaker on the selection criterion, 33% charged against 58% — and it folded
*better*: eight designed positions reproduced against shell 2's one, and the tightest fold in the
project. It also scored worse, 5 of 12 net favourable against 9 and 8. Total fragment interaction energy
therefore ranks shells in roughly the reverse of how well they fold, and should be replaced by a
geometric criterion or reweighted to neutralise the charge term. See
[the third worked example](RESULTS.md#third-worked-example-testing-the-selection-criterion-itself).

**Peptide strain is still not reported.** Relaxing a free peptide in vacuum measures collapse rather
than strain, so the term is disabled. Referencing it instead against the peptide's Boltz apo fold was
tried and does not work — vacuum still rewards the compact bound folds over an extended apo helix by
over 100 kcal/mol, and the apo prediction is only confident for sequences whose fold does not change on
binding. See `archive/peptide_strain/` for the attempt and the evidence. It needs implicit solvent to be
meaningful, which puts it in MM/GBSA territory rather than here.

**No experimental validation.** Nothing here has been synthesised or measured. The energies come from a
machine learning potential in the gas phase, with no solvent, no entropy and no desolvation. The six
MM/GBSA free energies are force-field estimates in implicit solvent, which is a better question than the
gas-phase one but still not a measurement.

**Six trajectories.** Everything in the dynamics section rests on six 20 ns runs: four designs and two
null controls, one ligand. That is enough to show that the static energies rank structures backwards,
that a designed arrangement beats a shuffle of itself by 9.2 kcal/mol, and that ligand release tracks
arrangement across two architectures. It is not enough to establish any metric as a predictor -- the best
one found, mean simultaneous side-chain engagement, fell from +0.90 to +0.83 when the sixth structure was
added, and got that structure's pair backwards. Nor is one design against one shuffle enough to put a
number on what the ordering search is worth in general.

**The enclosure–energy relationship is not established across sequences.** Within the glycine design
enclosure ordered the interaction energy perfectly across four structures. Across the three sequences it
did not: forcing twelve contacts raised enclosure substantially in all three but raised the interaction
energy by 30 kcal/mol, 4.9 and 0.9 respectively.

**Small numbers.** Four structures per sequence. Rank correlations on four points mean little, however
large the energy range.

**The reachability criterion is a necessary condition, not a sufficient one.** It counts degrees of
freedom against constraints and adds a steric term. A path that satisfies it may still fail for
reasons the count does not capture, and one three-fragment path that met its budget closed while
another did not.

---

## Credits and provenance

**The original repository is [lucia-71/peptidebuilder](https://github.com/lucia-71/peptidebuilder)**, and
the peptide design idea in it is the work of a masters student, developed under supervision. That
version contributed the amino-acid fragment library, pose selection, sequence assembly, the linker
estimate and the command-line pipeline — the whole peptide-specific conception — and the majority of
the repository's founding commits.

The fragment placement and scoring engine beneath it comes from
**[CafChemFragGrow](https://github.com/MauricioCafiero/CafChem)**, driven by
`notebooks/FragGrow_CafChem.ipynb` in that repository. `define_fragments`, `get_frag_coordinates`,
`get_binding_site_dims`, `grow_fragments`, `calc_frag_energy` and the viewers descend from it. That
version works with generic chemical fragments — water, cyclopropyl, acetylene, methanol, phenyl — and
the amino-acid library was written for this project.

The explicit-water machinery in `solvate.py` is adapted from the `solvation` class in
**[UMADock](https://github.com/MauricioCafiero/UMADock)**, modified here: the water count scales with
solute size rather than being fixed, the rejection distance is a settable 2.4 Å rather than 1.7, the
test covers all three atoms of a candidate water rather than only its oxygen, and the seed is settable.

External tools:

- **[FAIRChem UMA](https://github.com/FAIR-Chem/fairchem)** (Meta) — the interatomic potential behind
  every energy here
- **[Boltz-2](https://github.com/jwohlwend/boltz)** — co-folding, and the contact constraints the design
  is handed to it through
- **[ESM2](https://github.com/facebookresearch/esm)** (Meta) — linker substitution, via the
  `facebook/esm2_t12_35M_UR50D` checkpoint
- **[ASE](https://wiki.fysik.dtu.dk/ase/)** — geometry optimisation and structure handling
- **[RDKit](https://www.rdkit.org/)** — SMILES handling, conformer generation, bond perception
- **[pdbfixer](https://github.com/openmm/pdbfixer)** and **[OpenMM](https://openmm.org/)** — adding
  hydrogens to folded structures
- **[py3Dmol](https://github.com/3dmol/3Dmol.js)** — in-notebook visualisation

---

## Citation

If you use this work, please cite this repository and the original
[lucia-71/peptidebuilder](https://github.com/lucia-71/peptidebuilder), along with the UMA, Boltz-2 and
ESM2 papers for the models it depends on.

---

## Licence

MIT. See [LICENSE](LICENSE).
