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

---

## Contents

- [The problem this solves](#the-problem-this-solves)
- [End-to-end workflow](#end-to-end-workflow)
- [Installation](#installation)
- [Repository layout](#repository-layout)
- [What each file does](#what-each-file-does)
- [Where output goes](#where-output-goes)
- [Worked example: octinoxate](#worked-example-octinoxate)
- [Second worked example: a different shell around the same ligand](#second-worked-example-a-different-shell-around-the-same-ligand)
- [Third worked example: testing the selection criterion itself](#third-worked-example-testing-the-selection-criterion-itself)
- [What happens after this pipeline: dynamics and MM/GBSA](#what-happens-after-this-pipeline-dynamics-and-mmgbsa)
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

### Boltz-2, the one external dependency

Boltz-2 does the co-folding. It is heavy and GPU-bound, so it may live in its own environment rather
than this one — but nothing here assumes where. `code/boltz_env.py` resolves how to run it, first
match winning:

1. `--boltz-cmd`, or `$PEPTIDEBUILDER_BOLTZ_CMD` — a complete command, used as given
2. `--boltz-venv`, or `$PEPTIDEBUILDER_BOLTZ_VENV` — a virtualenv holding Boltz
3. `boltz` importable in this environment — the plain `pip install boltz` case
4. `boltz` on `$PATH`
5. an MPS wrapper beside the venv, if that layout happens to exist

If none match, the error lists these options rather than failing obscurely.

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
| `cif_to_md.py` | splits a Boltz complex into the protein PDB and ligand SDF the OpenMM MD pipeline wants, assigning ligand bond orders from the run's SMILES. See [dynamics and MM/GBSA](#what-happens-after-this-pipeline-dynamics-and-mmgbsa) |
| `md_contacts.py` | contacts, centroid separation and peptide Rg in windows across a trajectory: whether the designed pose survived, which one ΔG cannot tell you. Heavy atoms only, so the counts stay comparable with `fold_check*.csv` |
| `md_frames.py` | representative end-of-run frames as PDBs for figures — a ten-model ensemble superposed on the peptide, and the single frame nearest the window mean. `--window-ns` picks a different interval where the end of a run is unrepresentative |

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

---

## Worked example: octinoxate

Octinoxate — 2-ethylhexyl 4-methoxycinnamate — is a UV filter and the ligand this has been developed
against. The fragment library was screened against it, giving 165 accepted poses across ten side-chain
types.

Pose selection is greedy and order-dependent, so rotating which fragment picks first produces a family
of different shells. Ten starting fragments at one and two copies each gives twenty shells, from eight
to thirteen fragments, drawing on 48 distinct poses. The richest holds twelve fragments with a total
fragment interaction energy of −113.12 kcal/mol.

Sweeping that shell found **all 132 ordered pairs connectable** at some spacer count, and the ordering
search returned a path through all twelve poses:

```
RGGDGGKGGGGLGGGGKGIGEGWGGDGSGGEGS      33 residues, 21 spacers against a budget of 16
```

Two ESM2-substituted variants were generated from it, replacing the glycine spacers with probable
residues while leaving the designed positions untouched:

```
esm1   RLVDALKNKTPLSAAEKSILEDWLADKSKSESS
esm2   RKKDAKKRIRDLESIQKLILEKWGVDSSSNEDS
```

Each of the three sequences was co-folded four ways: unconstrained, and with the four, eight and twelve
tightest designed contacts supplied as constraints. Twelve structures in total.

### Results

| structure | enclosed | wrapped | engaged | centroid sep | interaction | ligand strain | sum | hints |
|---|---|---|---|---|---|---|---|---|
| `orig_control` | 0.495 | 0.75 | 15/20 | 20.9 Å | −15.83 | 14.66 | −1.17 | — |
| `orig_f4` | 0.595 | 0.85 | 17/20 | 7.1 Å | −26.32 | 11.64 | −14.68 | 3/4 |
| `orig_f8` | 0.515 | 0.65 | 13/20 | 8.4 Å | −22.30 | 42.06 | +19.76 | 1/8 |
| **`orig_f12`** | **0.96** | **1.00** | 20/20 | 4.2 Å | **−45.86** | 12.80 | **−33.07** | 0/12 |
| `esm1_control` | 0.455 | 0.60 | 12/20 | 8.4 Å | −12.22 | 7.20 | −5.02 | — |
| `esm1_f4` | 0.605 | 0.65 | 13/20 | 6.6 Å | −20.89 | 27.16 | **+6.27** | 0/4 |
| `esm1_f8` | 0.75 | 0.90 | 18/20 | 4.1 Å | −29.10 | 36.91 | **+7.81** | 1/8 |
| `esm1_f12` | 0.76 | 0.80 | 16/20 | 3.4 Å | −17.11 | 11.26 | −5.85 | 1/12 |
| `esm2_control` | 0.53 | 0.40 | 8/20 | 7.7 Å | −20.27 | 19.14 | −1.14 | — |
| `esm2_f4` | 0.46 | 0.50 | 10/20 | 11.0 Å | −39.90 | 25.91 | −13.99 | 0/4 |
| `esm2_f8` | 0.625 | 0.75 | 15/20 | 10.4 Å | −29.66 | 15.61 | −14.05 | 0/8 |
| `esm2_f12` | 0.905 | 1.00 | 20/20 | 7.1 Å | −21.14 | 9.86 | −11.28 | 0/12 |

Energies in kcal/mol. All three sequences are 33 residues.
\* `orig_f8`'s free-ligand relaxation did not converge within its step budget when this was first
scored. That no longer affects the number: ligand strain is now measured against a single reference
conformer rather than a per-structure relaxation, so no per-row convergence caveat applies. See
[what the metrics mean](#what-the-metrics-mean).

### What the example shows

**Forcing the full set of contacts produces the design that was wanted.** `orig_f12` encloses 96% of
the directions out of the ligand, contacts every one of its heavy atoms, binds at −45.86 kcal/mol, and
distorts the ligand *less* than the unconstrained fold of the same sequence. It is the only structure
in the set where the ligand is genuinely inside a shell rather than pressed against a surface.

**The constraints work, but not by being satisfied.** `orig_f12` honoured none of its twelve requested
contacts, and `esm1_f8` honoured one of eight. Judged on whether the requested residue landed against
the requested ligand atom, the experiment looks like a failure. Judged on whether the ligand ended up
enclosed, it is the most effective intervention available. Every f12 fold is the most enclosed member
of its own set.

**Partial constraint sets damage the ligand.** The four- and eight-contact folds repeatedly gain
interaction energy by bending the ligand rather than by surrounding it: 42.06, 36.91 and 27.16 kcal/mol
of ligand strain, against 9.86 to 12.80 for every twelve-contact fold. **Three of them end up net
unfavourable** — `orig_f8` at +19.76, `esm1_f8` at +7.81 and `esm1_f4` at +6.27. A peptide can satisfy
a few contacts by pulling the ligand toward whichever residues are nearby, but satisfying twelve at
once requires surrounding it, and a ligand can only be surrounded in a conformation it can actually
adopt.

**The ESM2 substitution locks the fold.** All four esm1 structures are the same shape to within 0.3–0.8 Å
CA RMSD, whether they were given zero, four, eight or twelve constraints. Backbone dihedrals explain it:
the glycine design is **0% helical**, while both substituted sequences are **77–81% helical** and their
helical content barely changes under constraint. ESM2 fills the linkers with alanine, leucine, lysine
and glutamate, and the rigid helix that results cannot be reshaped. Forcing twelve contacts gains the
glycine design 30 kcal/mol of interaction energy and the esm1 sequence 4.9.

**No common binding mode exists across the twelve.** Superposed on the ligand, the folds differ by a
median 18.1 Å CA RMSD. Each one wraps the ligand its own way.

**Geometry predicts energy only within a binding mechanism.** Ranking each sequence's four folds by
enclosure and again by interaction energy gives the identical order for the glycine design, one
swap for esm1, and a near-reversal for esm2. The reason is visible in the contacts. `orig_f12` binds
through 47 contacts, none of them charged, dominated by backbone carbonyl oxygens approaching the
ligand's ester carbons — 2.44 and 2.57 Å at `GLY14`, in the plane of the carbonyl (O···C=O 70°, not
the ~107° of an n→π* approach), so a dipolar contact that needs the whole wrap to accumulate.
`esm2_f4` reaches −39.90 kcal/mol from **eight** contacts, five of them charged: an N-terminal
arginine canted 30° against the methoxyphenyl ring at 5.1 Å with a glutamate carboxylate at the far
end. A terminal side chain reaching in from outside does not need to enclose anything, so enclosure
stops predicting anything. Forcing more contacts then makes `esm2` monotonically *worse* — −39.90 at
0.50 wrapped, −29.66 at 0.75, −21.14 at 1.00 — because the wrap is incompatible with the
arrangement that was doing the binding.

Two cautions follow. Charged-residue binding is real signal and worth having, but a solvent-exposed
guanidinium pays no desolvation penalty in a gas-phase interaction energy, so it is the mode most
likely to weaken in water; `orig_f12`'s buried, chargeless contact set is the one expected to
survive, and it is the strongest here anyway. And a single geometric criterion is not a sufficient
objective: wrapping and enclosure identify the best structure in this set, but they do so because
the best structure happens to bind by the mechanism they measure.

### Figures

All twelve structures, peptide in green, ligand in red. Each row is one sequence and each column one
rung of the constraint ladder, so reading left to right shows what forcing more contacts does to that
sequence. Figures are `wrapped` / `interaction` in kcal/mol.

**The glycine design.** Wrapping and binding strengthen together, ending in the best structure in the
set.

| unconstrained | 4 contacts | 8 contacts | 12 contacts |
|---|---|---|---|
| ![orig_control](runs/octinoxate/figures/orig_control.png) | ![orig_f4](runs/octinoxate/figures/orig_f4.png) | ![orig_f8](runs/octinoxate/figures/orig_f8.png) | ![orig_f12](runs/octinoxate/figures/orig_f12.png) |
| 0.75 / −15.83 | 0.85 / −26.32 | 0.65 / −22.30 | **1.00 / −45.86** |
| ligand against the outside, centroids 20.9 Å apart | | ligand bent to buy contacts: 42.06 kcal/mol of strain, net unfavourable | every ligand atom contacted, 47 pairs, none charged |

**The first ESM2 variant.** The helix barely changes shape across the whole ladder.

| unconstrained | 4 contacts | 8 contacts | 12 contacts |
|---|---|---|---|
| ![esm1_control](runs/octinoxate/figures/esm1_control.png) | ![esm1_f4](runs/octinoxate/figures/esm1_f4.png) | ![esm1_f8](runs/octinoxate/figures/esm1_f8.png) | ![esm1_f12](runs/octinoxate/figures/esm1_f12.png) |
| 0.60 / −12.22 | 0.65 / −20.89 | 0.90 / −29.10 | 0.80 / −17.11 |
| | | | 77–81% helical throughout, 0.3–0.8 Å CA RMSD between all four |

**The second ESM2 variant.** The counter-example: binding is strongest where wrapping is nearly
weakest, and gets worse as the wrap tightens.

| unconstrained | 4 contacts | 8 contacts | 12 contacts |
|---|---|---|---|
| ![esm2_control](runs/octinoxate/figures/esm2_control.png) | ![esm2_f4](runs/octinoxate/figures/esm2_f4.png) | ![esm2_f8](runs/octinoxate/figures/esm2_f8.png) | ![esm2_f12](runs/octinoxate/figures/esm2_f12.png) |
| 0.40 / −20.27 | **0.50 / −39.90** | 0.75 / −29.66 | 1.00 / −21.14 |
| | `ARG1` on the ring, 8 contacts, ligand largely outside | | fully wrapped and half the binding energy of `orig_f12` |

The pair worth looking at hardest is `orig_control` against `orig_f12` — same sequence, same ligand,
the only difference being twelve forced contacts, and the ligand goes from lying against the outside
of the peptide to sitting inside it.

`load_folds.pml` loads all twelve superposed on the ligand, with the designed shell in marine.
`style.pml` applies the representation used here to anything already loaded — sticks, black
background, green peptide, red ligand, residue labels on CA. Three further scripts —
`close_contact.pml`, `contact_orig_esm1.pml`, `contact_orig_f8.pml` — focus on individual short
contacts with their distances labelled.

---

## Second worked example: a different shell around the same ligand

Every number above comes from one designed arrangement. Pose selection is greedy and order-dependent,
so rotating which fragment picks first produces a family of shells, and the obvious question is
whether any of it generalises. This is the second-richest shell, at −109.15 kcal/mol against the
first's −113.12, sharing only five of its twelve poses:

```
phenylalanine:2, arginine:3, lysine:10, aspartic:21, glutamic:15, isoleucine:8,
leucine:12, serine:32, tryptophan:5, glutamic:10, leucine:13, serine:2
```

Its sweep found **all 132 ordered pairs connectable**, as the first shell's did, and the ordering
search returned a path through all twelve poses. The design is independently identical to the first on
every structural count:

| | shell 1 | shell 2 |
|---|---|---|
| sequence | `RGGDGGKGGGGLGGGGKGIGEGWGGDGSGGEGS` | `RGEGGEGGKGGFGGDGGIGLGGSGWGGGGLGGS` |
| residues | 33 | 33 |
| spacers used, against a budget of 16 | 21 | 21 |
| poses in the path | 12 / 12 | 12 / 12 |
| worst closure | 0.196 Å | 0.120 Å |
| total fragment interaction energy | −113.12 | −109.15 |

Two shells sharing five poses, ordered by independent exact searches, arriving at the same length and
the same overshoot against the same budget. The overshoot is therefore a property of twelve-fragment
shells around this ligand rather than an accident of one arrangement.

Two ESM2 linker variants were generated as before, and each of the three sequences co-folded four ways.

```
esm1   RSELAEAAKRGFLTDGGIDLVSSGWLTLRLAPS      net charge  0
esm2   RGERKEKSKLIFIIDSKIFLSFSIWRFLILLRS      net charge +5
```

### Results

| structure | enclosed | wrapped | engaged | centroid sep | interaction | ligand strain | sum | hints | designed positions |
|---|---|---|---|---|---|---|---|---|---|
| **`s2_orig_control`** | 0.60 | **1.00** | 20/20 | 7.7 Å | −25.18 | 9.23 | **−15.95** | — | 0/12 |
| `s2_orig_f4` | 0.625 | 0.85 | 17/20 | 2.9 Å | −17.48 | 15.87 | −1.61 | 0/4 | 0/12 |
| `s2_orig_f8` | 0.43 | 0.85 | 17/20 | 9.6 Å | −11.07 | 36.80 | **+25.73** | 1/8 | 0/12 |
| `s2_orig_f12` | 0.70 | 0.85 | 17/20 | 7.2 Å | −17.56 | 16.07 | −1.49 | 1/12 | 0/12 |
| `s2_esm1_control` | 0.375 | 0.55 | 11/20 | 11.2 Å | −9.70 | 11.64 | **+1.93** | — | 0/12 |
| `s2_esm1_f4` | 0.64 | 0.80 | 16/20 | 6.8 Å | −29.55 | 21.43 | −8.12 | 1/4 | 0/12 |
| `s2_esm1_f8` | 0.67 | 0.80 | 16/20 | 5.5 Å | −19.91 | 33.74 | **+13.83** | 1/8 | 0/12 |
| **`s2_esm1_f12`** | 0.675 | 0.90 | 18/20 | 5.9 Å | **−42.39** | 23.37 | **−19.03** | **3/12** | **1/12** |
| `s2_esm2_control` | 0.51 | 0.85 | 17/20 | 6.7 Å | *−74.89* | 9.05 | *−65.84* | — | 0/12 |
| `s2_esm2_f4` | 0.50 | 0.55 | 11/20 | 16.0 Å | −7.53 | 29.32 | **+21.80** | 0/4 | 0/12 |
| **`s2_esm2_f8`** | **0.86** | **0.95** | 19/20 | 9.5 Å | −28.29 | 23.30 | −4.99 | 0/8 | 0/12 |
| `s2_esm2_f12` | 0.61 | 0.75 | 15/20 | 5.9 Å | *−72.72* | 12.78 | *−59.94* | 2/12 | 0/12 |

Energies in kcal/mol; all three sequences are 33 residues. **The `esm2` energies are italicised
because they are artifacts, not measurements** — see below.

### What replicates, and what does not

**The eight-contact fold is usually the damaged one, with a consistent exception.** `s2_orig_f8` and
`s2_esm1_f8` carry the highest ligand strain of their ladders, at 36.80 and 33.74 kcal/mol, and both
end net unfavourable at +25.73 and +13.83. Across the six ladders in both shells, f8 is the
highest-strain rung in **four** — both `orig` and both `esm1` — at 33.7 to 42.1 kcal/mol.

**The two exceptions are both `esm2`, in both shells**, where f4 carries the higher strain instead:
25.91 against f8's 15.61 in shell 1, and 29.32 against 23.30 in shell 2. That is the same sequence
family that breaks the enclosure ladder and whose energies are charge artifacts, so it is the branch
that behaves differently on every measure rather than a random exception.

**Shell 3 replicates both halves of this**, which makes it the most reproducible energetic result here:
f8 is again the highest-strain rung for `orig` (55.77) and for `esm1` (61.67) — the two worst strains in
the project — and `esm2` is again the exception, with f4 at 15.16 just above f8's 14.70. Nine ladders
across three shells, and the rule and its exception hold in all nine.

A peptide can satisfy a few contacts by pulling the ligand toward whichever residues are nearby; eight
is apparently enough to demand serious distortion and too few to require the wrap that would relieve
it — except where the substitution has already locked the fold into a shape that four contacts fight
harder than eight.

**Forced contacts are a rescue mechanism, not an improvement.** This is the sharpest disagreement
between the shells, and it resolves rather than contradicts the first example:

| glycine design | control | twelve forced contacts | effect |
|---|---|---|---|
| shell 1 | −1.17 | **−33.07** | forcing gains 32 kcal/mol |
| shell 2 | **−15.95** | −1.49 | forcing loses 14 kcal/mol |
| shell 3 | +15.03 | +18.92 | forcing costs 4 kcal/mol; neither is favourable |

Shell 1's unconstrained control left the ligand pressed against the outside — 0.75 wrapped, centroid
20.9 Å — so the constraints had everything to gain. Shell 2's control already contacted all twenty
ligand heavy atoms at 7.7 Å against a radius of gyration of 13.3 Å, so they could only disturb it,
and the geometry says so: wrapping fell from 1.00 to 0.85. **Forcing the full contact set produces
the designed arrangement when the unconstrained fold has failed, and degrades it when the
unconstrained fold has already succeeded.** Read the control's wrapping first to know which case you
are in.

Shell 3 is the third case, and the one the rule does not cover: its control wraps 0.80 with the ligand
9.5 Å out, neither failed nor succeeded, and forcing moves the energy barely at all — from +15.03 to
+18.92, both unfavourable. What forcing *does* achieve there is entirely geometric: enclosure goes from
0.735 to 0.965 and engagement from 16/20 to 20/20, at 30 kcal/mol of extra ligand strain that the energy
then charges for. **So the rule holds for the energy and breaks for the geometry** — forcing reliably
improves the cavity and unreliably improves the score, which is the same disagreement that runs through
the whole of shell 3.

**A shell can succeed at the objective without realising the design.** `s2_orig_control` wraps every
ligand heavy atom while reproducing none of the twelve designed side-chain positions. The first shell
could not separate these, because its control failed at both at once.

**Designed positions are still essentially never reproduced: 1 of 144** across all twelve folds, at
ligand superposition RMSD 0.06–0.85 Å, so the comparison is sound. But the distances say something the
count hides. Forcing contacts moves side chains markedly closer to where the design asked for them —
the median distance to a designed position falls from 12–14.5 Å unconstrained to 7.6–9.8 Å forced, in
all three ladders — and the median distance to the nearest side chain **of any type** falls to 3.0–4.9
Å, reaching 3.02 Å for `s2_esm2_f8`. Something is arriving at the designed position; it is the wrong
residue. That relocates the failure from reachability, which the sweep shows is satisfiable for every
pair, to the sequence assembly that decides which residue lands where.

**The claim that every twelve-contact fold is the most enclosed of its set does not survive.** It
holds for `orig` at 0.70 and barely for `esm1` at 0.675 against f8's 0.67, but `esm2`'s best-enclosed
fold is **f8 at 0.86**, with f12 down at 0.61. `s2_esm2_f8` is the best-enclosed structure of the
whole set — 0.95 wrapped, 19 of 20 atoms engaged, 31 contacts — and it scores poorly. Enclosure and
energy disagree again, and this time enclosure is right about the structure and wrong about the
number.

### The charge artifact, quantified

The linker fill took the design from net charge −1 to **+5**, and `s2_esm2_control` posts −74.89
kcal/mol, nominally the strongest binding anywhere in this project. It should not be read as binding
at all. Across `esm2`'s four folds the interaction energy correlates **+0.93 with the distance between
the ligand and peptide centroids** and only **+0.18 with enclosure**, with the wrong sign:

| `esm2` fold | centroid sep | interaction | enclosed | wrapped | contacts |
|---|---|---|---|---|---|
| f12 | 5.9 Å | −72.72 | 0.61 | 0.75 | 17 |
| control | 6.7 Å | −74.89 | 0.51 | 0.85 | 18 |
| f8 | 9.5 Å | −28.29 | **0.86** | **0.95** | **31** |
| f4 | 16.0 Å | −7.53 | 0.50 | 0.55 | 20 |

The best-wrapped structure in the set scores 46 kcal/mol worse than a control that wraps less and
encloses far less, purely because the control's centroid sits 2.8 Å closer. That is long-range
electrostatics on a highly charged peptide, not an interface.

The control for this reading is the glycine design at net charge −1, where the same relationship is
absent: its f4 and f12 folds give −17.48 and −17.56 kcal/mol despite centroid separations of 2.9 and
7.2 Å. Grouped by sequence, the correlation between interaction energy and centroid separation runs
**+0.23** at charge −1, **+0.70** at charge 0 and **+0.93** at charge +5, while the correlation with
enclosure runs −0.56, −0.76 and +0.18. Pooling all twelve reports −0.14 for enclosure and +0.42 for
separation, hiding both.

So the distance-dominated energy is a property of the charged variant rather than of the scoring
method, which means the `orig` and `esm1` numbers can be read as chemistry and `esm2`'s cannot.
**Constrain the residue set when filling linkers, and report net charge beside any interaction
energy.** `correlate.py --group` computes these.

### Two binding modes, and what actually predicts them

Inspecting the twelve renders by eye splits them into two mechanisms, and neither the enclosure nor
the wrapping column reliably tells them apart.

**Encapsulation** — the ligand inside the peptide: `s2_esm1_f4`, `s2_esm1_f8`, `s2_esm1_f12` and
`s2_esm2_f8`.

**A groove on an elongated structure** — the ligand held against a surface channel:
`s2_orig_control`, `s2_orig_f4` and `s2_esm2_f12`.

**Radius of gyration separates them cleanly and the purpose-built metric does not.** The four
encapsulating folds span Rg 9.3–9.8 Å; the three groove binders span 11.8–14.9 Å. Nothing falls in
between. Enclosure separates the same two groups by 0.015 — 0.625 against 0.64 — which is far too
narrow to act on, and **wrapping inverts the verdict outright**: the highest wrapped value in the whole
set, 1.00, belongs to `s2_orig_control`, a groove binder. Wrapping counts ligand atoms in contact, and
a groove can contact all of them. When judging whether a design achieved encapsulation, read
`peptide_rg` first and treat wrapping as a contact census rather than a verdict.

**Helicity explains how each mode arises, and it is not simply "helical means elongated".** The groove
binders span 3% to 100% helical, so there are two unrelated routes to an extended structure:

| sequence | helical | Rg | what it is |
|---|---|---|---|
| `orig` (glycine design) | 3–6% | 7.7–13.3 Å | never helical; extended because a 64%-glycine chain has no fold |
| `esm1` | **42% in all four folds** | 9.3–10.3 Å | locked, and the locked shape is compact |
| `esm2` | 87–100% | 9.8–15.0 Å | a helical rod, except where the helix broke |

`s2_esm2_f12` is a groove binder because it is a rigid 100% helix at Rg 14.9 Å and the ligand can only
lie along it. `s2_orig_control` is a groove binder for the opposite reason — 3% helical, extended
because nothing folds it. And the one `esm2` fold that encapsulates, `s2_esm2_f8`, is the one where the
helix **broke**: 87% helical at Rg 9.8 Å, against 97–100% and ~15 Å for its three siblings.
Encapsulation by a substituted sequence required disrupting the structure the substitution created.

This is the same lock-in the first shell showed, where all four `esm1` structures held the same shape
regardless of constraint. Here `esm1` is locked at 42% helicity across all four folds and that
conformation happens to be compact, which is why every forced `esm1` fold encapsulates. The lock is
not itself good or bad; it decides the mode before any constraint is applied, and whether that helps
is luck.

**The glycine design's extended folds look like β sheets and are not.** `s2_orig_control` is **81%
β-basin**, the most extended backbone in the set, which is why it reads as pleated on inspection. But it
carries only **6 non-local backbone hydrogen bonds**, where a sheet needs a ladder of them between
paired strands. It is an extended, unpaired chain — the expected outcome for a sequence that is 64%
glycine, glycine being both a poor β-former and the residue with the most backbone freedom to give
away. `s2_orig_f4` is the same at 36% β and 6 non-local bonds.

The exception is informative: `s2_orig_f12` has the **most non-local hydrogen bonds of any fold in the
set**, 11, at the most compact Rg of 7.7 Å. Forcing twelve contacts is the only intervention that gave
the glycine design any tertiary structure at all, even though it cost binding energy. The `esm2` folds
sit at the opposite extreme — 0% β and zero non-local bonds, with 34 to 43 helical ones, a helix bonded
only to itself.

`check_fold.py` reports `helical_fraction`, `beta_fraction`, `helical_hbonds` and `nonlocal_hbonds`, so
all of this is measured rather than eyeballed. On the first shell's structures the helicity gives 81%
for `esm1_f12` and 6% for `orig_f12`, reproducing the figures quoted in that example.

### Figures

All twelve structures, peptide in green, ligand in red, laid out as before: one row per sequence, one
column per rung of the constraint ladder. Captions are `wrapped` / `interaction` in kcal/mol.

**The glycine design.** The unconstrained fold is the best of the four, and forcing contacts makes it
worse — the reverse of the first shell.

| unconstrained | 4 contacts | 8 contacts | 12 contacts |
|---|---|---|---|
| ![s2_orig_control](runs/octinoxate/figures_shell2/s2_orig_control.png) | ![s2_orig_f4](runs/octinoxate/figures_shell2/s2_orig_f4.png) | ![s2_orig_f8](runs/octinoxate/figures_shell2/s2_orig_f8.png) | ![s2_orig_f12](runs/octinoxate/figures_shell2/s2_orig_f12.png) |
| **1.00 / −25.18** | 0.85 / −17.48 | 0.85 / −11.07 | 0.85 / −17.56 |
| **groove**, Rg 13.3 Å, 3% helical — every ligand atom contacted, none of it enclosed | **groove**, Rg 11.8 Å | 36.80 kcal/mol of ligand strain, net unfavourable | most compact of all at Rg 7.7 Å, 36 contacts, still worse than the control |

**The first ESM2 variant.** The only ladder where geometry and energy agree throughout, and the only
fold in 144 to reproduce a designed position.

| unconstrained | 4 contacts | 8 contacts | 12 contacts |
|---|---|---|---|
| ![s2_esm1_control](runs/octinoxate/figures_shell2/s2_esm1_control.png) | ![s2_esm1_f4](runs/octinoxate/figures_shell2/s2_esm1_f4.png) | ![s2_esm1_f8](runs/octinoxate/figures_shell2/s2_esm1_f8.png) | ![s2_esm1_f12](runs/octinoxate/figures_shell2/s2_esm1_f12.png) |
| 0.55 / −9.70 | 0.80 / −29.55 | 0.80 / −19.91 | **0.90 / −42.39** |
| 42% helical, as all four are | **encapsulates**, Rg 9.6 Å | **encapsulates**, but 33.74 kcal/mol of strain, net unfavourable | **encapsulates**, Rg 9.3 Å — 3/12 hints honoured and the only fold in 144 to reproduce a designed position |

**The second ESM2 variant.** Net charge +5. The energies track how close the ligand sits, not how
well it is held, so read the geometry columns and ignore the numbers.

| unconstrained | 4 contacts | 8 contacts | 12 contacts |
|---|---|---|---|
| ![s2_esm2_control](runs/octinoxate/figures_shell2/s2_esm2_control.png) | ![s2_esm2_f4](runs/octinoxate/figures_shell2/s2_esm2_f4.png) | ![s2_esm2_f8](runs/octinoxate/figures_shell2/s2_esm2_f8.png) | ![s2_esm2_f12](runs/octinoxate/figures_shell2/s2_esm2_f12.png) |
| 0.85 / *−74.89* | 0.55 / −7.53 | **0.95 / −28.29** | 0.75 / *−72.72* |
| 100% helical rod, Rg 15.0 Å | ligand thrown 16.0 Å from the centroid | **encapsulates** — the one fold whose helix broke, 87% at Rg 9.8 Å; best geometry in the set, mediocre score | **groove** along a 100% helix at Rg 14.9 Å; fully forced, and the score barely moves from the control |

`load_folds.pml` in `figures_shell2/` loads all twelve superposed on the ligand with the designed
shell in marine, and `style.pml` beside it applies the representation used here.

---

## Third worked example: testing the selection criterion itself

The first two shells were chosen the same way — by total fragment interaction energy, −113.12 and
−109.15 kcal/mol, the two richest arrangements the sweep found. They also turned out to be the two most
charged, and that is not a coincidence: across all fragment poses, charged ones average **−10.79
kcal/mol against −4.71 aromatic and −2.94 aliphatic**, and across 36 folds the interaction energy
correlates **−0.56** with charged-residue fraction. The criterion that picks shells is substantially
measuring charge.

Which makes one experiment obvious. Take a shell that is much *weaker* on that criterion and much less
charged, and see what the folds do. If the criterion is sound, the geometry should get worse. If the
criterion is mostly counting charge, it might not.

Shell 3 is seeded on **tryptophan pose 2** and comes in at −60.40 kcal/mol — **47% weaker than shell
1** — with its composition shifted from charge to aromatics:

| | shell 1 | shell 2 | shell 3 |
|---|---|---|---|
| sequence | `RGGDGGKGGGGLGGGGKGIGEGWGGDGSGGEGS` | `RGEGGEGGKGGFGGDGGIGLGGSGWGGGGLGGS` | `YGGDGFGKGGGLGGGKGGGGLGGWGGEGSGGWGS` |
| residues | 33 | 33 | 34 |
| spacers used, against a budget of 16 | 21 | 21 | 22 |
| poses in the path | 12 / 12 | 12 / 12 | 12 / 12 |
| worst closure | 0.196 Å | 0.120 Å | 0.188 Å |
| **charged poses** | **7 (58%)** | 5 (42%) | **4 (33%)** |
| **aromatic poses** | 1 (8%) | 2 (17%) | **4 (33%)** |
| **total fragment interaction energy** | **−113.12** | −109.15 | **−60.40** |

The structural counts hold for a third time — twelve poses ordered by an independent exact search,
arriving at 33–34 residues and 21–22 spacers against the same budget of 16. Three shells now overshoot
that budget by the same margin, so the overshoot is a property of twelve-fragment shells around this
ligand, not of any one arrangement.

```
esm1   YLSDAFGKSSKLADLKDGSTLPSWASETSSLWKS      net charge  0
esm2   YINDSFIKINYLINKKKFKILKIWLLEISFLWSS      net charge +4
```

### Results

| structure | enclosed | wrapped | engaged | sep | Rg | helix | β | interaction | ligand strain | sum | hints | designed positions |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| `s3_orig_control` | 0.735 | 0.80 | 16/20 | 9.5 Å | 7.8 | 22% | 22% | +0.92 | 14.12 | **+15.03** | — | 2/12 |
| `s3_orig_f4` | 0.655 | 0.80 | 16/20 | 7.2 Å | 7.8 | 3% | 31% | −28.45 | 39.45 | **+11.01** | 1/4 | 1/12 |
| `s3_orig_f8` | 0.880 | 0.85 | 17/20 | 7.9 Å | 7.8 | 19% | 16% | −5.62 | 55.77 | **+50.15** | 1/8 | 0/12 |
| **`s3_orig_f12`** | **0.965** | **1.00** | **20/20** | **3.6 Å** | 8.1 | 19% | 16% | −25.13 | 44.05 | **+18.92** | 2/12 | 0/12 |
| `s3_esm1_control` | 0.755 | 0.85 | 17/20 | 5.0 Å | 9.4 | 59% | 6% | −26.88 | 13.68 | −13.21 | — | 1/12 |
| `s3_esm1_f4` | 0.785 | 0.70 | 14/20 | 6.5 Å | 9.5 | 72% | 9% | −19.30 | 51.17 | **+31.87** | 1/4 | 1/12 |
| `s3_esm1_f8` | 0.850 | **0.95** | 19/20 | 5.3 Å | 9.4 | 69% | 16% | +2.00 | 61.67 | **+63.67** | 1/8 | 1/12 |
| **`s3_esm1_f12`** | 0.880 | 0.90 | 18/20 | 4.3 Å | 9.3 | 72% | 3% | −25.14 | 11.65 | **−13.50** | 1/12 | 1/12 |
| `s3_esm2_control` | 0.360 | 0.45 | 9/20 | 13.6 Å | 15.0 | 100% | 0% | −10.78 | 13.13 | +2.35 | — | 1/12 |
| **`s3_esm2_f4`** | 0.795 | 0.90 | 18/20 | 7.0 Å | 10.3 | 91% | 0% | −34.59 | 15.16 | **−19.43** | 0/4 | 0/12 |
| `s3_esm2_f8` | 0.795 | 0.80 | 16/20 | 6.3 Å | 10.2 | 91% | 0% | −21.56 | 14.70 | −6.85 | 0/8 | 0/12 |
| `s3_esm2_f12` | 0.640 | 0.75 | 15/20 | 6.3 Å | 9.7 | 88% | 0% | −25.05 | 10.83 | −14.22 | 0/12 | 0/12 |

Energies in kcal/mol. Ligand strain is measured against the single global-minimum reference at
−557162.371 kcal/mol-equiv, as everywhere else in this README.

### The answer: geometry improved, energy got worse

**A 47% weaker shell folded better, not worse.**

| | shell 1 | shell 2 | shell 3 |
|---|---|---|---|
| designed positions reproduced, out of 144 | not measured | **1** | **8** |
| best fold | `orig_f12`, 0.960 enclosed | `s2_esm2_f8`, 0.860 | `s3_orig_f12`, 0.965 |
| best enclosed fraction | 0.960 | 0.860 | **0.965** |
| net favourable folds | **9 / 12** | 8 / 12 | **5 / 12** |
| best sum | −33.07 | −65.84 *(artifact)* | −19.43 |
| median ligand strain | 15.61 | 21.43 | 15.16 |

Shell 3 reproduces **eight** designed side-chain positions across its twelve folds where shell 2
reproduced **one**. (Shell 1's folds were never measured this way — `overlay.csv` holds its sequence-level
comparisons, not its twelve folds, so the column is blank rather than zero.)

More striking is what its best fold is. `s3_orig_f12` and shell 1's `orig_f12` are near-twins, arrived at
from shells 47% apart in fragment interaction energy:

| | `orig_f12` (shell 1) | `s3_orig_f12` (shell 3) |
|---|---|---|
| enclosed | 0.960 | 0.965 |
| wrapped | 1.00 | 1.00 |
| engaged | 20/20 | 20/20 |
| centroid separation | 4.2 Å | **3.6 Å** |
| contacts under cutoff | 47 | **67** |
| mean ligand distance | 3.49 Å | **3.39 Å** |
| closest approach | **2.44 Å** (a clash) | 2.65 Å |
| peptide Rg | 8.0 Å | 8.1 Å |

Two independent sweeps, two different seeds, two shells with only a fraction of their poses in common,
both producing a glycine-spaced design whose fully forced fold closes around the ligand at 0.96 enclosure
with every ligand atom engaged. Shell 3's is the tighter of the two on every count and it does not carry
shell 1's 2.44 Å clash. **That convergence is the most reproducible result in the project**, and it says
the fully forced glycine design is a robust target rather than one lucky fold.

By energy it is the worst: 5 of 12 net favourable against 9 and 8, and its best sum is −19.43 where
shell 1 managed −33.07. **The two families of metric disagree, and they disagree systematically** — the
shell that was selected for having less charge scores worse on a criterion that rewards charge, while
folding into better cavities. Shell 2's headline −65.84 is itself an artifact of net charge +5, so the
energetic ranking of all three shells rests on the term least worth trusting.

This is the strongest statement this repository can make about its own scoring: **total fragment
interaction energy is the wrong criterion for choosing which shell to fold.** It should be replaced by
something geometric — closure quality and pose diversity — or at minimum weighted to neutralise the
charge bias. What it should *not* be is read as a prediction of how well a shell will fold, because on
three shells it got the order backwards.

### The helix that folds instead of breaking

The `esm2` variant is the most interesting fold in the set, and it works by a mechanism neither earlier
shell used. Unconstrained it is a 100% helical rod at Rg 15.0 Å with the ligand flung 13.6 Å away and
only 7 contacts — the worst geometry in shell 3. Force contacts on it and it does not break: it **bends
double, keeping its helix intact, and clamps the ligand between the two arms.**

| | `s3_esm2_control` | `s3_esm2_f4` | `s3_esm2_f8` | `s3_esm2_f12` | `s2_esm2_f8` |
|---|---|---|---|---|---|
| Rg | 15.0 Å | 10.3 Å | 10.2 Å | 9.7 Å | 9.8 Å |
| helical fraction | 100% | **91%** | **91%** | **88%** | 87% |
| β fraction | 0% | 0% | 0% | 0% | **6%** |
| helical H-bonds | 42 | **45** | **46** | 38 | 34 |
| non-local H-bonds | 0 | **0** | **0** | **0** | 0 |
| CA end-to-end | 49.8 Å | **17.6 Å** | **18.1 Å** | **8.5 Å** | 11.4 Å |
| CA contour length | 125.2 Å | 125.4 Å | 125.5 Å | 125.4 Å | 121.2 Å |
| closest approach of the two halves | 8.7 Å | 6.4 Å | 6.3 Å | 5.5 Å | 5.9 Å |

Read the last three rows together. The contour length along the CA trace is **unchanged** at 125 Å while
the end-to-end distance collapses from 49.8 Å to 8.5–18 Å: the chain is not compacting, it is folding
back on itself. The two halves come within 5.5–6.4 Å, which is helix-packing contact. And the helical
H-bond count goes *up*, to 45 and 46 — the helix is not merely surviving the fold, it is better
hydrogen-bonded after it.

Contrast shell 2's `s2_esm2_f8`, which reached a similar Rg by the opposite route: its helicity fell to
87%, it picked up 6% β, and its helical H-bonds dropped from 36 to 34. That one broke its helix to bind.
Shell 3's keeps it.

What holds the hairpin shut is worth stating because of what it is not: **zero non-local backbone
H-bonds** in any of these folds. There is no β-sheet, no backbone clasp. The two helical arms are held
together by side-chain packing alone — which is exactly what a designed shell is, twelve side chains
arranged around a ligand, and the closest this project has come to Boltz building the thing that was
designed rather than something else that happens to contain it.

It is also the one sequence in shell 3 whose forced folds are all net favourable (−19.43, −6.85,
−14.22), at ligand strains of 10.8–15.2 kcal/mol, the lowest in the shell. **And it satisfies none of
the designed contact hints — 0/4, 0/8, 0/12.** It wraps the ligand well, holds it cheaply, and does so
without honouring a single hint it was given, which is a useful reminder that hint satisfaction is a
third metric that can disagree with the other two.

### Figures

All twelve structures, peptide in green, ligand in red, one row per sequence and one column per rung of
the constraint ladder. Captions are `wrapped` / `interaction` in kcal/mol.

**The glycine design.** Only the fully forced fold encapsulates — and it is the best-enclosed structure
in the project, at a ligand strain of 44.05 kcal/mol that makes it net unfavourable anyway.

| unconstrained | 4 contacts | 8 contacts | 12 contacts |
|---|---|---|---|
| ![s3_orig_control](runs/octinoxate/figures_shell3/s3_orig_control.png) | ![s3_orig_f4](runs/octinoxate/figures_shell3/s3_orig_f4.png) | ![s3_orig_f8](runs/octinoxate/figures_shell3/s3_orig_f8.png) | ![s3_orig_f12](runs/octinoxate/figures_shell3/s3_orig_f12.png) |
| 0.80 / +0.92 | 0.80 / −28.45 | 0.85 / −5.62 | **1.00 / −25.13** |
| 22% helix and 22% β at Rg 7.8 Å; reproduces 2 designed positions, the most of any control | ligand 7.2 Å out, helix collapsed to 3% | 59 contacts but 55.77 kcal/mol of strain, the second worst in the shell | **encapsulates** — 0.965 enclosed, 20/20 engaged, 3.6 Å separation, 67 contacts, no clash; a tighter twin of shell 1's `orig_f12`, at +18.92 net |

**The first ESM2 variant.** Two of the four encapsulate, and they are separated by 50 kcal/mol of ligand
strain — the clearest single illustration that enclosure and energy are independent.

| unconstrained | 4 contacts | 8 contacts | 12 contacts |
|---|---|---|---|
| ![s3_esm1_control](runs/octinoxate/figures_shell3/s3_esm1_control.png) | ![s3_esm1_f4](runs/octinoxate/figures_shell3/s3_esm1_f4.png) | ![s3_esm1_f8](runs/octinoxate/figures_shell3/s3_esm1_f8.png) | ![s3_esm1_f12](runs/octinoxate/figures_shell3/s3_esm1_f12.png) |
| 0.85 / −26.88 | 0.70 / −19.30 | **0.95 / +2.00** | **0.90 / −25.14** |
| 59% helical at Rg 9.4 Å, net −13.21 without any forcing | the only rung that loses ground, 14/20 engaged | **encapsulates** — 0.95 wrapped, 19/20 engaged, 55 contacts, but **61.67 kcal/mol of strain**, the worst in the shell and net +63.67 | **encapsulates** — 0.88 enclosed at 4.3 Å, and only 11.65 kcal/mol of strain: net **−13.50**, the one fold here that is both well enclosed and cheap |

**The second ESM2 variant.** Net charge +4. A folded helix wrapping the ligand: the rod bends double
without breaking, and the two arms close on the ligand.

| unconstrained | 4 contacts | 8 contacts | 12 contacts |
|---|---|---|---|
| ![s3_esm2_control](runs/octinoxate/figures_shell3/s3_esm2_control.png) | ![s3_esm2_f4](runs/octinoxate/figures_shell3/s3_esm2_f4.png) | ![s3_esm2_f8](runs/octinoxate/figures_shell3/s3_esm2_f8.png) | ![s3_esm2_f12](runs/octinoxate/figures_shell3/s3_esm2_f12.png) |
| 0.45 / −10.78 | **0.90 / −34.59** | 0.80 / −21.56 | 0.75 / −25.05 |
| 100% helical rod, Rg 15.0 Å, ligand 13.6 Å away with 7 contacts — worst geometry in the shell | **folded helix**, 91% helical at Rg 10.3 Å, end-to-end 49.8 → 17.6 Å, 36 contacts; best sum in the shell at −19.43 | **folded helix**, 91% helical, 33 contacts, net −6.85 | most folded of all, end-to-end 8.5 Å, 88% helical; net −14.22 |

`load_folds.pml` in `figures_shell3/` loads all twelve superposed on the ligand with the designed shell
in marine. It loads and aligns only — run `@style.pml` beside it for the representation used here,
which is what draws the residue labels.

---

## What happens after this pipeline: dynamics and MM/GBSA

Everything above is static. A fold is predicted, its geometry measured and its energy computed at that
one geometry, in the gas phase, with no solvent and no time. That is enough to rank and to reject, and
it is not enough to believe. **This repository is the filter; molecular dynamics is the verdict.**

The distinction matters because the filter's own metrics say so. Enclosure and wrapping are the design
objective, interaction energy is a proxy that has inverted the verdict more than once, and ligand strain
is a gas-phase quantity with no desolvation term. None of those can tell you whether a complex survives
in water. What they can do is choose which two or three structures out of thirty-six are worth an hour
of GPU each.

### The handoff

MD and MM/GBSA live in a separate project, **`~/python_mac/openmm`** (the `omd` CLI), for a hard reason:
GAFF2/OpenFF ligand parameterisation is not pip-installable, so that pipeline needs a conda environment
(`openmm-md`: openmm, openmmforcefields, openff-toolkit, mdtraj, parmed, AmberTools). This repository's
`.venv` has openmm but not the parameterisation stack, and cannot get it. Two environments called at
their boundary, as with Boltz.

`code/cif_to_md.py` is the adapter. Boltz writes one CIF holding peptide and ligand together, heavy
atoms only; `omd` wants a protein PDB and a ligand SDF apart. The ligand needs care: a PDB or CIF
records no bond orders, so the SDF is built by assigning bond orders from the run's own SMILES as a
template, then adding hydrogens with coordinates — the same route `uma_binding.protonate_ligand` takes.
Without the template step the aromatic ring and the ester come out as single bonds and the ligand is
parameterised as something it is not. **Check the printed formula before going further.**

```bash
python code/cif_to_md.py runs/octinoxate --structure orig_f12

OMD=~/miniforge3/envs/openmm-md/bin/omd
M=runs/octinoxate/md/orig_f12
$OMD prep-protein --pdb $M/orig_f12_protein.pdb --out $M/protein_fixed.pdb
$OMD prep-ligand  --sdf $M/orig_f12_ligand.sdf  --out $M/ligand_prepped.sdf
$OMD build --protein $M/protein_fixed.pdb --ligand $M/ligand_prepped.sdf \
           --out-dir $M/system --no-auto-cofactors --box-shape dodecahedron
$OMD run     --system $M/system/system.xml --topology $M/system/complex.pdb \
             --out-dir $M/prod --steps 10000000 --platform OpenCL
$OMD analyze --traj $M/prod/traj.dcd --topology $M/system/complex.pdb --out-dir $M/prod
$OMD mmgbsa  --protein $M/protein_fixed.pdb --ligand $M/ligand_prepped.sdf \
             --traj $M/prod/traj_wrapped.xtc --topology $M/prod/traj_wrapped.pdb \
             --out-dir $M/prod/mmgbsa --no-auto-cofactors --run
```

`analyze` is not optional: it is what writes the wrapped, solute-only `traj_wrapped.xtc` and matching
`traj_wrapped.pdb` that `mmgbsa` requires.

### Where to run the dynamics

Only `run` wants a GPU. Everything either side of it — parameterisation, MM/GBSA — is conda-only and
cheap on CPU, so the split is clean, and `code/modal_md.py` takes just that one leg to a rented A10G.
It imports `openmm_md.dynamics.run` rather than reimplementing the protocol, so a Modal trajectory and a
local one come out of the same tested code: same minimisation, same restrained equilibration, same
integrator, same timestep. That is not a nicety — these trajectories are compared as binding energies,
and a protocol difference would be indistinguishable from a result.

| | local, six-core Apple Silicon | Modal A10G |
|---|---|---|
| platform | OpenCL (FAIRChem has no MPS backend, but OpenMM reaches the GPU this way) | CUDA, passed explicitly |
| rate | 1.72 ms/step at 8,647 particles | 0.376 ms/step at 23,121 particles |
| 20 ns of the 23k system | ~14 hours | **63 minutes, $1.15** |
| equilibration | forced onto CPU — OpenCL in single precision cannot run the restrained phase | fine on CUDA |

Two things are worth knowing before renting anything. `dynamics.py`'s `auto` platform probe does not
list CUDA, so on a GPU host it silently selects CPU — renting a GPU to compute on a processor. Pass
`CUDA` explicitly, which takes the branch that raises rather than degrades. And warm up before timing:
OpenMM's CUDA platform compiles and autotunes lazily, so an un-warmed probe reads 2–5× too slow.

**Box size, not peptide size, sets the cost.** A compact fold and an extended one of the same sequence
length differ enormously: `orig_f12` at radius of gyration 8.0 Å solvates to 8,647 particles, while
`s2_esm2_control` at 15.0 Å — a 100% helical rod — needs 33,045 in a cubic box. Read `peptide_rg` before
estimating a run, and use `--box-shape dodecahedron`: a rhombic dodecahedron holds the same minimum
image distance in ~71% of the volume, which took that system from 33,045 particles to **23,121**, a 30%
saving for no change in physics.

### How long is long enough

`orig_f12`, the best structure in the first worked example, computed on the leading *n* ns of one 20 ns
trajectory:

| window | VDWAALS | EEL | EGB | ESURF | **ΔG bind** |
|---|---|---|---|---|---|
| 40 ps (pilot) | −38.28 | −14.88 | +32.28 | −4.48 | **−25.37 ± 0.50** |
| 5 ns | −31.23 | −10.81 | +24.67 | −4.14 | **−21.52 ± 0.29** |
| 10 ns | −23.46 | −8.17 | +17.89 | −3.12 | **−16.86 ± 0.37** |
| 15 ns | −20.43 | −6.08 | +14.60 | −2.70 | **−14.61 ± 0.38** |
| 20 ns | −19.75 | −5.32 | +14.01 | −2.64 | **−13.71 ± 0.04** |

The successive changes are 3.85, 4.66, 2.25, 0.90 — halving, so 20 ns is converged and anything shorter
is not. **A short window does not give a noisy answer, it gives a confidently wrong one:** the 40 ps
pilot's ±0.50 is the sampling error on forty correlated frames drawn from the starting pose, and it
overstates binding by 11.7 kcal/mol. Every window is tight and every short one is too negative, because
they are all still measuring the Boltz pose rather than the ensemble. **Run 20 ns; do not trust 5.**

What relaxes is the pose itself. Over the run the ligand leaves the cavity it was designed into —
heavy-atom contacts within 4 Å fall 71 → 19 and the centroid separation goes 4.7 → 9.9 Å by 6 ns — sits
on the surface for twelve nanoseconds, then re-inserts in the last two (5.07 Å, 31 contacts). It never
dissociates: 8 frames out of 20,000 have no contact at all. The designed 2.44 Å static clash never
recurs, the closest heavy-atom approach over the whole trajectory being 2.61 Å. So the enclosure was
real but not a minimum — which is exactly the kind of statement no static score can make.

### The charge artifact, tested

The static pipeline's worst known bias is that gas-phase interaction energy rewards net charge: across
36 folds it correlates −0.56 with charged-residue fraction, and charged fragment poses average −10.79
kcal/mol against −4.71 aromatic and −2.94 aliphatic. `EGB` is precisely the term that penalises exposed
charge, so MM/GBSA should take most of that advantage away. That is a falsifiable prediction about a
specific number, and the second 20 ns run was run to test it.

| | `orig_f12` | `s2_esm2_control` |
|---|---|---|
| peptide charge | −1 | **+5** |
| peptide Rg | 8.0 Å (compact) | 15.4 Å (100% helical rod) |
| UMA gas-phase interaction | −45.86 | **−74.89** |
| VDWAALS | −19.75 | −19.12 |
| EEL | −5.32 | −2.62 |
| EGB | +14.01 | +8.17 |
| ESURF | −2.64 | −2.67 |
| **MM/GBSA ΔG bind** | **−13.71 ± 0.04** | **−16.25 ± 0.08** |

**The prediction holds.** A 29.0 kcal/mol gas-phase advantage becomes a 2.5 kcal/mol advantage in
water: 91% of it evaporates. Anyone ranking these two on the static number would have believed the
+5 sequence was in a different class, and it is not.

The mechanism is not the one the artifact's name suggests. The +5 peptide's electrostatic attraction to
the ligand is *weaker*, not stronger (−2.62 against −5.32), which makes sense the moment one remembers
the ligand is neutral — there is no charge for a charge to pull on. What UMA was rewarding is long-range
polarisation of a neutral ligand by a charged rod, and GB screens that almost completely. The +5 peptide
then wins its small margin back through desolvation, not attraction: it buries less (+8.17 against
+14.01) because its ligand lies in a shallow surface groove rather than a closed cavity. The two
structures have essentially the same van der Waals contact (−19.1 against −19.8) reached two different
ways, and van der Waals dominates electrostatics in both — which is what the contact analysis of
`orig_f12` already said, 47 contacts and none of them charged.

That is also the flat answer to whether these two can be reconciled: **in water neither is
remarkable.** −13.7 and −16.3 kcal/mol are ordinary millimolar-to-micromolar affinities, not the −45 to
−75 the gas phase advertised.

### Which metric was right: four structures in water

Two shell-3 structures were then run to 20 ns for the same cost as one of the above, chosen because the
static pipeline ranked them at opposite extremes. `s3_orig_f12` is the tightest fold in the project and
the *worst* static score of the four; `s3_esm2_f4` is the folded helical hairpin.

| structure | ΔG bind | VDWAALS | EEL | EGB | UMA interaction | ligand strain | static sum | enclosed | contacts |
|---|---|---|---|---|---|---|---|---|---|
| **`s3_orig_f12`** | **−24.33 ± 0.08** | −34.41 | −7.99 | +22.74 | −25.13 | 44.05 | **+18.92** | **0.965** | **67** |
| **`s3_esm2_f4`** | **−21.08 ± 0.08** | −23.45 | −2.60 | +8.58 | −34.59 | 15.16 | −19.43 | 0.795 | 36 |
| `s2_esm2_control` | −16.25 ± 0.08 | −19.12 | −2.62 | +8.17 | −74.89 | 9.05 | **−65.84** | 0.510 | 18 |
| `orig_f12` | −13.71 ± 0.04 | −19.75 | −5.32 | +14.01 | −45.86 | 12.80 | −33.07 | 0.960 | 47 |

**The static score ranks these four in close to reverse order.** Spearman against ΔG is −0.80 for the
static sum, −0.80 for UMA interaction alone, +0.40 for enclosed fraction and +0.80 for MM/GBSA's own
van der Waals term. At n = 4 those coefficients are indicative rather than statistical — one swapped
pair moves them a long way — but the endpoints need no statistics: **the best static score is third of
four in water, and the worst static score is first.**

Two of those comparisons are worth isolating.

**`orig_f12` against `s3_orig_f12`: the same fold twice, ranked backwards.** These are geometric
near-twins from shells 47% apart in fragment energy — 0.960 against 0.965 enclosed, both 1.00 wrapped
and 20/20 engaged, Rg 8.0 against 8.1 Å. The static sum separates them by 52 kcal/mol in favour of
`orig_f12`; MM/GBSA separates them by 10.6 in favour of `s3_orig_f12`. **And the inversion is not the
strain term's fault**, which matters because single-trajectory MM/GBSA reports BOND, ANGLE and DIHED as
exactly zero — intramolecular terms cancel by construction, so ΔG contains no ligand strain at all and
comparing it to a static sum that does would be unfair. Comparing interaction alone:

| | `orig_f12` | `s3_orig_f12` | difference |
|---|---|---|---|
| UMA gas-phase interaction | −45.86 | −25.13 | **+20.73**, UMA prefers `orig_f12` |
| MM/GBSA gas-phase, VDW + EEL | −25.08 | −42.41 | **−17.33**, the force field prefers `s3_orig_f12` |
| of which VDWAALS | −19.75 | −34.41 | −14.66 |

The sign inverts on the interaction term by itself: a 38 kcal/mol swing in the difference between two
folds no geometric measure can tell apart. So the charge artifact is not the only problem with the
static energies — the interaction term misranks a pair where charge is not even the variable.

**And the geometry called it.** `s3_orig_f12` has 67 contacts against 47, and its VDWAALS is 14.7
kcal/mol stronger. The force field rewards buried contact area, which is exactly the quantity
`enclosed_fraction` and `contacts_under_cutoff` estimate, and it charges 8.7 more desolvation for the
privilege and still comes out 10.6 ahead. **Read the geometry, not the energy.**

**Does the designed pose survive?** `code/md_contacts.py` profiles a trajectory in windows; the answer
differs per structure and only one of the four decays:

| structure | contacts, first → last decile | centroid separation | frames with no contact |
|---|---|---|---|
| `s3_orig_f12` | 50.7 → **60.8** (mean 56.6) | 5.02 → **3.75 Å**, tightest 2.24 Å at 6–8 ns | **0 / 2000** |
| `s3_esm2_f4` | 19.8 → 20.6 (mean 20.1) | 7.67 → 5.82 Å | **0 / 2000** |
| `s2_esm2_control` | 22.0 → 22.4 (mean 23.4) | 8.61 → 11.38 Å | 1 / 2000 |
| `orig_f12` | 71.0 → 31.4 (mean 34.3) | 4.7 → 5.07 Å, out to **13.6 Å** | 8 / 20,000 |

`s3_orig_f12` is the first structure here whose designed enclosure **improves** in water rather than
merely surviving: it gains ten contacts and closes from 5.0 to 3.8 Å. Its twin does the opposite,
leaving the cavity at 6 ns, sitting on the surface for twelve nanoseconds and re-inserting only in the
last two. Neither ever releases the ligand, and no static clash recurs — closest heavy-atom approaches
over whole trajectories are 2.59 to 2.66 Å, against the 2.44 Å that `orig_f12`'s prediction contained.

**The helical hairpin holds.** `s3_esm2_f4` was the structural question mark: two helical arms clamped
on the ligand by side-chain packing with **zero non-local backbone H-bonds**, which looked like the sort
of thing water would prise open. It does not.

| | Boltz prediction | over 20 ns |
|---|---|---|
| helical fraction | 91% | 81% |
| β fraction | 0% | **0%**, never converts |
| CA end-to-end | 17.6 Å | **18.1 ± 2.9 Å** |
| closest approach of the two arms | 6.4 Å | **6.4 ± 0.5 Å** |
| peptide Rg | 10.3 Å | 10.4–10.8 Å |

End-to-end never drifts back toward the unfolded rod's 49.8 Å, and the arms hold at 6.4 Å with a
standard deviation of 0.5 Å. Side-chain packing alone keeps the hairpin shut for 20 ns at 300 K, which
is the closest this project has come to a designed tertiary arrangement being confirmed rather than
merely predicted.

### Figures: what 20 ns does to a designed pose

Peptide in green, ligand in red, as in the static figures above — but these are frames from a
trajectory, not predictions. Each is the medoid of a representative window: the single frame closest to
that window's mean structure, written by `code/md_frames.py`. **The two on the left are the same
designed fold, arrived at from different shells, and they are the whole argument of this section.**

| `s3_orig_f12` at 18.3 ns | `orig_f12` at 11.0 ns | `s3_esm2_f4` at 19.4 ns | `s2_esm2_control` at 18.9 ns |
|---|---|---|---|
| ![s3_orig_f12](runs/octinoxate/md/figures/s3_orig_f12_medoid.png) | ![orig_f12](runs/octinoxate/md/figures/orig_f12_plateau_medoid.png) | ![s3_esm2_f4](runs/octinoxate/md/figures/s3_esm2_f4_medoid.png) | ![s2_esm2_control](runs/octinoxate/md/figures/s2_esm2_control_medoid.png) |
| **ΔG −24.33** | **ΔG −13.71** | **ΔG −21.08** | **ΔG −16.25** |
| 0.965 enclosed as predicted, and it **closed further** — separation 5.02 → 3.75 Å, contacts 50.7 → 60.8, never released | 0.960 enclosed as predicted, and the ligand **left** at 6 ns; shown on the plateau it held for twelve, at 9–10 Å | the **helical hairpin**, still 81% helical with its arms at 6.4 ± 0.5 Å, holding the ligand on side-chain packing alone | the 100% helical **rod**, ligand in a surface groove, sliding along it 8.6 → 11.4 Å without letting go |

`orig_f12` is shown at 11.0 ns rather than at the end, because its last 2 ns are the only stretch where
the ligand is back in the cavity and a frame from there would imply the enclosure held. It did not.
That is the point of `md_frames.py --window-ns`: the end of a run is usually the settled part, and when
it is not, saying so in the caption is better than picking the flattering frame.

Two predictions at 0.960 and 0.965 enclosed, indistinguishable on every static geometric measure,
ranked 52 kcal/mol apart by the static energies in the wrong direction — and one holds its ligand while
the other does not. **No number computed at a single geometry could have told them apart.**

Results across all legs are in `runs/octinoxate/md/mmgbsa_summary.csv`. `code/md_frames.py` writes the
PDBs these are rendered from — the medoid above, and a ten-model ensemble superposed on the peptide for
showing the spread instead of one frame.

---

## What the metrics mean

Two of these are easy to confuse, and the fact that they disagree is the point.

**`enclosed`** — the fraction of directions out of the ligand's centre that run into peptide. Two
hundred directions are sampled, each treated as a ray with a 3 Å tolerance out to 12 Å. This asks
whether the ligand sits in a shell or against a face.

**`wrapped`** — the fraction of the ligand's heavy atoms with a peptide heavy atom within 4.5 Å. This
asks how much of the ligand is in contact.

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
[the third worked example](#third-worked-example-testing-the-selection-criterion-itself).

**Peptide strain is still not reported.** Relaxing a free peptide in vacuum measures collapse rather
than strain, so the term is disabled. Referencing it instead against the peptide's Boltz apo fold was
tried and does not work — vacuum still rewards the compact bound folds over an extended apo helix by
over 100 kcal/mol, and the apo prediction is only confident for sequences whose fold does not change on
binding. See `archive/peptide_strain/` for the attempt and the evidence. It needs implicit solvent to be
meaningful, which puts it in MM/GBSA territory rather than here.

**No experimental validation.** Nothing here has been synthesised or measured. The energies come from a
machine learning potential in the gas phase, with no solvent, no entropy and no desolvation.

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
