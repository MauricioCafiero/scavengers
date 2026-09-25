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
binding_energy.py          interaction energy and ligand strain for each folded complex
        │
        ├─ boltz/binding_*.csv
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

### Optional overrides

pdbfixer and ESM2 run in process. These only matter if an environment cannot install them, in which
case each falls back to calling another interpreter:

| flag | environment variable | what it provides |
|---|---|---|
| `--fixer-venv` | `$PEPTIDEBUILDER_FIXER_VENV` | pdbfixer / OpenMM, adds hydrogens to Boltz output |
| `--genmask-venv` | `$PEPTIDEBUILDER_GENMASK_VENV` | ESM2 / transformers, fills glycine linkers |

Both routes run identical code and return identical numbers.

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
| `make_figures.py` | flat filenames, a manifest, and a PyMOL loading script |

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
├── condense_<shell>.csv         reachability sweep: closure per pair, spacer count, CA choice
├── design.json                  chosen path, spacer counts, sequence, budget check
├── condense/                    condensed structures from the chain route
├── boltz/
│   ├── *.yaml, *.log            Boltz inputs and logs
│   ├── boltz_results_*/         folded complexes, as .cif
│   ├── structures/              every geometry a relaxation produced
│   ├── fold_check.csv           enclosure, wrapping, contacts, hints honoured
│   ├── binding_*.csv            interaction and strain per complex
│   └── partial_*.json           per-term results, saved as computed
└── figures/
    ├── *.cif                    every fold under a readable name
    ├── designed_shell.xyz       the arrangement that was asked for
    ├── manifest.csv             every structure joined to its numbers
    └── load_folds.pml           PyMOL script: loads all, superposed on the ligand
```

Scoring logs are written to `runs/*.log` and are kept, not ignored: they hold each relaxation's step
count, timing and convergence, which is the only record of how a number was arrived at.

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
| `orig_control` | 0.495 | 0.75 | 15/20 | 20.9 Å | −15.83 | 9.14 | −6.68 | — |
| `orig_f4` | 0.595 | 0.85 | 17/20 | 7.1 Å | −26.32 | 10.32 | −16.01 | 3/4 |
| `orig_f8` | 0.515 | 0.65 | 13/20 | 8.4 Å | −22.30 | 36.44 \* | +14.14 | 1/8 |
| **`orig_f12`** | **0.96** | **1.00** | 20/20 | 4.2 Å | **−45.86** | 8.95 | **−36.91** | 0/12 |
| `esm1_control` | 0.455 | 0.60 | 12/20 | 8.4 Å | −12.22 | 4.80 | −7.42 | — |
| `esm1_f4` | 0.605 | 0.65 | 13/20 | 6.6 Å | −20.89 | 20.66 | −0.23 | 0/4 |
| `esm1_f8` | 0.75 | 0.90 | 18/20 | 4.1 Å | −29.10 | 27.78 | −1.32 | 1/8 |
| `esm1_f12` | 0.76 | 0.80 | 16/20 | 3.4 Å | −17.11 | 7.13 | −9.98 | 1/12 |
| `esm2_control` | 0.53 | 0.40 | 8/20 | 7.7 Å | −20.27 | 10.89 | −9.38 | — |
| `esm2_f4` | 0.46 | 0.50 | 10/20 | 11.0 Å | −39.90 | 15.70 | −24.20 | 0/4 |
| `esm2_f8` | 0.625 | 0.75 | 15/20 | 10.4 Å | −29.66 | 8.85 | −20.81 | 0/8 |
| `esm2_f12` | 0.905 | 1.00 | 20/20 | 7.1 Å | −21.14 | 6.92 | −14.22 | 0/12 |

Energies in kcal/mol. All three sequences are 33 residues.
\* `orig_f8`'s ligand strain did not converge within its step budget, so 36.44 is a lower bound and
that row is worse than shown.

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
interaction energy by bending the ligand rather than by surrounding it: 36.44, 27.78 and 20.66 kcal/mol
of ligand strain, against 6.92 to 9.14 for every twelve-contact fold. Two of them end up net
unfavourable. A peptide can satisfy a few contacts by pulling the ligand toward whichever residues are
nearby, but satisfying twelve at once requires surrounding it, and a ligand can only be surrounded in a
conformation it can actually adopt.

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
| ligand against the outside, centroids 20.9 Å apart | | ligand bent to buy contacts: 36.44 kcal/mol of strain | every ligand atom contacted, 47 pairs, none charged |

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

**`interaction`** — `E(complex) − E(peptide) − E(ligand)`, all three at the complex geometry, so it is a
rigid interaction energy containing no strain. Hydrogens are relaxed first with heavy atoms fixed,
which leaves Boltz's predicted geometry untouched while removing the arbitrariness of where pdbfixer
placed them.

**`ligand strain`** — `E(ligand at its geometry in the complex) − E(ligand relaxed alone)`. The
conformational price the ligand pays. This is real signal rather than noise: it varied 4.34 kcal/mol
between two structures whose repeat measurements agree to 0.6.

**`hints`** — requested contacts satisfied, of those requested. Worth reporting and worth not trusting
as a measure of success: it is uncorrelated with everything else in the table.

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

**One shell, mostly.** Twenty shells were enumerated and one was swept in full. A second is in progress.
The twelve-structure comparison rests on a single designed arrangement.

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
