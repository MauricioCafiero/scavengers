# Handoff: fix the code, then run the second shell

Written 2026-09-25, ~21:00, for a fresh session with no prior context. Read this first.

Two jobs, in order: **fix four things in the code**, then **run the pipeline on a second designed
shell that is already on disk**. After that, the README gains the second shell's results and the
corrected setup instructions.

---

## 1. Where things stand

### Data on disk, all under `runs/octinoxate/`

* **165 poses** in `poses/`, with `energies.csv` listing every pose's UMA interaction energy and
  which twelve are selected. The fragment library was repaired earlier today (see §5), so these
  poses supersede any results predating 2026-09-25 07:41.
* **Shell 1 swept**: `condense_shell12.csv`, all 132 ordered pairs of the richest 12-pose shell at
  1–4 glycine spacers. Every pair is connectable at some spacer count.
* **The design**: `design.json` — path through all twelve poses,
  `RGGDGGKGGGGLGGGGKGIGEGWGGDGSGGEGS`, 33 residues, 21 spacers against a budget of 16.
* **Twelve folded complexes** in `boltz/boltz_results_*/`, also copied under readable names into
  `figures/`. Three sequences (the glycine design plus two ESM2 linker variants) × four constraint
  levels (unconstrained control, 4, 8 and 12 forced contacts).
* **Scoring**: `boltz/binding_is.csv` (six structures), `boltz/binding_esmf.csv` (three),
  `boltz/binding_esm2f.csv` (the last three). **Check this last file is complete** — when this was
  written, `esm2_f4` and `esm2_f8` were still running. If the job died, rerun only the missing ones:

  ```bash
  .venv/bin/python code/binding_energy.py runs/octinoxate \
      --terms interaction,strain_ligand --out binding_esm2f.csv \
      --structures esm2_f4,esm2_f8
  ```

* **Geometry for all twelve** in `boltz/fold_check.csv` — enclosure, wrapping, contacts, hints
  satisfied. Complete regardless of scoring.
* `figures/manifest.csv` joins every structure to its numbers; `figures/load_folds.pml` loads them
  all superposed on the ligand.

### Documents

* `README.md` — the public README. Complete: overview, workflow, installation, file-by-file
  description, the octinoxate worked example with its twelve-row table and twelve renders, metrics,
  traps, limitations and credit. Shell 2 is to be added to it as a second worked example (§4).
* `README_previous.md` — the old README, carrying trials and negative results, kept as an internal
  document.
* `PROGRESS.md` — the full working log, with a findings list and the README summary table at the top.

### What is *not* done

**Shell 2 ran overnight on 2026-09-26 and is all but complete.** Sweep (132/132 ordered pairs, 4752
rows), design, two ESM2 variants, twelve folds, geometry check, overlay and figures are all done; the
scoring of the twelve complexes was still finishing at the time of writing. §8 holds the queue of
shells to test next and the condition for starting them. §9 records what shell 2 showed.

Everything else listed here in earlier drafts is finished, so do not redo it:

* All four code fixes in §2 — Boltz resolver wired in, pdbfixer and ESM2 local, ligand naming local,
  README installation section rewritten. Verified, committed, pushed.
* All twelve structures scored. `figures/manifest.csv` carries twelve complete rows of enclosure,
  wrapping, contacts, interaction and ligand strain, and is the joined table to read.
* The README is written and is now `README.md`, with all twelve renders in place and captions
  checked against the manifest. The previous one is `README_previous.md`.
* The repository has been cleaned and the history rebuilt: no inherited commits or files, `origin`
  removed, `scavengers` the only remote. See §7 and `NEXT_STEPS.md`.

---

## 2. The four code fixes — 1–3 are DONE (2026-09-25, 22:22), 4 remains

**Do not redo fixes 1–3.** They are applied, verified, and in the working tree (uncommitted).
`grep -rn "python_mac/" code/*.py` returns nothing but `boltz_env.py`'s documented fallback and one
provenance comment in `solvate.py`. What was done, and the evidence:

| fix | what changed | verified by |
|---|---|---|
| 1. Boltz resolver | `boltz_hints.fold()` and `boltz_check.run_boltz()` now call `boltz_env.run_boltz()`; defaults read `$PEPTIDEBUILDER_BOLTZ_*` | `resolve_boltz()` falls through to the legacy MPS wrapper on this machine, so old behaviour is intact; explicit `--boltz-venv` honoured |
| 2a. pdbfixer local | `uma_binding.protonate_peptide()` imports `PDBFixer` in process, subprocess kept as fallback; `fixer_venv` now optional | `esm2_f4`: 271 heavy → **557 atoms**, exact match to the subprocess route that produced the scored CSVs |
| 2b. ESM2 local | `fill_linkers.fill()` runs `FILL_SCRIPT` in process via `exec` (same text, so the two routes cannot drift), subprocess kept as fallback | `RGYGGLGGEGFGSGDGWGKGI` → `RLYLALLREGFLSSDLWRKSI`, 11 linkers filled |
| 3. ligand naming | `ligand_atom_names()` defaults to local RDKit; `_ligand_atom_names_via_boltz()` kept for the affinity case only | identical to Boltz's names on **all 44 atoms** |

Environment, as actually installed (note: **`.venv` is a `uv` venv with no `pip` in it** — use
`uv pip install --python .venv/bin/python …`, *not* `.venv/bin/python -m pip`, which fails with
"No module named pip"):

```bash
uv pip install --python .venv/bin/python pdbfixer transformers   # openmm comes in with pdbfixer
```

`pdbfixer 1.12.0`, `openmm 8.6.1`, `transformers 5.17.0` — all on PyPI, no conda needed.
`transformers` downgrades `huggingface-hub` 2.0.0 → 1.33.0. **This is safe and was checked**:
fairchem's entire use of that library is one call to `hf_hub_download`, and `~/python_mac/neb/.venv`
already runs fairchem 2.21.0 on hub 1.24.0. A UMA water single point gave −2079.8644 eV before and
after, identical. Pre-install snapshot for rollback: `runs/_env/venv_snapshot_20260925_2151.txt`.

Remaining: **Fix 4 only** (README installation section, below).

---

### The original problem, for context

This repo hardwired absolute paths into another user's home directory, with no
fallback, no environment variable, and no usable error when they were absent. That made it
unrunnable for anyone else. Four paths, three of them introduced by earlier sessions of this
assistant:

```
~/python_mac/boltz_local          boltz_hints.py, boltz_check.py   (Boltz-2)
~/python_mac/pocket_assist/venv   uma_binding.py, binding_energy.py (pdbfixer/OpenMM)
~/python_mac/GenMaskFill/.venv    fill_linkers.py                   (ESM2/transformers)
~/python_mac/UMADock              docstring reference only          (attribution, fine as is)
```

The user's decision: **Boltz may stay a separate environment** — it is heavy and GPU-bound — but it
must be discoverable rather than hardwired. **pdbfixer and ESM2 should run locally**, in this repo's
own `.venv`.

### Fix 1 — wire in the Boltz resolver — **DONE**

`code/boltz_env.py` already exists and is tested to import, but **nothing uses it**. It resolves how
to run Boltz in this order: explicit `--boltz-cmd` or `$PEPTIDEBUILDER_BOLTZ_CMD`; then
`--boltz-venv` or `$PEPTIDEBUILDER_BOLTZ_VENV`; then `boltz` importable in this environment; then
`boltz` on `$PATH`; then the original developer's layout; then a clear error listing the options.

Replace the hardwired `BOLTZ_REPO` / `BOLTZ_VENV` defaults and the `subprocess.call` in:

* `code/boltz_hints.py` — `fold()`, currently builds `[python, wrapper, "predict", …]`
* `code/boltz_check.py` — `run_boltz()`, the same pattern

Use `boltz_env.run_boltz(yaml, out_dir, log_path, extra_args=…)` and `boltz_env.add_arguments(parser)`.
Keep passing `--use_potentials` through `extra_args` — see §5, it is essential.

### Fix 2 — pdbfixer and transformers in-process — **DONE**

```bash
uv pip install --python .venv/bin/python transformers pdbfixer   # openmm comes with pdbfixer
```

Note `uv pip`: this is a `uv` venv and has **no pip at all**, so both `.venv/bin/pip` and
`.venv/bin/python -m pip` fail.
**Do this only when no scoring job is running**: a scoring process holds 0.5–2.5 GB and this machine
has run with under 1 GB of free swap all day.

Then:

* `code/uma_binding.py` — `protonate_peptide()` writes a PDB, shells out to `$FIXER_VENV/bin/python`
  with `FIXER_SCRIPT`, and reads the result back. Replace with a direct `from pdbfixer import
  PDBFixer` in-process. `FIXER_SCRIPT` already contains exactly the code needed; it just runs
  elsewhere. Keep the subprocess path as a fallback when the import fails.
* `code/fill_linkers.py` — `fill()` shells out with `FILL_SCRIPT`. Same treatment with
  `transformers`. This is the easier of the two: the user's words, "that one requires one little
  piece of code."
* `code/binding_energy.py` imports `protonate_peptide` from `uma_binding`, so it inherits the fix
  for free. Check nothing else needs touching.

Verify by rescoring one already-scored structure and confirming the interaction energy reproduces:
`orig_control` should give **−15.8 ± 0.1 kcal/mol**.

### Fix 3 — ligand atom naming — **DONE**

`boltz_hints.ligand_atom_names()` runs a script *inside Boltz's environment* because Boltz names a
SMILES ligand's atoms `element + canonical rank` computed after its own `standardize()`, and the
contact constraints refer to those names. So even the offline export path currently needs a
reachable Boltz install.

Reimplement it with local RDKit — `standardize` is a thin wrapper, and the naming is
`AddHs` then `CanonicalRankAtoms` then `symbol.upper() + str(rank + 1)`. **Verify the names match
Boltz's before trusting it**: the current implementation is known-correct (checked today: 44 atoms,
zero element mismatches, composition C17H24O3), so compare output against it while a Boltz install
is still available. Keep the subprocess version as the fallback.

### Fix 4 — README installation section — **DONE**

The section described three external environments at absolute paths. It now describes one required
environment (this repo's `.venv`, with the `uv pip` note) and one optional external one, Boltz,
resolved at call time by `boltz_env.py` with `boltz_offline.py` for machines that have none. No
absolute paths remain anywhere in the README.

---

## 3. Then: run the second shell

### Before anything: ask, then keep the machine awake

**Ask before starting the run.** §6 is not a formality — it is why this context was reset. Say what
you intend to run and how long it will take, and wait. Do not begin because this document describes
the commands.

**Then arm a sleep assertion for the whole session, before the first job.** This machine sleeps
mid-conversation, not only mid-job, which stalls detached work and interrupts the session, so the
assertion has to cover the session rather than an individual calculation:

```sh
caffeinate -is -t 43200 &          # 12 h, dies on its own if forgotten
pmset -g assertions | grep -E 'PreventSystemSleep +[01]'   # 1 once armed
```

Arm one, not several, and give it a timeout so it cannot outlive the session indefinitely. Disarm
it when the work ends:

```sh
pgrep -f 'caffeinate -is'          # ignore the Bash tool's own short -t 300 wrapper
kill <pid>
```

The assertion armed on 2026-09-25 at 16:49 ran 12 hours and expired by itself around 04:49 on the
26th. Do not assume one is running — check with `pmset` before starting a long job.


The point is to stop resting every conclusion on one designed arrangement. **Everything must be
scored**, not a subset — the user was explicit about that.

### The shell

Pose selection is greedy and order-dependent, so rotating which fragment picks first gives a family
of shells. Twenty were enumerated (ten starting fragments × one and two copies); shell 1 was the
richest at −113.12 kcal/mol. Shell 2 is the second richest, at −109.15, and shares only five of its
twelve poses with shell 1:

```
phenylalanine:2,arginine:3,lysine:10,aspartic:21,glutamic:15,isoleucine:8,leucine:12,serine:32,tryptophan:5,glutamic:10,leucine:13,serine:2
```

### The commands

```bash
SHELL2="phenylalanine:2,arginine:3,lysine:10,aspartic:21,glutamic:15,isoleucine:8,leucine:12,serine:32,tryptophan:5,glutamic:10,leucine:13,serine:2"

# 1. reachability sweep. ~1 h, pure geometry, no UMA, safe alongside other work.
#    k=0 is skipped deliberately: it is settled, 0 successes in 810 attempts.
.venv/bin/python code/condense.py sweep runs/octinoxate \
    --poses "$SHELL2" --kmin 1 --kmax 4 --restarts 24 --out condense_shell2.csv

# 2. ordering and spacer counts
.venv/bin/python code/assign.py runs/octinoxate --poses "$SHELL2" \
    --pairs-csv condense_shell2.csv --save runs/octinoxate/design_shell2.json

# 3. twelve folds: three sequences x (control, f4, f8, f12). ~3 min each.
#    force:true AND --use_potentials are both required -- see §5.
#    Use design_test.py, or boltz_hints.py per fold as the twelve were done for shell 1.

# 4. geometry check -- free, and the primary measure
.venv/bin/python code/check_fold.py runs/octinoxate

# 5. score all twelve. THE BOTTLENECK: 30-70 min each on CPU, so ~10 h.
.venv/bin/python code/binding_energy.py runs/octinoxate \
    --terms interaction,strain_ligand --out binding_shell2.csv --structures <all twelve>
```

### On the GPU

The user has **about four hours of A100 on Colab**. Scoring is the bottleneck and
`binding_energy.py` already selects CUDA when `torch.cuda.is_available()`, so it runs there
unchanged — minutes per structure instead of tens of minutes. `code/boltz_offline.py` exists for
exactly this split: `export` writes every Boltz input into one portable folder with a `RUN.md`, and
`import` brings whatever comes back into the layout the rest of the pipeline expects. Offer this
before committing ten hours of laptop CPU.

### Naming

Shell 1's structures are `orig_*`, `esm1_*`, `esm2_*`. Shell 2's must not collide — use a distinct
prefix such as `s2_orig_*`. `make_figures.py` maps fold names to figure names and will need the new
prefix handling.

---

## 4. Then: add shell 2 to the README

The README is otherwise complete. Only one thing is wanted from shell 2:

* add its results as a **second worked example**, and say whether the shell-1 findings replicate —
  particularly whether f12 folds again enclose the ligand with low ligand strain while partial
  constraint sets distort it, and whether wrapping tracks interaction energy within a sequence or
  inverts as it did for `esm2`.

The existing worked example, its twelve-row table and its three figure ladders stay as they are.
That is the comparison shell 2 is being run to provide.

---

## 5. Traps that will cost you hours

**Boltz contact constraints.** `force: false` is **discarded** by the featurizer, not softened — a
fold with unforced constraints is an unconstrained fold. And `force: true` alone is only
conditioning; the guidance that acts on the structure is gated behind **`--use_potentials`**, off by
default. Both are required. Four folds were wasted today learning this.

**Boltz affinity is broken here.** Boltz 2.2.1 calls `AllChem.Descriptors.MolWt`, and RDKit
2026.03.6 no longer exposes `Descriptors` through `AllChem`, so requesting affinity fails at parse
time. `boltz_hints.py` omits it. Little lost — that head spans 1.65 kcal/mol across 37 sequences.
Note this also changes ligand atom naming, since `standardize()` is only called when affinity is
requested; the namer takes the same flag and they must agree.

**Hydrogens.** Boltz writes heavy atoms only. pdbfixer's placement is not deterministic and moved
one interaction energy by **1.2 kcal/mol** between runs. Hydrogen-only relaxation with heavy atoms
fixed brings reproducibility to 0.08 and is now the default — leave it on.

**Peptide strain is meaningless in vacuum**: 237.8 kcal/mol, unconverged, because a free peptide
collapses into a compact ball. Off by default (`--terms strain_peptide` to enable). Ligand strain is
fine and is real signal.

**Explicit-water cavity desolvation is unaffordable on CPU** — ~420 atoms, still descending
1.6 kcal/mol per step at step 102 of 200. For same-length designs it largely cancels anyway.

**FAIRChem has no MPS backend** — `cpu` or `cuda` only. And it refuses the UMA 1.0 checkpoint on
2.22+; use `uma-s-1p2p1`.

**Short contacts are not automatically defects.** Of three sub-2.6 Å contacts found today, one was
an n→π\* interaction (70° out of the carbonyl plane), one a CH···O hydrogen bond (17° off the C–H
vector), and only the third — backbone carbonyl O to ligand ether O, no donor available — was
unphysical. Inspect the geometry before concluding.

**Enclosure and wrapping are different measures and they disagree.** `orig_control` contacts 75% of
the ligand's atoms while enclosing only 49% of the directions out of it, centroids 20.9 Å apart —
the ligand is stuck to the outside. **Wrapping and enclosure are the design objective; hint
satisfaction and interaction energy are proxies that have both inverted the verdict.** Report
enclosure beside any energy.

---

## 6. How the user wants to be worked with

These are corrections from today, and the reason the context is being reset.

**Ask before starting or stopping any calculation.** Not after. This was said explicitly and then
violated repeatedly — a long run launched before being mentioned, jobs reniced rather than queried, a
sweep started unprompted, and a three-hour scoring job killed on a one-word reply that turned out to
be about something else.

**A rhetorical question is not an instruction.** "Now you realize memory is tight?" and "Stop" in
the middle of a discussion are criticism, not commands to kill things. Ask which thing is meant.

**Do the thing asked, not an adjacent thing.** A request to reorder a table means reorder the table,
not rewrite the generator that produces it. A request for clarification means answer, not launch.

**Command output is not visible to the user.** Tables printed inside a tool call go nowhere. Put
results in the reply text. Several messages were wasted before this was understood.

**Fix the generator, not the generated file.** A hand-edit to `load_folds.pml` was undone the next
time `make_figures.py` ran, reinstating two display problems the user had already objected to.

**Keep a monitor armed on long jobs and re-arm on expiry** (the cap is 30 minutes; structures take
about an hour, so it will expire between results). A lapsed monitor meant results were only reported
when asked for.

**Do not over-generalise from three data points.** Two claimed patterns were announced and killed by
the next measurement today.

**READMEs and repo documents are yours to write.** The `own-prose` skill applies to papers and
publications only; its wording was narrowed today after it caused a README to be refused.

---

## 7. Provenance, for the README and any publication

* The original repository is **[lucia-71/peptidebuilder](https://github.com/lucia-71/peptidebuilder)**,
  the work of the user's masters student under his supervision. It contributed the amino-acid
  fragment library, pose selection, sequence assembly, the linker estimate and the CLI — the whole
  peptide-specific conception — and the majority of the founding commits. She is to be credited at
  the top of the README and for most of the initial work in any publication.
* The placement and scoring engine beneath it descends from the user's own
  **[CafChemFragGrow](https://github.com/MauricioCafiero/CafChem)** (`notebooks/FragGrow_CafChem.ipynb`),
  which works with generic fragments — water, cyclopropyl, acetylene, methanol, phenyl. The
  amino-acid library was written for this project, which is also where its two defects originated.
* `code/solvate.py` is adapted from **[UMADock](https://github.com/MauricioCafiero/UMADock)**'s
  `solvation` class, modified here in four ways, all documented in that file's docstring.
* The older `frag_grow` repository referred to in conversation is **not on this machine** and not on
  the user's GitHub; the only trace is the filename `notebooks/frag_grow.ipynb`.

---

## 8. The next two shells: `tryptophan/2` then `isoleucine/2`

Decided 2026-09-26 with the user. **Do not launch either until the shell-2 run is fully wrapped** —
see the gate at the end of this section.

### Why these two, and not the next-richest

The twenty enumerated shells are **18 distinct**, and nine of those are strict subsets: each
`copies=1` shell is literally the first pass of its `copies=2` counterpart, and the `lysine` and
`aspartic` walks converge to identical selections. So the family offers only nine independent
arrangements plus nested subsets — considerably less variety than "twenty shells" implies.

More important, **the selection criterion is biased by the same artifact that makes the ESM2
energies untrustworthy.** Mean UMA interaction energy per pose, over all 165:

| class | n | mean | best |
|---|---|---|---|
| charged (R, K, D, E) | 67 | −10.79 | −23.29 |
| aromatic (F, W, Y) | 28 | −4.71 | −8.94 |
| neutral aliphatic | 70 | −2.94 | −6.01 |

Charged fragments score 2.3x better than aromatics and 3.7x better than aliphatics on a gas-phase
energy with no desolvation term. Across the twenty runs, IE per fragment correlates **−0.56 with
charged fraction** and **+0.63 with aromatic fraction**. Greedy selection ranks shells by total IE,
so the ranking rewards charge — and shells 1 and 2, the two richest, are the two most
charge-dominated arrangements in the family (58% and 42% charged; shell 1 has a single aromatic among
twelve fragments).

That is backwards for this ligand. Octinoxate is a methoxyphenyl chromophore conjugated to an ester
on a branched C8 chain; a real binder grips the ring and the tail. Both shells tested largely ignore
both. **Tryptophan pose 2, at −8.94, is the best non-charged pose in the whole library, and neither
shell uses it** — both take tryptophan pose 5 at −5.28, because when arginine or phenylalanine picks
first the tryptophan is crowded out.

### Test 1: `tryptophan/2` — composition, at matched size

```
tryptophan:2,tyrosine:0,phenylalanine:4,lysine:1,aspartic:21,glutamic:8,leucine:5,serine:4,tryptophan:5,lysine:8,leucine:13,serine:32
```

`n=12`, total IE −60.40, budget 16, 33% charged, 33% aromatic (F, W, W, Y), overlap 5/12 with shells
1+2 combined. Holds fragment count and spacer budget identical to both shells already run, so
composition moves and size does not. It secures tryptophan's best pose, the untested aromatic site.
Its other aromatics are compromised by mutual clashes — tyrosine lands at rank 4/10 and
phenylalanine at 8/11 — so this is a less charged shell, not an idealised aromatic interface.

It also tests something neither shell could: **whether total fragment IE predicts anything.** Shells
1 and 2 differ by 3.5% (−113.12 against −109.15), far too close to tell. At −60.40 this shell is 47%
weaker. If it folds and wraps comparably, the pipeline is ranking shells by the artifact rather than
the chemistry, which would be the most consequential methodological result here so far.

### Test 2: `isoleucine/2` — maximum independence, lowest charge

```
isoleucine:18,leucine:2,serine:4,tryptophan:5,tyrosine:0,phenylalanine:7,lysine:1,aspartic:18,glutamic:8,isoleucine:2,leucine:5,serine:32,tryptophan:0
```

`n=13`, total IE −57.87, budget 18, **23% charged** (the lowest in the family), 31% aromatic, overlap
only **2/13** — the most genuinely independent arrangement available. Note that size moves too
(13 fragments, budget 18), so it confounds composition with size; run it second, after
`tryptophan/2` has isolated composition at matched size.

### A cheap controlled extra, if wanted

`arginine/1` is a strict subset of shell 1 — the same eight poses as its first pass:

```
arginine:3,lysine:5,aspartic:21,glutamic:9,isoleucine:4,leucine:13,serine:31,tryptophan:5
```

`n=8`, IE −89.73, budget 10. Varies **only** fragment count and spacer budget, with zero
compositional change. Shells 1 and 2 both needed 21 spacers against a budget of 16 and both produced
33-residue peptides that are 64% glycine; that overshoot has never been varied, and it is a leading
candidate explanation for why designed positions are not reproduced.

### The gate: do not start until shell 2 is wrapped

All of these must be true first:

1. `runs/octinoxate/boltz/binding_shell2.csv` and `binding_shell2_forced.csv` between them hold all
   twelve `s2_*` structures, and no scoring process is running.
2. `runs/octinoxate/figures_shell2/manifest.csv` has twelve complete rows — rerun
   `make_figures.py runs/octinoxate --fig-dir figures_shell2 --match s2_ --design-sequence <s2 seq>`
   after the last structure lands, because the pass that ran mid-run only saw three.
3. Shell 2 is written into the README as the second worked example (§4), including the findings in §9.
4. The working tree is committed.

### Naming for these runs

Use `s3_*` for `tryptophan/2` and `s4_*` for `isoleucine/2`, with `condense_shell3.csv` /
`condense_shell4.csv`, `design_shell3.json` / `design_shell4.json`,
`fold_check_shell3.csv` / `_shell4.csv`, `binding_shell3.csv` / `_shell4.csv`, and
`--fig-dir figures_shell3` / `figures_shell4`. `code/run_shell2_forced.sh` is the working template —
it has the correct `--boltz-args=--use_potentials` form. **Do not copy `code/run_shell2.sh`**: it
passes that flag space-separated, which argparse rejects, and it cost nine folds.

The enumeration is reproducible with
`/private/tmp/.../scratchpad/enumerate_shells.py` — if that scratchpad is gone, it is 60 lines
reimplementing `peptide_builder.select_poses` over `energies.csv` and `poses/*.xyz`, walking
`FRAGS[start:] + FRAGS[:start]` for each start at `copies` 1 and 2. Worth moving into `code/` if a
third shell family is ever enumerated.

---

## 9. What shell 2 showed

Run 2026-09-26. Shell is `phenylalanine/2`, the second richest at −109.15 kcal/mol, sharing five of
its twelve poses with shell 1. Sweep found **all 132 ordered pairs connectable**, 4752 rows, as shell
1's did. Energies below are complete for the controls and partial for the forced folds at the time of
writing; finish the table from `binding_shell2.csv` and `binding_shell2_forced.csv`.

### The design replicates shell 1's shape exactly

| | shell 1 | shell 2 |
|---|---|---|
| sequence | `RGGDGGKGGGGLGGGGKGIGEGWGGDGSGGEGS` | `RGEGGEGGKGGFGGDGGIGLGGSGWGGGGLGGS` |
| residues | 33 | 33 |
| spacers / budget | 21 / 16 | 21 / 16 |
| poses in path | 12 / 12 | 12 / 12 |
| worst closure | 0.196 A | 0.120 A |

Two shells sharing only five poses, ordered by independent Held-Karp searches, arrive at the same
length, the same spacer count and the same overshoot. The single-shell caveat in the README's
Limitations now has a second data point.

### Findings

* **Shell 2's unconstrained control is a far better structure than shell 1's.** `s2_orig_control`
  wraps 20/20 ligand atoms with no constraints at all, encloses 0.60, sits 7.7 A from the peptide
  centroid against an Rg of 13.3, and binds at −25.18 kcal/mol. Shell 1's `orig_control` managed
  0.75 wrapped, 0.495 enclosed, 20.9 A separation and −15.83 — the ligand stuck to the outside. So
  shell 1's poor control was a property of that arrangement, not of the method.
* **Designed positions are still not reproduced: 1 of 144** across all twelve folds, at sound ligand
  superposition (RMSD 0.06–0.85 A). Historical baseline was 10/396. But the *distances* improve
  markedly under constraint — median to the designed position falls from ~12–14.5 A unconstrained to
  ~7.6–9.8 A forced, in all three ladders — and the median to the nearest side chain **of any type**
  falls to **3.0–4.9 A**. So the constraints put side chains in roughly the right places with the
  wrong identities there. That points at the ordering and sequence-assembly step, not at
  reachability, which the sweep shows is satisfiable for every pair.
* **`s2_orig_control` wraps completely while reproducing none of the design.** Shell 1 could not
  separate these, because its control failed at both. Wrapping success and design realisation are
  independent.
* **The f8 rung is the damaged one in four of the six ladders**, and `esm2` is the exception in both
  shells. Recomputed against a shared reference conformer (see the note at the end of this section),
  f8 carries the highest ligand strain of both `orig` and both `esm1` ladders at 33.7-42.1 kcal/mol,
  and drives the sum net-unfavourable in three: `orig_f8` +19.76 and `esm1_f8` +7.81 in shell 1,
  `s2_orig_f8` +25.73 and `s2_esm1_f8` +13.83 in shell 2. For `esm2` it is f4 that strains the ligand
  more, in both shells (25.91 against 15.61; 29.32 against 23.30). Shell 2's `orig_f8` is also worst
  of its ladder on enclosure
  (0.43). The eight-contact set specifically buys contacts by bending the ligand.
* **Ranked by sum (interaction + ligand strain), f12 wins in both shells' `esm1` ladders** — shell 1:
  −5.02, +6.27, +7.81, **−5.85**; shell 2: +1.93, −8.12, +13.83, **−19.03**. This is a better
  headline than the enclosure-based claim, because it survives cases where raw interaction energy
  misleads: shell 1's `esm1_f8` has the strongest interaction of its ladder (−29.10) yet is nearly
  net-zero once strain is counted.
* **But for the glycine design the two shells reverse, and this is the most consequential result of
  the run.** Forcing twelve contacts *gains* shell 1 thirty kcal/mol and *costs* shell 2 ten:

  | glycine design | control sum | f12 sum | effect of forcing |
  |---|---|---|---|
  | shell 1 | −1.17 | **−33.07** | gains 32 kcal/mol |
  | shell 2 | **−15.95** | −1.49 | loses 14 kcal/mol |

  It tracks how good the unconstrained fold already was. Shell 1's control had the ligand on the
  outside (0.75 wrapped, centroid 20.9 A), so the constraints had everything to gain; shell 2's
  already wrapped 20/20 at 7.7 A, so they could only disturb it, and the geometry agrees — wrapping
  fell 1.00 to 0.85, interaction −25.18 to −17.56.

  So **forced contacts are a rescue mechanism, not an improvement**: they help when the unconstrained
  fold has failed and hurt when it has already succeeded. The README's "forcing the full set of
  contacts produces the design that was wanted" is true of shell 1's starting point and not in
  general. Rewrite that claim conditionally rather than restating it, and report the control's
  wrapping beside any f12 result so the reader can see which regime they are in.
* **The README's claim that every f12 fold is the most enclosed of its set fails here.** It holds for
  `orig` (0.70) and barely for `esm1` (0.675 against f8's 0.67), but `esm2`'s best-enclosed fold is
  **f8 at 0.86** with f12 at 0.61. `s2_esm2_f8` is the best-enclosed structure in the shell-2 set
  (0.95 wrapped, 19/20 engaged, 31 contacts). Qualify that claim rather than restating it.
* **The ESM2 charge artifact is worse here and now quantified.** The fill took the design from charge
  −1 to **+5** (`RGERKEKSKLIFIIDSKIFLSFSIWRFLILLRS`), and `s2_esm2_control` posts **−74.89
  kcal/mol** — nominally the strongest binding in the project. **Exclude the whole `esm2` set from
  any energy-based conclusion**, on this evidence: across its four folds, interaction energy
  correlates **+0.93 with ligand-to-peptide centroid separation** and only **+0.18 with enclosure**,
  with the wrong sign.

  | fold | centroid | interaction | enclosed | wrapped | contacts |
  |---|---|---|---|---|---|
  | f12 | 5.9 A | −72.72 | 0.61 | 0.75 | 17 |
  | control | 6.7 A | −74.89 | 0.51 | 0.85 | 18 |
  | f8 | 9.5 A | −28.29 | **0.86** | **0.95** | **31** |
  | f4 | 16.0 A | −7.53 | 0.50 | 0.55 | 20 |

  `s2_esm2_f8` has the best geometry in the entire shell-2 set and scores 46 kcal/mol *worse* than the
  control, purely because the control's centroid sits 2.8 A closer. That is long-range electrostatics
  on a +5 peptide, not a binding interface. For contrast, `esm1` — net charge 0 — moves 33 kcal/mol
  across comparable geometry changes and its enclosure ladder is monotonic. n=4, and the control/f12
  pair is effectively tied (0.8 A, 2.2 kcal/mol), but f8 and f4 are unambiguous.

  Constrain the residue set during linker filling, and report net charge beside any interaction
  energy. An earlier reading of this data — that the energy "barely responds" to geometry, from
  control-versus-f12 alone — was wrong: it is not inert, it tracks distance.
* **`esm1` enclosure is monotonic in constraint** (0.375, 0.64, 0.67, 0.675) while its interaction
  energy is not (f8 dips). Geometry and energy diverge inside a single ladder — as do enclosure and
  wrapping in the `orig` ladder, where the control is best wrapped and f12 best enclosed.

### Process notes worth keeping

* `boltz/structures/` was created for the first time by this run. The twelve shell-1 complexes were
  all scored before commit `d28ac25` added the geometry saving, so none of their relaxed structures
  were ever written — they were never written, not lost.
* Scoring runs **40–50 min per structure** on this laptop, not the 28 min the first structure
  suggested. Twelve complexes is most of a night. The A100 route via `boltz_offline.py` is the real
  answer for any further shell.
* A scoring process holds 0.5–2.0 GB, matching the README. `top`'s `MEM` column is not resident size
  and will mislead you; use `ps -o rss`. Output through `tee` is block-buffered, so BFGS steps appear
  in bursts of ~100 — absence of recent log lines is not a stall. Check `ps -o time` against `etime`.

---

## 10. Ligand strain was redefined on 2026-09-26. Read this before quoting any strain number.

Every ligand strain in §9 and in both README worked examples was recomputed. The old definition
relaxed each complex's own bound pose to get the free-ligand energy, which measured every structure
against a **different** local minimum. Three separate defects fed into it, found while validating GPU
scoring on Modal:

1. **The default `fmax 0.1` stops well short of a minimum.** Tightening to 0.01 moved individual
   strains by +1.5 to +7.7 kcal/mol; a sweep on one structure gave 6.96, 8.08, 8.56, 8.61 at fmax
   0.10, 0.05, 0.02, 0.01, converging only by 0.02.
2. **A generous step cap lets the ligand change conformer rather than relax.** The largest shifts came
   with up to 1.0 A of heavy-atom drift -- a rotor flipping. Raising the cap from 75 to 1000 was a
   mistake made while investigating this; it confounds the cutoff with the search.
3. **The bound-state hydrogens came from `pdbfixer`, which is not deterministic here.** Three
   protonations of one CIF on this Mac gave hydrogen positions differing by up to 2.3 A, RMSD 0.78 A,
   with heavy atoms untouched. That alone was worth 3.6 kcal/mol of apparent strain. The same call on
   a Linux container was stable to 0.01 kcal/mol across runs.

### The definition now

`strain_i = E(bound_i) - E_ref`, one shared `E_ref` per run:

* `code/ligand_reference.py` embeds 20 conformers with ETKDG, ranks them with MMFF94, and relaxes the
  five lowest with UMA at fmax 0.01. The lowest is the reference, written to `ligand_reference.json`
  and `boltz/structures/ligand_global_reference.xyz`. For octinoxate it is conformer 18 at
  **-557162.371 kcal/mol-equivalent**, with a 0.996 kcal/mol spread across the five. Note that ETKDG
  returned duplicates (conformers 1 and 11 were identical), so the five kept were fewer than five
  distinct geometries -- raise `--n-conformers` if a strain ever comes back negative, which would mean
  a bound pose below the reference and hence a missed global minimum. None did.
* `code/strain_global.py` recomputes the term for every fold from the ligand's heavy atoms in the
  Boltz CIF plus RDKit hydrogens relaxed with the heavy atoms fixed. **No complex relaxation**, so it
  costs seconds per structure and works on folds whose complex geometries were never saved -- which is
  how shell 1's twelve were corrected despite predating the geometry-saving commit.
* For the 24 folds already scored, `E(bound_i)` was recovered with **no computation at all**: it is the
  step-0 energy of the free-ligand BFGS block in each scoring log, in eV, immediately after the line
  `relaxing free ligand`. Results in `boltz/strain_global.csv`.

### What changed, and what did not

Interaction energy is **unchanged and was never affected** -- it compares the complex against its own
parts at one fixed geometry, so hydrogen-placement error largely cancels. That is why the CPU-vs-GPU
gap was 0.9 kcal/mol on interaction and 4.6 on strain.

Surviving unchanged: f12 wins by sum in `S1 orig`, `S1 esm1` and `S2 esm1`; shell 2's `orig` control
still beats its f12, so forced contacts rescue rather than improve; `esm2`'s energies still track
centroid distance rather than wrap.

Changed: the "five for five" f8 claim was wrong on two counts -- it only counted four ladders, and on
all six f8 leads in four with `esm2` the exception in both shells. Three more structures are now
net-unfavourable (`esm1_f4` +6.27, `esm1_f8` +7.81, `s2_esm1_control` +1.93). `S1 esm2`'s best rung by
sum nominally flips f4 to f8, but at -13.99 against -14.05 they are tied and neither should be called
the winner. Every absolute strain rose by +1.3 to +9.1.

### For shell 3 onwards

Run `ligand_reference.py` once per ligand, then `strain_global.py` after folding. There is no need for
`binding_energy.py` to compute `strain_ligand` at all any more -- dropping that term also removes the
per-structure free-ligand relaxation, which was a meaningful slice of its runtime. `binding_shell*.csv`
keeps the interaction energy; strain comes from `strain_global.csv`.

`boltz/strain_tight_shell2.csv` is superseded. It holds the intermediate experiment -- shell 2's
strains at fmax 0.01 with a 1000-step cap, against per-structure references -- and is kept only as the
record of how the defect was characterised. Do not quote from it.
