# Working rules for this repo

**Start of every session: read `NEXT_STEPS.md`'s "Active right now" block and its "Resume
checklist" first, then `YOUNG.md` (untracked, Young HPC lane) and `RACC.md` (untracked, racc lane)
before touching either cluster.** Live jobs, drains and the next batch to submit are recorded
there, not in this file.

**Use the scripts that already exist; parameterise, don't reimplement.** The established flow for a
docked Modal leg is: `run_dock_pose_md.sh BUILD_ONLY=1` → `modal push` → `probe` → `produce` →
`md_window_modal.sh` → whole-run `omd mmgbsa` → `md_contacts.py` → `md_frames.py` → `pair_contacts.py`,
each invoked directly. The rule bans *reimplementing* a tool's steps, not automation: drivers that
invoke these scripts are fine (`code/racc_feed.sh`, `code/racc_drain.sh`), and orchestration is not
forbidden here. Reimplementing a driver instead is how four purpose-written drivers (2026-10-01,
499 lines) silently dropped three of `run_dock_pose_md.sh`'s steps across nine legs.

## Dynamics runs

**Use `code/run_dock_pose_md.sh`. Do not write a new driver.** It already handles every step
correctly: stage guards so a re-run resumes, MM/GBSA invoked from inside the leg directory, the
leading convergence windows, `caffeinate`, `md_contacts.py` and `md_frames.py` at the end. If it
hardcodes something you need to vary — it sets `N=s3_orig_f12_dock1` — that is a parameter to add,
not a reason to start over. **Never "fix" a failed leg by nudging the ligand**: reproduce the
failure locally with the real tool first — a tight 1.69 Å contact that looked like a broken build
was a deeply enclosed pose whose resubmission produced the campaign's strongest result.

**`omd build` needs `--box-shape dodecahedron` explicitly.** It defaults to a cube and the
co-folded baselines in `runs/octinoxate/md/` are dodecahedral. Same physics and no effect on
MM/GBSA (it strips to the solute), but ~30% more waters — invisible on compact solutes, real on
elongated ones (one control cost $0.41 and 24 minutes extra on one leg). Switching to dodecahedral
*restores* comparability with the baselines, so it is always safe to fix.

**`run_md20.sh` is not in this repo.** It lives in the sibling BoltzGen repo at
`~/python_mac/boltzgen_local/md/run_md20.sh`, which is where `racc_drain.sh` calls it from. Only
`run_dock_pose_md.sh` is in `code/`; a change to the co-folded tail lands in a different repo and
git history.

**Figures are a pipeline step, with their own environment.** `code/render_md_figures.py` renders a
leg's `frame_medoid.pdb` headlessly at the end of the tail, so every leg gets its figure as it
lands. It needs **PyMOL, which is not in the MD env** — its own conda env because the
`pymol-open-source` PyPI wheel links to libraries that cannot load here:
```sh
~/miniforge3/bin/mamba create -y -n pymol-render -c conda-forge python=3.12 pymol-open-source
~/miniforge3/envs/pymol-render/bin/pip install pillow   # alpha-channel cropping
```
mdtraj writes trajectory PDBs with **no HELIX/SHEET records**, so PyMOL loads every residue as coil
and a cartoon draws a helix flat — the script runs `dss` first to assign secondary structure from
geometry. The committed figures are **transparent RGBA**: why `figures_shell*/style.pml` says
`bg_color black` while its own output reads white.

**`md_contacts.py` and `md_frames.py` are pipeline steps, not extras.** Contacts give residence and
separation, which a scalar ΔG cannot and which separated every design from every null; frames give
the end-of-run ensemble. A leg without them is incomplete.

**Record input → production ligand RMSD on every leg** with `code/ligand_slide.py <prod dirs...>`
(`--csv` to write the table; it superposes on protein Cα first, so the number means ligand movement
relative to its peptide, not box translation). Equilibration restrains protein heavy atoms only
(`openmm_md/dynamics.py:144`), so the ligand slides 1.2–15.1 Å before production. **Do not "fix"
this by restraining the ligand** — every run equilibrated it freely and a restrained leg would not
be comparable. Record the number; it says whether a cell is worth reading. **A second molecule must
pass `SMILES` to `boltzgen_local/tools/metrics.sh`** — it defaults to the octinoxate analogue;
perceive it from the run's own geometry rather than typing it:
`boltz_check.ligand_smiles('runs/<mol>')` reads `ligand.xyz`.

## Modal

**Launch every `modal run` through `code/modal_run.sh` — it is not optional.** The client stays
attached for the whole leg even with `--detach`, and when the Mac sleeps the client drops and Modal
cancels the in-flight leg at full price. The script `exec`s the run under `caffeinate -i -s`, which
holds idle *and* system-sleep assertions for exactly the client's lifetime (a `-i`-only guard lost
one leg to a 16-minute maintenance dark sleep). There is no warm restart: `omd run` writes
`checkpoint.chk` but nothing loads it. The guard protects a *launch*, not a day of work: after any
sleep-inhibitor release (e.g. an overnight laptop sleep), verify one is armed with
`pgrep -f caffeinate` before launching anything.

**One Modal leg at a time** — verify every stage of it (probe rate, fetched files, MM/GBSA row)
before launching the next.

**When a run cancels, pull both sides' logs before relaunching:** `modal billing report --for
today --show-resources` and the `pmset -g log` sleep/wake timeline (this Mac does normal 45 s–2 min
dark wakes for maintenance/keepalive that look like noise but kill an unguarded client). Diagnose,
then relaunch — a guessed relaunch into the same failure pays for the leg twice.

**Leg progress is unobservable mid-run.** The volume's `energy.csv` is a stale early flush and the
container prints nothing between start and finish; any percent-complete figure is a guess. What
*is* observable: accrued spend per app from the billing report, and `produce`'s own `minutes`/`usd`
on completion.

**Preempted time is billed** — a preempted leg restarts from step zero, so it pays twice.

**Cost model, measured on A10G:** $0.45 fixed per leg plus $0.041 per 1,000 particles; the fixed
half is ~53% of a typical leg, so atom-count savings deliver less than the atom ratio suggests.
Rate scales as roughly N^0.77 (an A10G number — **do not apply it to racc**, which is much faster).
Read the particle count off the build before quoting a time.

**Modal returns a strided trajectory** (stride 10, 2,000 frames): MM/GBSA on it takes ~2 minutes
where an unstrided local leg takes ~8, with SE ~0.07 against ~0.04 — same 20 ns, uniform subsample,
so the mean is unbiased. The full `.dcd` stays on the volume if more frames are ever wanted.

## Reading ARC cluster (racc)

The second compute lane, for legs Modal would bill. **Everything cluster-side lives in the
untracked `RACC.md` — access recipe, measured H100 rate model, queue discipline (`racc_feed.sh` is
the only queue filler; `SBATCH_OPTS --time=est+1h30 --cpus-per-task=8` on every submission; 3 jobs
always in queue), the drain loop and fetch rules, scratch/dcd policy, and how a campaign is driven.
Read it before submitting anything.** Repo-side notes:

- Build locally (`run_dock_pose_md.sh BUILD_ONLY=1`, or `boltzgen_local/md/prep.sh` for a co-fold),
  run remotely; stage by tar-over-ssh into `~/remote_work/openmm/runs/<mol>/md/<leg>/` mirroring
  the repo layout (rsync and scp -r fail through the tunnel).
- A leg that dies on the cluster is invisible from the queue — sacct stays clean while the drain
  waits forever; the checks are in RACC.md.
- The drain's tails and `code/record_leg.py` feed the folds_md record pipeline below; the
  reporting cadence (report every leg, commit every ~5) is under "Reporting results".

## folds_md records & the whole-table scripts

**The cofold tables are a three-step chain.** `code/cofold_tables.py --write` splices Tables 1 and 2
into `runs/folds_md/cofold_binding_tables.md` itself and **refuses if any published cell would
change** (only cells in newly added columns and the GNINA-pick row label may change). Without
`--write` it only prints a `missing / pending` report and writes table blocks to `$OUT`. The
markdown is the source of truth; `code/render_cofold_pdf.py` rebuilds `cofold_binding_tables.pdf`
from it — **with the venv or conda python, not bare `python3`, which has no matplotlib**; skip that
step and the PDF rots invisibly. `cofold_tables.py` resolves each leg's residence from
`runs/folds_md/record_<leg>.log` first, prose tables second — doc-scanning alone once harvested a
static fold-geometry percentage and printed it as MD residence.

**`record_leg.py`'s own `audit: ... -> OK` is not the done-signal.** It checks windows, summary
rows, frames and scatter, and nothing else — run the four-line audit under "Audit the recording
set" against the leg. It has twice mis-resolved pair tables (now it globs `pair_contacts*.csv`,
including the bare base-lane file, and takes the design from the parent's own row rather than the
table's first row); **after any pair backfill, check each new row's `design` against its run-1
row's**.

**A fold lane needs a `runs/folds_md/dock_legs.txt` entry per docked leg before it is recorded.**
That file is where FOLDS.md's pose label comes from; a missing entry now writes `UNREGISTERED` and
warns (older ones defaulted to `GNINA`).

**Whole-table scripts rewrite their tables: name every structure, not just the new one.**
`vina_redock.py` (dock_summary.csv, dock_poses.csv), `make_gnina_bundle.py` (`reference.csv`),
`md_stability.py` (`protein_stability.csv`), `ligand_slide.py --csv`
(`runs/{octinoxate,oxybenzone}/md/ligand_slide.csv`) and `pair_contacts.py --csv`
(`runs/<mol>/md/pair_contacts_<lane>_<arm>.csv`) all write mode `w` over only the structures named
— naming only the new ones silently drops every other row. None warns (except `md_stability.py`'s
shrink guard). The rules:

- Pass **every lane's globs** to `md_stability.py` and `ligand_slide.py` — the BoltzGen baselines
  live in `~/python_mac/boltzgen_local/md/*/prod_20ns`, so a run globbing only `runs/<mol>/md/*`
  drops them; oxybenzone has `prod_L1_modal` legs too, and `boltzgen_local` has no `prod_L1_modal`.
  Note zsh aborts the whole launch if one glob matches nothing.
- `md_stability.py` takes the leg list from argv — a **no-argument run writes nothing**, and a
  single-leg pass *overwrites* `protein_stability.csv` with one row; it refuses to write fewer rows
  than the file holds unless `--allow-shrink`.
- `pair_contacts.py`: **`--run-dir runs/<mol>` is required and has no default** — it chooses whose
  design json defines the slots; `--design design.json` resolves against `--run-dir`. A default
  would score one molecule's legs against another molecule's slots.
- **The route that avoids the damage:** `vina_redock.py --out` docks to a scratch directory and
  rows are appended to the real tables; verify with `git diff --numstat` that the change is
  additions-only. `make_gnina_bundle.py` has no `--out` (hardcoded `out_root`), so back
  `reference.csv` up, build only the new bundles, restore, and append. Re-docking everything is
  *not* the fix — it rewrites files belonging to done and committed legs (a re-dock with seed 42 is
  byte-identical but still rewrites; `--all` silently pulled in an eighth structure never in the
  docking matrix). **Run `git diff --numstat` on every leg's commit** — additions-only is the only
  real proof no whole-table rewrite shrank a table.

**`vina_redock.py` must run with `~/python_mac/dock_assist/dock-env/bin/python`**, never
peptidebuilder's `.venv` — it locates Vina through `import dockstring`, installed only in `dock-env`
(the `DockError`'s advice to "pass --vina-bin" is a dead end; the script does not expose that
flag). `dock-env` carries rdkit and numpy, so the script runs there unchanged.

## MM/GBSA housekeeping

**Strip the scatter — `run_dock_pose_md.sh` does not.** After it runs `omd mmgbsa` for the whole
run and for each window, strip with
`find <leg> \( -name 'reference.frc' -o -name '_MMPBSA_*' -o -name 'traj.nc' \) -delete`
(`traj.nc` is MMPBSA.py's NetCDF conversion of the kept `traj_wrapped.xtc`, regenerates in
seconds; use `find -delete`, never a multi-path `rm` whose zsh globs abort the whole command when
one pattern matches nothing). Window `traj_wrapped.xtc` slices are not kept either — a slice is two
seconds of work off the kept `prod_20ns`. **Run MM/GBSA from inside the leg directory**: `cd` into
the leg in the same command that launches it, never trust the session's cwd — the scratch (~300 MB
a leg; a 1.2 GB `reference.frc` on big systems) lands wherever the shell is, and it has landed in
the repo root. Verify the strip is complete before staging a leg for commit.

**Keep what the September legs keep:** no `reference.frc`/`_MMPBSA_*` anywhere;
`FINAL_RESULTS_MMPBSA.dat`, the prmtops, `energy.csv`, `traj_wrapped.*`, `frame_medoid.pdb` and
`frames_last.pdb` stay. Window directories are `first_Nns`; window rows come from the
`Differences (Complex - Receptor - Ligand)` block of the `.dat`, not the complex section (a naive
`VDWAALS`/`EEL` grep returns the complex's −hundreds).

**Every result goes into its own molecule's `runs/<mol>/md/mmgbsa_summary.csv` as it is computed**,
not in a batch at the end. Pass **both** of `collect_mmgbsa_rows.py`'s paths —
`--md-root runs/<mol>/md --outdir runs/<mol>` — they do different jobs (where legs are *read* vs
which summary is *written*); `--outdir` is required and has no default. Nothing reads these files —
`make_gnina_bundle.py` carries the same numbers in a hardcoded `REFERENCE` dict — so when a
peptide's docked legs land, extend its `REFERENCE` entry with the flat `vina_p1_*`/`gnina_p9_*`
fields.

**Audit the recording set per leg; do not tick it off from memory.** One query over the finished
leg:

```
leg dir: <n>/4 FINAL_RESULTS_MMPBSA.dat, md_contacts.csv, prod_20ns/{frames_last,frame_medoid}.pdb,
         and no reference.frc / _MMPBSA_* / traj.nc anywhere
tables : 4 rows in runs/<mol>/md/mmgbsa_summary.csv, 1 in protein_stability.csv,
         1 in runs/<mol>/md/ligand_slide.csv, 1 in a pair_contacts_<lane>_<arm>.csv
```

Match `pair_contacts_*.csv` on its **`structure`** column, not its first column — `leg` holds the
full production path and matching that reports every row missing. The two whole-table rewrites key
differently: `protein_stability.csv`'s `leg` is the full production path, `ligand_slide.csv`'s the
**bare leg name** — grepping both for `/<leg>/` reports missing rows that exist. **A leg is not
done until `protein_stability.csv` has its row.**

## GNINA rescoring (the `gnina_local` repo) — retrieving results

GNINA CNNaffinity rescoring runs OUTSIDE this repo, through the user's **public** GitHub repo
**`github.com/MauricioCafiero/gnina_local`** — not a local directory, not `dock_assist`, and must
not be cloned. Flow: `make_gnina_bundle.py --system <mol>` writes `runs/<mol>/gnina/` (per
structure `receptor.pdb`, `receptor_h.pdb`, `poses_with_reference.sdf`, `ref_pose.sdf`, plus
`reference.csv`), the bundle is committed to scavengers, the user runs GNINA on the gnina_local
side, and the results are committed back into `gnina_local` — so the recurring task is
*retrieving* them. Follow this verbatim:

- It is a **public repo the user has told you to look in**: navigate it directly with the GitHub
  API / raw URLs. Do not ask where it is, and do not treat exploring it as a banned "search" —
  that rule is about the local filesystem only. **Do not clone it**; `gh` is not installed, use raw
  URLs.
- Layout: `rescoring/<molecule>/results_<mode>/`, `<mode>` ∈ {`receptor_h_score_only`,
  `receptor_score_only`, `receptor_minimize`}. Each dir holds `gnina_scores.csv` (the result),
  `gnina_version.txt`, and per-structure `<name>.log` + `<name>_scored.sdf`.
- If unsure which files/modes to pull, mirror the previous molecule's `runs/<prev>/gnina/`
  (octinoxate, the template, keeps all three modes, each a subdir `runs/<mol>/gnina/results_<mode>/`).
  Pull every structure's `*_scored.sdf` into those subdirs — the script below needs the SDFs, not
  just the CSV:
  `curl -sf https://raw.githubusercontent.com/MauricioCafiero/gnina_local/main/rescoring/<mol>/results_<mode>/<struct>_scored.sdf -o runs/<mol>/gnina/results_<mode>/<struct>_scored.sdf`
- **Then USE the existing script — do NOT hand-roll the ranking or the pick.** From
  `runs/<mol>/gnina/`, run `collect_gnina.py` (it lives at `runs/octinoxate/gnina/collect_gnina.py`)
  once per mode: `python <path>/collect_gnina.py results_<mode> --reference reference.csv`. It
  reads the scored SDFs, (re)writes `results_<mode>/gnina_scores.csv`, and prints where GNINA ranks
  the predicted pose 0 by each score. `gnina_scores.csv` IS this script's output, not a file to
  parse by hand.
- **Primary mode is `receptor_score_only`** — heavy-atom receptor, `--score_only` (CNN scores the
  exact pose, no minimisation); the other two are second passes. Columns include `CNNscore`,
  `CNNaffinity`, `minimizedAffinity`.
- **The GNINA-pick pose** (queued as a docked MD leg AFTER the Vina-p1 legs) is the docked pose
  (pose ≥ 1) with the **max `CNNaffinity`, consistent across all three passes**, corroborated by
  `minimizedAffinity` (argmin). It becomes `<struct>_dock<K>` via `run_dock_pose_md.sh
  STRUCT=<struct> POSE=<K>`; if the pick is pose 1 it coincides with the Vina-p1 leg, so no new
  leg. Feeds the `gnina_p*` REFERENCE fields.

## Downstream analysis repos (`pepgraph`, `morph`) — pipeline steps, not side projects

Two sibling repos consume `runs/<mol>/folds/`. **A new fold lane is not finished until both have
been extended and re-run**, the same way a leg is not done until `protein_stability.csv` has its
row. Neither is discoverable from this repo, so they are listed here.

**`~/python_mac/pepgraph`** — inverse-folding perplexity per lane: given a backbone, how much does
its geometry want the sequence it was given (21 = no information), plus `z` against 20 shuffles of
the same sequence (`z ≤ 0` = the backbone carries no residue-order information). Add the lane to
`LANES` in `scripts/compare_lanes.py`, then
`.venv/bin/python scripts/compare_lanes.py --runs ~/python_mac/peptidebuilder/runs --lanes <lane>`.
**Pass `--lanes`** — a full run re-scores all 202 structures and, because the shuffles come off one
RNG stream consumed in structure order, re-draws every other lane's z and moves it by a few tenths.
Its `RESULTS.md` and `outputs/inverse/*.json` are **gitignored** — a rerun that rewrites them
cannot be recovered from git; copy them aside first.

**`~/python_mac/morph`** — apo → co-fold morph movies, MM strain (ff14SB + GBn2) and
charged/hydrophobic surface-patch changes. `$PY code/run_batch.py` with the conda `openmm-md`
python; `--set matched` (the default) runs the model-matched pairs `boltzapo→boltz` /
`of3apo→openfold3` / `rf3apo→rf3`, so a ΔE is the ligand's effect inside one model. A new apo lane
means a new pair in `MATCHED_PAIRS`. `--set legacy` reruns the superseded all-apo→one-Boltz
analysis into `results/legacy/`; nothing should be concluded from it, and ESMFold exists only
there.

**ESMFold is retired from downstream analysis (2026-10-11)** on pepgraph's evidence: last lane at
mean perplexity 20.4, the no-information limit, owning most of the campaign's zero-information
structures. Its folds and docks stay committed; do not build anything new on them, and do not
bundle its docks for GNINA.

## Reporting results

**Lead with MM/GBSA free energy and ligand retention, grouped by pipeline.** That is the
comparison that matters here. Pose-selection questions, scorer correlations and protocol
diagnostics are methodology and belong in their own labelled section.

**Report every leg as it lands; commit every ~five.** Two different cadences — collapsing them is
a mistake in either direction. The report is the scientific output and is owed immediately: it
needs nothing committed, because `FINAL_RESULTS_MMPBSA.dat` and the leg's tables are already on
disk when the tail finishes. **Report off the MM/GBSA, not the record** (waiter `--on mmgbsa`) —
the dG, windows and release (from `tail.log`) do not wait for `record_leg.py`; slide and the tally
follow when the record lands. The commit is bookkeeping: batch it (~5 legs, or a natural boundary
— a legfile finishing, a campaign ending, stopping for the day) and run `git diff --numstat` on
each.

**MM/GBSA standard errors are not uncertainty on a comparison.** They measure frames within one
trajectory (0.04–0.10 kcal/mol); run-to-run variation is far larger. Never write "N standard
errors" as though it bore on reproducibility.

**"Residence" and "release" are different measures, and only release is shape-independent.** In
`md_contacts.py`, `residence_pct` is the fraction of frames whose ligand-COM-to-peptide-COM
distance is ≤ 10 Å, while `released_pct` is the fraction with zero atom–atom contacts under 4 Å.
A ligand held firmly at one end of an elongated peptide is never released yet scores low on
residence, because the far half of the rod drags the peptide's COM away — **a COM-based residence
number penalises a fold for being elongated**, and the unforced co-folds are rods. Quote release
(plus episode count, longest episode and the late share) when retention is the claim; use
residence only within one shell, where the geometry is comparable.

## Paths

Scripts derive the repo root from their own location (`REPO=${0:A:h:h}`) and read `OMD_ENV`,
`BOLTZGEN_ROOT` and `DOCK_ASSIST`, defaulting to `~/miniforge3/envs/openmm-md`,
`~/python_mac/boltzgen_local` and `~/python_mac/dock_assist`. The MD environment **must be conda**:
`omd build` needs OpenFF/GAFF2 and `omd mmgbsa` needs AmberTools, and `.venv` has no `mdtraj`.