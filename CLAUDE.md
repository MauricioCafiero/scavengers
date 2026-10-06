# Working rules for this repo

Written 2026-10-02 after a day in which every rule below was learned by breaking it.

**Use the scripts that already exist; parameterise, don't reimplement.** The established flow for a
docked Modal leg is: `run_dock_pose_md.sh BUILD_ONLY=1` → `modal push` → `probe` → `produce` →
`md_window_modal.sh` → whole-run `omd mmgbsa` → `md_contacts.py` → `md_frames.py` → `pair_contacts.py`,
each invoked directly. On 2026-10-04 an orchestration driver for this was drafted twice and rejected
twice — the tools each hold one step of the tacit knowledge, and gluing them into a new driver is where
steps silently drop.

## Dynamics runs

**Use `code/run_dock_pose_md.sh`. Do not write a new driver.** It is the one MD driver and it already
handles every step correctly: stage guards so a re-run resumes, MM/GBSA invoked from inside the leg
directory, the leading convergence windows, `caffeinate`, `md_contacts.py` and `md_frames.py` at the end.
If it hardcodes something you need to vary — it sets `N=s3_orig_f12_dock1` — that is a parameter to add,
not a reason to start over. On 2026-10-01 four purpose-written drivers (499 lines) were written instead,
and they silently dropped three of its steps across nine legs.

**`omd build` needs `--box-shape dodecahedron` explicitly.** It defaults to a cube, and six of the seven
co-folded baselines in `runs/octinoxate/md/` are dodecahedral. Same physics and no effect on MM/GBSA,
which strips to the solute — but ~30% more waters. The waste is invisible on compact solutes ($0.09 a
leg at 7,384 atoms) and real on elongated ones: `s2_esm2_control` built to 33,057 atoms cubic against
23,142 dodecahedral, costing $0.41 and 24 minutes on one leg. Switching to dodecahedral *restores*
comparability with the baselines, so it is always safe to fix.

**`md_contacts.py` and `md_frames.py` are pipeline steps, not extras.** Contacts give residence and
separation, which a scalar ΔG cannot and which separated every design from every null. Frames give the
end-of-run ensemble, needed because the predicted pose and the pose after 20 ns are different
structures. A leg without them is incomplete; compare against `s3_orig_f12_dock1`, which has the full
file set.

**Record input → production ligand RMSD on every leg.** Equilibration restrains protein heavy atoms only
(`openmm_md/dynamics.py:144`), so the ligand slides 1.2 to 15.1 Å before production starts depending on
the run. **Do not "fix" this by restraining the ligand** — every run in the project equilibrated it
freely and a restrained leg would not be comparable. Record the number; it is what says whether a cell
is worth reading.

## Modal

**Launch every `modal run` through `code/modal_run.sh` — it is not optional.** The client stays
attached for the whole leg even with `--detach`, and when the Mac sleeps the client drops and Modal
cancels the in-flight leg at full price. The script `exec`s the run under `caffeinate -i`, which holds
an idle-sleep assertion for exactly the client's lifetime, so the guard is self-cleaning and covers the
whole leg (2026-10-04: two cancelled legs, $1.06, for skipping this). There is no warm restart: `omd run`
writes `checkpoint.chk` but nothing loads it. The guard protects a *launch*, not a day of work: any
instruction to release sleep inhibitors is scoped to the work that motivated it — before resuming any
work after a release (e.g. an overnight laptop sleep), verify one is armed with `pgrep -f caffeinate`
before launching anything.

**One Modal leg at a time.** Run one leg, verify every stage of it (probe rate, fetched files, MM/GBSA
row) before launching the next — a repeated error across parallel legs doubles the loss, and 2026-10-04
cost $1.06 the exact way single-leg checking exists to prevent.

**When a run cancels, pull both sides' logs before relaunching.** The service side (`modal billing
report --for today --show-resources`) and the machine side (`pmset -g log` sleep/wake timeline — this
Mac also does normal 45 s–2 min dark wakes for maintenance and TCP keepalive, which look like noise but
are what kill an unguarded client). On 2026-10-04 billing alone left the cause "unexplained" for an
hour; the power log dated the cancellations to the exact sleep events. Diagnose, then relaunch — a
guessed relaunch into the same failure pays for the leg twice.

**Leg progress is unobservable mid-run.** The volume's `energy.csv` is a stale early flush and the
container prints nothing between start and finish. Any percent-complete figure is a guess; say so rather
than offering scaled arithmetic as data. What *is* observable: `modal billing report --for today
--show-resources` gives accrued spend per app, and `produce` prints its own `minutes` and `usd` on
completion.

**Preempted time is billed.** A preempted leg restarts from step zero — `omd run` writes a checkpoint
but nothing loads it — so it pays twice. One of nine legs on 2026-10-01 cost $1.90 against ~$0.88 for a
clean one.

**Cost model, measured on A10G:** $0.45 fixed per leg plus $0.041 per 1,000 particles; the fixed half is
~53% of a typical leg, so atom-count savings deliver less than the atom ratio suggests. Rate scales as
roughly N^0.77, not linearly — 0.264 ms/step at 14,105 particles, 0.362 at 23,142, 0.507 at 33,057. Read
the particle count off the build before quoting a time.

**Modal returns a strided trajectory** (stride 10, 2,000 frames). MM/GBSA on it takes ~2 minutes where
an unstrided local leg takes ~8, and the standard error is ~0.07 against ~0.04. Same 20 ns, uniform
subsample, so the mean is unbiased. The full `.dcd` stays on the volume if more frames are ever wanted.

## Reading ARC cluster (racc)

The second compute lane, for legs Modal would bill. RACC.md (untracked) holds the access recipe
(two user-typed ControlMaster commands), paths, quotas and the account facts; the workflow here:

**Stage with tar, drive with the socket.** Every command is `ssh -o ControlPath=~/.ssh/cm-racc
racc.rdg.ac.uk '<cmd>'` and must be wrapped in a login shell (`bash -lc "module load anaconda;
source activate openmm-md; ..."` — non-interactive ssh has no module init). Send legs with
tar-over-ssh into `~/remote_work/openmm/runs/octinoxate/md/<leg>/` mirroring the repo layout; rsync
and scp -r both fail through the tunnel. Build locally (`run_dock_pose_md.sh BUILD_ONLY=1`), run
there: `submit_gpu.sh racc_run.py` from inside an activated env (sbatch exports it; a
`python3.12 -> python` symlink in the env covers the script's hardcoded interpreter)

**Platform is OpenCL, and the first run's wall time is kernel compile.** The conda-forge OpenMM
build has no CUDA platform (user declined the CUDA-variant env 2026-10-04); OpenCL works on the
H100 NVL nodes. A smoke run's total wall time is JIT-dominated — 100k steps cost 6:47 while the
live rate is 0.45 ms/step — so **measure the rate from the live `energy.csv`, never from a smoke's
wall time**, and expect the real leg to be several times faster than the smoke suggested. Measured
across six oxybenzone null legs 2026-10-06: **0.29 ms/step at 5,425 particles** (10M steps, 47:57
wall including equilibration), so a compact 20 ns leg is 46–48 min and a 12,800-particle one 56–66.
The 0.45 figure came from one earlier leg and is pessimistic; quote whichever matches the particle
count, and read that count off the build.

**energy.csv flushes continuously here, unlike on Modal** — mid-run progress and ms/step are
readable any time (row cadence 500 steps) without contacting produce. Use it to revise estimates
within the first 20 minutes instead of waiting for completion to find out.

Waits go in background watchers that poll sacct (`--format=State --noheader | tail -1`), never
foreground sleep loops. Partition is `gpuscavenger`: free, preemptible, and — like on Modal — a
preempted leg restarts from zero, so a long leg is not protected; check `sacct` before assuming one
is still the run it was.

**Large trajectories are written to scratch from the first byte.** The dcd is 1.8–5 GB per leg (it
carries every particle) and the home quota filled mid-write on 2026-10-04, killing both legs at
once — so `racc_run.py` points `--out-dir` at a scratch directory through a symlink in the leg dir
(`prod_20ns -> /scratch5/gaussian/io927423/dcd/<leg>_prod`): the dcd lands on scratch from frame
one, home never carries it, and every path reading `<leg>/prod_20ns` stays unchanged. A move-at-end
was tried first and is wrong: it saves nothing from the kill that happens mid-write.

**Strip to the solute on the cluster; fetch only the wrapped solute, never the full dcd (racc mirrors
Modal).** After `omd run`, `racc_run.py` runs the same `omd analyze` the local tail runs to write
`traj_wrapped.{xtc,pdb}` (solute-only, ~100 MB) beside the full dcd on scratch; the full `traj.dcd`
stays on scratch. Fetch `dcd/<leg>_prod/traj_wrapped.{xtc,pdb}` + `energy.csv` + the small analyze
outputs — **not** the 5 GB dcd. The Mac has 8 GB RAM and mdtraj loads a trajectory whole, so fetching
the full box and stripping locally OOM'd the analyze at ~4–5 GB. With only the wrapped solute local,
run_dock_pose_md.sh's whole-run analyze guard skips and its window loop slices the wrapped solute
(no full-box load) — the Modal legs always worked this way, the racc lane just hadn't.

**Whole-table scripts rewrite their table: name every structure, not just the new one.**
`vina_redock.py` and `make_gnina_bundle.py` write `dock_summary.csv`, `dock_poses.csv` and
`reference.csv` wholesale, so a run naming only the new structures drops every other row — the same trap
`md_stability.py` has. On 2026-10-06 docking two nulls cut `dock_summary.csv` from five structures to
two. The repair is `git checkout` the tables and dock only the new structures, **not** re-docking
everything: a re-dock is reproducible (seed 42 gave byte-identical rows) but it rewrites files belonging
to legs whose GNINA and MD are already done and committed. `--all` is not the fix either — it silently
pulled in an eighth structure that had never been in the docking matrix. Name the targets.

## MM/GBSA housekeeping

**`run_dock_pose_md.sh` does not strip its own scatter — you must.** It runs `omd mmgbsa` for the whole
run and for each window and leaves every `reference.frc` and `_MMPBSA_*` in place; only
`md_window_modal.sh` carries the `rm`. On 2026-10-06 four docked legs tailed through the driver were
holding **24 GB** (5.0–7.2 GB each, against 294–427 MB for two co-folded legs tailed through
`md_window_modal.sh`). Strip with `find <leg> \( -name 'reference.frc' -o -name '_MMPBSA_*' \) -delete`
— `find -delete` rather than a multi-path `rm`, whose zsh globs abort the whole command when one pattern
matches nothing — and re-check `FINAL_RESULTS_MMPBSA.dat` afterwards. **Window `traj_wrapped.xtc` slices
are not kept either**: the committed legs hold only `first_Nns/mmgbsa/` plus `traj_wrapped.pdb`, because
a slice is two seconds of work off `prod_20ns`, which is kept. Leaving them in put 250 MB a leg into the
staging area.

Run it from inside the leg directory: MMPBSA.py scatters `reference.frc` plus a `_MMPBSA_*` set into the
working directory, ~300 MB a leg, and `reference.frc` reached 1.9 GB once when run from the repo root.
**The shell's cwd at launch decides where the scratch lands** — on 2026-10-04 leg 2's whole-run
MM/GBSA ran with the shell at the repo root and scattered a 272 MB `reference.frc` there even though
`--out-dir` pointed into the leg. `cd` into the leg directory in the same command that launches it, do
not trust the session's cwd. Before staging a leg for commit, verify the strip is complete with
`find <leg> repo-root -name 'reference.frc' -o -name '_MMPBSA_*'` — a zsh glob that matches nothing
silently aborts a multi-path `rm` and leaves the old files in place.
**Keep what the September legs keep:** no `reference.frc` and no `_MMPBSA_*` anywhere;
`FINAL_RESULTS_MMPBSA.dat`, the prmtops, `energy.csv`, `traj_wrapped.*`, `frame_medoid.pdb` and
`frames_last.pdb` stay. Window directories are named `first_Nns`; window legs get `first_Nns` rows in
the summary CSV with their full component set (parse the `Differences (Complex - Receptor - Ligand)`
block, not the complex section — a naive grep of `VDWAALS`/`EEL` returns the complex's −hundreds).

Every result goes in `runs/octinoxate/md/mmgbsa_summary.csv` **as it is computed**, not in a batch at the
end. Nothing reads that file — `make_gnina_bundle.py` carries the same numbers in a hardcoded `REFERENCE`
dict — so a new row has to be added in both places or the next GNINA bundle is built on stale values.
When a peptide's docked legs land, extend its `REFERENCE` entry with the flat `vina_p1_*`/`gnina_p9_*`
fields, as both shuffle peers now have.

A leg is not done until `protein_stability.csv` has its row: `md_stability.py` takes the leg list from
argv, so a **no-argument run writes nothing** — pass the full glob
(`runs/octinoxate/md/*/prod_20ns runs/octinoxate/md/*/prod_L1_modal ~/python_mac/boltzgen_local/md/*/prod_20ns`),
never a single leg (a single-leg pass *overwrites* the file with one row). zsh also aborts the whole
launch if one glob in the list matches nothing (`boltzgen_local` has no `prod_L1_modal` legs).

## GNINA rescoring (the `gnina_local` repo) — retrieving results

GNINA CNNaffinity rescoring runs OUTSIDE this repo, on a GPU, through the user's **public** GitHub repo
**`github.com/MauricioCafiero/gnina_local`** — it is not a local directory, not `dock_assist`, and must
not be cloned. The flow: `make_gnina_bundle.py --system <mol>` writes `runs/<mol>/gnina/` (per structure:
`receptor.pdb`, `receptor_h.pdb`, `poses_with_reference.sdf`, `ref_pose.sdf`, plus `reference.csv`) and
that bundle is committed to scavengers; the user runs GNINA against it on the gnina_local side; the
**results are committed back into `gnina_local`**, not here. So the recurring task is *retrieving* those
results.

**How to retrieve (this step has repeatedly cost me turns — follow it verbatim):**
- gnina_local is a **public repo the user has told you to look in**: navigate it directly with the
  GitHub API / raw URLs. Do NOT ask where it is, and do NOT treat exploring it as a banned "search" —
  the never-search-the-disk rule is about the *local* filesystem only, not a public repo you were pointed at.
- **Do not clone it.** `gh` is not installed; use raw URLs.
- Layout: `rescoring/<molecule>/results_<mode>/`, `<mode>` ∈
  {`receptor_h_score_only`, `receptor_score_only`, `receptor_minimize`}. Each dir holds
  `gnina_scores.csv` (the result), `gnina_version.txt`, and per-structure `<name>.log` + `<name>_scored.sdf`.
- **If unsure which files/modes to pull, mirror the previous molecule's `runs/<prev>/gnina/`.** Octinoxate
  (the template) keeps **all three modes**, each a subdir `runs/<mol>/gnina/results_<mode>/` holding the
  per-structure `*_scored.sdf` files and a `gnina_scores.csv`. Pull every structure's `*_scored.sdf` into
  those subdirs (the script below needs the SDFs, not just the CSV):
  `curl -sf https://raw.githubusercontent.com/MauricioCafiero/gnina_local/main/rescoring/<mol>/results_<mode>/<struct>_scored.sdf -o runs/<mol>/gnina/results_<mode>/<struct>_scored.sdf`
- **Then USE the existing script — do NOT hand-roll the ranking or the pick** (project rule 1: use the
  scripts that exist, don't reimplement). From `runs/<mol>/gnina/`, run `collect_gnina.py` (it lives at
  `runs/octinoxate/gnina/collect_gnina.py`) once per mode:
  `python <path>/collect_gnina.py results_<mode> --reference reference.csv`. It reads the scored SDFs,
  (re)writes `results_<mode>/gnina_scores.csv`, and prints where GNINA ranks the predicted pose 0 by each
  score. `gnina_scores.csv` IS this script's output, not a file to parse by hand.
- **Primary mode is `receptor_score_only`** — heavy-atom receptor, `--score_only` (CNN scores the exact
  pose, no minimisation), per octinoxate's gnina README; `receptor_h_score_only` and `receptor_minimize`
  are second passes. Columns include `CNNscore`, `CNNaffinity`, `minimizedAffinity`.
- **The GNINA-pick pose** (queued as a docked MD leg AFTER the Vina-p1 legs) is the docked pose (pose >= 1)
  with the **max `CNNaffinity`, consistent across all three passes**, corroborated by `minimizedAffinity`
  (argmin). It becomes `<struct>_dock<K>` via `run_dock_pose_md.sh STRUCT=<struct> POSE=<K>`; if the pick is
  pose 1 it coincides with the Vina-p1 leg, so no new leg. This feeds the `gnina_p*` REFERENCE fields.

## Reporting results

**Lead with MM/GBSA free energy and ligand retention, grouped by pipeline.** That is the comparison that
matters here. Pose-selection questions, scorer correlations and protocol diagnostics are secondary and
belong in their own section, labelled as methodology. On 2026-10-01 an invented "which starting pose won
this row" metric was reported as a finding and buried the pipeline comparison entirely, which inverted
the conclusion: peptidebuilder's four designs average −18.84 kcal/mol at 84.1% retention with one
release, against BoltzGen's two at −15.76, 70.5% and eighteen.

**MM/GBSA standard errors are not uncertainty on a comparison.** They measure frames within one
trajectory (0.04–0.10 kcal/mol). Run-to-run variation is far larger. Never write "this difference is N
standard errors" as though it bore on reproducibility.

## Paths

Scripts derive the repo root from their own location (`REPO=${0:A:h:h}`) and read `OMD_ENV`,
`BOLTZGEN_ROOT` and `DOCK_ASSIST`, defaulting to `~/miniforge3/envs/openmm-md`,
`~/python_mac/boltzgen_local` and `~/python_mac/dock_assist`. The MD environment **must be conda**:
`omd build` needs OpenFF/GAFF2 and `omd mmgbsa` needs AmberTools, and `.venv` has no `mdtraj`.
