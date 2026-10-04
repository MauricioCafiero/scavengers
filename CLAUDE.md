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

## MM/GBSA housekeeping

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
