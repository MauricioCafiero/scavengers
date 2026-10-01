# Working rules for this repo

Written 2026-10-02 after a day in which every rule below was learned by breaking it.

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
**Keep what the September legs keep:** no `reference.frc` and no `_MMPBSA_*` anywhere;
`FINAL_RESULTS_MMPBSA.dat`, the prmtops, `energy.csv`, `traj_wrapped.*`, `frame_medoid.pdb` and
`frames_last.pdb` stay. Window directories are named `first_Nns`.

Every result goes in `runs/octinoxate/md/mmgbsa_summary.csv` **as it is computed**, not in a batch at the
end. Nothing reads that file — `make_gnina_bundle.py` carries the same numbers in a hardcoded `REFERENCE`
dict — so a new row has to be added in both places or the next GNINA bundle is built on stale values.

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
