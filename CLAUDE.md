# Working rules for this repo

Written 2026-10-02 after a day in which every rule below was learned by breaking it.

**Use the scripts that already exist; parameterise, don't reimplement.** The established flow for a
docked Modal leg is: `run_dock_pose_md.sh BUILD_ONLY=1` → `modal push` → `probe` → `produce` →
`md_window_modal.sh` → whole-run `omd mmgbsa` → `md_contacts.py` → `md_frames.py` → `pair_contacts.py`,
each invoked directly. The rule is about *reimplementing* a tool's steps, not about automation: a
driver that invokes these scripts is fine, and `code/racc_feed.sh` (queue top-up) and
`code/racc_drain.sh` (fetch + tail) are exactly that — they call `run_dock_pose_md.sh` and
`run_md20.sh` rather than inlining their stages, so each step's tacit knowledge stays where it lives.

An earlier version of this file said such a driver had been "drafted twice and rejected twice" on
2026-10-04. **That was wrong** — the user never rejected it, and says a long-context session was
misreading them at the time (corrected 2026-10-07). Do not treat orchestration as forbidden here.

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

**`run_md20.sh` is not in this repo.** It lives in the sibling BoltzGen repo at
`~/python_mac/boltzgen_local/md/run_md20.sh`, which is where `racc_drain.sh` calls it from
(`$BG/md/run_md20.sh`). Only `run_dock_pose_md.sh` is in `code/`. A change to the co-folded tail
therefore lands in a different repo and a different git history — noted 2026-10-14 after an earlier
version of this file implied both were local.

**Figures are a pipeline step too, and they need their own environment.** `code/render_md_figures.py`
renders a leg's `frame_medoid.pdb` to a PNG headlessly and runs at the end of the tail, so every leg
gets its figure as it lands instead of by hand months later. It needs **PyMOL, which is not in the MD
env** — it has its own conda env because the `pymol-open-source` wheel on PyPI links to libraries
that do not exist on this machine and cannot be made to load:
```sh
~/miniforge3/bin/mamba create -y -n pymol-render -c conda-forge python=3.12 pymol-open-source
~/miniforge3/envs/pymol-render/bin/pip install pillow   # alpha-channel cropping
```
Two things about these structures that are not obvious. mdtraj writes the trajectory PDBs with **no
HELIX/SHEET records**, so PyMOL loads every residue as coil and an automatic cartoon draws a helix
flat — the script runs `dss` first to assign secondary structure from geometry, without which every
helical leg in the campaign renders as a meaningless squiggle. And the committed figures are
**transparent RGBA**, not white-backed, which is why `figures_shell*/style.pml` says `bg_color black`
while its own output reads white: there is no background in them at all.

**`md_contacts.py` and `md_frames.py` are pipeline steps, not extras.** Contacts give residence and
separation, which a scalar ΔG cannot and which separated every design from every null. Frames give the
end-of-run ensemble, needed because the predicted pose and the pose after 20 ns are different
structures. A leg without them is incomplete; compare against `s3_orig_f12_dock1`, which has the full
file set.

**Record input → production ligand RMSD on every leg** with `code/ligand_slide.py <prod dirs...>`
(`--csv` to write the table; it superposes on protein Cα first, so the number means the ligand moved
relative to its peptide rather than the box translating). Equilibration restrains protein heavy atoms
only (`openmm_md/dynamics.py:144`), so the ligand slides 1.2 to 15.1 Å before production starts
depending on the run. **Do not "fix" this by restraining the ligand** — every run in the project
equilibrated it freely and a restrained leg would not be comparable. Record the number; it is what
says whether a cell is worth reading. Oxybenzone's twelve BoltzGen legs came in at 1.05–3.15 Å, the
tightest spread in the project.

**A second molecule must pass `SMILES` to `boltzgen_local/tools/metrics.sh`.** It defaults to the
octinoxate analogue, so `lig_geom.py` would otherwise measure every predicted bond against the wrong
reference and report plausible-looking errors. Perceive it from the run's own geometry rather than
typing it: `boltz_check.ligand_smiles('runs/<mol>')` reads `ligand.xyz`.

## Modal

**Launch every `modal run` through `code/modal_run.sh` — it is not optional.** The client stays
attached for the whole leg even with `--detach`, and when the Mac sleeps the client drops and Modal
cancels the in-flight leg at full price. The script `exec`s the run under `caffeinate -i -s`, which
holds idle and system-sleep assertions for exactly the client's lifetime, so the guard is self-cleaning
and covers the whole leg (2026-10-04: two cancelled legs, $1.06, for skipping this; 2026-10-10: a
16-minute maintenance/dark sleep dropped a `-i`-only produce client mid-leg — the client's disconnect
handler cancelled the input itself — which is why `-s` was added; valid on AC, and all legs run on AC).
There is no warm restart: `omd run`
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

**Do not use Modal's N^0.77 scaling here — it is an A10G number and badly overestimates the H100.**
Measured on racc: 0.236 ms/step at 5,838 particles (`bgox31_3`, 2026-10-07, from two `energy.csv`
samples 180 s apart), against wall times of 45:56–47:57 at 5,425–5,761 particles and 56:28–65:49 at
12,826–12,841. Going 2.3× in particle count adds only 10–18 minutes of wall, so the effective
exponent is nearer **0.4**: these systems are far too small to saturate an H100 NVL, and 7–9 minutes
of every leg's wall is fixed JIT-plus-equilibration overhead that does not scale with N at all.
Applying the A10G exponent to an 18,650-particle leg predicted 1.6–1.8 h where ~70–80 min is right.
For a quick figure: **wall ≈ 8 min + 10M × 0.236 ms × (N/5838)^0.4**.

**The 0.4 exponent holds only to ~13k particles. Above that use 0.64** — measured 2026-10-18 on
`orig_control` at 19,617 particles: 0.512 ms/step net of overhead (9,945,000 steps in 1:32:51 from
`sacct` start to `energy.csv` mtime, less ~8 min JIT/equilibration), which is 2.17× the 5,838-particle
rate for 3.36× the particles → exponent 0.64 across that range. The 0.4 model predicted 71 min for
that leg against ~94 min actual, a 1.3× under-estimate, because it was calibrated only on
5.4k–12.8k-particle legs. So: **≤13k particles, `(N/5838)^0.4`; above, `(N/5838)^0.64`.** This
matters for sizing a `--time` request — unforced folds are rods and run 18k–23k particles, where
the old model silently buys a limit that cannot fit the batch (`bzfree_oct_2`'s three 23k legs came
out at 5.14 h against the 5:20 limit it had been given, 11 minutes of slack).

**A pending job's TimeLimit can be RAISED, a running one's cannot** (measured 2026-10-18:
`scontrol update JobId=<pending> TimeLimit=06:40:00` succeeded, the same call on the running job
returned "Access/permission denied"). So a mis-sized limit is recoverable while the job waits and
unrecoverable once it starts — which is the real reason to size it correctly at submission.

**energy.csv flushes continuously here, unlike on Modal** — mid-run progress and ms/step are
readable any time (row cadence 500 steps) without contacting produce. Use it to revise estimates
within the first 20 minutes instead of waiting for completion to find out.

**Two watchers do the queue work; neither reimplements a tool.** `code/racc_feed.sh <legfile>` keeps
the queue topped up from an explicit ordered list of `<local dir>:<leg name>` lines, staging each leg
by tar before submitting; `BATCH=1` instead submits the whole list as ONE job. `code/racc_drain.sh
<legfile>` watches for completions, fetches each one's wrapped solute and runs that leg's tail,
dispatching to `run_dock_pose_md.sh` for a `<leg>_dock<N>` and `run_md20.sh` for a co-fold, then
strips the MM/GBSA scatter. Launch both detached (`start_new_session`), and give every ssh a
`ConnectTimeout` — a wedged ControlMaster with the network gone under it blocked one poll for 40
minutes on 2026-10-07, where a timeout costs seconds.

**The drain loop as reworked (2026-10-16).** `racc_drain.sh <legfile>` eager-fetches each poll's
completions rather than checking then fetching, so a poll lag never costs analysis time. Each
completion's tail (`run_dock_pose_md.sh` continuation or `md_window_modal.sh`-style record) launches
**detached under the global tail lock** `lockf -k /tmp/peptidebuilder_tail.lock` — one heavy tail at
a time, ~10–18 min per tail — with the drain never blocked on it; markers `.drained`, `.tail_launched`
and `.tail_pid` in `runs/folds_md/` make the bookkeeping safe against reconnects and restarts (a
re-launched drain resumes rather than refetches). Recording goes through
`code/record_leg.py <mol> <leg>`, which appends to all the molecule's tables and to FOLDS.md; the
per-leg log `runs/folds_md/record_<leg>.log` ends with
`<leg>: dG ... | windows ... | residence ... | slide ... | Ca ... | pairs ...` and an
`audit: windows 4/4, summary rows 4, frames True, scatter 0 -> OK` line, which is the
leg-is-done signal. A campaign is a chain of legfiles `runs/folds_md/dk*.legs` (lines
`<abs local dir>:<leg name>`, `#` comments), one drain per legfile; a committed batch is every
recorded leg since the last commit, verified additions-only.

**The waiter is `code/watch_legs.py`, not a fresh scratchpad script.** It takes `mol/leg` specs,
`--legfile <racc legfile>` (repeatable) or `--any`, waits for the next leg to finish, prints its
dG + windows + signed ddG-against-run-1 for a `_r2` leg + contacts, and **exits** -- because the
harness re-invokes on a background job's exit, not on its output, so a waiter that reports and
keeps looping is mute however much it prints. Re-arm it after every report. `--on record` (the
default) keys on the `audit:` line, which `record_leg.py` writes last; `--on mmgbsa` keys on the
whole-run `FINAL_RESULTS_MMPBSA.dat` and fires ~16 min earlier, while the window slices are still
computing. **Run `--selftest` before trusting it**: it resolves every trigger path and says which
already exist, which is the only cheap check against the two bugs that hide longest -- a path that
can never exist, and a signal that fires early. Both of those shipped in throwaway versions
(2026-10-19, 2026-10-20); the second swallowed three legs' reports because its check was
`runs/<mol>/<leg>/...`, missing the `md/`, while the python inside it built the path correctly.

**Report every leg as it lands; commit every five.** These are two different cadences and
collapsing them into one is a mistake in either direction. The report is the scientific output and
is owed immediately — it needs nothing committed, because `FINAL_RESULTS_MMPBSA.dat` and the
leg's tables are already on disk when the tail finishes. The commit is bookkeeping, and one per leg
turns a 26-leg campaign into 26 pushes of near-identical table deltas. So: read the number out and
say it the moment the tail is done, let the recorded legs accumulate, and commit in batches of
about five (or at a natural boundary -- a legfile finishing, a campaign ending, or stopping for the
day). Re-stated by the user 2026-10-20 after three replicate legs went out as three separate
pushes. Host/access facts and the measured
H100 rate model live in `RACC.md`, which is untracked on purpose.

**`gpuscavenger` is `MaxSubmitJobsPU=3`, `MaxJobsPU=1`.** Only three jobs may exist and only one
runs, so **spare slots cannot add throughput** — a second job holding different legs just waits its
turn, exactly as those legs would have inside the first. The only thing idle slots can usefully hold
is a standby copy for failover, and with `Requeue=1` already set on these jobs even that is close to
redundant. The real consequence is that one-leg-per-slot needs the socket re-touched every ~50
minutes to submit the next, and that socket needs a password and TFA only the user can give. For any
unattended stretch, bank the remainder in one job: `racc_run.py` takes `LEGS` as a comma-separated
list and runs them in sequence, guarding each on its own `traj_wrapped.xtc` on scratch so a
resubmission after preemption skips what is already done. The partition's **DEFAULTTIME is
1-00:00:00 against a 5-day cap**, so an unspecified request is 24 h.

**Every submission passes `SBATCH_OPTS="--time=<estimate + 1 h 30 min> --cpus-per-task=8"`, and a
batch is sized to ≤6 h of estimate.** Those are the campaign's settings, not suggestions:
`racc_feed.sh` *refuses* to run (exit 2) unless `SBATCH_OPTS` carries both, because an optional
step is a skipped one. The estimate is the sum of `8 min + 39.3 min × (N/5838)^e` over the
batch's legs (e = 0.4 up to ~13k particles, 0.64 above — see the rate note above), N read off
each build; the ≤6 h split is why `dk5_*`/`dk6_*` are many small legfiles
rather than one — a preempted leg restarts from zero, so a smaller batch bounds the loss. A 24 h
default request scores worse on JOBSIZE priority, cannot backfill into a short gap, and sits
exposed to preemption for hours it will never use; 8 cores lift JOBSIZE priority and must come
through `SBATCH_OPTS`, since `gpu.slurm`'s own `--cpus-per-task=1` beats `SBATCH_TRES_PER_TASK`.
Each legfile header carries its batch's values (`runs/folds_md/bzfree_*.legs`).

**This was learned by getting it wrong on 2026-10-18**, when five batches went out with
`SBATCH_OPTS` empty — the convention was in `racc_feed.sh`'s own parameter comment, which had been
read and not applied, and is in neither RACC.md nor this file until now. The after-the-fact repair
is `scontrol update JobId=<id> TimeLimit=<hh:mm:ss> CPUsPerTask=8 MinCPUsNode=8`, and it took only
on the PENDING jobs: a running job keeps the cores it started with, so the one already executing
stayed at 2. Getting it right at submission is the only reliable route.

**Do not repair the core count with `NumCPUs` alone — it wedges the job.** `scontrol update
NumCPUs=8` is accepted silently but raises only the job's total CPU count, leaving `CPUs/Task=1` and
`MinCPUsNode=1`. The request is then self-inconsistent (8 CPUs on one node, one task of one CPU,
one-CPU-per-node minimum) and Slurm parks the job at **`Reason=BadConstraints`** forever — it never
starts and the reason never clears. On 2026-10-18 that silently killed two of the three queue slots;
because `gpuscavenger` is `MaxSubmitJobsPU=3`, the two dead jobs also blocked the third batch's
feeder from ever submitting, stalling 12 of 23 legs behind a job that was never going to run. The
fix is `scontrol update JobId=<id> CPUsPerTask=8 MinCPUsNode=8`, after which the reason becomes
`QOSMaxJobsPerUserLimit` — which is the *healthy* pending reason here, meaning only "another of your
jobs is running". **Read a pending job's `Reason` as data, not decoration:** `BadConstraints` on a
job whose peers run means the request is broken, not that it is queued.

**A leg that dies on the cluster is invisible from the queue — grep the job log for it.**
`racc_run.py` catches a failed leg, prints `!!! <leg> FAILED: <exception>` and moves to the next one,
which is right for a batch but means the job stays `RUNNING`, `sacct` shows nothing wrong, and the
drain waits forever for a `traj_wrapped.xtc` that will never appear. Nothing writes a marker, so a
watcher keyed only on completions is silent on exactly the case it most needs to report. Check with
`grep -h '^!!!' ~/remote_work/openmm/*_run.out`, and the real stack trace is in
`~/remote_work/openmm/gpu_out.txt` — the `CalledProcessError` in the job log carries only the command
line. Found 2026-10-19 when `ox1_esm1_control` NaN'd in restrained equilibration and the batch
carried on through its other five legs.

**Try a plain resubmission before concluding a leg's build is broken.** That failure looked like a
steric clash — 1.69 Å between a ligand and a peptide heavy atom, against 2.19 Å on its docked
siblings — and the obvious fix was a harder minimisation. Both readings were wrong: `minimize_max_iter`
is already `0` (uncapped, runs to the 10 kJ/mol/nm tolerance), so there was no cap to raise, and the
*same system through the identical code path* (`omd run --platform CPU --steps 2000` from the leg
directory) minimises that contact to 2.6 Å and equilibrates clean. The resubmitted leg then produced
**the strongest result in the campaign** (−26.09, zero release frames, 46 mean contacts): the tight
contact was a deeply enclosed pose, not a broken one. So reproduce the failure locally with the real
tool before changing anything, and never "fix" it by nudging the ligand — that would have destroyed
the result. A one-leg retry legfile (`bzfree_oxb_retry.legs`) resubmits with the protocol untouched.

**A fetch can wedge with its data already complete; check the payload before refetching.** The drain
streams `ssh ... | tar xzf -`, and if the ssh end dies the `tar` blocks on an EOF that never comes —
it sat 40 minutes on 2026-10-19 holding `drain_bzfree_oxb_1` behind it, with every file already
written. Verify before assuming loss: `mdtraj` frame count on `traj_wrapped.xtc` (expect 20,000 at
20 ns) and `wc -l` on `energy.csv` (expect 20,001). Both were complete, so killing the `tar` was
enough — the drain re-polled, found the leg, tailed it locally and recorded it with no network at
all. A mid-tail leg transiently holds ~1.5 GB of MM/GBSA scatter (a 1.2 GB `reference.frc` plus the
`_MMPBSA_*` mdcrds) which the drain's `find -delete` strips on completion, so a leg directory in the
gigabytes during its tail is normal and is not a bad fetch.

**But a wedged fetch is sometimes genuinely TRUNCATED, and then killing the `tar` is not enough —
the partial file must be deleted.** Four fetches wedged on 2026-10-20; three held complete payloads
and one (`orig_f12_dock1_r2`) held 16,317,952 bytes of a 31,758,108-byte `traj_wrapped.xtc`, with
`energy.csv` and `traj_wrapped.pdb` never written at all. **Compare the local size against the
cluster's** (`ls -la` on both) rather than only counting frames — mdtraj cannot even open a
half-written xtc without its `.pdb` topology, so the frame check reports a confusing
`No such file: traj_wrapped.pdb` instead of a clean "truncated". The danger is specific and silent:
`racc_drain.sh`'s `remote_phase` treats *the mere presence* of `prod_20ns/traj_wrapped.xtc` as "leg
finished, never re-fetch, tail it", so a truncated file left in place yields a tail and an MM/GBSA
number computed over half a trajectory, which looks entirely plausible. Delete the partial and let
the drain re-poll; it costs one ~30 MB re-download and no compute. For the same reason, when staging
a fetch by hand, move `traj_wrapped.xtc` into place **last**, after the small files.

**Submission is socket-free once batched; draining is not.** A network loss stops nothing on the
cluster — the legs keep computing and their wrapped solutes wait on scratch indefinitely — but the
drain cannot fetch until the user rebuilds the two ControlMasters. The worst case is delayed
analysis, never lost compute; say so rather than implying the run is at risk.

Waits go in background watchers that poll sacct (`--format=State --noheader | tail -1`), never
foreground sleep loops. **A completion is better judged by the leg's own output than by a job
state**: a leg running inside a multi-leg batch has no `<leg>_run` job for `sacct --name=` to find,
so key on `traj_wrapped.xtc` and consult sacct only to tell "not started" from "its job died". Partition is `gpuscavenger`: free, preemptible, and — like on Modal — a
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

**The cofold tables are a three-step chain; the splice is now scripted and the PDF still is not.**
`code/cofold_tables.py --write` splices Tables 1 and 2 into `runs/folds_md/cofold_binding_tables.md`
itself and **refuses if any published cell would change**, allowing only cells in a newly added
column and the GNINA-pick row label, which carries one pose tag per column. That check is the
"cross-check the co-fold cells against the published values" the docstring used to ask for in prose,
and the splice was a manual paste until 2026-10-19 — which is why the pocket-guided lane's first six
legs sat in the CSVs with no row in the markdown at all. Without `--write` the script only prints a
`missing / pending` report and writes table *blocks* to `$OUT` (default `/tmp/tables.md`). The
markdown is still the source of truth, and `code/render_cofold_pdf.py` rebuilds
`cofold_binding_tables.pdf` from it — **with the venv or conda python, not the bare `python3`, which
has no matplotlib**. Skip that third step and the PDF rots invisibly: the committed one went four days stale across a restructure and a whole
campaign, and had no committed generator at all until 2026-10-19 (its metadata said matplotlib, and
that was the only clue). Also note `cofold_tables.py` resolves each leg's residence from
`runs/folds_md/record_<leg>.log` first and only then from the prose tables — scanning the docs alone
once harvested a *static fold-geometry* percentage (`RESULTS.md`'s leg-keyed geometry rows) and
printed it as MD residence.

**`record_leg.py` can extend a lane's tables but could not start one, and said so misleadingly.**
It finds the `pair_contacts_<lane>_<arm>.csv` already holding a leg's parent structure and re-runs
over it; with no such file it printed `no pair table (no design json for this arm)`, which is false —
`pair_contacts.design_for()` resolves the json by peptide name and every file exists. The real
situation is a lane with no table yet. That message reads as "nothing to do", so all five of the
pocket-guided lane's first legs recorded with no pair row while the audit line still said `OK`.
It now starts the lane's table — but **that file is keyed on the design arm, not the peptide**, so a
"new" lane table is often one another peptide of the same shell already owns (`s3_esm2_f4` and
`s3_orig_f12` are both `design_shell3`). Starting it from an empty row list rewrote it with only the
new leg, and the `before <= set(after)` assertion passed vacuously because `before` was empty: five
rows went that way across three files before a `git diff --numstat` caught it on 2026-10-20. It now
loads the file's rows when it exists. **Run `git diff --numstat` on every leg's commit** — an
additions-only diff is the only thing that actually proves no whole-table rewrite shrank a table.
The lesson generalises: **`record_leg.py`'s own `audit: ... -> OK`
is not the done-signal** — it checks windows, summary rows, frames and scatter, and nothing else.
Run the four-line audit in "Audit the recording set" against the leg, which is what caught this.

**A fold lane also needs a `runs/folds_md/dock_legs.txt` entry per docked leg before it is recorded.**
That file is where FOLDS.md's pose label comes from, and a missing entry used to default to `GNINA` —
so the guided lane's two Vina-p1 legs were published as `p1 GNINA`. A missing entry now writes
`UNREGISTERED` and warns. It deliberately does not abort: the MM/GBSA rows are appended earlier in
the same `record_leg.py` run, so exiting there would mean a re-run duplicates them.

**Whole-table scripts rewrite their table: name every structure, not just the new one.**
`vina_redock.py` writes `dock_summary.csv` and `dock_poses.csv`; `make_gnina_bundle.py` writes
`reference.csv` — one script each, not both, whatever an earlier version of this file said. All three
are rewritten wholesale over only the structures named, so a run naming only the new ones drops every
other row — the same trap `md_stability.py` has. On 2026-10-06 docking two nulls cut
`dock_summary.csv` from five structures to two.

**`ligand_slide.py --csv` is in this family as well**, and worse, it had no table at all for the
octinoxate lane until 2026-10-14 — ten legs' slides had been read off stdout, quoted in `FOLDS.md`
and never persisted. It writes with mode `w` over only the legs named, and its leg list has the same
lane problem as `md_stability.py`: the BoltzGen baselines live in
`~/python_mac/boltzgen_local/md/<leg>/prod_20ns`, so a run globbing only `runs/<mol>/md/*` silently
dropped `bgox31_2`–`bgox31_5`. Pass every lane, same globs as `md_stability.py`, and **write the
table on every leg** — the number is what says whether a cell is worth reading, and back-filling it
retired a counterexample this file had relied on for weeks (`shuffle_control_esm0_o3cof`, 10.62 Å:
production never started from the OF3 pose). Tables now at `runs/{octinoxate,oxybenzone}/md/ligand_slide.csv`.

**`pair_contacts.py --csv` is in this family too**, and was missing from this list until 2026-10-14,
when writing `bg33_2_o3cof` alone silently dropped `bg33_1_o3cof`'s row from
`pair_contacts_o3cof_boltzgen.csv`. It has no `--allow-shrink` guard, so nothing warns you. Name
every leg on every run:
`$OMD/bin/python code/pair_contacts.py <leg dirs...> --run-dir runs/<mol> --design design.json --csv runs/<mol>/md/pair_contacts_<lane>_<arm>.csv`
— `--design` resolves against `--run-dir`, so it takes `design.json`, not a full path. **`--run-dir`
is required and has no default** (changed 2026-10-14, with `collect_mmgbsa_rows.py`'s `--outdir`): it
chooses which molecule's design json defines the slots, and a default would score one molecule's legs
against another's slots and print plausible-looking counts. The script also refuses a run whose legs
are not under `runs/<that molecule>/`.

**The route that avoids the damage rather than repairing it:** `vina_redock.py` has an `--out`, so
dock the new structures into a scratch directory and append the rows to the real tables, which leaves
committed legs' pose files untouched; verify with `git diff --numstat` that the change is
additions-only. `make_gnina_bundle.py` has no `--out` and its `out_root` is hardcoded, so back
`reference.csv` up first, build only the new bundles, then restore and append — it does not rebuild
other structures' bundle directories that way, which matters once their GNINA results are retrieved.
Both verified on 2026-10-07 adding four BoltzGen legs to oxybenzone's eleven-structure tables.

Re-docking everything is *not* the fix: a re-dock is reproducible (seed 42 gives byte-identical rows)
but it rewrites files belonging to legs whose GNINA and MD are already done and committed. `--all` is
not the fix either — it silently pulled in an eighth structure that had never been in the docking
matrix. Name the targets.

`md_stability.py` is the worst of the three because its one file is shared across **every** ligand
plus the BoltzGen legs (62 rows as of 2026-10-07), and the glob in older notes omitted
`runs/oxybenzone` entirely. It now refuses to write fewer rows than the file holds unless
`--allow-shrink` is passed. Pass every lane: `runs/octinoxate/md/*/prod_20ns`,
`runs/octinoxate/md/*/prod_L1_modal`, `runs/oxybenzone/md/*/prod_20ns`,
`runs/oxybenzone/md/*/prod_L1_modal` (oxybenzone has Modal legs since 2026-10-10 — omitting this
glob silently drops them from the rewritten table),
`~/python_mac/boltzgen_local/md/*/prod_20ns` — and not `boltzgen_local`'s `prod_L1_modal`, which does
not exist and whose empty glob aborts the whole zsh launch.

**`vina_redock.py` must be run with `~/python_mac/dock_assist/dock-env/bin/python`**, never
peptidebuilder's `.venv`: it locates Vina through `import dockstring`, which is installed only in
`dock-env`, and the `DockError`'s advice to "pass --vina-bin" is a dead end because `vina_redock.py`
does not expose that flag. `dock-env` carries rdkit and numpy, so the script runs there unchanged.

## MM/GBSA housekeeping

**`run_dock_pose_md.sh` does not strip its own scatter — you must.** It runs `omd mmgbsa` for the whole
run and for each window and leaves every `reference.frc` and `_MMPBSA_*` in place; only
`md_window_modal.sh` carries the `rm`. On 2026-10-06 four docked legs tailed through the driver were
holding **24 GB** (5.0–7.2 GB each, against 294–427 MB for two co-folded legs tailed through
`md_window_modal.sh`). Strip with
`find <leg> \( -name 'reference.frc' -o -name '_MMPBSA_*' -o -name 'traj.nc' \) -delete`
— **`traj.nc` belongs in that list and was missing from it until 2026-10-07**, by which point 89 of
them held 6.0 GB across the two campaigns (one per window plus one per whole run, ~70 MB each). It is
MMPBSA.py's NetCDF conversion of `traj_wrapped.xtc`, which is kept, so it regenerates in seconds and
is tracked by neither repo; `code/racc_drain.sh` now strips it automatically
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

Every result goes in **its own molecule's** `runs/<mol>/md/mmgbsa_summary.csv` **as it is computed**,
not in a batch at the end. One file per molecule: oxybenzone has had its own summary since its legs
landed. Pass **both** of `collect_mmgbsa_rows.py`'s paths —
`--md-root runs/<mol>/md --outdir runs/<mol>` — because they do different jobs: `--md-root` says
where legs are *read* and `--outdir` says which summary is *written*, so pointing only `--md-root` at
another molecule reads the right legs and appends to the wrong table. That put 16 oxybenzone rows into
octinoxate's file on 2026-10-14.

**`--outdir` is required and has no default** (changed 2026-10-14 at the user's instruction: a script
that points at one molecule or the other should error rather than guess), and it refuses an `--outdir`
whose summary does not exist, listing the molecules that have one. Still worth the check afterwards:
`grep -c '^<leg>,' runs/<mol>/md/mmgbsa_summary.csv`.

Nothing reads these files — `make_gnina_bundle.py` carries the same numbers in a hardcoded `REFERENCE`
dict — so a new row has to be added in both places or the next GNINA bundle is built on stale values.
When a peptide's docked legs land, extend its `REFERENCE` entry with the flat `vina_p1_*`/`gnina_p9_*`
fields, as both shuffle peers now have.

**Audit the recording set per leg; do not tick it off from memory.** Working from memory is how
`ligand_slide.py`'s `--csv` got dropped on ten consecutive legs on 2026-10-14 — each slide was read
off stdout, quoted in `FOLDS.md`, and never written to a table, which left a 10.62 Å cell being
cited as a counterexample for weeks. One query over the finished leg catches it:

```
leg dir: <n>/4 FINAL_RESULTS_MMPBSA.dat, md_contacts.csv, prod_20ns/{frames_last,frame_medoid}.pdb,
         and no reference.frc / _MMPBSA_* / traj.nc anywhere
tables : 4 rows in runs/<mol>/md/mmgbsa_summary.csv, 1 in protein_stability.csv,
         1 in runs/<mol>/md/ligand_slide.csv, 1 in a pair_contacts_<lane>_<arm>.csv
```

Match `pair_contacts_*.csv` on its **`structure`** column, not its first column — `leg` holds the
full production path, and matching that reports every row missing. **The two whole-table rewrites key
differently from each other**, which is its own false-alarm source: `protein_stability.csv`'s `leg` is
the full production path (`runs/<mol>/md/<leg>/prod_20ns`) while `ligand_slide.csv`'s is the **bare
leg name**. Grepping both for `/<leg>/` reports the slide row missing on a leg that has it.

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
Perplexity is deterministic; z is not. Its `RESULTS.md` and `outputs/inverse/*.json` are
**gitignored**, so a rerun that rewrites them cannot be recovered from git — copy them aside first.

**`~/python_mac/morph`** — apo → co-fold morph movies, MM strain (ff14SB + GBn2) and
charged/hydrophobic surface-patch changes. `$PY code/run_batch.py` with the conda `openmm-md`
python; `--set matched` (the default) runs the model-matched pairs
`boltzapo→boltz` / `of3apo→openfold3` / `rf3apo→rf3`, so a ΔE is the ligand's effect inside one
model. A new apo lane means a new pair in `MATCHED_PAIRS`. `--set legacy` reruns the superseded
all-apo-lanes→one-Boltz-co-fold analysis into `results/legacy/`; nothing should be concluded from
it, and ESMFold exists only there.

**ESMFold is retired from downstream analysis (2026-10-11)** on pepgraph's evidence: last of six
lanes at mean perplexity 20.4 — the no-information limit — owning 7 of the campaign's 11
zero-information structures, 6 of those with exactly zero residue contacts. Its folds and docks stay
committed as the lane that exposed the apo collapse; do not build anything new on them, and do not
bundle its docks for GNINA.

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

**"Residence" and "release" are different measures, and only release is shape-independent.** In
`md_contacts.py`, `residence_pct` is the fraction of frames whose **ligand-COM-to-peptide-COM
distance is ≤ 10 Å** (line 101), while `released_pct` is the fraction with **zero atom–atom contacts
under 4 Å**. A ligand held firmly at one end of an elongated peptide is never released and yet
scores low on residence, because the far half of the rod drags the peptide's centre of mass away
from it. So **a COM-based residence number penalises a fold for being elongated**, and the unforced
co-folds are rods: `s3_esm2_control` came in at 29.4% residence with 0.15% release — permanent
atomic contact — against `ox2_esm1_control_dock1`'s 21.2% residence with 31.35% release, which
really was letting go. Quote release (plus episode count, longest episode and the late share) when
retention is the claim; use residence only within one shell, where the geometry is comparable.
Learned 2026-10-19 after several legs had been contrasted on residence alone.

## Paths

Scripts derive the repo root from their own location (`REPO=${0:A:h:h}`) and read `OMD_ENV`,
`BOLTZGEN_ROOT` and `DOCK_ASSIST`, defaulting to `~/miniforge3/envs/openmm-md`,
`~/python_mac/boltzgen_local` and `~/python_mac/dock_assist`. The MD environment **must be conda**:
`omd build` needs OpenFF/GAFF2 and `omd mmgbsa` needs AmberTools, and `.venv` has no `mdtraj`.
