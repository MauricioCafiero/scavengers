# HANDOFF — 2026-10-15 21:35

Supersedes the 2026-10-15 07:00 handoff, whose entire to-do list is now done. Read this, then
`CLAUDE.md`, then `FOLDS.md`.

---

## 1. Nothing is running on this Mac

The user disconnected for the night and asked for everything stopped. Verified clear: no
`caffeinate`, no `racc_feed.sh`, no `racc_drain.sh`, no `feed_chain.sh`, no tails. **Do not assume a
watcher is armed — there isn't one.**

**On racc, three jobs survive the disconnect and keep computing** (a network loss stops nothing
there; their output waits on scratch indefinitely):

| job | state at 21:30 | legs |
|---|---|---|
| `fold_dk5_oct_2_run` 521673 | RUNNING | 5 |
| `fold_dk5_oxb_2_run` 521739 | PENDING | 5 |
| `fold_dk6_oct_1_run` 521742 | PENDING | 3 |

**Ten batches were never submitted** because the feed chain was stopped with them queued:
`dk6_oct_2`, `dk6_oct_3`, `dk5_oct_3..6`, `dk5_oxb_3..8`. Their `.legs` files are committed and
their systems are built; they need only a chain launch.

### Resuming in the morning

1. The user rebuilds the two ControlMasters (their password and TFA; nothing can be done remotely
   without them). Verify with
   `ssh -o ControlPath=~/.ssh/cm-racc -o ConnectTimeout=15 racc.rdg.ac.uk true`.
2. Arm a session `caffeinate` before launching anything (`caffeinate -i -s -t 43200 &`).
3. Relaunch the feed chain over the **unsubmitted** batches, alternating molecules, detached:
   ```sh
   ~/miniforge3/bin/python -c "import subprocess; subprocess.Popen(['zsh','runs/folds_md/feed_chain.sh',
     'dk6_oct_2:octinoxate:6:00:00','dk5_oxb_3:oxybenzone:6:00:00','dk6_oct_3:octinoxate:2:30:00',
     'dk5_oxb_4:oxybenzone:6:00:00','dk5_oct_3:octinoxate:6:00:00','dk5_oxb_5:oxybenzone:6:00:00',
     'dk5_oct_4:octinoxate:6:00:00','dk5_oxb_6:oxybenzone:6:00:00','dk5_oct_5:octinoxate:6:00:00',
     'dk5_oxb_7:oxybenzone:6:00:00','dk5_oct_6:octinoxate:6:00:00','dk5_oxb_8:oxybenzone:6:00:00'],
     stdout=open('runs/folds_md/feed_chain_dk.log','ab'), stderr=subprocess.STDOUT,
     stdin=subprocess.DEVNULL, start_new_session=True, cwd='$PWD')"
   ```
   It launches a detached `racc_drain.sh` per batch itself.
4. **The three live jobs need their drains relaunched by hand** — their old ones were killed:
   ```sh
   for b in dk5_oct_2 dk5_oxb_2 dk6_oct_1; do mol=octinoxate; [[ $b == *oxb* ]] && mol=oxybenzone
     ~/miniforge3/bin/python -c "import subprocess,os; subprocess.Popen(['zsh','code/racc_drain.sh','runs/folds_md/$b.legs'],
       stdout=open('runs/folds_md/drain_$b.log','ab'), stderr=subprocess.STDOUT, stdin=subprocess.DEVNULL,
       start_new_session=True, env={**os.environ,'SYSNAME':'$mol'})"; done
   ```

---

## 2. The docked campaign — 10 of 76 runs done

76 docked legs, one per OF3/RF3 co-fold structure's **Vina pose 1**, plus its
**`receptor_score_only` max-`CNNaffinity` GNINA pick** where that is not pose 1. List with
corroboration flags: `runs/folds_md/dock_legs.txt`. All 76 systems are **already built**
(`runs/folds_md/dock_build.log`), so a leg only needs submitting.

Rows land in FOLDS.md under "Fold-lane docked legs"; 10 are there. Each row carries ΔG with a
convergence marker, residence, released frames, slide, and the **Δ against that structure's
cofold-pose leg**, which is the comparison the campaign exists for.

**What the first ten say.** Seven are weaker than their cofold pose (+0.87 to +1.74), three
stronger. The three gains are informative: `bg33_4_rf3cof_dock1/dock4` (−2.46, −2.92) and
`bgox31_3_o3cof_dock1` (**−5.17**, the largest), and in each case the cofold placement was poor —
`bgox31_3_o3cof` is the weakest OF3 cell in the oxybenzone set at −5.82. **Docking rescues a bad
co-fold pose and slightly degrades a good one** is the working reading; it is ten runs old, has not
been tested against a null, and every mechanistic reading in this campaign so far has been
falsified by a later leg. Hold it loosely.

`bg33_4_rf3cof_dock4` is one of only **4 of 35 picks meeting the full criterion** (same pose in all
three GNINA passes *and* `minimizedAffinity` argmin) and it is the best docked cell so far at
−19.06. RESULTS.md treats a pick on `CNNaffinity` alone as unsupported; 31 of the 35 are that.

### Recording a finished run

```sh
~/miniforge3/envs/openmm-md/bin/python code/record_leg.py <mol> <leg>
```
New this session. It drives the existing table scripts for one leg — `collect_mmgbsa_rows.py`,
`ligand_slide.py`, `md_stability.py`, `pair_contacts.py` — **refusing any whole-table rewrite that
would drop an existing row**, appends the FOLDS.md row, and prints the audit line. It is idempotent,
so re-running it is safe. It takes 2–4 minutes (the whole-table rewrites re-read every leg's
trajectory) — run it backgrounded.

**Find what still needs recording** (do not trust `tail -1` on the drain logs — they are several
files and not in time order):
```sh
grep -h 'tail exit' runs/folds_md/drain_dk*.log | sed 's/.*for \([^ ]*\) .*/\1/' | sort -u |
  while read l; do grep -q "| \`$l\` |" FOLDS.md || echo $l; done
```

---

## 3. Conventions set today — follow these

* **Commit once per 5 MD runs**, not per run. The user asked for this explicitly after the
  per-run rate (the convention for the previous several days) became too noisy. Each run is still
  recorded and audited the moment it lands; only the commit waits.
* **Report each run in 2–3 lines**: ΔG, retention, Δ vs the cofold pose, and anything genuinely odd.
  The user confirmed this verbosity.
* **Local-only ignores go in `.git/info/exclude`, never `.gitignore`** — a committed `.gitignore`
  advertises what exists locally and is being hidden. The `*_movie.mp4` rule is there now.
* **Every racc job ≤6 h.** `gpuscavenger` preempts, and a requeue restarts the in-flight leg from
  zero while the job holds the single run slot. Proof, 2026-10-15: the 14.5 h `dk_oct_1` produced
  **1 leg of 10 in ~7 h** across repeated preemptions while `dk5_oct_1` (4:16) and `dk5_oxb_1`
  (3:24) both completed in the same window. It was cancelled and its 9 survivors re-split into
  `dk6_oct_1..3`. The user had said this before it was proven; it was not acted on fully, and that
  cost the day's throughput.
* **One MM/GBSA tail at a time.** `code/run_dock_pose_md.sh` now re-execs under
  `lockf -k /tmp/peptidebuilder_tail.lock` (builds skip it), so any number of drains may fetch in
  parallel and their tails queue. This replaced hand-suspending drains with SIGSTOP. Concurrent
  tails OOM this 8 GB Mac.

---

## 4. The co-fold campaign is complete — all 29 legs

Both cofold lanes now have all 21 structures per molecule through 20 ns and MM/GBSA, with every
table row written. **RF3 beats OF3 in 15 of 21 cells** (octinoxate 7/10, oxybenzone 8/11), but four
of those wins are against OF3 cells that slid >7 Å in equilibration and so never tested their
predicted pose — check `runs/<mol>/md/ligand_slide.csv` before reading any cell.

The pre-registered `bg33` conflict **resolved against the molecule split**: `bg33_1` lost by 2.47
and `bg33_2` tied, both against clean comparators, and `bg33_4` lost by 0.98. The one regular
pattern is BoltzGen-specific: every `bgox31` RF3 cell beats OF3 *and* falls short of its BoltzGen
baseline.

**Presentation tables**: `runs/folds_md/cofold_binding_tables.{md,pdf}`, two tables (one per
molecule), ΔG with retention and a convergence marker per cell. The user asked these be left as
they are — the Boltz column is a constrained co-fold for the `_f4/_f8/_f12` peptides while OF3/RF3
are unconstrained, and **there is no unconstrained Boltz MD to swap in** because dynamics were never
run on the poor unconstrained folds. Do not "fix" this; it was raised and closed.

---

## 5. Uncommitted / untracked, deliberately

* `code/render_md_movie.py` — new, headless-PyMOL trajectory → mp4. The user asked it stay
  uncommitted. Renders a leg's `traj_wrapped.{xtc,pdb}`: `intra_fit` on peptide CAs to kill the
  tumble, one camera for the whole movie, buffer −8 (tested: fills ~70% of frame, no clipping).
  `ox2_orig_f8_movie.mp4` (20 s, 9.8 MB) lives beside that leg and is excluded locally.
* `CH3PR1_pres_2026.pptx`, `Claude outputs/` — the user's, not mine. Leave alone.
* A shell render (ligand + 12 placed side-chain fragments, pre-backbone) was made from
  `runs/octinoxate/figures/designed_shell.xyz` — ligand is atoms 1–44, fragments 45–190. It was
  scratch-only; nothing committed. No shell image from 2026-10-14 was ever found.

The 2.18 GB `peptidebuilder-backup-20261006-063634.bundle` in `~/python_mac/` was **deleted** at the
user's request after verification: the only thing it held that current `main` lacks was 120
`first_Nns/traj_wrapped.xtc` window slices, which the repo deliberately does not keep. It had been
made as a safety net before those slices were stripped from history — a 2.18 GB backup of the files
the strip existed to remove.
