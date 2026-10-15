# HANDOFF — 2026-10-15 07:00

Written to hand the fold-cofold MD campaign to a fresh session. Read this first, then
`CLAUDE.md`, then `FOLDS.md`.

---

## 1. STOP — read before touching anything

**Processes are suspended on purpose. Do not kill them, do not restart them blindly.**

Four `racc_drain.sh` instances and one `feed_chain.sh` are detached at PPID 1. Two drains and two
tail process-trees are in `T` (SIGSTOP) state because **concurrent MM/GBSA tails OOM this 8 GB Mac**
— three were running at once at 06:37 with 282 MB free. The user's instruction is **serial, not
parallel: one tail at a time.**

| pid | what | state |
|---|---|---|
| 10439 | drain `oxb_2.legs` | running (its legs are done) |
| 30961 | drain `oct_3.legs` | running |
| 43870 | drain `oxb_3.legs` | **SIGSTOP** |
| 27492 | drain `oct_4.legs` | **SIGSTOP** |
| 27426 | tail tree (oct_3) | resumed 06:53, should be running |
| 27737 | tail tree (oxb_3) | **SIGSTOP** |
| 28106 | session `caffeinate -i -s -t 43200` | runs to ~18:35, not tied to a drain |

**To resume a suspended tree**, walk its descendants and SIGCONT each — `pkill -P` is not enough,
`mmpbsa_py_energy` is a grandchild:

```python
import subprocess, os, signal
out = subprocess.run(['ps','-ax','-o','pid=,ppid='], capture_output=True, text=True).stdout
kids = {}
for l in out.strip().split('\n'):
    p, pp = l.split(); kids.setdefault(int(pp), []).append(int(p))
def tree(r):
    s=[r]; st=[r]
    while st:
        for c in kids.get(st.pop(), []): s.append(c); st.append(c)
    return s
for pid in tree(ROOT_PID):
    try: os.kill(pid, signal.SIGCONT)   # SIGSTOP to suspend
    except ProcessLookupError: pass
```

**Check before resuming anything:** `pgrep -f 'omd mmgbsa' | wc -l` must be 0, and
`vm_stat` free+inactive should be comfortably above ~1 GB.

**The structural fix, not yet done:** `racc_drain.sh` launches its tail with no mutual exclusion, so
N drains means up to N concurrent MM/GBSA. It needs a lockfile around the tail invocation. **It could
not be edited** because four instances are live and zsh re-reads a running script by byte offset —
editing it mid-run corrupts execution. Do it in a quiet window, or write the lock into
`code/run_dock_pose_md.sh` and `~/python_mac/boltzgen_local/md/run_md20.sh` instead (same hazard:
check `pgrep` first).

---

## 2. Where the campaign is

29-leg fold-cofold MD campaign on racc (`gpuscavenger`, free, preemptible). **All 8 OF3 legs are
done and committed.** RF3 is the remaining lane.

**Jobs:** `oxb_2` and `oct_3` COMPLETED; `oxb_3` RUNNING (started ~05:26, 6 h limit); `oct_4`
submitted 06:36.

**Legs recorded in FOLDS.md: 10 RF3** — `orig_f12`, `s2_esm2_control`, `s3_esm2_f4`, `s3_orig_f12`
(octinoxate); `ox1_esm1_f8`, `ox1_orig_f12`, `ox2_esm1_f8`, `ox2_orig_f8`, `ox3_orig_f8`,
`ox2_shuffle_esm0` (oxybenzone).

**Fetched/tailing but NOT yet written up — do these next:**

| leg | status |
|---|---|
| `ox2_shuffle_esm0_rf3cof` | **tables NOT yet updated** — ΔG −11.4357 ±0.039, residence 32.6%, 38 released frames / 2 episodes / longest 370 ps / **0% late**, slide 1.71 Å, Cα 5.9 final / 7.0 max, Rg 9.8 → 10.7, contacts rise to 44.8, separation 10.0 Å. CSV rows for mmgbsa/slide/stability ARE written; only FOLDS.md tables + commit are missing. |
| `shuffle_control_esm0_rf3cof` | tail running (tree 27426) |
| `shuffle_control_rf3cof` | fetched, tail pending |
| `bg33_1_rf3cof` | fetched, tail pending |
| `bg33_2_rf3cof` | fetched, tail pending |
| `ox2_shuffle_rf3cof` | tail suspended (tree 27737) |

Then `oxb_3`'s remaining 4 and `oct_4`'s 2.

### Per-leg recording set (all of it, every leg)

```sh
PY=~/miniforge3/envs/openmm-md/bin/python
# 1. four windows + residence/release from <leg>/tail.log and md_contacts.csv
# 2. rows  (BOTH flags required now, no defaults)
~/miniforge3/bin/python code/collect_mmgbsa_rows.py <leg> --md-root runs/<mol>/md --outdir runs/<mol>
# 3. slide + stability + pairs: whole-table rewrites, NAME EVERY LEG, build lists in python not zsh
#    globs: runs/{octinoxate,oxybenzone}/md/*/{prod_20ns,prod_L1_modal}
#            + ~/python_mac/boltzgen_local/md/{bg33,bgox31}_*/prod_20ns   <-- omitting this drops rows
$PY code/pair_contacts.py <legs sharing a design...> --run-dir runs/<mol> --design <design*.json> --csv runs/<mol>/md/pair_contacts_rf3cof_<arm>.csv
# 4. audit: 4/4 windows, 4 summary rows, 1 stability, 1 slide, 1 pair row, frames present, 0 scatter
```

Designs: `orig_f12`/`bg33_*` → `design.json`; `s2_*` → `design_shell2.json`; `s3_*`,
`shuffle_control*` → `design_shell3.json`; `ox1_*` → `design_shell1.json`; `ox2_*` →
`design_shell2.json`; `ox3_*` → `design_shell3.json`. **No design json exists for the `bgox31` arm**
— its pair columns stay dashed.

---

## 3. The scientific state — read this before writing any conclusion

**Every mechanistic reading built during collection was falsified by a later leg.** In order:
shrinking margins (killed by `s3_orig_f12`), static separation predicts outcome (killed by
`s2_esm2_control`), floor-dependence and regression-to-the-mean (both killed by `ox2_orig_f8`). The
user's assessment, which the file now adopts: **patterns are not yet predictable.**

**One reading survives: RF3 beats OF3 on 4/4 octinoxate and 2/5 oxybenzone** — a molecule split, now
deconfounded from floor. Hold it loosely: it is the only survivor of five readings, which is about
what one of five does by chance.

**An open conflict, pre-registered in FOLDS.md and unresolved:** the two `bg33` RF3 folds have the
worst static reads in the campaign (`bg33_1` 0.30 enclosed / 7 of 20 engaged; `bg33_2` 0.345 / **4**)
against OF3 counterparts that each gained 3.6. The static read says they should be the weakest RF3
cells yet; the molecule split says they should still beat −14.90 and −11.58. **Both cannot hold.**
These two legs are the most informative thing left in the campaign — do not bury them.

**Also unresolved:** FOLDS.md and `$ODLOC/publications/peptidebuilder/paper_outline.md` §3.5 disagree
about which predictor is better. pepgraph says OF3 (px 10.7 vs 12.2, paired +3.6σ, 21 sequences);
dynamics say RF3 (4/4 octinoxate, 8 pairs, no replicates). pepgraph is the better-powered side. Do
not restate pepgraph's tables in FOLDS.md — cross-reference §3.5, one authoritative copy.

**Data-quality traps that bit, now flagged in-file:** three re-fold cells and one baseline slid so
far during equilibration that they are not tests of the predicted pose —
`ox2_shuffle_esm0_o3cof` 13.14 Å, `shuffle_control_esm0_o3cof` 10.62 Å, `shuffle_control_o3cof`
7.59 Å, and baseline `s3_esm2_f4` **7.81 Å** (so every Δ for that peptide is compromised, marked ⚠).
Check `runs/<mol>/md/ligand_slide.csv` before reading any cell. Also: residence below ~60% is
uninterpretable without its separation trace — four cells have low residence purely because the
ligand parks just outside the 10 Å cutoff with no real detachment.

---

## 4. Figures — done, working, one cosmetic question open

`code/render_md_figures.py`, headless PyMOL from its own conda env:

```sh
~/miniforge3/bin/mamba create -y -n pymol-render -c conda-forge python=3.12 pymol-open-source
~/miniforge3/envs/pymol-render/bin/pip install pillow
```

84 medoid PNGs rendered (42 + 42). Defaults: `violetpurple` peptide, `skyblue` ligand, `cartoon
automatic` after **`dss`** (the PDBs carry no HELIX/SHEET records, so without it every helix renders
as coil), 40% cartoon transparency, orient+zoom on the ligand, buffer 9, transparent RGBA, 1600×1200,
skip-existing, medoids only. **Wired into both tails**, so new legs render themselves.

**Open, cosmetic only:** ~¾ of the median frame is empty. Batch crop is a no-op because half the set
touches a canvas edge; tightening the buffer makes it strictly worse (clips fold); `--crop each`
breaks the 4-across markdown grids by force-scaling. The user liked "buffer 12 then centre/crop",
which is viable only if no image touches an edge at 12 — **unverified across all 83, and a previous
agent ran the buffer search away to 24.** Cap any search at ~14 and special-case outlier folds
instead. Buffer 9 is perfectly serviceable; leaving it is a fine answer.

**8 leftover `_buf*.png` comparison files** in `runs/octinoxate/figures_md/` can be deleted.

**Not started:** embedding the 84 images into FOLDS.md / RESULTS.md. RESULTS.md has an established
convention to copy — 4-across markdown table, image row, then a numbers row, then a caption row
(see RESULTS.md ~line 752).

---

## 5. Cost warning

**22% of allocation in 30 minutes**, caused mainly by a subagent reading many images across six
review rounds. If delegating visual work again: have the agent verify by file size and exit status,
read at most one or two images, and send candidates to the user to judge rather than inspecting them
all itself. Do not spawn background shells that poll other background shells — one agent accumulated
seven simultaneous waiters.

---

## 6. Uncommitted

`code/run_dock_pose_md.sh` (tail render block) and `~/python_mac/boltzgen_local/md/run_md20.sh`
(same, **different repo**) are modified and uncommitted. Four pre-existing modified files from before
this session: `code/build_cofolds.sh`, `code/dock_folds.py`, `code/racc_drain.sh`, `code/racc_feed.sh`.
Everything else from the last two days is committed.

Attribution for commits:
```
Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01ALr2nsS4EbqzD3LX3xTB4x
```
