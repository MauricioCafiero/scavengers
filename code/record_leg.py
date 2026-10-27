"""Record one finished MD leg into every table, then audit it (CLAUDE.md, "Audit the recording set").

Orchestrates the existing scripts; it does not reimplement them:
  collect_mmgbsa_rows.py  -> runs/<mol>/md/mmgbsa_summary.csv   (appends this leg's 4 rows)
  ligand_slide.py         -> runs/<mol>/md/ligand_slide.csv     (whole-table rewrite, every lane)
  md_stability.py         -> runs/octinoxate/md/protein_stability.csv (whole-table rewrite, every lane)
  pair_contacts.py        -> the pair_contacts_*.csv that already holds the leg's parent structure
                             (<pep>_o3cof / <pep>_rf3cof), re-run over every leg already in it plus
                             this one, with that file's own design json.
The two whole-table rewrites go to a scratch copy first and are refused if any existing row would be
lost. Prints a one-paragraph summary of the leg for the write-up.

    ~/miniforge3/envs/openmm-md/bin/python code/record_leg.py <mol> <leg>
"""
import csv
import glob
import os
import re
import shutil
import subprocess
import sys
import tempfile

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BG = os.path.expanduser("~/python_mac/boltzgen_local/md")
sys.path.insert(0, os.path.join(REPO, "code"))
os.chdir(REPO)


def dg_total(f):
    block = open(f).read().split("Differences")[1]
    return float(next(l for l in block.split("\n") if l.startswith("DELTA TOTAL")).split()[2])


def rows_by(path, key):
    return {r[key] for r in csv.DictReader(open(path))} if os.path.exists(path) else set()


def upsert(path, key, new):
    """Put this leg's row(s) into a table: replace in place if the key is there, else append.

    Every other row is copied through untouched. Until 2026-10-27 recording re-ran ligand_slide.py,
    md_stability.py and pair_contacts.py over EVERY leg in every lane and rewrote the whole file --
    ~8 min of re-reading ~150 trajectories to add one row whose neighbours cannot change (fixed
    trajectory, deterministic code), and the source of every dropped-row incident in CLAUDE.md.
    A full rebuild is still available by running those scripts by hand over the lane globs.
    """
    rows = list(csv.DictReader(open(path))) if os.path.exists(path) else []
    fields = list(rows[0].keys()) if rows else list(new[0].keys())
    by_key = {r[key]: r for r in new}
    out = [by_key.pop(r[key], r) for r in rows] + list(by_key.values())
    tmp = tempfile.mktemp(suffix=".csv", dir=os.path.dirname(path))
    with open(tmp, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        w.writerows(out)
    assert rows_by(path, key) <= rows_by(tmp, key), f"{path} lost rows"
    shutil.move(tmp, path)


def one_leg_rows(run):
    """Run a table writer on this leg alone into a temp csv and return its rows."""
    tmp = tempfile.mktemp(suffix=".csv")
    run(tmp)
    rows = list(csv.DictReader(open(tmp)))
    os.remove(tmp)
    return rows


def main(mol, leg):
    import ligand_slide
    import md_stability
    import pair_contacts

    d = f"runs/{mol}/md/{leg}"
    prod = f"{d}/prod_20ns"
    py = sys.executable
    # a replicate (<leg>_r2, _r3, ...) files its MM/GBSA rows in mmgbsa_replicates.csv, never the
    # canonical summary, which holds one row set per leg. This tested `endswith("_r2")` until
    # 2026-10-21, so the first THIRD run (orig_f12_dock1_r3) put four rows into octinoxate's
    # mmgbsa_summary.csv -- and replicate_tally.py, which reads run 3 from the replicates table,
    # reported the triple as still pending while the canonical table silently carried a second
    # value for a published leg. Match any run number.
    summary = "mmgbsa_replicates.csv" if re.search(r"_r\d+$", leg) else "mmgbsa_summary.csv"

    subprocess.run([os.path.expanduser("~/miniforge3/bin/python"), "code/collect_mmgbsa_rows.py", leg,
                    "--md-root", f"runs/{mol}/md", "--outdir", f"runs/{mol}", "--summary", summary], check=True,
                   stdout=subprocess.DEVNULL)

    upsert(f"runs/{mol}/md/ligand_slide.csv", "leg",
           one_leg_rows(lambda out: ligand_slide.main([prod, "--csv", out])))
    st_row = md_stability.stats(prod)
    st_row["peptide"], st_row["pose"] = md_stability.label(prod)
    upsert("runs/octinoxate/md/protein_stability.csv", "leg", [st_row])

    # strip the run suffix first, then the pose suffix: _r\d+ not _r2, so a third run resolves to the
    # same parent structure its earlier runs did rather than to a name no pair table has ever held
    parent = re.sub(r"_dock\d+$", "", re.sub(r"_r\d+$", "", leg))
    pair = None
    for f in glob.glob(f"runs/{mol}/md/pair_contacts_*.csv"):
        rs = list(csv.DictReader(open(f)))
        if any(r["structure"] == parent for r in rs):
            pair = (f, rs)
            break
    if pair is None:
        # First leg of a new co-fold lane: there is no table to extend, so start one. Until
        # 2026-10-19 this fell through to "no design json for this arm", which is simply wrong --
        # pair_contacts.design_for() resolves the json by peptide name and the files all exist. The
        # message was read as "nothing to do" and all five legs of the pocket-guided lane landed
        # with no pair row at all. Name the file the way the lane's others are named:
        # pair_contacts_<lane>_<arm>.csv, lane off the parent's suffix, arm off the design json.
        lane = parent.rsplit("_", 1)[-1]
        if lane.endswith("cof"):
            design = pair_contacts.design_for(parent)
            arm = re.sub(r"^design(_)?", "", design[:-len(".json")]) or "design"
            f = f"runs/{mol}/md/pair_contacts_{lane}_{arm}.csv"
            # The file may ALREADY EXIST and hold another peptide: the name keys on the design arm,
            # not the peptide, and two peptides of the same shell share an arm (s3_esm2_f4 and
            # s3_orig_f12 are both design_shell3). Starting it from [] then rewrote it with only the
            # new leg -- and `before <= after` passed vacuously because before was empty. That cost
            # five rows across three files on 2026-10-19/20 before a `git diff --numstat` caught it.
            pair = (f, list(csv.DictReader(open(f))) if os.path.exists(f) else [])
    pair_line = "no pair table (not a co-fold lane leg)"
    if pair:
        f, rs = pair
        legs = [r["leg"] for r in rs if r["structure"] != leg] + [prod]
        design = rs[0]["design"] if rs else pair_contacts.design_for(parent)
        before = {r["structure"] for r in rs}
    # design_for() falls through to the unnumbered "design.json" for any name it does not
    # recognise. That file exists for octinoxate and NOT for oxybenzone, so a BoltzGen structure
    # on a co-fold lane (bgox31_4_o3cof_dock7_r2, 2026-10-21) crashed the whole recording here
    # with FileNotFoundError -- after the MM/GBSA rows had already been appended, leaving the leg
    # half-recorded. The comment in design_for() documents this exact fall-through being fixed for
    # ox*_ names; generator-arm names still reach it.
    #
    # Skipping is also the right answer on the merits, not just a way to avoid the crash: the pair
    # measure is entirely slot-derived, a BoltzGen generator arm has no design slots, and its
    # numbers were retracted on 2026-10-18 for that reason. Borrowing another molecule's design
    # would compute a plausible-looking number that means nothing.
    #
    # The open question above -- whether the bg33 rows that DO exist (octinoxate, where design.json
    # resolves and so nothing crashed) should be skipped too -- was answered by the user on
    # 2026-10-22: "No generator arms should have pair rows." So the skip is by ARM, not by whether
    # the design file happens to resolve, and those 20 rows were struck from the committed tables
    # that day. The authoritative test is pair_contacts.is_generator_arm().
    if pair and pair_contacts.is_generator_arm(parent):
        pair_line = (f"no pair table ({parent} is a generator arm -- no designed slots, "
                     f"so the pair measure does not apply)")
        pair = None
    if pair and not os.path.exists(os.path.join(f"runs/{mol}", design)):
        pair_line = (f"no pair table (no {design} for {mol}; "
                     f"{parent} has no design arm -- slot measures do not apply)")
        pair = None
    if pair:
        upsert(f, "structure", one_leg_rows(lambda out: subprocess.run(
            [py, "code/pair_contacts.py", prod, "--run-dir", f"runs/{mol}", "--design", design, "--csv", out],
            check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)))
        after = {r["structure"]: r for r in csv.DictReader(open(f))}
        assert before <= set(after), f"pair table lost rows: {before - set(after)}"
        r = after[leg]
        pair_line = f"pairs {r['pairs_realised']} (band {r['band_position_pct']}%), {r['dead_slots']} dead slots"

    # summary + audit
    win = []
    for n in (5, 10, 15):
        fs = glob.glob(f"{d}/first_{n}ns/**/FINAL_RESULTS_MMPBSA.dat", recursive=True)
        win.append(dg_total(fs[0]) if fs else None)
    # "/last_" as well as "/first_": a trailing diagnostic window is not the whole run, and because
    # glob order put last_10ns first, this leg reported its back-half value (-18.72) as its headline
    # dG in place of -15.03 on 2026-10-21. Exclude every window directory, not the leading ones.
    whole = [f for f in glob.glob(f"{d}/**/FINAL_RESULTS_MMPBSA.dat", recursive=True)
             if "/first_" not in f and "/last_" not in f]
    win.append(dg_total(whole[0]) if whole else None)
    log = open(f"{d}/tail.log").read()
    res = re.search(r"residence within 10 A: ([\d.]+)%", log).group(1)
    rel = re.search(r"release .*?: (\d+) frames \(([\d.]+)%\) in (\d+) episodes, longest ([\d.]+) ps.*?(\d+)% of it in the second half", log)
    slide = next(r["ligand_slide_A"] for r in csv.DictReader(open(f"runs/{mol}/md/ligand_slide.csv")) if r["leg"] == leg)
    st = next(r for r in csv.DictReader(open("runs/octinoxate/md/protein_stability.csv")) if r["leg"] == prod)
    scatter = [p for p in glob.glob(f"{d}/**/*", recursive=True)
               if os.path.basename(p) in ("reference.frc", "traj.nc") or os.path.basename(p).startswith("_MMPBSA_")]
    n_sum = sum(1 for r in csv.DictReader(open(f"runs/{mol}/md/{summary}")) if r["structure"] == leg)
    frames = all(os.path.exists(f"{prod}/{x}") for x in ("frames_last.pdb", "frame_medoid.pdb"))
    print(f"{leg}: dG {win[3]:.2f} | windows " + " / ".join(f"{w:.2f}" for w in win)
          + f" | residence {res}% | released {rel.group(1)} ({rel.group(2)}%) in {rel.group(3)}, longest "
          f"{rel.group(4)} ps, {rel.group(5)}% late | slide {slide} A | Ca {st['ca_rmsd_final_A']}/"
          f"{st['ca_rmsd_max_A']} | {pair_line}")
    ok = len(whole) + sum(w is not None for w in win[:3]) == 4 and n_sum == 4 and frames and not scatter
    # FOLDS.md: one row per docked fold leg, appended under its table (header written once by hand)
    m = re.match(r"(.+_(?:o3p?|rf3)cof)_dock(\d+)$", leg)
    if m and ok:
        why = next((l.split(None, 3)[3].strip() for l in open("runs/folds_md/dock_legs.txt")
                    if l.split()[:3] == [mol, m.group(1), m.group(2)]), None)
        # An unregistered leg used to default to "GNINA", which silently mislabelled the pocket-guided
        # lane's two Vina p1 legs on 2026-10-19 -- the lane had no dock_legs.txt entries at all. A
        # missing entry is unknown provenance, not a GNINA pick. Do not exit here: the MM/GBSA rows
        # were already appended above, so a re-run after fixing the file would duplicate them. Write
        # a row whose pose cell is visibly wrong instead, and say so on stderr.
        if why is None:
            print(f"WARNING {leg}: no entry in runs/folds_md/dock_legs.txt for {mol} {m.group(1)} "
                  f"pose {m.group(2)} -- pose cell written as UNREGISTERED. Add "
                  f"'<mol> <struct> <pose> vina|gnina [un]corroborated' and edit the cell by hand.",
                  file=sys.stderr)
        pose = (f"p{m.group(2)} UNREGISTERED" if why is None else
                f"p{m.group(2)} " + ("Vina" if why == "vina"
                                     else "GNINA" + (" ✓" if why.endswith(" corroborated") else "")))
        pdg = next(float(r["dg_bind"]) for r in csv.DictReader(open(f"runs/{mol}/md/mmgbsa_summary.csv"))
                   if r["structure"] == m.group(1) and not r["leg"].startswith("first_"))
        drift = win[3] - win[1]
        mark = "•" if abs(drift) < 0.5 else ("↓" if drift < 0 else "↑")
        row = (f"| `{leg}` | {mol[:3]} | {pose} | **{win[3]:.2f}** {mark} | {res}% | {rel.group(1)} | {slide} | "
               f"{pdg:.2f} | {win[3] - pdg:+.2f} |\n").replace("-", "−").replace("`" + leg.replace("-", "−"), "`" + leg)
        txt = open("FOLDS.md").read()
        hdr = "| leg | mol | pose | ΔG 20 ns | residence | released | slide Å | cofold-pose ΔG | Δ vs cofold |\n|---|---|---|---|---|---|---|---|---|\n"
        i = txt.index(hdr) + len(hdr)
        while txt.startswith("| `", i):
            i = txt.index("\n", i) + 1
        if f"| `{leg}` |" not in txt:
            open("FOLDS.md", "w").write(txt[:i] + row + txt[i:])
    print(f"audit: windows {sum(w is not None for w in win)}/4, summary rows {n_sum}, frames {frames}, "
          f"scatter {len(scatter)} -> {'OK' if ok else 'INCOMPLETE'}")


if __name__ == "__main__":
    main(*sys.argv[1:3])
