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


def rewrite_checked(path, key, run):
    """Run a whole-table writer into a temp copy; install it only if no existing row is lost."""
    tmp = tempfile.mktemp(suffix=".csv")
    run(tmp)
    lost = rows_by(path, key) - rows_by(tmp, key)
    if lost:
        sys.exit(f"REFUSED {path}: would drop {sorted(lost)}")
    shutil.move(tmp, path)


def main(mol, leg):
    import ligand_slide
    import md_stability
    import pair_contacts

    d = f"runs/{mol}/md/{leg}"
    prod = f"{d}/prod_20ns"
    py = sys.executable
    # a replicate (<leg>_r2) files its MM/GBSA rows in mmgbsa_replicates.csv, not the canonical summary
    summary = "mmgbsa_replicates.csv" if leg.endswith("_r2") else "mmgbsa_summary.csv"

    subprocess.run([os.path.expanduser("~/miniforge3/bin/python"), "code/collect_mmgbsa_rows.py", leg,
                    "--md-root", f"runs/{mol}/md", "--outdir", f"runs/{mol}", "--summary", summary], check=True,
                   stdout=subprocess.DEVNULL)

    lanes = lambda m: sorted(glob.glob(f"runs/{m}/md/*/prod_20ns") + glob.glob(f"runs/{m}/md/*/prod_L1_modal"))
    bgpat = "bg33_*" if mol == "octinoxate" else "bgox31_*"
    rewrite_checked(f"runs/{mol}/md/ligand_slide.csv", "leg",
                    lambda out: ligand_slide.main(lanes(mol) + sorted(glob.glob(f"{BG}/{bgpat}/prod_20ns"))
                                                  + ["--csv", out]))
    stab_legs = lanes("octinoxate") + lanes("oxybenzone") + sorted(glob.glob(f"{BG}/*/prod_20ns"))
    md_stability.main(stab_legs)   # refuses to shrink on its own

    parent = re.sub(r"_dock\d+$", "", re.sub(r"_r2$", "", leg))
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
            pair = (f"runs/{mol}/md/pair_contacts_{lane}_{arm}.csv", [])
    pair_line = "no pair table (not a co-fold lane leg)"
    if pair:
        f, rs = pair
        legs = [r["leg"] for r in rs if r["structure"] != leg] + [prod]
        design = rs[0]["design"] if rs else pair_contacts.design_for(parent)
        before = {r["structure"] for r in rs}
        subprocess.run([py, "code/pair_contacts.py", *legs, "--run-dir", f"runs/{mol}", "--design", design,
                        "--csv", f], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        after = {r["structure"]: r for r in csv.DictReader(open(f))}
        assert before <= set(after), f"pair table lost rows: {before - set(after)}"
        r = after[leg]
        pair_line = f"pairs {r['pairs_realised']} (band {r['band_position_pct']}%), {r['dead_slots']} dead slots"

    # summary + audit
    win = []
    for n in (5, 10, 15):
        fs = glob.glob(f"{d}/first_{n}ns/**/FINAL_RESULTS_MMPBSA.dat", recursive=True)
        win.append(dg_total(fs[0]) if fs else None)
    whole = [f for f in glob.glob(f"{d}/**/FINAL_RESULTS_MMPBSA.dat", recursive=True) if "/first_" not in f]
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
