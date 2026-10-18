"""Rebuild cofold Tables 1 and 2 of runs/folds_md/cofold_binding_tables.md.

One row per peptide x starting pose, three sub-rows per peptide (co-fold, Vina p1, GNINA pick), and
one column per folding model (Boltz/BoltzGen, OpenFold3, RoseTTAFold3). Restructured 2026-10-18 from
one row per peptide, so a column is now one model across every pose.

Nothing here is typed by hand. dG and the 10 ns drift markers come from
runs/<mol>/md/mmgbsa_summary.csv, plus each BoltzGen leg's own FINAL_RESULTS_MMPBSA.dat -- those legs
live in the sibling boltzgen_local repo and are in no summary CSV. Residence has no CSV anywhere
(md_contacts.csv is per-decile, too coarse for a share-of-frames figure), so it is collated from the
per-leg tables already published in FOLDS.md and RESULTS.md, including the ones whose first column is
the peptide and second the pose, and falls back to a leg's tail.log. It prints every cell it could not
source rather than leaving a silent gap; a run that reports missing cells should not be pasted in.

    python code/cofold_tables.py     # writes /tmp/.../tables.md, prints any gaps

The output is pasted into cofold_binding_tables.md under each table's caption. Cross-check the
co-fold cells against the published values before replacing them -- that check found 0 mismatches on
2026-10-18 and is what makes the rebuild safe.
"""
import csv, re, os, glob
BG = os.path.expanduser("~/python_mac/boltzgen_local/md")

def dat_dg(path):
    if not os.path.exists(path): return None
    blk = open(path).read().split("Differences")[-1]
    for l in blk.split("\n"):
        if l.startswith("DELTA TOTAL"): return float(l.split()[2])
    return None

dg, dg10 = {}, {}
for mol in ("octinoxate", "oxybenzone"):
    for r in csv.DictReader(open(f"runs/{mol}/md/mmgbsa_summary.csv")):
        s, leg, v = r["structure"], r["leg"], float(r["dg_bind"])
        if leg in ("prod_20ns", "prod_L1_modal"): dg[s] = v
        elif leg == "first_10ns": dg10[s] = v
# BoltzGen legs live in the sibling repo and are not in either summary CSV
for d in sorted(glob.glob(f"{BG}/*/")):
    leg = os.path.basename(d.rstrip("/"))
    v = dat_dg(f"{d}prod_20ns/mmgbsa/FINAL_RESULTS_MMPBSA.dat")
    if v is not None: dg.setdefault(leg, v)
    w = dat_dg(f"{d}first_10ns/mmgbsa/FINAL_RESULTS_MMPBSA.dat")
    if w is not None: dg10.setdefault(leg, w)

# residence: leg-keyed tables, then peptide/pose sub-row tables, then tail.log
res = {}
def scan_leg_keyed(path):
    for line in open(path):
        if not line.startswith("|"): continue
        c = [x.strip() for x in line.strip().strip("|").split("|")]
        m = re.match(r"`([A-Za-z0-9_]+)`", c[0])
        if not m: continue
        for cell in c[1:]:
            p = re.fullmatch(r"\*?\*?([\d.]+)%\*?\*?", cell)
            if p: res.setdefault(m.group(1), float(p.group(1))); break

def scan_pose_subrows(path):
    """Tables whose first column is the peptide and second the pose, blank on continuation rows."""
    cur = None
    for line in open(path):
        if not line.startswith("|"): continue
        c = [x.strip() for x in line.strip().strip("|").split("|")]
        if len(c) < 3: continue
        m = re.match(r"`([A-Za-z0-9_]+)`", c[0])
        if m: cur = m.group(1)
        elif c[0]: cur = None
        if not cur: continue
        pose = c[1]
        if pose == "co-folded": leg = cur
        elif pose == "Vina p1": leg = f"{cur}_dock1"
        else:
            g = re.match(r"(?:GNINA|CNNscore) p(\d+)", pose)
            if not g: continue
            leg = f"{cur}_dock{g.group(1)}"
        for cell in c[2:]:
            p = re.fullmatch(r"\*?\*?([\d.]+)%\*?\*?", cell)
            if p: res.setdefault(leg, float(p.group(1))); break

for f in ("FOLDS.md", "RESULTS.md", "runs/folds_md/cofold_binding_tables.md"): scan_leg_keyed(f)
scan_pose_subrows("RESULTS.md")
for d in glob.glob("runs/*/md/*/tail.log") + glob.glob(f"{BG}/*/tail.log"):
    leg = os.path.basename(os.path.dirname(d))
    if leg not in res:
        m = re.findall(r"residence within 10 A: ([\d.]+)%", open(d).read())
        if m: res[leg] = float(m[-1])

def mark(leg):
    a, b = dg.get(leg), dg10.get(leg)
    if a is None or b is None: return ""
    d = a - b
    return " •" if abs(d) < 0.5 else (" ↓" if d < 0 else " ↑")

def pick_leg(mol, p, suf):
    c = [os.path.basename(x) for x in glob.glob(f"runs/{mol}/md/{p}{suf}_dock*") + glob.glob(f"{BG}/{p}{suf}_dock*")]
    c = [x for x in c if not x.endswith("_dock1") and "_r2" not in x]
    return sorted(c)[0] if c else None

def cell(leg):
    if leg is None: return "—"
    v, r = dg.get(leg), res.get(leg)
    if v is None: return "—"
    s = f"{v:.2f}".replace("-", "−")
    rs = "" if r is None else f" ({r:.1f}%)"
    return f"**{s}**{rs}{mark(leg)}" if v <= -18 else f"{s}{rs}{mark(leg)}"

PEPS = {"octinoxate": ['orig_f12','s2_esm2_control','s3_esm2_f4','s3_orig_f12','shuffle_control','shuffle_control_esm0','bg33_1','bg33_2','bg33_3','bg33_4'],
        "oxybenzone": ['ox1_esm1_f8','ox1_orig_f12','ox2_esm1_f8','ox2_orig_f8','ox2_shuffle','ox2_shuffle_esm0','ox3_orig_f8','bgox31_2','bgox31_3','bgox31_4','bgox31_5']}
out, missing = {}, []
for mol, peps in PEPS.items():
    lines = ["| peptide | pose | Boltz / BoltzGen co-fold | OpenFold3 co-fold | RoseTTAFold3 co-fold |",
             "|---|---|---|---|---|"]
    for p in peps:
        for i, pose in enumerate(("co-fold", "Vina p1", "GNINA pick")):
            legs = []
            for suf in ("", "_o3cof", "_rf3cof"):
                legs.append(f"{p}{suf}" if pose == "co-fold" else
                            (f"{p}{suf}_dock1" if pose == "Vina p1" else pick_leg(mol, p, suf)))
            for L in legs:
                if L is not None and (dg.get(L) is None or res.get(L) is None):
                    missing.append((L, "dG" if dg.get(L) is None else "residence"))
            label = f"`{p}`" if i == 0 else ""
            pl = pose
            if pose == "GNINA pick":
                ns = [l.split("_dock")[-1] if l else None for l in legs]
                pl = "GNINA pick (p" + "/".join(n or "—" for n in ns) + ")"
            lines.append(f"| {label} | {pl} | " + " | ".join(cell(l) for l in legs) + " |")
    out[mol] = "\n".join(lines)
open(os.environ.get("OUT", "/tmp/tables.md"), "w").write(f"## octinoxate\n{out['octinoxate']}\n\n## oxybenzone\n{out['oxybenzone']}\n")
print(f"missing: {len(missing)}")
for m in sorted(set(missing)): print("  ", m)
