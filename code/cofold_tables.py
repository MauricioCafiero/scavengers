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
import csv, re, os, glob, sys
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

# residence: the leg's own record log first, then leg-keyed tables, then peptide/pose sub-row
# tables, then tail.log. Every source below uses setdefault, so first write wins and the order
# here is the priority order.
#
# The record log has to come first because the doc scans are too permissive to trust: they take
# the first `N%` cell in any leg-keyed row of FOLDS.md / RESULTS.md, and RESULTS.md's *static*
# fold-geometry tables are keyed by leg too. On 2026-10-19 that harvested `s3_esm2_control`'s
# 100% from a fold-geometry row (RESULTS.md:1159) and printed it as MD residence beside the MD
# dG, against the leg's real 29.4%. record_leg.py writes `residence NN.N%` per leg from
# md_contacts.py's own profile, so it is the pipeline's answer rather than a guess at which
# column of which table meant retention.
res = {}
for _f in glob.glob("runs/folds_md/record_*.log"):
    _leg = re.sub(r"^record_|\.log$", "", os.path.basename(_f))
    _m = re.findall(r"residence ([\d.]+)%", open(_f, errors="replace").read())
    if _m: res.setdefault(_leg, float(_m[-1]))
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

# A GNINA-pick cell can be empty for two unrelated reasons, and rendering both the same way was
# wrong (caught 2026-10-18). runs/folds_md/dock_legs.txt is the campaign's own record of which pose
# each structure's pick was, and its convention is the one CLAUDE.md states: a structure gets a
# `gnina` line only when the pick differs from pose 1, because "if the pick is pose 1 it coincides
# with the Vina-p1 leg, so no new leg". So a `vina`-only entry means the pick WAS pose 1 -- not that
# GNINA never ran. The scored poses themselves live in the gnina_local repo and only some were
# pulled into runs/<mol>/gnina/, so the local scores file is not evidence of absence; this function
# must consult dock_legs.txt first or it will report 44 structures as unscored that were all scored.
def gnina_pick(mol, struct):
    """(pose, status): 'ok' = its own leg, 'is_p1' = the Vina-p1 leg, 'unknown' = no record."""
    for line in open("runs/folds_md/dock_legs.txt"):
        w = line.split()
        if len(w) >= 4 and w[0] == mol and w[1] == struct and w[3] == "gnina":
            return int(w[2]), "ok"
    for line in open("runs/folds_md/dock_legs.txt"):
        w = line.split()
        if len(w) >= 4 and w[0] == mol and w[1] == struct and w[3] == "vina":
            return 1, "is_p1"
    f = f"runs/{mol}/gnina/results_receptor_score_only/gnina_scores.csv"
    if os.path.exists(f):
        rs = [r for r in csv.DictReader(open(f))
              if r["structure"] == struct and int(r["pose"]) >= 1 and r.get("CNNaffinity")]
        if rs: return int(max(rs, key=lambda r: float(r["CNNaffinity"]))["pose"]), "ok"
    return None, "unknown"


def pick_leg(mol, p, suf):
    c = [os.path.basename(x) for x in glob.glob(f"runs/{mol}/md/{p}{suf}_dock*") + glob.glob(f"{BG}/{p}{suf}_dock*")]
    c = [x for x in c if not x.endswith("_dock1") and "_r2" not in x]
    return sorted(c)[0] if c else None

def cell(leg, kind=None):
    if leg is None: return kind or "—"
    v, r = dg.get(leg), res.get(leg)
    if v is None:
        return "*running*" if os.path.exists(os.path.join("runs", MOL, "md", leg, "system", "system.xml")) else "—"
    s = f"{v:.2f}".replace("-", "−")
    rs = "" if r is None else f" ({r:.1f}%)"
    return f"**{s}**{rs}{mark(leg)}" if v <= -18 else f"{s}{rs}{mark(leg)}"

# The Boltz column must be UNFORCED to be comparable with the OF3 and RF3 columns, which were run
# with no constraints at all. That is the whole reason the unforced campaign exists: the published
# version of these tables set hint-forced Boltz legs (_f4/_f8/_f12) against unforced OF3/RF3, so a
# row mixed "which model" with "forced or not". Each forced peptide therefore maps to the hint-free
# Boltz fold of the SAME sequence (its *_control twin); a peptide that was already unforced -- the
# shuffles, s2_esm2_control, and every BoltzGen bg* design, which the generator makes without hints
# -- maps to itself and its published numbers stand.
UNFORCED = {'orig_f12': 'orig_control', 's3_esm2_f4': 's3_esm2_control',
            's3_orig_f12': 's3_orig_control', 'ox1_esm1_f8': 'ox1_esm1_control',
            'ox1_orig_f12': 'ox1_orig_control', 'ox2_esm1_f8': 'ox2_esm1_control',
            'ox2_orig_f8': 'ox2_orig_control', 'ox3_orig_f8': 'ox3_orig_control'}

# The pocket-guided OpenFold3 lane (NEXT_STEPS.md Part 2) covers only the nine non-null, non-BoltzGen
# rows of these two tables -- one guided fold per peptide that has a designed shell to take a pocket
# from. The nulls deliberately get none (a shuffled sequence has no designed pocket to restrain to,
# and two ways of manufacturing one were considered and declined on 2026-10-19), and the BoltzGen
# designs are made without hints by construction. So a blank `_o3pcof` cell on those rows means "not
# attempted, by design", which is why those legs must not be built into `missing` below.
GUIDED = {'orig_f12', 's2_esm2_control', 's3_esm2_f4', 's3_orig_f12',
          'ox1_esm1_f8', 'ox1_orig_f12', 'ox2_esm1_f8', 'ox2_orig_f8', 'ox3_orig_f8'}
SUFFIXES = ("", "_o3cof", "_rf3cof", "_o3pcof")

PEPS = {"octinoxate":['orig_f12','s2_esm2_control','s3_esm2_f4','s3_orig_f12','shuffle_control','shuffle_control_esm0','bg33_1','bg33_2','bg33_3','bg33_4'],
        "oxybenzone": ['ox1_esm1_f8','ox1_orig_f12','ox2_esm1_f8','ox2_orig_f8','ox2_shuffle','ox2_shuffle_esm0','ox3_orig_f8','bgox31_2','bgox31_3','bgox31_4','bgox31_5']}
out, missing, pending, unscored = {}, [], [], []
for mol, peps in PEPS.items():
    MOL = mol
    lines = ["| peptide | pose | Boltz / BoltzGen co-fold | OpenFold3 co-fold | RoseTTAFold3 co-fold "
             "| OpenFold3 pocket-guided |",
             "|---|---|---|---|---|---|"]
    for p in peps:
        for i, pose in enumerate(("co-fold", "Vina p1", "GNINA pick")):
            legs = []
            for suf in SUFFIXES:
                # the Boltz column (suf == "") takes the unforced twin; OF3/RF3 were already unforced
                base = UNFORCED.get(p, p) if suf == "" else p
                if suf == "_o3pcof" and p not in GUIDED:
                    legs.append(None)
                    continue
                legs.append(f"{base}{suf}" if pose == "co-fold" else
                            (f"{base}{suf}_dock1" if pose == "Vina p1" else pick_leg(mol, base, suf)))
            for L in legs:
                if pose == "GNINA pick" and L is None: continue
                if L is None or dg.get(L) is not None and res.get(L) is not None: continue
                # a staged-but-unfinished leg of the running unforced campaign is pending, not missing
                if os.path.exists(f"runs/{mol}/md/{L}/system/system.xml") and dg.get(L) is None:
                    pending.append(L)
                else:
                    missing.append((L, "dG" if dg.get(L) is None else "residence"))
            # Row label is the SEQUENCE, so drop the _fN fold-version suffix on a remapped
            # peptide: the row's Boltz cell is <stem>_control and its OF3/RF3 cells are
            # <stem>_fN_{o3,rf3}cof, all three unforced folds of one sequence. The _fN survives in
            # the OF3/RF3 leg names only because they were named after the forced Boltz structure
            # whose sequence they refold -- no constraint was applied to them.
            stem = re.sub(r"_f\d+$", "", p) if p in UNFORCED else p
            label = f"`{stem}`" if i == 0 else ""
            pl = pose
            if pose == "GNINA pick":
                tags = []
                for suf, l in zip(SUFFIXES, legs):
                    # a lane the peptide was never folded in has no pick to report and is not
                    # "unscored" -- without this the nulls and bg* designs would all land in the
                    # never-GNINA-scored list purely for having no pocket-guided fold
                    if suf == "_o3pcof" and p not in GUIDED: tags.append("—"); continue
                    struct = (UNFORCED.get(p, p) if suf == "" else p) + suf
                    pk, st = gnina_pick(mol, struct)
                    if l: tags.append(f"p{l.split('_dock')[-1]}")
                    elif st == "is_p1": tags.append("p1")
                    else: tags.append("?"); unscored.append(struct)
                pl = "GNINA pick (" + "/".join(tags) + ")"
            kinds = [None] * len(SUFFIXES)
            if pose == "GNINA pick":
                kinds = []
                for suf, l in zip(SUFFIXES, legs):
                    if suf == "_o3pcof" and p not in GUIDED: kinds.append("—"); continue
                    struct = (UNFORCED.get(p, p) if suf == "" else p) + suf
                    pk, st = gnina_pick(mol, struct)
                    kinds.append("*no record*" if st == "unknown" else "see Vina p1")
            lines.append(f"| {label} | {pl} | " + " | ".join(cell(l, k) for l, k in zip(legs, kinds)) + " |")
    out[mol] = "\n".join(lines)
open(os.environ.get("OUT", "/tmp/tables.md"), "w").write(f"## octinoxate\n{out['octinoxate']}\n\n## oxybenzone\n{out['oxybenzone']}\n")


def splice(path=os.path.join("runs", "folds_md", "cofold_binding_tables.md")):
    """Replace Tables 1 and 2 in the markdown in place, refusing if a published cell would change.

    Both halves of this were manual until 2026-10-19: the docstring said to paste the blocks under
    each caption and to "cross-check the co-fold cells against the published values before replacing
    them", which is the check that found 0 mismatches on 2026-10-18. An optional check is a skipped
    check (CLAUDE.md), and a manual paste is a step that gets forgotten -- the pocket-guided lane's
    first six legs sat in the CSVs with no row in this markdown at all. So the paste and the check
    are one function here, and `--write` is how the markdown gets updated from now on.

    What may legitimately differ, and nothing else: cells in a newly added column, and the
    GNINA-pick row label, which carries one pose tag per column and so gains a slot when a column
    is added. Every other published cell must match exactly, compared cell-by-cell after stripping
    -- a whitespace-sensitive diff of whole rows reports every row as changed and says nothing.
    """
    txt = open(path).read()
    lines = txt.split("\n")
    heads = [i for i, l in enumerate(lines) if l.startswith("| peptide | pose |")]
    if len(heads) != 2:
        raise SystemExit(f"splice: expected 2 'peptide | pose' tables in {path}, found {len(heads)}")
    changed, newcols = [], 0
    for start, mol in sorted(zip(heads, ["octinoxate", "oxybenzone"]), reverse=True):
        end = start
        while end < len(lines) and lines[end].startswith("|"):
            end += 1
        old, new = lines[start:end], out[mol].split("\n")
        if len(old) != len(new):
            raise SystemExit(f"splice: {mol} has {len(old)} published rows, generator made {len(new)}")
        cells = lambda l: [c.strip() for c in l.strip().strip("|").split("|")]
        newcols = len(cells(new[0])) - len(cells(old[0]))
        for a, b in zip(old, new):
            ca, cb = cells(a), cells(b)
            if cb[1].startswith("GNINA pick"):     # label gains one pose tag per added column
                ca, cb = ca[:1] + ca[2:], cb[:1] + cb[2:]
            for j, (x, y) in enumerate(zip(ca, cb[:len(ca)])):
                if x != y and not set(x) <= set("- "):
                    changed.append(f"{mol} row {a[:40]!r} col {j}: published {x!r} -> {y!r}")
        lines[start:end] = new
    if changed:
        raise SystemExit("splice REFUSED, published cells would change:\n  " + "\n  ".join(changed))
    open(path, "w").write("\n".join(lines))
    print(f"spliced Tables 1 and 2 into {path}" + (f" (+{newcols} column)" if newcols else ""))


if "--write" in sys.argv:
    splice()
print(f"structures never GNINA-scored: {len(sorted(set(unscored)))} -> {sorted(set(unscored))}")
print(f"missing: {len(missing)}   pending (legs still running): {len(sorted(set(pending)))}")
for q in sorted(set(pending)): print("   pending", q)
for m in sorted(set(missing)): print("  ", m)
