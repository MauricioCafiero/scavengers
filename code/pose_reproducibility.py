"""Pose-reproducibility screen at campaign scale.

Recomputes the 2026-09-30 screen (n=8, rho -0.743 vs residence): Vina re-docks the ligand
into a structure's own fold; `best_rmsd_to_boltz` (dock_summary.csv) is how closely any
docked pose reproduces that structure's predicted pose. Correlated against the residence
the record flow keeps for the structure's BASE leg (same name, no _dock<N> suffix).

Residence is never home-made: parsed from the same md_contacts output (leg tail.log) that
code/record_leg.py parses, then from its record log, then from FOLDS.md's residence tables;
last resort it is recomputed with md_contacts.profile (the same function that produced all
of them). A leg with none of those and no trajectory is skipped and reported.

Writes (do not overwrite anything else):
  runs/<mol>/md/pose_reproducibility.csv         first analysis: structure | base leg residence
  runs/<mol>/md/pose_reproducibility_docked.csv  second analysis: <struct>_dock<K> leg vs pose K

    .venv/bin/python code/pose_reproducibility.py [--csv]
"""
import argparse
import csv
import glob
import os
import re
import sys

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BG = os.path.expanduser("~/python_mac/boltzgen_local/md")
BOLTZGEN_PREFIX = ("bg33_", "bgox31_")
sys.path.insert(0, os.path.join(REPO, "code"))
os.chdir(REPO)

DOCK_SUFFIX = re.compile(r"_dock(\d+)$")
# the exact two lines record_leg.py parses out of md_contacts' output
RES_RE = re.compile(r"residence within 10 A: ([\d.]+)%")
REL_RE = re.compile(r"release .*?: (\d+) frames \(([\d.]+)%\) in (\d+) episodes, longest ([\d.]+) ps")
REC_RES_RE = re.compile(r"residence ([\d.]+)% \| released (\d+) \(([\d.]+)%\) in (\d+)")
LANE_RE = re.compile(r"_(o3cof|rf3cof|o3apo|rf3apo|esmapo|bzapo)$")


def prod_dir(leg_dir):
    for p in (f"{leg_dir}/prod_20ns", f"{leg_dir}/prod_L1_modal"):
        if os.path.exists(p):
            return p
    return None


def leg_dir_of(mol, name):
    for d in (f"runs/{mol}/md/{name}", f"{BG}/{name}"):
        if os.path.isdir(d) and prod_dir(d):
            return d
    return None


def dg_total(f):
    """record_leg.py's own parser of FINAL_RESULTS_MMPBSA.dat (Differences block)."""
    block = open(f).read().split("Differences")[1]
    return float(next(l for l in block.split("\n") if l.startswith("DELTA TOTAL")).split()[2])


def dg_for(mol, struct):
    """Whole-run 20 ns row from that molecule's mmgbsa_summary (record field 9, dg_bind);
    falls back to the leg's own MMPBSA dat (BoltzGen legs never had rows collected)."""
    f = f"runs/{mol}/md/mmgbsa_summary.csv"
    if os.path.exists(f):
        for r in csv.DictReader(open(f)):
            if r["structure"] == struct and not r["leg"].startswith("first_"):
                return float(r["dg_bind"]), float(r["std_err"])
    for d in (f"runs/{mol}/md/{struct}", f"{BG}/{struct}"):
        f = f"{d}/prod_20ns/mmgbsa/FINAL_RESULTS_MMPBSA.dat"
        if os.path.exists(f):
            return dg_total(f), None
    return None, None


def from_md_contacts_text(text):
    m = RES_RE.search(text)
    if not m:
        return None
    out = {"residence_pct": float(m.group(1)), "source": "tail.log"}
    r = REL_RE.search(text)
    if r:
        out.update(released_frames=int(r.group(1)), released_pct=float(r.group(2)),
                   episodes=int(r.group(3)), longest_ps=float(r.group(4)))
    return out


def folds_tables():
    """Every FOLDS.md table whose header carries a residence column -> {leg: dict}."""
    out = {}
    lines = open("FOLDS.md").read().split("\n")
    i = 0
    while i < len(lines):
        if lines[i].startswith("|") and "|---" in (lines[i + 1] if i + 1 < len(lines) else ""):
            cols = [c.strip() for c in lines[i].strip("|").split("|")]
            if "residence" in cols:
                k = 2
                while i + k < len(lines) and lines[i + k].startswith("|"):
                    cells = [c.strip() for c in lines[i + k].strip("|").split("|")]
                    if len(cells) == len(cols):
                        leg = cells[cols.index("leg")].strip("`")
                        v = cells[cols.index("residence")].strip("*`")
                        vm = re.match(r"([\d.]+)%", v)
                        if vm and leg:
                            r = {"residence_pct": float(vm.group(1)), "source": "FOLDS.md"}
                            j = next((c for c in cols if c.startswith("released")), None)
                            if j:
                                rm = re.match(r"(\d+) \(([\d.]+)%\)", cells[cols.index(j)].strip("`*"))
                                if rm:
                                    r.update(released_frames=int(rm.group(1)),
                                             released_pct=float(rm.group(2)))
                            if "episodes" in cols:
                                em = re.match(r"(\d+)", cells[cols.index("episodes")].strip("*`"))
                                if em:
                                    r["episodes"] = int(em.group(1))
                            out[leg] = r
                        k += 1
                i += k
                continue
        i += 1
    return out


def residence_of(mol, leg, leg_dir, folds):
    """Same quantity record_leg records, from the first source that holds it."""
    f = f"{leg_dir}/tail.log"
    if os.path.exists(f):
        r = from_md_contacts_text(open(f, errors="ignore").read())
        if r:
            return r
    f = f"runs/folds_md/record_{leg}.log"
    if os.path.exists(f):
        m = REC_RES_RE.search(open(f, errors="ignore").read())
        if m:
            r = {"residence_pct": float(m.group(1)), "released_frames": int(m.group(2)),
                 "released_pct": float(m.group(3)), "episodes": int(m.group(4)),
                 "source": "record log"}
            return r
    if leg in folds:
        return dict(folds[leg])
    # last resort: recompute with the function that produced every recorded value above
    import md_contacts
    p = prod_dir(leg_dir)
    if p and os.path.exists(f"{p}/traj_wrapped.xtc"):
        try:
            prof = md_contacts.profile(p)
            rel = prof["release"]
            return {"residence_pct": prof["residence_pct"], "source": "recomputed",
                    "released_frames": rel["released_frames"], "released_pct": rel["released_pct"],
                    "episodes": rel["episodes"], "longest_ps": rel["longest_episode_ps"]}
        except Exception as e:
            print(f"  recompute failed for {leg}: {e}")
    return None


def lane_tag(struct):
    lane = re.search(LANE_RE, struct)
    if struct.startswith(BOLTZGEN_PREFIX):
        kind = "boltzgen"
    elif "shuffle" in struct:
        kind = "shuffle-null"
    elif "esm2" in struct:
        kind = "esm2-control"
    else:
        kind = "design"
    return f"{kind}|{lane.group(1) if lane else 'bare-cofold'}"


def ligand_cip(struct):
    """CIP label of the stereocentre the fold itself carries (RESULTS.md: 3 of the original 8
    folds hold the mirror-image ligand). Oxybenzone has none; '' if the SDF is unreadable."""
    from rdkit import Chem
    for d in (f"runs/octinoxate/md/{struct}", f"{BG}/{struct}",
              f"runs/oxybenzone/md/{struct}"):
        f = f"{d}/{struct}_ligand.sdf"
        if not os.path.exists(f):
            continue
        try:
            m = Chem.MolFromMolFile(f, removeHs=False)
            Chem.AssignStereochemistryFrom3D(m)
            Chem.AssignStereochemistry(m, cleanIt=False, force=True)
            c = [a.GetProp("_CIPCode") for a in m.GetAtoms()
                 if a.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED and a.HasProp("_CIPCode")]
            return "+".join(c) if c else "none"
        except Exception:
            return ""
    return ""


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--csv", action="store_true", help="write the two output tables")
    args = ap.parse_args()

    folds = folds_tables()
    per_mol, docked = {}, {}
    skipped = []
    for mol in ("octinoxate", "oxybenzone"):
        summary = list(csv.DictReader(open(f"runs/{mol}/dock/dock_summary.csv")))
        poses = {}
        for r in csv.DictReader(open(f"runs/{mol}/dock/dock_poses.csv")):
            poses.setdefault(r["structure"], {})[int(r["pose"])] = r
        rows, drows = [], []
        for s in summary:
            u = s["structure"]
            base_dir = leg_dir_of(mol, u)
            if not base_dir:
                skipped.append((mol, u, "no base leg (fold never had its own MD)"))
                continue
            rec = residence_of(mol, u, base_dir, folds)
            dg, dgse = dg_for(mol, u)
            if not rec:
                skipped.append((mol, u, "no recorded residence and no trajectory to recompute"))
            else:
                r = {"structure": u, "mol": mol, "fold_lane": lane_tag(u),
                     "ligand_cip": ligand_cip(u) if mol == "octinoxate" else "achiral",
                     "best_rmsd_to_boltz": float(s["best_rmsd_to_boltz"]),
                     "best_rmsd_pose_rank": int(s["best_rmsd_pose_rank"]),
                     "vina_best_score": float(s["vina_best_score"]),
                     "residence_pct": rec["residence_pct"], "residence_source": rec["source"],
                     "released_frames": rec.get("released_frames"),
                     "episodes": rec.get("episodes"), "longest_ps": rec.get("longest_ps"),
                     "dg_bind": dg, "dg_std_err": dgse}
                rows.append(r)
            # second analysis: this structure's docked-pose legs, paired pose K by name
            for leg_dir in sorted(glob.glob(f"runs/{mol}/md/{u}_dock*")):
                if not os.path.isdir(leg_dir):  # stray *_docks_driver.log / *_run.log files
                    continue
                leg = os.path.basename(leg_dir)
                m = DOCK_SUFFIX.search(leg)  # anchored at $: separates _r2 replicates
                if not m:
                    skipped.append((mol, leg, "r2 replicate of a paired pose (excluded)"))
                    continue
                elif prod_dir(leg_dir) is None:
                    skipped.append((mol, leg, "docked leg staged but not run yet"))
                    continue
                K = int(m.group(1))
                p = poses.get(u, {}).get(K)
                if p is None:
                    skipped.append((mol, leg, f"pose {K} missing from dock_poses.csv"))
                    continue
                drec = residence_of(mol, leg, leg_dir, folds)
                if not drec:
                    skipped.append((mol, leg, "docked leg has no residence"))
                    continue
                ddg, _ = dg_for(mol, leg)
                drows.append({"leg": leg, "structure": u, "mol": mol, "pose": K,
                              "ligand_cip": ligand_cip(u) if mol == "octinoxate" else "achiral",
                              "rmsd_to_boltz": float(p["rmsd_to_boltz"]),
                              "vina_score": p["vina_score"],
                              "residence_pct": drec["residence_pct"],
                              "residence_source": drec["source"],
                              "released_frames": drec.get("released_frames"),
                              "episodes": drec.get("episodes"),
                              "dg_bind": ddg})
        per_mol[mol], docked[mol] = rows, drows

    if args.csv:
        for mol, rows in per_mol.items():
            out = f"runs/{mol}/md/pose_reproducibility.csv"
            cols = list(rows[0])
            with open(out, "w", newline="") as f:
                w = csv.DictWriter(f, fieldnames=cols)
                w.writeheader()
                w.writerows(rows)
            print(f"wrote {out} ({len(rows)} rows)")
        for mol, rows in docked.items():
            if not rows:
                continue
            out = f"runs/{mol}/md/pose_reproducibility_docked.csv"
            with open(out, "w", newline="") as f:
                w = csv.DictWriter(f, fieldnames=list(rows[0]))
                w.writeheader()
                w.writerows(rows)
            print(f"wrote {out} ({len(rows)} rows)")

    report(per_mol, docked, skipped)


def corr(rows, x, y):
    import numpy as np
    from scipy.stats import pearsonr, spearmanr
    pts = [(float(r[x]), float(r[y])) for r in rows if r.get(y) is not None]
    if len(pts) < 3:
        return None
    xs, ys = (np.array(a) for a in zip(*pts))
    sr, sp = spearmanr(xs, ys)
    pr, pp = pearsonr(xs, ys)
    return len(pts), sr, sp, pr, pp


def fmt(c):
    return f"n={c[0]:<3} rho={c[1]:+.3f} (p={c[2]:.3g})   r={c[3]:+.3f} (p={c[4]:.3g})" if c else "n/a"


def report(per_mol, docked, skipped):
    print("\n=== FIRST ANALYSIS: structure best_rmsd_to_boltz vs its own base leg (docked poses excluded)")
    pool = []
    for mol, rows in per_mol.items():
        pool += rows
        print(f"\n{mol}  (n={len(rows)})")
        for x, lab in (("best_rmsd_to_boltz", "repro-RMSD"), ("vina_best_score", "vina best")):
            for y, ylab in (("residence_pct", "residence"), ("episodes", "episodes"),
                            ("dg_bind", "dG")):
                print(f"   {lab:>12} vs {ylab:<10} {fmt(corr(rows, x, y))}")
        if mol == "octinoxate":
            for c in ("S", "R"):
                sub = [r for r in rows if "ligand_cip" in r and r["ligand_cip"].startswith(c + "+") or r["ligand_cip"] == c]
                if len(sub) >= 4:
                    print(f"   {c}-fold only     " + fmt(corr(sub, 'best_rmsd_to_boltz', 'residence_pct')))
    print(f"\npooled  (n={len(pool)})")
    for y in ("residence_pct", "episodes", "dg_bind"):
        print(f"   repro-RMSD  vs {y:<10} {fmt(corr(pool, 'best_rmsd_to_boltz', y))}")

    print("\n=== SECOND ANALYSIS: docked-pose legs, pose-K RMSD vs THAT DOCKED LEG's residence")
    dpooled = [r for rows in docked.values() for r in rows]
    for rows, label in [(docked.get("octinoxate"), "octinoxate"), (docked.get("oxybenzone"), "oxybenzone"),
                        (dpooled, "pooled")]:
        rows = [r for r in rows if r["vina_score"]]
        if not rows:
            continue
        print(f"\n{label}  (n={len(rows)})")
        for y in ("residence_pct", "episodes", "dg_bind"):
            print(f"   pose RMSD   vs {y:<10} {fmt(corr(rows, 'rmsd_to_boltz', y))}")
            print(f"   vina score  vs {y:<10} {fmt(corr(rows, 'vina_score', y))}")

    print("\n=== SKIPPED")
    counts = {}
    for mol, leg, why in skipped:
        counts.setdefault(why.split(" (")[0], []).append(f"{mol}/{leg}")
    for why, legs in sorted(counts.items(), key=lambda kv: -len(kv[1])):
        print(f"  {len(legs)}  {why}: {', '.join(legs) if len(legs) <= 12 else legs[:12] + ['...']}")


if __name__ == "__main__":
    main()