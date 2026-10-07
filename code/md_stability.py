"""Peptide structural stability over a trajectory: C-alpha RMSD drift and radius of gyration.

MM/GBSA and residence say whether the ligand stays; this says whether the peptide held its fold
while doing it. Two signals, both already computed for every leg:

- protein_ca_rmsd_nm in rmsd.csv, written by `omd analyze`: C-alpha RMSD to the starting frame
  (aligned). A leg that keeps its fold drifts a few tenths of a nm; one that unravels climbs.
- peptide radius of gyration, the same per-decile measure md_contacts.py reports: a compact
  fold that unravels grows its Rg even while contacts persist.

Sampling differs between runs (1 ps local, 10 ps on Modal). RMSD is computed frame-by-frame
against the aligned starting structure, so final/mean values are sampling-insensitive; max is
not -- finer sampling reports more faithfully -- and is quoted as an upper-bound figure.

Pose labels are derived from the leg name: `<peptide>_dock1` is the Vina rank-1 pose, `<peptide>_dockN`
for N>1 is the GNINA pose of that rank, and a bare `<peptide>` leg is the co-folded pose.

Usage:
    ~/miniforge3/envs/openmm-md/bin/python code/md_stability.py runs/octinoxate/md/*/prod_*
"""
import csv
import os
import re

LIG = {"LIG", "UNK", "UNL"}
GNINA = {9, 6, 5}  # the ranks GNINA's CNN has made top pose on the systems simulated here


def label(leg_dir):
    """Peptide name and pose label from the leg's parent directory name."""
    parent = os.path.basename(os.path.dirname(os.path.normpath(leg_dir)))
    m = re.match(r"^(.+)_dock(\d+)$", parent)
    if m:
        pose = f"Vina p{m.group(2)}" if int(m.group(2)) == 1 else f"GNINA p{m.group(2)}"
        return m.group(1), pose
    return parent, "co-folded"


def stats(leg_dir):
    """RMSD and Rg stability figures for one production directory."""
    import mdtraj as md
    import numpy as np

    xtc = os.path.join(leg_dir, "traj_wrapped.xtc")
    pdb = os.path.join(leg_dir, "traj_wrapped.pdb")
    rmsd_csv = os.path.join(leg_dir, "rmsd.csv")
    if not (os.path.exists(xtc) and os.path.exists(pdb)):
        return None

    t = md.load(xtc, top=pdb)
    top = t.topology
    if os.path.exists(rmsd_csv):
        ca = np.genfromtxt(rmsd_csv, delimiter=",", names=True)["protein_ca_rmsd_nm"] * 10.0
    else:
        # Modal legs were never put through `omd analyze`, so rmsd.csv is absent: compute the
        # same measure here -- C-alpha RMSD to the starting frame, aligned (mdtraj rmsd does
        # the least-squares superposition).
        ca_idx = t.topology.select("element C and name CA")
        ca = md.rmsd(t.atom_slice(ca_idx), t.atom_slice(ca_idx), 0) * 10.0

    pep_all = [a.index for a in top.atoms if a.residue.name not in LIG]
    rg = md.compute_rg(t.atom_slice(pep_all)) * 10.0  # nm -> A
    k = max(1, len(rg) // 10)

    return {"leg": leg_dir,
            "ca_rmsd_final_A": round(float(ca[-1]), 1),
            "ca_rmsd_max_A": round(float(ca.max()), 1),
            "ca_rmsd_mean_last_quarter_A": round(float(ca[len(ca) * 3 // 4 :].mean()), 1),
            "rg_first_A": round(float(rg[:k].mean()), 1),
            "rg_last_A": round(float(rg[-k:].mean()), 1)}


def main(argv):
    # Only the standard production directories; skips stray files and extra lengths (prod_25ns)
    legs = sorted({p for p in argv
                   if os.path.basename(os.path.normpath(p)) in ("prod_20ns", "prod_L1_modal")})
    out = []
    for leg_dir in legs:
        r = stats(leg_dir)
        if r is None:
            print(f"skipping {leg_dir}: no traj_wrapped.xtc/pdb")
            continue
        pep, pose = label(leg_dir)
        r["peptide"], r["pose"] = pep, pose
        out.append(r)
        print(f"{pep:<24}{pose:<12}{r['ca_rmsd_final_A']:6.1f}{r['ca_rmsd_max_A']:6.1f}"
              f"{r['ca_rmsd_mean_last_quarter_A']:8.1f}{r['rg_first_A']:7.1f}"
              f" ->{r['rg_last_A']:7.1f}")

    if out:
        # fixed location, derived from this script's own position (code/ -> repo root). One shared
        # file for every ligand -- octinoxate, oxybenzone and the BoltzGen legs all land here.
        repo = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        csv_path = os.path.join(repo, "runs", "octinoxate", "md", "protein_stability.csv")

        # This file is rewritten wholesale, so a run naming fewer legs than it already holds silently
        # deletes the rest -- the trap CLAUDE.md records for the dock tables, and a single-leg run
        # here reduces 62 rows to 1. Refuse to shrink it unless that is asked for explicitly. The
        # glob to pass is every lane: runs/octinoxate/md/*/prod_20ns, .../prod_L1_modal,
        # runs/oxybenzone/md/*/prod_20ns and ~/python_mac/boltzgen_local/md/*/prod_20ns (boltzgen has
        # no prod_L1_modal, and zsh aborts the whole launch if one glob matches nothing).
        if os.path.exists(csv_path) and "--allow-shrink" not in argv:
            with open(csv_path) as fh:
                had = max(0, sum(1 for _ in fh) - 1)
            if len(out) < had:
                print(f"REFUSING to write {len(out)} rows over a file holding {had}: this script "
                      f"rewrites the whole table, so the missing {had - len(out)} would be lost.\n"
                      f"Pass every lane's glob, or --allow-shrink if the shrink is intended.")
                return 1
        with open(csv_path, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=["peptide", "pose", "leg", "ca_rmsd_final_A",
                                               "ca_rmsd_max_A", "ca_rmsd_mean_last_quarter_A",
                                               "rg_first_A", "rg_last_A"])
            w.writeheader()
            w.writerows(out)
        print(f"\nwrote {csv_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(__import__("sys").argv[1:]))