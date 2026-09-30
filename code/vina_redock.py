#!/usr/bin/env python3
"""Redock the ligand into its own Boltz-folded peptide with AutoDock Vina.

Every peptide-octinoxate complex in this project was produced either by the
shell search (fragments placed around a fixed ligand) or by Boltz co-folding
(peptide and ligand placed jointly). Neither is a docking calculation, so the
ligand's position has never been checked by an independent pose search. This
script does that check: it takes the receptor and reference ligand that the MD
runs already prepared, and asks Vina to re-find the ligand's place in a rigid
copy of the same peptide.

What the numbers mean, and what they do not: the peptide fold was predicted
*with* the ligand inside it, so the pocket was carved around the ligand. A low
RMSD therefore says Vina agrees with Boltz about where the ligand sits in that
fold -- self-consistency of the placement. It does not say the peptide would
bind the ligand de novo, and the receptor is rigid here, so no induced fit is
explored.

Inputs come from runs/<system>/md/<structure>/ as prepared for the dynamics:

    <structure>_protein.pdb    receptor, heavy atoms, ligand already stripped
    <structure>_ligand.sdf     reference ligand in the Boltz pose

The ligand is converted straight from that SDF rather than re-embedded from
SMILES. The built-in ligand is an analogue of octinoxate (C17H24O3, one CH2
short of the real thing), and a SMILES round-trip risks both losing its
stereocentre and quietly docking a different molecule than MD scored.

Requires dock_assist's environment, which carries Open Babel and the Vina
binary vendored in dockstring:

    ~/python_mac/dock_assist/dock-env/bin/python code/vina_redock.py orig_f12
"""

import argparse
import csv
import os
import re
import shutil
import sys
from types import SimpleNamespace

import numpy as np
from rdkit import Chem

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
DOCK_ASSIST = os.path.expanduser("~/python_mac/dock_assist")

# dock_assist owns the docking machinery; import it rather than reimplement it.
sys.path.insert(0, os.path.join(DOCK_ASSIST, "code"))
# obabel lives in dock-env/bin, which is not on PATH just because we run its
# python; vina_dock.require() looks the binary up with shutil.which.
os.environ["PATH"] = os.path.join(DOCK_ASSIST, "dock-env", "bin") + os.pathsep + os.environ["PATH"]

from vina_dock import (DockError, build_receptor_pdbqt, convert, find_vina_bin,
                       parse_vina_log, run_vina)

# The project's own wrapping/enclosure definitions, so a docked pose is measured by exactly the
# same cutoffs as a Boltz fold in fold_check.csv rather than by a second set invented here.
from check_fold import geometry  # noqa: E402

# Cubic box edge. The peptides are small enough that 22 A around the ligand may
# cover most of the chain; frac_receptor_in_box in the CSV records how much, so
# a "site" dock that was really a whole-peptide search is visible as such.
BOX_SIZE = 22.0
EXHAUSTIVENESS = 16
NUM_MODES = 9
SEED = 42


def heavy_atom_coords(pdb_path):
    """Coordinates of ATOM/HETATM records, skipping hydrogens."""
    coords = []
    with open(pdb_path) as fh:
        for line in fh:
            if not line.startswith(("ATOM", "HETATM")):
                continue
            if line[76:78].strip() == "H":
                continue
            coords.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
    return np.array(coords)


def pdbqt_models(path):
    """Heavy-atom coordinates of each MODEL in a PDBQT, in file order.

    PDBQT from Open Babel is united-atom: non-polar hydrogens are already gone,
    so every ATOM record is a heavy atom.
    """
    models, cur = [], []
    for line in open(path):
        if line.startswith("MODEL"):
            cur = []
        elif line.startswith("ENDMDL"):
            models.append(np.array(cur))
            cur = []
        elif line.startswith(("ATOM", "HETATM")):
            cur.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
    if cur:
        models.append(np.array(cur))
    return models


def ref_order_poses(poses_pdbqt, lig_pdbqt, ref_coords):
    """Docked pose coordinates, reindexed into the reference SDF's atom order.

    Vina preserves the atom order of the input ligand PDBQT in its output poses, and that PDBQT was
    built from the reference SDF without touching coordinates -- so the PDBQT-to-reference mapping
    is recoverable exactly, by coordinate identity, with no chemical perception step. This matters
    because Open Babel's pdbqt->sdf round trip loses valence information, and a pose read back that
    way will not substructure-match the reference at all.
    """
    lig = pdbqt_models(lig_pdbqt)[0]
    if lig.shape != ref_coords.shape:
        raise DockError(f"ligand PDBQT has {lig.shape[0]} heavy atoms, "
                        f"reference SDF has {ref_coords.shape[0]}")
    dist = np.linalg.norm(lig[:, None, :] - ref_coords[None, :, :], axis=2)
    perm = dist.argmin(axis=1)
    if dist.min(axis=1).max() > 0.01 or len(set(perm.tolist())) != len(perm):
        raise DockError("could not map ligand PDBQT atoms onto the reference SDF "
                        "by coordinate identity; did the SDF get re-embedded?")
    out = []
    for pose in pdbqt_models(poses_pdbqt):
        ordered = np.empty_like(ref_coords)
        ordered[perm] = pose
        out.append(ordered)
    return out


def pose_rmsds(poses_pdbqt, lig_pdbqt, ref_mol):
    """Symmetry-corrected in-place RMSD of each docked pose against the reference.

    Deliberately avoids RDKit's CalcRMS/GetBestRMS on an SDF converted back from
    PDBQT. Open Babel's pdbqt->sdf loses valence information (carbons come back
    as bracketed [C] with no implicit hydrogens), so the docked molecule will not
    substructure-match the reference and every RMSD comes back as a failure.

    Instead: Vina preserves the atom order of the input ligand PDBQT in its
    output poses, and that PDBQT was built from the reference SDF without
    touching coordinates. So the PDBQT->reference atom mapping is recoverable
    exactly, by coordinate identity, and pose coordinates can be read straight
    out of the PDBQT with no perception step at all.

    RMSD is in place, never superimposed: the poses and the reference share a
    coordinate frame, and aligning them would discard the displacement being
    measured. The minimum is taken over the reference's automorphisms so that a
    symmetry-equivalent orientation (here the para-substituted ring flip) is not
    counted as a difference.
    """
    ref = ref_mol.GetConformer().GetPositions()
    poses = ref_order_poses(poses_pdbqt, lig_pdbqt, ref)
    autos = ref_mol.GetSubstructMatches(ref_mol, uniquify=False, useChirality=False)
    if not autos:
        autos = [tuple(range(ref_mol.GetNumAtoms()))]

    rmsds, centroid_shifts = [], []
    ref_center = ref.mean(axis=0)
    for ordered in poses:
        rmsds.append(min(float(np.sqrt(((ordered[list(a)] - ref) ** 2).sum(axis=1).mean()))
                         for a in autos))
        centroid_shifts.append(float(np.linalg.norm(ordered.mean(axis=0) - ref_center)))
    return rmsds, centroid_shifts


def redock(structure, md_dir, out_root, verbose=False):
    """Dock one structure's ligand back into its own receptor.

    Returns (summary row, per-pose rows). The per-pose rows include pose 0, which is the Boltz
    reference pose measured with the same receptor -- the baseline every docked pose is compared to.
    """
    src_pdb = os.path.join(md_dir, f"{structure}_protein.pdb")
    src_sdf = os.path.join(md_dir, f"{structure}_ligand.sdf")
    for p in (src_pdb, src_sdf):
        if not os.path.exists(p):
            raise DockError(f"missing input: {p}")

    # Everything for this structure lands in its own folder. The receptor and
    # reference ligand are copied in, so the PDBQT and pose files Vina and
    # obabel write beside them never touch the MD directories.
    work = os.path.join(out_root, structure)
    os.makedirs(work, exist_ok=True)
    rec_pdb = os.path.join(work, f"{structure}_protein.pdb")
    ref_sdf = os.path.join(work, "ligand_ref.sdf")
    shutil.copy2(src_pdb, rec_pdb)
    shutil.copy2(src_sdf, ref_sdf)

    rec_pdbqt = os.path.join(work, f"{structure}_protein.pdbqt")
    lig_pdbqt = os.path.join(work, "ligand.pdbqt")
    poses_pdbqt = os.path.join(work, "poses.pdbqt")
    poses_sdf = os.path.join(work, "poses.sdf")
    log_path = os.path.join(work, "vina.log")

    build_receptor_pdbqt(rec_pdb, rec_pdbqt)
    convert(["obabel", "-isdf", ref_sdf, "-opdbqt", "-O", lig_pdbqt],
            f"reference ligand SDF->PDBQT ({structure})")

    ref_mol = Chem.MolFromMolFile(ref_sdf, removeHs=True)
    if ref_mol is None:
        raise DockError(f"could not read reference ligand: {ref_sdf}")
    ref_coords = ref_mol.GetConformer().GetPositions()
    center = ref_coords.mean(axis=0)

    rec_coords = heavy_atom_coords(rec_pdb)
    half = BOX_SIZE / 2.0
    inside = np.all(np.abs(rec_coords - center) <= half, axis=1)
    frac_in_box = float(inside.mean())

    args = SimpleNamespace(exhaustiveness=EXHAUSTIVENESS, num_modes=NUM_MODES,
                           cpu=0, seed=SEED)
    run_vina(find_vina_bin(None), rec_pdbqt, lig_pdbqt, center,
             [BOX_SIZE] * 3, args, poses_pdbqt, log_path)
    best_score, n_modes = parse_vina_log(log_path)
    if best_score is None:
        raise DockError(f"no score parsed from {log_path}")

    convert(["obabel", "-ipdbqt", poses_pdbqt, "-osdf", "-O", poses_sdf],
            f"poses PDBQT->SDF ({structure})")

    rmsds, centroid_shifts = pose_rmsds(poses_pdbqt, lig_pdbqt, ref_mol)
    scores = []
    row_re = re.compile(r"^\s*\d+\s+(-?\d+\.\d+)\s+")
    with open(log_path) as fh:
        for line in fh:
            m = row_re.match(line)
            if m:
                scores.append(float(m.group(1)))

    # Pose 1 is Vina's best-scoring pose; the best-RMSD pose may be a different
    # one, and the gap between them is the interesting part -- it says whether
    # Vina found the Boltz pose at all but ranked something else above it.
    best_rmsd_idx, best_rmsd = min(enumerate(rmsds), key=lambda t: t[1])
    pose1_rmsd = rmsds[0]
    pose1_centroid = centroid_shifts[0]

    # Wrapping and enclosure, on the same footing as the Boltz folds. The reference row is the
    # baseline: it is the geometry check_fold.py already reports for this structure, recomputed here
    # from the MD receptor so that pose and reference differ only in where the ligand sits. A Vina
    # score cannot say whether a pose keeps the ligand wrapped, and that is the design objective --
    # a pose can score better while sitting in a surface groove with half the ligand in solvent.
    ref_geom = geometry(rec_coords, ref_coords)
    pose_geoms = [geometry(rec_coords, pose) for pose in pdbqt_models(poses_pdbqt)]
    pose_rows = []
    for i, g in enumerate(pose_geoms):
        pose_rows.append({
            "structure": structure,
            "pose": i + 1,
            "vina_score": scores[i] if i < len(scores) else None,
            "rmsd_to_boltz": round(rmsds[i], 3),
            "ligand_centroid_shift": round(centroid_shifts[i], 3),
            "wrapped_fraction": g["wrapped_fraction"],
            "enclosed_fraction": g["enclosed_fraction"],
            "engaged": g["engaged"],
            "centroid_separation": g["centroid_separation"],
            "mean_ligand_distance": g["mean_ligand_distance"],
        })
    pose_rows.append({
        "structure": structure,
        "pose": 0,  # 0 = the Boltz reference pose, for comparison
        "vina_score": None,
        "rmsd_to_boltz": 0.0,
        "ligand_centroid_shift": 0.0,
        "wrapped_fraction": ref_geom["wrapped_fraction"],
        "enclosed_fraction": ref_geom["enclosed_fraction"],
        "engaged": ref_geom["engaged"],
        "centroid_separation": ref_geom["centroid_separation"],
        "mean_ligand_distance": ref_geom["mean_ligand_distance"],
    })
    g1 = pose_geoms[0]
    gbest = pose_geoms[best_rmsd_idx]

    return {
        "structure": structure,
        "source": os.path.relpath(md_dir, REPO) if md_dir.startswith(REPO) else md_dir,
        "n_receptor_heavy_atoms": len(rec_coords),
        "n_ligand_heavy_atoms": ref_mol.GetNumAtoms(),
        "box_center_x": round(float(center[0]), 3),
        "box_center_y": round(float(center[1]), 3),
        "box_center_z": round(float(center[2]), 3),
        "box_size": BOX_SIZE,
        "frac_receptor_in_box": round(frac_in_box, 3),
        "n_modes": n_modes,
        "vina_best_score": best_score,
        "vina_worst_score": scores[-1] if scores else None,
        "pose1_rmsd_to_boltz": round(pose1_rmsd, 3),
        "pose1_centroid_shift": round(pose1_centroid, 3),
        "best_rmsd_to_boltz": round(best_rmsd, 3),
        "best_rmsd_pose_rank": best_rmsd_idx + 1,
        "best_rmsd_pose_score": (scores[best_rmsd_idx] if best_rmsd_idx < len(scores) else None),
        "peptide_rg": ref_geom["peptide_rg"],
        "boltz_wrapped_fraction": ref_geom["wrapped_fraction"],
        "boltz_enclosed_fraction": ref_geom["enclosed_fraction"],
        "boltz_engaged": ref_geom["engaged"],
        "boltz_centroid_separation": ref_geom["centroid_separation"],
        "boltz_mean_ligand_distance": ref_geom["mean_ligand_distance"],
        "pose1_wrapped_fraction": g1["wrapped_fraction"],
        "pose1_enclosed_fraction": g1["enclosed_fraction"],
        "pose1_engaged": g1["engaged"],
        "pose1_centroid_separation": g1["centroid_separation"],
        "pose1_mean_ligand_distance": g1["mean_ligand_distance"],
        "bestrmsd_wrapped_fraction": gbest["wrapped_fraction"],
        "bestrmsd_enclosed_fraction": gbest["enclosed_fraction"],
        "poses_sdf": os.path.relpath(poses_sdf, REPO),
    }, pose_rows


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("structures", nargs="*",
                    help="MD run names under runs/<system>/md/ (e.g. orig_f12)")
    ap.add_argument("--system", default="octinoxate", help="ligand system (default: octinoxate)")
    ap.add_argument("--md-root", action="append", default=None, metavar="DIR",
                    help="extra directory of prepared MD runs to take structures from, repeatable. "
                         "Use it to dock structures prepared outside this repository -- BoltzGen's "
                         "live in ~/python_mac/boltzgen_local/md -- into the same tables, so they "
                         "compare directly. runs/<system>/md is always searched first")
    ap.add_argument("--all", action="store_true",
                    help="dock every md/ subdirectory that has the protein+ligand pair")
    ap.add_argument("--out", default=None, help="output root (default: runs/<system>/dock)")
    ap.add_argument("--verbose", action="store_true", help="show Vina and obabel output")
    a = ap.parse_args()

    roots = [os.path.join(REPO, "runs", a.system, "md")]
    roots += [os.path.expanduser(r) for r in (a.md_root or [])]
    out_root = a.out or os.path.join(REPO, "runs", a.system, "dock")
    os.makedirs(out_root, exist_ok=True)

    def prepared(root, d):
        sub = os.path.join(root, d)
        return (os.path.isdir(sub)
                and os.path.exists(os.path.join(sub, f"{d}_protein.pdb"))
                and os.path.exists(os.path.join(sub, f"{d}_ligand.sdf")))

    # Where each structure's inputs live. First root wins, so a name present in both is taken from
    # this repository.
    where = {}
    for root in roots:
        if not os.path.isdir(root):
            ap.error(f"not a directory: {root}")
        for d in sorted(os.listdir(root)):
            if prepared(root, d) and d not in where:
                where[d] = root

    names = list(a.structures)
    if a.all:
        names = sorted(set(names) | set(where))
    if not names:
        ap.error("name at least one structure, or pass --all")
    missing = [n for n in names if n not in where]
    if missing:
        ap.error(f"no prepared protein+ligand pair found for: {', '.join(missing)}\n"
                 f"searched: {', '.join(roots)}")

    rows, all_poses, failed = [], [], []
    for name in names:
        print(f"\n=== {name} ===", flush=True)
        try:
            row, pose_rows = redock(name, os.path.join(where[name], name), out_root, a.verbose)
        except Exception as e:
            print(f"  FAILED: {type(e).__name__}: {e}", flush=True)
            failed.append((name, f"{type(e).__name__}: {e}"))
            continue
        rows.append(row)
        all_poses.extend(pose_rows)
        print(f"  score {row['vina_best_score']:.2f} kcal/mol   "
              f"pose1 RMSD {row['pose1_rmsd_to_boltz']} A   "
              f"best RMSD {row['best_rmsd_to_boltz']} A (rank {row['best_rmsd_pose_rank']})",
              flush=True)
        print(f"  wrapped  Boltz {row['boltz_wrapped_fraction']:.2f} -> pose1 "
              f"{row['pose1_wrapped_fraction']:.2f}    enclosed  "
              f"{row['boltz_enclosed_fraction']:.2f} -> {row['pose1_enclosed_fraction']:.2f}"
              f"    centroid sep {row['boltz_centroid_separation']:.1f} -> "
              f"{row['pose1_centroid_separation']:.1f} A", flush=True)

    if rows:
        csv_path = os.path.join(out_root, "dock_summary.csv")
        with open(csv_path, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
            w.writeheader()
            w.writerows(rows)
        print(f"\nwrote {csv_path}  ({len(rows)} structures)")
    if all_poses:
        pose_csv = os.path.join(out_root, "dock_poses.csv")
        with open(pose_csv, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(all_poses[0].keys()))
            w.writeheader()
            w.writerows(sorted(all_poses, key=lambda r: (r["structure"], r["pose"])))
        print(f"wrote {pose_csv}  ({len(all_poses)} rows, pose 0 = Boltz reference)")
    if failed:
        print(f"\n{len(failed)} failed:")
        for n, why in failed:
            print(f"  {n}: {why}")
    return 1 if failed and not rows else 0


if __name__ == "__main__":
    sys.exit(main())
