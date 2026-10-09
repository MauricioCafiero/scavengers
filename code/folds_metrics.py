"""ESMFold/OpenFold3-vs-co-fold metrics for the fold comparison in FOLDS.md.

Reads the per-molecule manifests (`runs/<mol>/folds/manifest.csv`), each peptide's fold and the
co-fold it is compared against (the `protein_fixed.pdb` the dynamics leg was built from), and
writes one row per peptide to `runs/<mol>/folds/metrics/<lane>.csv`:

    helical/beta/hbonds   check_fold.secondary_structure, same dict format, same cutoffs
    rg                    peptide CA radius of gyration, both folds
    ca_rmsd               fold CAs superposed on the co-fold CAs (kabsch; same sequence)
    mean_plddt            fold's B-factor column, the model's per-residue confidence

Lanes:
    esmfold     `runs/<mol>/folds/esmfold/<peptide>.pdb`   — ligand-free; wrapped/enclosed/engaged
                are N/A because there is no ligand in the structure.
    openfold3   `runs/<mol>/folds/openfold3/<peptide>.pdb` — cofolded WITH the ligand (converted
                from OF3's best-ranked cif), so `check_fold.geometry()` runs unmodified and the
                wrapped/enclosed/engaged columns read straight against fold_check.csv and
                dock_poses.csv.

Usage:  python code/folds_metrics.py [esmfold|openfold3]   # default esmfold
"""
import csv
import glob
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__))))
from check_fold import geometry, secondary_structure  # noqa: E402

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BG_MD = os.path.expanduser("~/python_mac/boltzgen_local/md")


def pdb_atoms(path):
    """Heavy-atom list of dicts in check_fold's format, plus CA coords and B-factors."""
    atoms = []
    for line in open(path):
        if not line.startswith(("ATOM", "HETATM")):
            continue
        name = line[12:16].strip()
        if name.startswith("H") or line[77:79].strip() == "H":
            continue
        atoms.append({"name": name, "resseq": int(line[22:26]), "chain": line[21],
                      "resname": line[17:20].strip(),
                      "xyz": [float(line[30:38]), float(line[38:46]), float(line[46:54])],
                      "b": float(line[60:66])})
    return atoms


def split_pep_lig(atoms):
    """(peptide atoms, ligand atoms) by chain: chains holding many residues are peptide, the rest
    (OF3's ligand on chain Z) are ligand. Falls back to element HETATM if only one chain."""
    chains = {}
    for a in atoms:
        chains.setdefault(a["chain"], []).append(a)
    if len(chains) == 1:
        return atoms, []
    pep_key = max(chains, key=lambda c: len({a["resseq"] for a in chains[c]}))
    pep = chains[pep_key]
    lig = [a for c, aa in chains.items() if c != pep_key for a in aa]
    return pep, lig


def cofold_reference_atoms(legdir):
    """Peptide and LIG heavy atoms from the leg's own system/complex.pdb -- the structure the
    dynamics leg was actually built from, read through the same geometry() the fold_check tables
    used. Peptide = ATOM records (chain A + hydrogens excluded), ligand = HETATM resname LIG
    (HOH/NA/CL ride along as other HETATM chains and are dropped)."""
    path = os.path.join(legdir, "system", "complex.pdb")
    if not os.path.exists(path):
        return None, None
    pep, lig = [], []
    for a in pdb_atoms(path):
        if a["chain"] == "A":
            pep.append(a)
        elif a["resname"] == "LIG":
            lig.append(a)
    return pep, lig


def geom_of(pep_at, lig_at, which=""):
    ph = np.array([a["xyz"] for a in pep_at])
    lh = np.array([a["xyz"] for a in lig_at])
    g = geometry(ph, lh)
    base = {"enclosed": g.get("enclosed_fraction"),
            "wrapped": g.get("wrapped_fraction"),
            "engaged": g.get("engaged"), "lig_heavy": g.get("ligand_heavy_atoms"),
            "cent_sep": g.get("centroid_separation")}
    return {f"{k}{which}": v for k, v in base.items()}


def bfac_to_plddt(atoms):
    vals = [a["b"] for a in atoms if a["name"] == "CA"]
    if not vals:
        return ""
    # ESMFold writes pLDDT 0-1, OF3/RF3 write 0-100; normalise to a 0-100 scale.
    m = float(np.mean(vals))
    return round(m * 100 if m <= 1.0 else m, 1)


def ca_of(atoms):
    idx = [i for i, a in enumerate(atoms) if a["name"] == "CA"]
    return np.array([atoms[i]["xyz"] for i in idx]), idx


def rg(ca):
    c = ca - ca.mean(axis=0)
    return round(float(np.sqrt((c ** 2).sum(axis=1).mean())), 2)


def kabsch_rmsd(P, Q):
    """RMSD after optimal superposition of P onto Q (both (n,3), same order)."""
    Pc, Qc = P - P.mean(0), Q - Q.mean(0)
    U, S, Vt = np.linalg.svd(Pc.T @ Qc)
    d = np.sign(np.linalg.det(U @ Vt))
    R = U @ np.diag([1, 1, d]) @ Vt
    diff = (Pc @ R) - Qc
    return round(float(np.sqrt((diff ** 2).sum(1).mean())), 2)


def seq_of(atoms):
    ca_order = sorted({a["resseq"] for a in atoms if a["name"] == "CA"})
    return ca_order


def main():
    lane = sys.argv[1] if len(sys.argv) > 1 else "esmfold"
    total = 0
    for mol in ("octinoxate", "oxybenzone"):
        man = f"{REPO}/runs/{mol}/folds/manifest.csv"
        if not os.path.exists(man):
            continue
        out_rows = []
        for row in csv.DictReader(open(man)):
            pep = row["peptide"]
            fold_path = f"{REPO}/runs/{mol}/folds/{lane}/{pep}.pdb"
            if not os.path.exists(fold_path):
                continue  # folds as they land; re-run the script to refill
            src = (f"{REPO}/runs/{mol}/md/{pep}/protein_fixed.pdb"
                   if os.path.isdir(f"{REPO}/runs/{mol}/md/{pep}")
                   else f"{BG_MD}/{pep}/protein_fixed.pdb")
            fold_atoms, cof = pdb_atoms(fold_path), pdb_atoms(src)
            ss_e, ss_c = secondary_structure(fold_atoms), secondary_structure(cof)
            fCA, _ = ca_of(fold_atoms)
            cCA, _ = ca_of(cof)
            rmsd = kabsch_rmsd(fCA, cCA) if len(fCA) == len(cCA) else ""
            base = {"peptide": pep, "kind": row["kind"],
                    "helix_esm": ss_e["helical_fraction"], "helix_cofold": ss_c["helical_fraction"],
                    "beta_esm": ss_e["beta_fraction"], "beta_cofold": ss_c["beta_fraction"],
                    "nlocal_esm": ss_e["helical_hbonds"], "nnonlocal_esm": ss_e["nonlocal_hbonds"],
                    "rg_esm": rg(fCA), "rg_cofold": rg(cCA),
                    "ca_rmsd": rmsd, "mean_plddt": bfac_to_plddt(fold_atoms)}
            if lane != "esmfold":
                # cofold lanes carry the ligand: check_fold.geometry() unmodified, same cutoffs
                # as fold_check.csv and dock_poses.csv read -- on the fold, and on the leg's own
                # complex.pdb (the Boltz/BoltzGen co-fold) as the reference column set.
                pep_at, lig_at = split_pep_lig(fold_atoms)
                if lig_at:
                    base.update(geom_of(pep_at, lig_at))
                rpep, rlig = cofold_reference_atoms(os.path.dirname(src))
                if rpep and rlig:
                    base.update(geom_of(rpep, rlig, "ref"))
            out_rows.append(base)
            total += 1
        if out_rows:
            path = f"{REPO}/runs/{mol}/folds/metrics/{lane}.csv"
            os.makedirs(os.path.dirname(path), exist_ok=True)
            with open(path, "w") as f:
                w = csv.DictWriter(f, list(out_rows[0]))
                w.writeheader(); w.writerows(out_rows)
            for r in out_rows:
                line = (f"{mol[:3]} {r['peptide']:22s} helix {r['helix_esm']:>5} vs {r['helix_cofold']:>5}"
                        f" | Rg {r['rg_esm']:>5} vs {r['rg_cofold']:>5}"
                        f" | RMSD {r['ca_rmsd']:>5} | plddt {r['mean_plddt']}")
                if lane != "esmfold":
                    line += (f" | enc {r['enclosed']}/{r['enclosedref']} wrap {r['wrapped']}/{r['wrappedref']}"
                             f" eng {r['engaged']}/{r['engagedref']} d {r['cent_sep']}/{r['cent_sepref']}")
                print(line)
    print(f"{total} rows written -> runs/*/folds/metrics/{lane}.csv")


if __name__ == "__main__":
    main()