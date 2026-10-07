#!/usr/bin/env python3
"""How far the ligand moved between the built pose and the first production frame.

Equilibration restrains protein heavy atoms only (`openmm_md/dynamics.py:144`), so the ligand is free
the whole time and arrives at production somewhere other than where the build put it -- between 1.2
and 15.1 A away across this project's legs. CLAUDE.md requires this number on every leg, because it
is what says whether a cell is attributable to the starting pose it is labelled with. A leg whose
ligand slid 10 A is not reporting the docked pose's affinity; it is reporting wherever the ligand
drifted to.

Do NOT "fix" a large slide by restraining the ligand through equilibration: every run in this project
equilibrated it freely, and a restrained leg would not be comparable to any of them.

The two structures have different atom sets -- `system/complex.pdb` is the full solvated box, while
`prod_*/traj_wrapped.pdb` is the stripped solute -- so they are matched on the solute only: protein
C-alphas are superposed (the fold itself drifts, and that drift is md_stability.py's business, not
this script's), then the ligand heavy-atom RMSD is read off. Superposing first is what makes the
number mean "the ligand moved relative to its peptide" rather than "the whole complex translated in
the box".

    python code/ligand_slide.py runs/oxybenzone/md/*/prod_20ns ~/python_mac/boltzgen_local/md/*/prod_20ns
    python code/ligand_slide.py --csv runs/oxybenzone/md/ligand_slide.csv <legs...>
"""
import argparse
import csv
import os

LIG = {"LIG", "UNK", "UNL"}


def slide(prod_dir):
    """Ligand heavy-atom RMSD, built pose -> first production frame, in angstrom."""
    import mdtraj as md
    import numpy as np

    leg = os.path.dirname(os.path.normpath(prod_dir))
    built = os.path.join(leg, "system", "complex.pdb")
    top = os.path.join(prod_dir, "traj_wrapped.pdb")
    xtc = os.path.join(prod_dir, "traj_wrapped.xtc")
    for f in (built, top, xtc):
        if not os.path.exists(f):
            return None, f"missing {os.path.relpath(f, leg)}"

    a = md.load(built)
    b = md.load(xtc, top=top)[0]

    def sel(t):
        lig = t.topology.select(" or ".join(f"resname {r}" for r in LIG))
        lig = [i for i in lig if t.topology.atom(i).element.symbol != "H"]
        ca = t.topology.select("name CA and protein")
        return np.asarray(lig), np.asarray(ca)

    lig_a, ca_a = sel(a)
    lig_b, ca_b = sel(b)
    if len(lig_a) == 0 or len(lig_b) == 0:
        return None, "no ligand found (resname not in LIG/UNK/UNL)"
    if len(lig_a) != len(lig_b) or len(ca_a) != len(ca_b):
        return None, (f"atom-count mismatch: ligand {len(lig_a)} vs {len(lig_b)}, "
                      f"CA {len(ca_a)} vs {len(ca_b)}")

    # Kabsch on the C-alphas, applied to everything, then ligand RMSD. nm throughout in mdtraj.
    P = a.xyz[0][ca_a]
    Q = b.xyz[0][ca_b]
    Pc, Qc = P - P.mean(0), Q - Q.mean(0)
    V, _, W = np.linalg.svd(Pc.T @ Qc)
    d = np.sign(np.linalg.det(V @ W))
    R = V @ np.diag([1.0, 1.0, d]) @ W
    moved = (b.xyz[0] - Q.mean(0)) @ R.T + P.mean(0)
    diff = a.xyz[0][lig_a] - moved[lig_b]
    return float(np.sqrt((diff ** 2).sum(axis=1).mean()) * 10.0), None


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("legs", nargs="+", help="production directories with traj_wrapped.{xtc,pdb}")
    ap.add_argument("--csv", help="write the table here as well as printing it")
    a = ap.parse_args(argv)

    rows = []
    for prod in sorted(a.legs):
        leg = os.path.basename(os.path.dirname(os.path.normpath(prod)))
        rmsd, why = slide(prod)
        if rmsd is None:
            print(f"{leg:24s}   --      ({why})")
            continue
        flag = "" if rmsd < 2.5 else ("  <- slid" if rmsd < 6 else "  <- SLID FAR")
        print(f"{leg:24s} {rmsd:6.2f} A{flag}")
        rows.append({"leg": leg, "ligand_slide_A": round(rmsd, 2)})

    if a.csv and rows:
        with open(a.csv, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=["leg", "ligand_slide_A"])
            w.writeheader()
            w.writerows(rows)
        print(f"\nwrote {a.csv}  ({len(rows)} legs)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
