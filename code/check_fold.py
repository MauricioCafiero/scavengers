"""Did the fold actually bind the ligand, and did it honour the contact hints?

Two questions that have to be asked of every fold before any energy is quoted, because the answer
is not always the obvious one:

1. **Is the ligand in contact at all?** `NEXT_STEPS.md` records Boltz giving `RYGLGEFSDWKI` a
   confident pIC50 of 6.96 with the ligand 6.6 A away. A score without a contact distance beside it
   is not evidence of binding.
2. **Were the contact hints satisfied?** A hinted fold that ignored its hints has not tested the
   designed arrangement -- it has just folded the sequence. Soft hints (`force: false`) are a
   suggestion, and Boltz is free to decline.

The hints are read back out of the yaml that was used, so this checks what was actually asked for
rather than what was intended.

Usage:
    python code/check_fold.py runs/octinoxate
    python code/check_fold.py runs/octinoxate --contact 4.0
"""
import argparse
import csv
import glob
import os
import re
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from uma_binding import parse_cif  # noqa: E402

HINT_RE = re.compile(r"token1: \[A, (\d+)\]\s*\n\s*token2: \[B, (\w+)\]\s*\n\s*"
                     r"max_distance: ([\d.]+)\s*\n\s*force: (\w+)")


def read_hints(yaml_path):
    """[(residue index, ligand atom name, max distance, forced)] from a written yaml."""
    if not os.path.exists(yaml_path):
        return []
    return [(int(r), a, float(d), f == "true")
            for r, a, d, f in HINT_RE.findall(open(yaml_path).read())]


def check(cif_path, boltz_dir, contact=4.0, wrap=4.5):
    name = os.path.basename(cif_path).replace("_model_0.cif", "")
    pep, lig = parse_cif(cif_path)
    P = np.array([a["xyz"] for a in pep])
    L = np.array([a["xyz"] for a in lig])
    d = np.linalg.norm(P[:, None, :] - L[None, :, :], axis=-1)
    row = {"name": name, "closest_approach": round(float(d.min()), 2),
           "contacts_under_cutoff": int((d < contact).sum()),
           "bound": int(d.min() < contact)}

    # How much of the ligand the peptide actually wraps. This is the design objective, and it is
    # not implied by the others: the forced-contact folds satisfied almost none of their hints while
    # engaging 100% of the ligand, and the ESM2 variants matched the original on summed energy while
    # engaging 60% and 40% against its 75%. Per ligand heavy atom, is there a peptide heavy atom
    # within `wrap` angstroms.
    ph = np.array([a["xyz"] for a in pep if a["element"] != "H"])
    lh = np.array([a["xyz"] for a in lig if a["element"] != "H"])
    per_ligand = np.linalg.norm(lh[:, None, :] - ph[None, :, :], axis=-1).min(axis=1)
    # Enclosure: is the peptide *around* the ligand, or pressed against one face of it? Wrapping
    # alone does not tell you -- the unconstrained glycine fold contacts 75% of the ligand's atoms
    # while enclosing only 49% of the directions out of it, with the two centroids 20.9 A apart.
    # Sample outward directions from the ligand centre and ask which ones run into peptide.
    n_dir = 200
    idx = np.arange(n_dir) + 0.5
    polar = np.arccos(1 - 2 * idx / n_dir)
    azim = np.pi * (1 + 5 ** 0.5) * idx
    dirs = np.c_[np.cos(azim) * np.sin(polar), np.sin(azim) * np.sin(polar), np.cos(polar)]
    centre = lh.mean(axis=0)
    v = ph - centre
    along = v @ dirs.T
    across = np.sqrt(np.maximum((v ** 2).sum(1)[:, None] - along ** 2, 0.0))
    reached = ((along > 0) & (along < 12.0) & (across < 3.0)).any(axis=0)
    row["enclosed_fraction"] = round(float(reached.mean()), 3)
    row["centroid_separation"] = round(float(np.linalg.norm(centre - ph.mean(axis=0))), 1)
    row["peptide_rg"] = round(float(np.sqrt(((ph - ph.mean(axis=0)) ** 2).sum(1).mean())), 1)

    row["ligand_heavy_atoms"] = len(lh)
    row["engaged"] = int((per_ligand < wrap).sum())
    row["wrapped_fraction"] = round(float((per_ligand < wrap).mean()), 3)
    row["mean_ligand_distance"] = round(float(per_ligand.mean()), 2)

    hints = read_hints(os.path.join(boltz_dir, name + ".yaml"))
    by_atom = {a["name"]: np.array(a["xyz"]) for a in lig}
    by_res = {}
    for a in pep:
        by_res.setdefault(int(a["resseq"]), []).append(a["xyz"])
    details = []
    satisfied = 0
    # A constraint with force:false never reaches the model -- the featurizer skips it -- so a
    # fold whose yaml carries only unforced hints is an unconstrained fold, and reporting "0/12
    # satisfied" for it invents a failure that was never attempted.
    applied = [h for h in hints if h[3]]
    if hints and not applied:
        row["hints_applied"] = 0
        hints = []
    else:
        row["hints_applied"] = len(applied)
        hints = applied
    for res, atom, want, forced in hints:
        if res not in by_res or atom not in by_atom:
            continue
        got = float(np.linalg.norm(np.array(by_res[res]) - by_atom[atom], axis=1).min())
        ok = got <= want
        satisfied += ok
        details.append((res, atom, want, got, ok, forced))
    row["n_hints"] = len(details)
    row["hints_satisfied"] = satisfied
    row["hint_rate"] = round(satisfied / len(details), 3) if details else ""
    row["worst_hint_miss"] = round(max((g - w for _, _, w, g, ok, _ in details if not ok),
                                      default=0.0), 2)
    row["forced"] = int(any(f for *_, f in details)) if details else ""
    return row, details


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("outdir")
    parser.add_argument("--contact", type=float, default=4.0,
                        help="heavy-atom separation that counts as a contact, A (default 4.0)")
    parser.add_argument("--wrap", type=float, default=4.5,
                        help="a ligand atom counts as engaged within this, A (default 4.5)")
    parser.add_argument("--quiet", action="store_true", help="table only, no per-hint lines")
    parser.add_argument("--out", default="fold_check.csv")
    args = parser.parse_args(argv)

    boltz_dir = os.path.join(args.outdir, "boltz")
    cifs = sorted(glob.glob(os.path.join(boltz_dir, "boltz_results_*", "predictions", "*",
                                         "*_model_0.cif")))
    if not cifs:
        sys.exit(f"no folded complexes under {boltz_dir}")

    rows = []
    for cif in cifs:
        row, details = check(cif, boltz_dir, args.contact, args.wrap)
        rows.append(row)
        if details and not args.quiet:
            print(f"\n{row['name']}")
            for res, atom, want, got, ok, forced in details:
                print(f"  residue {res:>3} -> {atom:<4} asked <= {want:4.1f} A, got {got:6.2f} A"
                      f"  {'ok' if ok else 'MISSED'}{'  (forced)' if forced else ''}")

    print(f"\n{'structure':<34} {'enclosed':>9} {'wrapped':>8} {'engaged':>8} {'cent sep':>9} "
          f"{'closest':>8} {'hints':>13}")
    for r in sorted(rows, key=lambda x: -x["enclosed_fraction"]):
        hint = (f"{r['hints_satisfied']}/{r['n_hints']}" if r["n_hints"]
                else ("none applied" if r.get("hints_applied") == 0 else "-"))
        print(f"{r['name'][:34]:<34} {r['enclosed_fraction']:>9.2f} {r['wrapped_fraction']:>8.2f} "
              f"{str(r['engaged']) + '/' + str(r['ligand_heavy_atoms']):>8} "
              f"{r['centroid_separation']:>9.1f} {r['closest_approach']:>8.2f} {hint:>13}")
    print("\n'enclosed' = fraction of directions out of the ligand that run into peptide: is it in a"
          "\nshell or against one face. 'wrapped' = fraction of ligand atoms with peptide within "
          "4.5 A.\nThey disagree usefully -- a fold can contact three quarters of the ligand while "
          "leaving it\noutside, which the centroid separation shows.")

    out = os.path.join(boltz_dir, args.out)
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print(f"\nresults in {out}")

    ignored = [r for r in rows if r["n_hints"] and r["hints_satisfied"] == 0]
    if ignored:
        print(f"\n{len(ignored)} fold(s) with enforced contacts satisfied none of them. Check the "
              f"fold log says\n'--use_potentials': without it a forced contact is conditioning "
              f"only, and the model may decline.")
    unapplied = [r for r in rows if r.get("hints_applied") == 0]
    if unapplied:
        print(f"\n{len(unapplied)} fold(s) carried hints that Boltz discarded (force:false) -- "
              f"treat these as unconstrained controls, not as failures.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
