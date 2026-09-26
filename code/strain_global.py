"""Recompute ligand strain for every fold against the one global reference.

`strain_i = E(bound_i) - E_ref`, where `E_ref` comes from `ligand_reference.py` and is the same number
for every structure in the run. That replaces the old definition, which relaxed each complex's own
bound pose and so measured every structure against a different local minimum.

**The bound state is the ligand as it exists in the complex, hydrogens included.** Strain is the energy
released on going from that state to the totally free one, so `E(bound_i)` must be the ligand with the
hydrogens the complex relaxation gave it -- positioned in the peptide's field. Re-placing those
hydrogens on the isolated ligand and relaxing them changes the initial state to one the ligand never
occupies, and relaxes away part of the very energy being measured. It was tried, and it lowered every
one of 24 strains by 0.2 to 6.1 kcal/mol. Do not do it.

That is also why the peptide's influence on those hydrogens is not double-counted against the
interaction energy: interaction is evaluated at one fixed geometry and says nothing about what is
released on reaching the free state. They are different legs of the same cycle.

So `E(bound_i)` comes from one of two places, in this order:

1. `boltz/structures/<name>_ligand_bound.xyz`, written by `binding_energy.py` -- the bound ligand with
   its complex-relaxed hydrogens. A single point on it, no relaxation.
2. failing that, the scoring log: the step-0 energy of the `relaxing free ligand` BFGS block is exactly
   that same quantity, which is how folds scored before geometry saving existed can still be corrected
   without re-running anything.

Results go to a csv of their own; nothing existing is modified.

Usage:
    python code/strain_global.py runs/octinoxate
    python code/strain_global.py runs/octinoxate --match s2_ --out strain_global_shell2.csv
"""
import argparse
import csv
import glob
import json
import os
import sys

import numpy as np


def bound_energy(name, boltz_dir, calc):
    """E(ligand in the complex, with its complex-relaxed hydrogens), in kcal/mol-equivalent.

    Two sources, in order. Both are the same physical state -- the bound ligand exactly as the complex
    relaxation left it -- so they are interchangeable and neither involves touching a hydrogen.

    1. `structures/<name>_ligand_bound.xyz`, written by `binding_energy.py`. A single point.
    2. the scoring log: the step-0 energy of that fold's `relaxing free ligand` BFGS block is the
       energy of the bound ligand, because step 0 is evaluated before the optimiser moves anything.
       This is how folds scored before geometry saving existed are handled.
    """
    import glob as _glob
    import re as _re
    from binding_energy import energy as _energy, EV_TO_KCAL as _EV

    path = os.path.join(boltz_dir, "structures", f"{name}_ligand_bound.xyz")
    if os.path.exists(path):
        lines = open(path).read().splitlines()
        n = int(lines[0])
        sym, xyz = [], []
        for ln in lines[2:2 + n]:
            f = ln.split()
            sym.append(f[0])
            xyz.append([float(f[1]), float(f[2]), float(f[3])])
        return _energy(sym, np.array(xyz), 0, calc) * _EV, "xyz"

    # A structure scored more than once appears in more than one log. Those are repeat measurements --
    # they differ by the pdbfixer hydrogen spread, 0.21 kcal/mol for the one structure scored twice
    # here -- so take the most recent and say that others existed, rather than silently depending on
    # filename order.
    logs = os.path.join(os.path.dirname(boltz_dir), "uma_logs", "score_*.log")
    hits = []
    for log in sorted(_glob.glob(logs), key=os.path.getmtime):
        lines = open(log, errors="replace").read().splitlines()
        cur = None
        for j, ln in enumerate(lines):
            m = _re.match(r"^([A-Za-z0-9_]+)\s*$", ln)
            if m and not ln.startswith(("BFGS", "Step")):
                cur = m.group(1)
            if "relaxing free ligand" in ln and cur == name:
                for k in range(j + 1, min(j + 5, len(lines))):
                    b = _re.match(r"BFGS:\s+0\s+\S+\s+(-?\d+\.\d+)", lines[k])
                    if b:
                        hits.append((float(b.group(1)) * _EV, os.path.basename(log)))
                        break
    if not hits:
        return None, None
    e, src = hits[-1]
    if len(hits) > 1:
        spread = max(h[0] for h in hits) - min(h[0] for h in hits)
        print(f"    note: {name} was scored {len(hits)} times, spread {spread:.3f} kcal/mol; "
              f"using the most recent ({src})", flush=True)
    return e, src


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0],
                               formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("outdir")
    p.add_argument("--match", help="only folds whose name contains this")
    p.add_argument("--reference", default="ligand_reference.json")
    p.add_argument("--out", default="strain_global.csv")
    p.add_argument("--model", default="uma-s-1p2p1")
    args = p.parse_args(argv)

    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from boltz_check import ligand_smiles
    from binding_energy import EV_TO_KCAL  # noqa: F401  (bound_energy imports its own)
    from fairchem.core import pretrained_mlip, FAIRChemCalculator

    ref_path = os.path.join(args.outdir, args.reference)
    if not os.path.exists(ref_path):
        sys.exit(f"no {ref_path}; run code/ligand_reference.py first")
    ref = json.load(open(ref_path))
    e_ref = float(ref["reference_kcal"])
    print(f"reference: {e_ref:.3f} kcal/mol-equivalent "
          f"(conformer {ref['reference_conformer']} of {ref['n_conformers']}, "
          f"spread {ref['spread_kcal']} across {ref['keep']} relaxed)")

    # The reference is a property of the ligand and the potential, so it is computed once and reused
    # for every shell of the same ligand -- which is the point. That also makes it the kind of file
    # that outlives the assumptions behind it, so check them: a reference from a different molecule or
    # a different model is not a smaller error than a missing one, it is a silently wrong strain for
    # every structure.
    smiles_now = ligand_smiles(args.outdir)
    if ref.get("smiles") != smiles_now:
        sys.exit(f"reference was computed for a different ligand:\n"
                 f"  reference: {ref.get('smiles')}\n"
                 f"  this run:  {smiles_now}\n"
                 f"delete {os.path.basename(ref_path)} and re-run ligand_reference.py")
    if ref.get("model") != args.model:
        sys.exit(f"reference was computed with model {ref.get('model')!r} but this run uses "
                 f"{args.model!r}. Energies from different potentials are not comparable; re-run "
                 f"ligand_reference.py with --model {args.model}")

    boltz_dir = os.path.join(args.outdir, "boltz")
    cifs = sorted(glob.glob(os.path.join(boltz_dir, "boltz_results_*", "predictions", "*",
                                         "*_model_0.cif")))
    if args.match:
        cifs = [c for c in cifs if args.match in os.path.basename(c)]
    if not cifs:
        sys.exit("no folded complexes found")
    print(f"{len(cifs)} folds\n")

    calc = FAIRChemCalculator(pretrained_mlip.get_predict_unit(args.model, device="cpu"),
                             task_name="omol")
    rows, missing = [], []
    for i, cif in enumerate(cifs, 1):
        name = os.path.basename(cif).replace("_model_0.cif", "")
        e_bound, src = bound_energy(name, boltz_dir, calc)
        if e_bound is None:
            missing.append(name)
            print(f"[{i}/{len(cifs)}] {name:<38} SKIPPED: no bound ligand geometry and none in the "
                  f"logs", flush=True)
            continue
        strain = e_bound - e_ref
        rows.append({"name": name, "e_bound_kcal": f"{e_bound:.3f}",
                     "e_reference_kcal": f"{e_ref:.3f}",
                     "strain_global_kcal": f"{strain:.3f}", "source": src})
        print(f"[{i}/{len(cifs)}] {name:<38} strain {strain:>8.3f}   [{src}]", flush=True)

    if missing:
        print(f"\n{len(missing)} fold(s) had no bound ligand: {', '.join(missing)}")
        print("Score them first -- binding_energy.py writes structures/<name>_ligand_bound.xyz.")
    if not rows:
        sys.exit("nothing to write")

    out = os.path.join(boltz_dir, args.out)
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    vals = [float(r["strain_global_kcal"]) for r in rows]
    print(f"\nstrain range {min(vals):.2f} to {max(vals):.2f} kcal/mol")
    if min(vals) < -0.5:
        print(f"WARNING: a negative strain means a bound pose below the reference, so the conformer "
              f"search missed the global minimum. Re-run ligand_reference.py with more conformers.")
    print(f"wrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
