"""How the designed pose holds up over a trajectory: contacts, separation, radius of gyration.

MM/GBSA gives one number for a whole run. This gives the shape of the run behind it, which is what
says whether that number describes the designed pose or something the pose decayed into. On
`orig_f12` the two answers are different: ΔG converges to −13.71 kcal/mol, but the ligand leaves the
cavity at about 6 ns, sits on the surface for twelve, and re-inserts at the end. A single ΔG cannot
show that, and the decision it implies -- believe the enclosure or not -- depends on seeing it.

Everything is measured on **heavy atoms only**. Including hydrogens roughly quadruples the contact
count and makes runs incomparable with the static `contacts_under_cutoff` in `fold_check*.csv`, which
is heavy-atom. The centroid separation uses all peptide atoms, matching `check_fold.py`.

Usage:
    python code/md_contacts.py runs/octinoxate/md/orig_f12/prod_20ns
    python code/md_contacts.py runs/octinoxate/md/*/prod* --deciles 10 --csv out.csv

Needs mdtraj, so run it with the openmm-md interpreter:
    ~/miniforge3/envs/openmm-md/bin/python code/md_contacts.py ...
"""
import argparse
import csv
import os
import sys


def profile(leg_dir, cutoff=4.0, deciles=10, resample_ps=10.0):
    """Per-decile contacts, centroid separation and peptide Rg for one production directory."""
    import mdtraj as md
    import numpy as np

    xtc = os.path.join(leg_dir, "traj_wrapped.xtc")
    pdb = os.path.join(leg_dir, "traj_wrapped.pdb")
    if not (os.path.exists(xtc) and os.path.exists(pdb)):
        return None
    t = md.load(xtc, top=pdb)
    top = t.topology

    lig_res = {"LIG", "UNK", "UNL"}
    lig = [a.index for a in top.atoms if a.residue.name in lig_res and a.element.symbol != "H"]
    pep_all = [a.index for a in top.atoms if a.residue.name not in lig_res]
    pep = [i for i in pep_all if top.atom(i).element.symbol != "H"]
    if not lig:
        raise ValueError(f"no ligand residue in {pdb}; looked for {sorted(lig_res)}")

    sep = np.linalg.norm(t.xyz[:, lig, :].mean(1) - t.xyz[:, pep_all, :].mean(1), axis=1) * 10
    d = md.compute_distances(t, np.array([(i, j) for i in lig for j in pep])) * 10
    nc = (d < cutoff).sum(1)
    rg = md.compute_rg(t.atom_slice(pep_all)) * 10

    # Release statistics, the measures that actually separate these trajectories. Mean contact
    # count does not: a ligand that wanders out and back averages to a respectable number while
    # spending a fifth of the run detached. What matters is whether it ever lets go, for how long,
    # and whether the episodes cluster late -- progressive loss rather than thermal flicker.
    #
    # These are computed on a COMMON time spacing, because they are otherwise not comparable
    # between runs. Sampled every 1 ps, orig_f12 shows 8 released frames; every 10 ps it shows 0,
    # because all eight were isolated single frames. Finer sampling catches more brief excursions,
    # so comparing a 1 ps run against a 10 ps one overstates the finely sampled one's release rate.
    dt_ps = float(t.time[1] - t.time[0]) if t.n_frames > 1 else 0.0
    stride = max(1, int(round(resample_ps / dt_ps))) if dt_ps else 1
    rel = (nc == 0)[::stride]
    runs, run = [], 0
    for x in rel:
        if x:
            run += 1
        elif run:
            runs.append(run)
            run = 0
    if run:
        runs.append(run)
    late = float(rel[len(rel) // 2:].sum() / rel.sum()) if rel.sum() else 0.0
    release = {"sample_ps": round(dt_ps * stride, 1),
               "released_frames": int(rel.sum()),
               "released_pct": round(100 * float(rel.mean()), 2),
               "longest_episode_ps": round(max(runs) * dt_ps * stride, 1) if runs else 0.0,
               "time_released_ns": round(rel.sum() * stride * dt_ps / 1000, 3),
               "episodes": len(runs),
               "share_in_second_half": round(late, 2)}

    edges = np.linspace(0, t.n_frames, deciles + 1).astype(int)
    rows = []
    for a, b in zip(edges[:-1], edges[1:]):
        if b <= a:
            continue
        rows.append({"ns": round(float(t.time[b - 1]) / 1000, 2),
                     "contacts": round(float(nc[a:b].mean()), 1),
                     "separation_A": round(float(sep[a:b].mean()), 2),
                     "peptide_rg_A": round(float(rg[a:b].mean()), 2)})
    return {"leg": leg_dir, "frames": int(t.n_frames),
            "n_ligand_heavy": len(lig), "n_peptide_heavy": len(pep),
            "contacts_mean": round(float(nc.mean()), 1),
            "contacts_first": rows[0]["contacts"], "contacts_last": rows[-1]["contacts"],
            "separation_first_A": rows[0]["separation_A"], "separation_last_A": rows[-1]["separation_A"],
            "separation_max_A": round(float(sep.max()), 2),
            "separation_mean_A": round(float(sep.mean()), 2),
            "closest_approach_A": round(float(d.min()), 2),
            # A frame with no contact at all is the only unambiguous sign of release. One or two out
            # of thousands is thermal; a run of them is dissociation.
            "zero_contact_frames": int((nc == 0).sum()),
            # fraction of the run with the ligand centroid within 10 A of the peptide centroid
            "residence_pct": round(100 * float((sep <= 10.0).mean()), 1),
            "release": release,
            "rows": rows}


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0],
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("legs", nargs="+", help="production directories holding traj_wrapped.xtc/.pdb")
    p.add_argument("--cutoff", type=float, default=4.0, help="contact cutoff, A (default 4.0)")
    p.add_argument("--deciles", type=int, default=10, help="windows across the run (default 10)")
    p.add_argument("--resample-ps", type=float, default=10.0,
                   help="common frame spacing for the release statistics, ps (default 10.0). Runs "
                        "saved at different intervals are not comparable without this")
    p.add_argument("--csv", help="also write the per-window rows here")
    args = p.parse_args(argv)

    out, summaries = [], []
    for leg in args.legs:
        r = profile(leg, args.cutoff, args.deciles, args.resample_ps)
        if r is None:
            print(f"skipping {leg}: no traj_wrapped.xtc/.pdb -- run `omd analyze` or "
                  f"`modal_md.py strip` first")
            continue
        summaries.append(r)
        print(f"\n=== {leg}")
        print(f"{r['frames']} frames, {r['n_peptide_heavy']} peptide / {r['n_ligand_heavy']} "
              f"ligand heavy atoms")
        print(f"{'ns':>7} {'contacts':>9} {'sep A':>7} {'pep Rg':>7}")
        for w in r["rows"]:
            print(f"{w['ns']:7.1f} {w['contacts']:9.1f} {w['separation_A']:7.2f} "
                  f"{w['peptide_rg_A']:7.2f}")
            out.append(dict(leg=leg, **w))
        rl = r["release"]
        print(f"whole run: contacts {r['contacts_first']} -> {r['contacts_last']} "
              f"(mean {r['contacts_mean']}), separation {r['separation_first_A']} -> "
              f"{r['separation_last_A']} A (max {r['separation_max_A']}), "
              f"closest heavy-atom approach {r['closest_approach_A']} A")
        print(f"  residence within 10 A: {r['residence_pct']}% of the run")
        print(f"  release (at {rl['sample_ps']} ps spacing): {rl['released_frames']} frames "
              f"({rl['released_pct']}%) in {rl['episodes']} episodes, longest "
              f"{rl['longest_episode_ps']} ps, {rl['time_released_ns']} ns detached, "
              f"{rl['share_in_second_half']:.0%} of it in the second half")

    if len(summaries) > 1:
        print(f"\n{'leg':<44}{'contacts':>16}{'mean sep':>10}{'resid%':>8}"
              f"{'rel%':>7}{'longest ps':>12}{'late':>7}")
        for r in summaries:
            rl = r["release"]
            contacts = f"{r['contacts_first']} -> {r['contacts_last']}"
            print(f"{os.path.relpath(r['leg']):<44}{contacts:>16}"
                  f"{r['separation_mean_A']:>10.2f}{r['residence_pct']:>8.1f}"
                  f"{rl['released_pct']:>7.2f}{rl['longest_episode_ps']:>12.1f}"
                  f"{rl['share_in_second_half']:>7.0%}")

    if args.csv and out:
        with open(args.csv, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(out[0]))
            w.writeheader()
            w.writerows(out)
        print(f"\nwrote {args.csv}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
