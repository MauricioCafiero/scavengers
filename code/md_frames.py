"""Pull representative frames from the end of a trajectory, as PDB files for figures.

`make_figures.py` does this for the static Boltz folds. This is the dynamics equivalent, and the
reason it is needed is that the Boltz pose and the pose after 20 ns are not the same structure -- on
`orig_f12` the ligand leaves the designed cavity at 6 ns, and a figure showing only the prediction
would misrepresent what the simulation found.

Frames come from the **end** of the run, because that is the part that has stopped depending on the
starting pose. Two files are written per leg:

* `frames_last.pdb`   -- a multi-model PDB, evenly spaced through the final window. PyMOL loads the
  models as states, so `set all_states, on` shows the ensemble and the spread of the ligand in it.
* `frame_medoid.pdb`  -- the single frame closest to the window's mean structure, for a clean render
  where an ensemble would be mush. This is a real frame, not an average: averaging coordinates over
  an ensemble produces bond lengths that no force field would allow.

Everything is superposed on the **peptide** heavy atoms, so the ligand's movement relative to its
binding site is what the picture shows, rather than the whole complex tumbling.

Usage:
    ~/miniforge3/envs/openmm-md/bin/python code/md_frames.py runs/octinoxate/md/*/prod*
    ... --last-fraction 0.1 --n-frames 10
"""
import argparse
import os


def extract(leg_dir, last_fraction=0.10, n_frames=10, window_ns=None, tag=None):
    import mdtraj as md
    import numpy as np

    xtc = os.path.join(leg_dir, "traj_wrapped.xtc")
    pdb = os.path.join(leg_dir, "traj_wrapped.pdb")
    if not (os.path.exists(xtc) and os.path.exists(pdb)):
        return None
    t = md.load(xtc, top=pdb)
    top = t.topology

    lig_res = {"LIG", "UNK", "UNL"}
    pep_heavy = [a.index for a in top.atoms
                 if a.residue.name not in lig_res and a.element.symbol != "H"]
    if not pep_heavy:
        raise ValueError(f"no peptide heavy atoms in {pdb}")

    if window_ns:
        lo, hi = window_ns
        keep = [i for i, ps in enumerate(t.time) if lo * 1000 <= ps <= hi * 1000]
        if not keep:
            raise ValueError(f"no frames in {lo}-{hi} ns; run spans "
                             f"{t.time[0]/1000:.2f}-{t.time[-1]/1000:.2f} ns")
        window = t[keep]
    else:
        start = max(0, int(t.n_frames * (1.0 - last_fraction)))
        window = t[start:]
    # Superpose on the peptide so the ligand's motion in the site is what varies between models.
    window.superpose(window, frame=0, atom_indices=pep_heavy)

    # Medoid: the frame closest to the window's mean coordinates. Chosen over an averaged structure
    # because an average is not a physical geometry.
    mean_xyz = window.xyz[:, pep_heavy, :].mean(axis=0)
    dev = np.sqrt(((window.xyz[:, pep_heavy, :] - mean_xyz) ** 2).sum(-1).mean(-1))
    medoid = int(dev.argmin())

    picks = np.unique(np.linspace(0, window.n_frames - 1, n_frames).astype(int))
    suffix = f"_{tag}" if tag else ""
    ens_path = os.path.join(leg_dir, f"frames_last{suffix}.pdb")
    med_path = os.path.join(leg_dir, f"frame_medoid{suffix}.pdb")
    window[picks].save(ens_path)
    window[medoid].save(med_path)

    return {"leg": leg_dir, "n_frames_total": int(t.n_frames),
            "window_from_ns": round(float(window.time[0]) / 1000, 2),
            "window_to_ns": round(float(window.time[-1]) / 1000, 2),
            "models_written": len(picks),
            "medoid_ns": round(float(window.time[medoid]) / 1000, 2),
            "medoid_rmsd_to_mean_A": round(float(dev[medoid]) * 10, 2),
            "spread_A": round(float(dev.mean()) * 10, 2),
            "ensemble": ens_path, "medoid": med_path}


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0],
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("legs", nargs="+", help="production directories holding traj_wrapped.xtc/.pdb")
    p.add_argument("--last-fraction", type=float, default=0.10,
                   help="fraction of the run to draw from, measured from the end (default 0.10)")
    p.add_argument("--n-frames", type=int, default=10,
                   help="models in the ensemble PDB (default 10)")
    # The end of a run is usually the settled part, but not always: orig_f12's ligand re-inserts in
    # its last 2 ns after twelve on the surface, so its final window misrepresents the run. Name the
    # window explicitly in that case, and tag the output so both sets can coexist.
    p.add_argument("--window-ns", nargs=2, type=float, metavar=("FROM", "TO"),
                   help="take frames from this interval instead of the trailing fraction")
    p.add_argument("--tag", help="suffix for the output filenames, e.g. --tag plateau")
    args = p.parse_args(argv)

    done = []
    for leg in args.legs:
        r = extract(leg, args.last_fraction, args.n_frames, args.window_ns, args.tag)
        if r is None:
            print(f"skipping {leg}: no traj_wrapped.xtc/.pdb")
            continue
        done.append(r)
        print(f"\n=== {leg}")
        print(f"  {r['n_frames_total']} frames total; window {r['window_from_ns']}-"
              f"{r['window_to_ns']} ns")
        print(f"  {r['ensemble']}  ({r['models_written']} models)")
        print(f"  {r['medoid']}  (frame at {r['medoid_ns']} ns, "
              f"{r['medoid_rmsd_to_mean_A']} A from the window mean; "
              f"window spread {r['spread_A']} A)")

    if done:
        print("\nin PyMOL, for the ensemble:")
        print("  load frames_last.pdb, ens")
        print("  set all_states, on")
        print("  hide everything; show cartoon, polymer; show sticks, not polymer")
        print("  color grey70, polymer; color red, not polymer")
        print("\nor for a single clean frame:")
        print("  load frame_medoid.pdb, best")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
