"""Render an MD trajectory to an mp4 for presentations, headless, from the pymol-render env.

Same look as code/render_md_figures.py (violetpurple peptide, skyblue ligand, dss before cartoon,
40% cartoon transparency) so a movie and the committed stills read as one set. Two differences a
movie needs:

  * `intra_fit` on the peptide CAs first. Without it the whole complex tumbles and rotates through
    the movie, which hides the only motion worth watching -- the ligand against the peptide.
  * one camera for the whole movie, set from the fitted ensemble, rather than re-orienting per frame.

The background is opaque (mp4 has no alpha), white by default.

    ~/miniforge3/envs/pymol-render/bin/python code/render_md_movie.py <prod dir> [--out x.mp4]
"""
import argparse
import glob
import os
import subprocess
import sys
import tempfile


def main():
    p = argparse.ArgumentParser()
    p.add_argument("prod_dir", help="production dir holding traj_wrapped.{xtc,pdb}")
    p.add_argument("--out", help="output mp4 (default: <prod dir>/../<leg>_movie.mp4)")
    p.add_argument("--frames", type=int, default=400, help="frames in the movie (default 400)")
    p.add_argument("--fps", type=int, default=20)
    p.add_argument("--width", type=int, default=1280)
    p.add_argument("--height", type=int, default=960)
    # Negative: cmd.zoom fits a bounding SPHERE, which leaves a compact solute filling only ~35% of
    # a frame at buffer 0. -8 fills ~70% with no frame touching any edge across the trajectory
    # (measured over 13 states of ox2_orig_f8); -10 already clips. The zoom spans every state, so
    # the camera accommodates the whole movie and nothing walks out of shot.
    p.add_argument("--buffer", type=float, default=-8.0, help="zoom padding, A (negative = tighter)")
    p.add_argument("--carbon-color", default="violetpurple")
    p.add_argument("--ligand-carbon-color", default="skyblue")
    p.add_argument("--background", default="white")
    a = p.parse_args()

    pdb = os.path.join(a.prod_dir, "traj_wrapped.pdb")
    xtc = os.path.join(a.prod_dir, "traj_wrapped.xtc")
    for f in (pdb, xtc):
        if not os.path.exists(f):
            sys.exit(f"missing {f}")
    leg = os.path.basename(os.path.dirname(os.path.abspath(a.prod_dir)))
    out = a.out or os.path.join(os.path.dirname(os.path.abspath(a.prod_dir)), f"{leg}_movie.mp4")

    import imageio
    import imageio_ffmpeg  # noqa: F401  (provides the bundled ffmpeg imageio calls)
    import pymol
    pymol.finish_launching(["pymol", "-cq"])
    from pymol import cmd

    cmd.load(pdb, "m")
    n_total = _frame_count(xtc)
    interval = max(1, round(n_total / a.frames))
    cmd.load_traj(xtc, "m", state=1, interval=interval)
    n = cmd.count_states("m")
    print(f"{leg}: {n_total} frames -> {n} states (every {interval})", flush=True)

    cmd.intra_fit("polymer and name CA")       # kill the tumble; see module docstring
    cmd.dss("m", state=1)                      # PDBs carry no HELIX/SHEET: assign from geometry
    cmd.hide("everything")
    cmd.show("cartoon", "polymer")
    cmd.show("sticks", "organic")
    cmd.color(a.carbon_color, "polymer and elem C")
    cmd.color(a.ligand_carbon_color, "organic and elem C")
    cmd.set("cartoon_transparency", 0.4)
    cmd.set("ray_opaque_background", 1)
    cmd.bg_color(a.background)
    cmd.orient("m")                            # one camera, from the fitted ensemble
    cmd.zoom("m", a.buffer, complete=1)

    cmd.mset(f"1 -{n}")                        # one movie frame per loaded state
    cmd.set("ray_trace_frames", 1)             # mpng has no ray= argument; this is the switch
    tmp = tempfile.mkdtemp(prefix="pymovie_")
    cmd.mpng(os.path.join(tmp, "f"), width=a.width, height=a.height)
    pngs = sorted(glob.glob(os.path.join(tmp, "f*.png")))
    print(f"rendered {len(pngs)} frames -> {out}", flush=True)
    with imageio.get_writer(out, fps=a.fps, codec="libx264", quality=8,
                            macro_block_size=8, ffmpeg_log_level="error") as w:
        for f in pngs:
            w.append_data(imageio.imread(f))
    for f in pngs:
        os.remove(f)
    os.rmdir(tmp)
    print(f"done: {out} ({os.path.getsize(out) / 1e6:.1f} MB, {len(pngs) / a.fps:.1f} s)")


def _frame_count(xtc):
    """Frame count without loading coordinates: mdtraj is not in this env, so read the header count
    the cheap way -- PyMOL reports it while loading, but we need it before. gmx is absent too, so
    fall back to a one-pass xdrfile-free scan: the magic number 1995 starts every xtc frame header."""
    import struct
    n = 0
    with open(xtc, "rb") as fh:
        data = fh.read(8)
        if struct.unpack(">i", data[:4])[0] != 1995:
            return 0
        natoms = struct.unpack(">i", data[4:8])[0]
        fh.seek(0)
        while True:
            head = fh.read(92)                      # xtc frame header through the compressed size
            if len(head) < 92:
                break
            size = struct.unpack(">i", head[88:92])[0]
            fh.seek((size + 3) // 4 * 4, os.SEEK_CUR)
            n += 1
    del natoms
    return n


if __name__ == "__main__":
    main()
