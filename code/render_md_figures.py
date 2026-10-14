"""Render MD solute structures to PNG with headless PyMOL, in the project's medoid-figure style.

The input is `frame_medoid.pdb` under `runs/<mol>/md/<leg>/prod_*/`, written by `md_frames.py` as
part of every leg: solute only, peptide plus ligand, and -- unlike the Boltz co-folds
`make_figures.py` gathers -- carrying hydrogens. The medoid is the whole figure set; the
`frames_last.pdb` sitting beside it is skipped even when named, unless `--include-last`. Medoids
exist only for legs whose MD has finished and been tailed, so the set grows as the campaign lands
and a re-run is expected to skip what it has already rendered.

**The style** follows the hand-made images in `runs/octinoxate/md/figures/`: transparent
background, the peptide as a cartoon with no side chains, the ligand as sticks coloured by element,
and no labels. It departs from them in three deliberate places, each revertible by one flag:

* the cartoon is `automatic`, not `tube`, because a tube draws a helical peptide as an unreadable
  squiggle (`--cartoon tube`);
* the cartoon is slightly transparent, so backbone passing in front of the ligand does not hide it
  (`--cartoon-transparency 0`);
* the ligand's carbons contrast with the peptide's instead of matching them, so the ligand is
  obvious at a glance (`--ligand-carbon-color teal` for the single-colour look).

The black-background all-sticks style of `figures_shell*/style.pml` is the *fold*-figure style and
is not offered here -- `--style` takes a .pml of your own.

**The camera.** The ligand is the subject of these pictures, and both halves of the view are built
around it. `orient` on the ligand puts the ligand's own smallest principal axis down the barrel of
the camera, which is the same thing as showing the molecule face-on rather than edge-on. The zoom
is taken on the ligand too, with a buffer wide enough that the surrounding fold stays in frame as
context, so the ligand reads centred and large rather than as one feature of a whole-complex shot
(`--zoom all` frames the complex instead, if that is ever wanted). The whole recipe is fixed, so
one leg's medoid is directly comparable with another's.

Usage:
    # one leg
    $PYMOL_PY code/render_md_figures.py runs/octinoxate/md/s3_orig_f12_dock1/prod_20ns/frame_medoid.pdb

    # a whole lane, named explicitly (the repo rule: name the targets, do not scan)
    $PYMOL_PY code/render_md_figures.py runs/octinoxate/md/*/prod_20ns/frame_medoid.pdb

where $PYMOL_PY is the pymol-render env's python (see --pymol). Outputs land in
`runs/<mol>/figures_md/<leg>_medoid.png` unless --out-dir says otherwise, and an existing PNG is
left alone unless --force is passed, so a re-run after an interrupted batch resumes.
"""
import argparse
import os
import re
import shutil
import subprocess
import sys

# Default interpreter for the headless render: a conda env holding conda-forge's
# pymol-open-source. The pip wheel on PyPI (3.2.0a0) links libpng/freetype/GLEW/netcdf through an
# rpath baked to its builder's own mamba prefix and cannot be loaded here, so conda is not a
# preference in this one case -- see this module's entry in the README notes.
DEFAULT_PYMOL = os.path.expanduser("~/miniforge3/envs/pymol-render/bin/pymol")

RENDER_SETUP = """
set antialias, 2
set ray_trace_mode, 0
"""

# Every PNG the user has committed -- the fold figures in figures_shell*/ and the medoid figures in
# md/figures/ alike -- is RGBA with a transparent background. That is PyMOL's default
# (`ray_opaque_background` off), and it is why `figures_shell3/style.pml` says `bg_color black`
# while the files it produced read as white: there is no background in them at all, so a viewer
# supplies one. Matching the existing figures means rendering transparent, which is the default.
BACKGROUND = {
    "transparent": "set ray_opaque_background, 0\n",
    "black": "bg_color black\nset ray_opaque_background, 1\n",
    "white": "bg_color white\nset ray_opaque_background, 1\n",
}

# The style of `runs/octinoxate/md/figures/*_medoid.png`, read off those images rather than
# guessed. It was never written down as a .pml -- the user made them at the GUI prompt -- so this
# is its first written form.
#
# `dss` first: mdtraj writes these PDBs with no HELIX or SHEET records, so PyMOL loads every
# residue as loop and an automatic cartoon would draw a helical peptide as a flat ribbon. `dss`
# assigns secondary structure from the geometry, which is what makes `--cartoon automatic` able to
# show a helix at all. It is harmless under the default tube.
#
# `automatic` rather than the hand-made set's `tube`, decided by rendering both: a helical leg
# (bg33_1_o3cof, helix fraction 1.0) drawn as a tube is an unreadable zigzag squiggle that gives no
# hint it is a helix at all, where automatic shows it as two helices with the ligand lying in the
# groove between them. On a coil peptide, which most of these are, automatic falls back to a loop
# trace that looks very much like the tube anyway, so nothing is lost on the structures the
# hand-made set covered. `--cartoon tube` restores the hand-made look exactly.
#
# `cartoon_loop_radius` up from its 0.2 default: automatic's loop is thinner than a tube and the
# coil legs looked flimsy beside the existing figures at the stock value.
#
# Side chains stay hidden, so the only sticks in the picture are the ligand.
#
# Colour: `color <carbon>`, then the ligand's carbons over the top, then `util.cnc`, which
# recolours everything except carbon by element. That is what gives the figures their red oxygens
# and white hydrogens. The hand-made set varies the carbon colour between legs only because PyMOL
# auto-assigns a new object colour on each load -- here it is fixed, so legs are comparable.
#
# `violetpurple` peptide, `skyblue` ligand, settled by rendering the candidates rather than
# reasoning about them. Three findings drove it.
#
# 1. The ligand's own red oxygens have to stay readable against its carbons. Orange fails outright
#    -- red sits next to it in hue and the oxygens vanish -- magenta passes, and blue is the widest
#    hue gap available, so the oxygens separate most strongly there. That is why the ligand is blue
#    and the peptide is not.
# 2. The cartoon is 40% transparent, so wherever it crosses the ligand the two colours blend, and
#    that blend must not invent a third substance. The violet family over blue blends to violet,
#    which stays in family. An earlier teal/yellow pair blended to green and a yellow-green/magenta
#    pair blended to brown -- both read as a third material.
# 3. Loudness here is driven by SATURATION, not darkness, and a loud cartoon inverts the point of
#    these figures by making the peptide the subject instead of the ligand. `purple` (0.75 0 0.75)
#    is one of the darkest candidates and among the loudest, because its green channel is zero.
#    The calm shades are the ones carrying some green. So do not "quieten" this by going darker
#    along the zero-green line -- that brings the problem back in a darker costume.
#
# `skyblue` not raw `blue`, which is near-black in a cartoon and kills the shading.
#
# Shades from the comparison ladder that are NOT in PyMOL's palette, kept so they are not lost.
# Pass one with --set-color, e.g.  --set-color amethyst=0.65,0.32,0.72 --carbon-color amethyst
# `amethyst` was the runner-up to violetpurple and is one step lighter.
EXTRA_COLORS = {
    "orchid": (0.80, 0.42, 0.85),    # one step below `violet`, still fairly pastel
    "amethyst": (0.65, 0.32, 0.72),  # runner-up: solid, muted, a touch lighter than violetpurple
}
#
# No semicolons: PyMOL splits commands on one even inside a comment.
MD_PRESET = """# MD medoid style -- see code/render_md_figures.py
hide everything
dss
show cartoon, polymer
cartoon {cartoon}
set cartoon_loop_radius, 0.3
set cartoon_transparency, {cartoon_transparency}
show sticks, not polymer

color {carbon}
color {ligand_carbon}, not polymer
util.cnc
"""


def structure_label(path):
    """`runs/oxybenzone/md/ox1_f12_dock1/prod_20ns/frame_medoid.pdb` -> `ox1_f12_dock1_medoid`.

    The leg directory is the parent of the production directory, which is what the campaign calls a
    leg everywhere else. A production directory other than `prod_20ns` (the Modal legs'
    `prod_L1_modal`) is folded into the name as well, because the same leg can hold both and they
    would otherwise collide onto one filename.
    """
    path = os.path.abspath(path)
    base = os.path.basename(path)
    stem = os.path.splitext(base)[0]
    kind = {"frame_medoid": "medoid", "frames_last": "last"}.get(stem, stem)
    prod = os.path.basename(os.path.dirname(path))
    leg = os.path.basename(os.path.dirname(os.path.dirname(path)))
    tag = "" if prod in ("prod_20ns", "") else "_" + re.sub(r"^prod_", "", prod)
    return f"{leg}{tag}_{kind}"


def default_out_dir(path):
    """`runs/<mol>/figures_md`, read off the structure's own path, or None if it is not under runs/."""
    parts = os.path.abspath(path).split(os.sep)
    if "runs" in parts:
        i = parts.index("runs")
        if i + 1 < len(parts):
            return os.path.join(os.sep.join(parts[:i + 2]), "figures_md")
    return None


def load_style(args):
    """The representation block to apply, and a human-readable note of where it came from."""
    if args.style:
        with open(args.style) as fh:
            return fh.read(), args.style
    return (MD_PRESET.format(carbon=args.carbon_color,
                             ligand_carbon=args.ligand_carbon_color or args.carbon_color,
                             cartoon=args.cartoon,
                             cartoon_transparency=args.cartoon_transparency),
            f"md preset (carbons {args.carbon_color}, ligand "
            f"{args.ligand_carbon_color or args.carbon_color}, cartoon {args.cartoon} at "
            f"transparency {args.cartoon_transparency})")


def crop_pngs(paths, mode, margin_frac):
    """Trim the transparent surround off the rendered PNGs, in place.

    The renders come out with a lot of empty alpha around the molecule, because PyMOL's `zoom`
    fits a *sphere* around the selection: an elongated complex leaves the corners of a 4:3 frame
    empty, and the ligand-centred zoom used here leaves more.

    `mode="batch"` (the default) crops every image in the run to ONE box -- the union of their
    content boxes -- and that choice is the whole point. Every leg is rendered with the same camera
    on the same ligand, so angstroms-per-pixel is already identical across the set; cropping each
    image to its own content would keep that true but leave the images different sizes, and
    anything that then displays them at a common width (a contact sheet, a slide, a README table)
    silently rescales them and a small peptide ends up looking as big as a large one. One shared
    box keeps both the scale and the canvas identical, so the images stay honestly comparable
    however they are viewed. `mode="each"` is the per-image crop for when that does not matter.

    The box is computed over the images of *this run*, so a full campaign render should be one
    command -- which is how the bulk run is written anyway. A later top-up run (skip-existing)
    computes its own box; `--force` over the whole set re-renders and re-crops to one box again.
    """
    from PIL import Image

    boxes = {}
    for p in paths:
        with Image.open(p) as im:
            if im.mode != "RGBA":
                print(f"crop: {os.path.basename(p)} has no alpha channel, left alone")
                continue
            box = im.getchannel("A").getbbox()
            size = im.size
        if box:
            boxes[p] = (box, size)
    if not boxes:
        return

    if mode == "batch":
        lefts, tops, rights, bottoms = zip(*(b for b, _ in boxes.values()))
        shared = (min(lefts), min(tops), max(rights), max(bottoms))
        boxes = {p: (shared, size) for p, (_, size) in boxes.items()}

    for p, (box, size) in boxes.items():
        left, top, right, bottom = box
        margin = round(margin_frac * max(right - left, bottom - top))
        left, top = max(0, left - margin), max(0, top - margin)
        right, bottom = min(size[0], right + margin), min(size[1], bottom + margin)
        with Image.open(p) as im:
            im.crop((left, top, right, bottom)).save(p)
    w, h = right - left, bottom - top
    print(f"cropped {len(boxes)} image(s), {mode} box -> {w}x{h}"
          + (" (shared)" if mode == "batch" else ""))


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("structures", nargs="+",
                        help="pdb/cif files to render, named explicitly")
    parser.add_argument("--out-dir",
                        help="where the PNGs go. Default: runs/<mol>/figures_md beside each input, "
                             "so a mixed octinoxate/oxybenzone list still sorts itself")
    parser.add_argument("--carbon-color", default="violetpurple",
                        help="carbon colour for the peptide (everything else is by element)")
    parser.add_argument("--ligand-carbon-color", default="skyblue",
                        help="carbon colour for the ligand. It contrasts with the peptide by "
                             "default so the ligand reads at a glance. Pass the same value as "
                             "--carbon-color for the hand-made figures' single-colour look")
    parser.add_argument("--cartoon", choices=("automatic", "tube", "loop", "putty"),
                        default="automatic",
                        help="peptide cartoon. automatic (default) draws helices and sheets as "
                             "such -- `dss` runs first, so it works despite mdtraj writing no "
                             "HELIX records. tube is the hand-made figures' look, but it renders a "
                             "helix as an unreadable squiggle")
    parser.add_argument("--cartoon-transparency", type=float, default=0.4,
                        help="0 is the hand-made look. The default lets the ligand show through "
                             "backbone that passes in front of it, which is the whole point of "
                             "these pictures -- pass 0 for an opaque peptide")
    parser.add_argument("--set-color", action="append", default=[], metavar="NAME=R,G,B",
                        help="define a colour before the style runs, e.g. "
                             "--set-color orchid=0.80,0.42,0.85 . Repeatable. PyMOL's named "
                             "palette has gaps -- notably between `violet` and `violetpurple` -- "
                             "and this fills them without hardcoding anything")
    parser.add_argument("--style",
                        help="a .pml of your own to apply instead of the built-in style")
    parser.add_argument("--pymol", default=DEFAULT_PYMOL, help="pymol executable")
    parser.add_argument("--width", type=int, default=1600)
    parser.add_argument("--height", type=int, default=1200)
    parser.add_argument("--dpi", type=int, default=300)
    parser.add_argument("--orient", default="not polymer",
                        help="selection whose principal axes set the camera. The default, the "
                             "ligand, puts its smallest axis down the camera barrel, i.e. shows it "
                             "face-on. `polymer` frames the fold instead and buries the ligand")
    parser.add_argument("--zoom", default="not polymer",
                        help="selection to frame after orienting. The default frames the ligand, "
                             "so it sits centred and large with the fold around it as context. "
                             "`all` instead frames the whole complex, which centres the fold and "
                             "shrinks the ligand")
    parser.add_argument("--buffer", type=float, default=9.0,
                        help="angstroms of padding around the zoom selection. Sets how much "
                             "peptide context surrounds the ligand -- smaller is a tighter "
                             "close-up, larger backs off towards the whole fold")
    parser.add_argument("--hide-hydrogens", action="store_true",
                        help="hide hydrogens, which these MD structures carry and the co-folds "
                             "do not")
    parser.add_argument("--crop", choices=("batch", "each", "none"), default="batch",
                        help="trim the empty alpha surround. batch (default) crops every image in "
                             "the run to one shared box, so scale AND canvas stay identical "
                             "across legs. each crops every image to its own content. none keeps "
                             "the full render")
    parser.add_argument("--crop-margin", type=float, default=0.03,
                        help="margin left around the content when cropping, as a fraction of the "
                             "box's longer side")
    parser.add_argument("--include-last", action="store_true",
                        help="also render any frames_last.pdb named on the command line. They are "
                             "skipped by default -- the medoid is the figure set")
    parser.add_argument("--background", choices=tuple(BACKGROUND), default="transparent",
                        help="transparent (default) matches every PNG already committed here. "
                             "black or white bake the colour in")
    parser.add_argument("--no-labels", action="store_true",
                        help="drop any `label` lines a --style file brings (the built-in style has "
                             "none)")
    parser.add_argument("--name-suffix", default="",
                        help="appended to each PNG's name, for rendering a variant without "
                             "overwriting the standard one")
    parser.add_argument("--force", action="store_true",
                        help="re-render structures whose PNG already exists (default: skip)")
    parser.add_argument("--keep-pml", action="store_true",
                        help="keep the generated .pml beside the PNGs instead of in a temp file")
    parser.add_argument("--dry-run", action="store_true",
                        help="write the .pml and print the plan without calling pymol")
    args = parser.parse_args(argv)

    style, style_src = load_style(args)
    if args.no_labels:
        style = "\n".join(l for l in style.splitlines()
                          if not l.strip().startswith("label "))

    suffix = args.name_suffix

    jobs = []
    for path in args.structures:
        if not os.path.exists(path):
            sys.exit(f"no such structure: {path}")
        # The medoid is the figure set. `frames_last.pdb` sits beside every one of them and a
        # convenient glob picks it up, so it is skipped here rather than silently doubling the job.
        if os.path.basename(path) == "frames_last.pdb" and not args.include_last:
            print(f"skip (not a medoid)  {path}")
            continue
        out_dir = args.out_dir or default_out_dir(path)
        if not out_dir:
            sys.exit(f"{path} is not under a runs/<mol>/ tree, so pass --out-dir")
        png = os.path.join(out_dir, structure_label(path) + suffix + ".png")
        if os.path.exists(png) and not args.force:
            print(f"skip (exists)  {png}")
            continue
        jobs.append((os.path.abspath(path), out_dir, png))

    if not jobs:
        print("nothing to render")
        return 0
    for _, out_dir, _ in jobs:
        os.makedirs(out_dir, exist_ok=True)

    # One pymol process for the whole batch: startup and the first ray cost more than the renders
    # do. `delete all` between structures rather than `reinitialize`, which would also throw away
    # the settings and the label state that style.pml just established.
    lines = [f"# Generated by code/render_md_figures.py -- representation from {style_src}",
             "# No semicolons anywhere: PyMOL splits commands on one even inside a comment.",
             RENDER_SETUP]
    for spec in args.set_color:
        name, _, rgb = spec.partition("=")
        if not rgb:
            sys.exit(f"--set-color wants NAME=R,G,B, got {spec!r}")
        lines.append(f"set_color {name.strip()}, [{rgb.strip()}]")
    for path, _, png in jobs:
        obj = structure_label(path)
        lines += [
            "",
            f"# --- {obj} " + "-" * max(0, 70 - len(obj)),
            "delete all",
            f"load {path}, {obj}",
            style,
            # after the style block, so it wins over the style's own bg_color
            BACKGROUND[args.background],
        ]
        if args.hide_hydrogens:
            lines.append("hide sticks, hydro")
        lines += [
            f"orient {args.orient}",
            f"zoom {args.zoom}, {args.buffer}",
            f'print "RENDER {obj} atoms", cmd.count_atoms("all"), '
            f'"polymer", cmd.count_atoms("polymer"), "ligand", cmd.count_atoms("not polymer")',
            f"png {png}, width={args.width}, height={args.height}, dpi={args.dpi}, ray=1",
        ]
    script = "\n".join(lines) + "\n"

    pml_path = (os.path.join(jobs[0][1], "render_md_figures.pml") if args.keep_pml
                else os.path.join(os.path.dirname(jobs[0][2]), ".render_md_figures.pml"))
    with open(pml_path, "w") as fh:
        fh.write(script)
    print(f"{len(jobs)} structure(s), style from {style_src} -> {pml_path}")

    if args.dry_run:
        print("dry run, not calling pymol")
        return 0
    if not (os.path.exists(args.pymol) or shutil.which(args.pymol)):
        sys.exit(f"no pymol at {args.pymol} -- create the env or pass --pymol")

    # -cq: no GUI, no splash. Ray tracing is PyMOL's own software renderer and needs no display.
    proc = subprocess.run([args.pymol, "-cq", pml_path], capture_output=True, text=True)
    sys.stdout.write(proc.stdout)
    if proc.stderr.strip():
        sys.stderr.write(proc.stderr)
    if not args.keep_pml:
        os.remove(pml_path)
    if proc.returncode != 0:
        sys.exit(f"pymol exited {proc.returncode}")

    # pymol exits 0 even when a load silently failed, so the PNGs are the check, by size: a
    # successful ray of one of these structures is tens of KB, an empty scene a couple of hundred
    # bytes.
    missing = 0
    for _, _, png in jobs:
        if not os.path.exists(png):
            print(f"MISSING  {png}")
            missing += 1
    if missing:
        sys.exit(f"{missing} structure(s) produced no PNG")

    if args.crop != "none":
        crop_pngs([png for _, _, png in jobs], args.crop, args.crop_margin)

    for _, _, png in jobs:
        size = os.path.getsize(png)
        flag = "  <-- suspiciously small" if size < 5000 else ""
        print(f"{size:>9,} B  {png}{flag}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
