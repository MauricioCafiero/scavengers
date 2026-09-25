"""Fragment condensation: grow a real backbone between two placed fragment poses (method A).

The pipeline's `build_sequence` chains fragment *centroids* and turns gaps into glycine with
`round(d / 3.8) - 1`, so a design can demand a backbone path no peptide can follow. This module
asks the question the other way round: hold two fragments exactly where the search put them, graft
a CA onto each capping carbon (`attach_geom`), and try to grow an ideal trans backbone from one CA
to the other through k glycine spacers. What comes out is not a score but a measurement.

Stage 1, here, is geometric and needs no potential. The backbone is built by NeRF from ideal bond
lengths and angles with omega fixed at 180 degrees, leaving as free parameters

    t1      rotation of residue A's backbone tripod about the CB-CA bond
    t2      the CA_A-C_A dihedral (psi of residue A, shifted by the CB reference)
    phi_i, psi_i   for each of the k spacer glycines

so 2 + 2k parameters. The closure condition is that residue B's amide nitrogen land on the circle
of positions compatible with B's fixed CA and CB -- a 111 degree cone of radius 1.458 A -- which is
two conditions, the remaining freedom on that circle being residue B's own phi. k = 0 is therefore
exactly determined and often has no solution; k >= 1 leaves a solution manifold. The residual, in
angstroms, is the reachability criterion: measured, not assumed.

Usage:
    python code/condense.py sweep runs/octinoxate                 # ordered pairs, k = 0..4
    python code/condense.py pair runs/octinoxate arginine 3 lysine 5 --k 2
"""
import argparse
import csv
import itertools
import os
import sys

import numpy as np
from scipy.optimize import minimize

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from attach_geom import pose_attachment  # noqa: E402
from peptide_builder import define_fragments, load_ligand, load_poses  # noqa: E402

# ideal trans peptide geometry (angstroms, degrees)
B_N_CA, B_CA_C, B_C_N, B_C_O, B_CA_CB = 1.458, 1.525, 1.329, 1.231, 1.530
B_N_H, B_CA_HA = 1.010, 1.090
A_N_CA_C, A_CA_C_N, A_C_N_CA = 111.2, 116.2, 121.7
A_CA_C_O, A_CB_CA_C, A_CB_CA_N, A_C_N_H = 120.8, 110.5, 110.5, 119.0
OMEGA = 180.0


def cross3(u, v):
    """Cross product of two 3-vectors. np.cross routes through moveaxis and dominates the
    closure optimisation's runtime; this does not."""
    return np.array([u[1] * v[2] - u[2] * v[1],
                     u[2] * v[0] - u[0] * v[2],
                     u[0] * v[1] - u[1] * v[0]])


def unit(v):
    """v normalised, without np.linalg.norm's overhead."""
    return v / np.sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2])


def place_atom(a, b, c, bond, angle_deg, dihedral_rad):
    """NeRF: the position bonded to c, at `angle_deg` to b-c, with dihedral a-b-c-d."""
    theta = np.radians(angle_deg)
    bc = unit(c - b)
    n = unit(cross3(b - a, bc))
    m = cross3(n, bc)
    sin_t = np.sin(theta)
    return c + bond * (-np.cos(theta) * bc
                       + sin_t * np.cos(dihedral_rad) * m
                       + sin_t * np.sin(dihedral_rad) * n)


def grow_backbone(cg_a, cb_a, ca_a, k, params):
    """Backbone atoms from residue A's CA up to residue B's amide nitrogen.

    Returns (chain, n_b). `chain` is a list of (label, xyz) in bonded order, starting at the
    newly grafted CA of residue A: CA_A, C_A, N_1, then (CA_i, C_i, N_i+1) per spacer, so
    3 + 3k atoms. Every consecutive pair is bonded, which is what the steric term relies on.
    """
    t1, t2 = params[0], params[1]
    chain = [("CA", np.asarray(ca_a, dtype=float))]
    c_prev = place_atom(cg_a, cb_a, ca_a, B_CA_C, A_CB_CA_C, t1)     # C of residue A
    chain.append(("C", c_prev))
    ca_prev = ca_a
    n_next = place_atom(cb_a, ca_a, c_prev, B_C_N, A_CA_C_N, t2)     # N of the next residue
    chain.append(("N", n_next))
    for i in range(k):
        phi, psi = params[2 + 2 * i], params[3 + 2 * i]
        ca_i = place_atom(ca_prev, c_prev, n_next, B_N_CA, A_C_N_CA, np.radians(OMEGA))
        c_i = place_atom(c_prev, n_next, ca_i, B_CA_C, A_N_CA_C, phi)
        n_i = place_atom(n_next, ca_i, c_i, B_C_N, A_CA_C_N, psi)
        chain += [("CA", ca_i), ("C", c_i), ("N", n_i)]
        ca_prev, c_prev, n_next = ca_i, c_i, n_i
    return chain, chain[-1][1]


def closure_error(chain, ca_b, cb_b):
    """How far the grown chain is from joining residue B properly, in angstroms.

    Residue B's CA is fixed by the pose, so the chain has to arrive at it: the trans peptide
    unit places a CA from the last three backbone atoms, and that position must be residue B's,
    with the N-CA-CB angle the residue's own geometry requires. Both parts go to zero together,
    and the angular part is reported as the arc it corresponds to at the N-CA bond length, so the
    total is a distance.

    Four conditions against 2 + 2k free dihedrals: k = 0 is generically unsolvable, k = 1 has
    isolated solutions, k >= 2 leaves a manifold.
    """
    ca_k, c_k, n_b = chain[-3][1], chain[-2][1], chain[-1][1]
    ca_pred = place_atom(ca_k, c_k, n_b, B_N_CA, A_C_N_CA, np.radians(OMEGA))
    diff = ca_pred - ca_b
    d_pos = float(np.sqrt(diff @ diff))
    v1, v2 = n_b - ca_pred, cb_b - ca_pred
    cos = float(v1 @ v2) / float(np.sqrt((v1 @ v1) * (v2 @ v2)))
    d_ang = B_N_CA * abs(np.arccos(np.clip(cos, -1.0, 1.0)) - np.radians(A_CB_CA_N))
    return float(np.hypot(d_pos, d_ang))


def closure_starts(n_par, restarts, rng, grid=6):
    """Starting dihedrals, structured first then random.

    t1 and t2 aim the whole chain, so a grid over those two with the spacers fully extended
    (all trans, the longest reach a backbone has) finds the far solutions that random starts in
    a 10-dimensional space regularly miss. Random starts then cover the folded-back ones.
    """
    extended = np.full(n_par, np.pi)
    for t1 in np.linspace(-np.pi, np.pi, grid, endpoint=False):
        for t2 in np.linspace(-np.pi, np.pi, grid, endpoint=False):
            start = extended.copy()
            start[0], start[1] = t1, t2
            yield start
    for _ in range(restarts):
        yield rng.uniform(-np.pi, np.pi, n_par)


CLASH_CUTOFF = 2.9      # heavy-atom contact the new backbone may not come inside of
SELF_CUTOFF = 3.0       # same, between backbone atoms five or more bonds apart
SELF_SKIP = 4


def build_obstacles(att_a, ca_a, att_b, ca_b, k, ligand_coords=None, cutoff=CLASH_CUTOFF):
    """Heavy atoms the grown backbone must stay clear of, with the pairs to ignore.

    The chain's own neighbourhood is legitimately close to the attachment: C_A is 1-3 to CB_A,
    N_B is bonded to CA_B, and so on, so each obstacle carries the chain positions it does not
    apply to. Chain order is [CA_A, C_A, N_1, (CA_i, C_i, N_i+1) ...], length 3 + 3k.
    """
    n_chain = 3 + 3 * k
    coords, skips = [], []

    def obstacle(xyz, skip):
        coords.append(np.asarray(xyz, dtype=float))
        skips.append(set(skip))

    obstacle(att_a["cb"], {0, 1, 2})             # bonded to CA_A, 1-3 to C_A, 1-4 to N_1
    obstacle(att_a["cg"], {0, 1})
    obstacle(ca_b, {n_chain - 1, n_chain - 2})   # bonded to N_B, 1-3 to the C before it
    obstacle(att_b["cb"], {n_chain - 1, n_chain - 2})
    obstacle(att_b["cg"], {n_chain - 1})
    skip_atoms = {"a": {att_a["site"]["c"], att_a["site"]["heavy"]},
                  "b": {att_b["site"]["c"], att_b["site"]["heavy"]}}
    for key, att in (("a", att_a), ("b", att_b)):
        for i, (sym, xyz) in enumerate(zip(att["frag_symbols"], att["frag_coords"])):
            if sym != "H" and i not in skip_atoms[key]:
                obstacle(xyz, set())
    if ligand_coords is not None:
        for xyz in ligand_coords:
            obstacle(xyz, set())

    coords = np.array(coords)
    mask = np.zeros((n_chain, len(coords)), dtype=bool)
    for j, skip in enumerate(skips):
        for i in skip:
            mask[i, j] = True
    return {"coords": coords, "mask": mask, "cutoff": cutoff}


def steric_penalty(chain_xyz, obstacles, self_cutoff=SELF_CUTOFF, skip=SELF_SKIP):
    """Sum of squared overlaps, in A^2: backbone against the fixed atoms, and against itself."""
    d = np.linalg.norm(chain_xyz[:, None, :] - obstacles["coords"][None, :, :], axis=-1)
    over = np.clip(obstacles["cutoff"] - d, 0.0, None)
    over[obstacles["mask"]] = 0.0
    total = float((over ** 2).sum())
    n = len(chain_xyz)
    ds = np.linalg.norm(chain_xyz[:, None, :] - chain_xyz[None, :, :], axis=-1)
    sep = np.abs(np.arange(n)[:, None] - np.arange(n)[None, :])
    over_self = np.clip(self_cutoff - ds, 0.0, None)
    over_self[sep < skip] = 0.0
    return total + 0.5 * float((over_self ** 2).sum())


def optimize_closure(cg_a, cb_a, ca_a, ca_b, cb_b, k, restarts=24, seed=0, tol=0.02,
                     obstacles=None, w_steric=1.0, w_closure=50.0, n_refine=6):
    """Minimise the closure error over the free dihedrals, optionally avoiding obstacles.

    Without `obstacles` this is pure reachability: can an ideal backbone bridge the two
    attachment points at all, sterics ignored. With them the chain also has to stay clear of
    both fragments, of the ligand and of itself, which is the criterion a design can use. The
    steric pass starts from the best few geometric solutions rather than from scratch, because
    the steric penalty is the expensive part of the objective and the geometric solutions are
    where the clash-free ones live. Closure carries much the larger weight: the molecule has to
    be bonded, whereas a clash is a cost. A case with a low closure and a penalty left over is
    one where no clash-free connection exists, which is the answer, not a failure.

    Returns (closure error in A, params, chain, penalty in A^2).
    """
    rng = np.random.default_rng(seed)
    n_par = 2 + 2 * k

    def geometry(p):
        chain, _ = grow_backbone(cg_a, cb_a, ca_a, k, p)
        return closure_error(chain, ca_b, cb_b) ** 2

    solutions = []
    for p0 in closure_starts(n_par, restarts, rng):
        res = minimize(geometry, p0, method="L-BFGS-B")
        solutions.append((float(res.fun), res.x))
        if obstacles is None and res.fun < tol ** 2:
            break
    solutions.sort(key=lambda s: s[0])

    def result(params):
        chain, _ = grow_backbone(cg_a, cb_a, ca_a, k, params)
        xyz = np.array([c for _, c in chain])
        pen = 0.0 if obstacles is None else steric_penalty(xyz, obstacles)
        return closure_error(chain, ca_b, cb_b), params, chain, pen

    if obstacles is None:
        return result(solutions[0][1])

    def combined(p):
        chain, _ = grow_backbone(cg_a, cb_a, ca_a, k, p)
        xyz = np.array([c for _, c in chain])
        return (w_closure * closure_error(chain, ca_b, cb_b) ** 2
                + w_steric * steric_penalty(xyz, obstacles))

    best = (np.inf, solutions[0][1])
    for _, p0 in solutions[:n_refine]:
        res = minimize(combined, p0, method="L-BFGS-B")
        if res.fun < best[0]:
            best = (float(res.fun), res.x)
    return result(best[1])


def angle_between(v1, v2):
    """Angle in degrees between two vectors."""
    c = float(v1 @ v2) / float(np.linalg.norm(v1) * np.linalg.norm(v2))
    return float(np.degrees(np.arccos(np.clip(c, -1.0, 1.0))))


def min_distance(points, others):
    """Closest approach between two coordinate sets, inf when either is empty."""
    if len(points) == 0 or len(others) == 0:
        return np.inf
    d = np.linalg.norm(np.asarray(points)[:, None, :] - np.asarray(others)[None, :, :], axis=-1)
    return float(d.min())


def pair_scan(att_a, att_b, k_values, fixed_coords, restarts=24, seed=0, steric=False):
    """Every CA-candidate combination and spacer count for one ordered pair of poses.

    `fixed_coords` are the atoms the new backbone must not run into (the two fragments and the
    ligand), given as {"fragments": array, "ligand": array}. With `steric` the optimisation
    itself avoids them; without it they are only measured afterwards.
    """
    rows = []
    cos_att = float(att_a["attach_vector"] @ att_b["attach_vector"])
    d_cb = float(np.linalg.norm(att_a["cb"] - att_b["cb"]))
    # what a clash against the fragments means: their heavy atoms other than the attachment
    # carbon and its neighbour, which the new backbone is bonded to and legitimately close to
    probe = []
    for att in (att_a, att_b):
        skip = {att["site"]["c"], att["site"]["heavy"]}
        probe += [xyz for i, (sym, xyz) in enumerate(zip(att["frag_symbols"], att["frag_coords"]))
                  if sym != "H" and i not in skip]
    probe = np.array(probe)
    for (h_a, ca_a), (h_b, ca_b) in itertools.product(att_a["ca"], att_b["ca"]):
        d_ca = float(np.linalg.norm(ca_a - ca_b))
        # how far behind its own fragment each CA's partner sits: the backbone leaves CA away
        # from CB, so a small angle here means the chain has to turn back on itself
        theta_a = angle_between(ca_b - ca_a, att_a["cb"] - ca_a)
        theta_b = angle_between(ca_a - ca_b, att_b["cb"] - ca_b)
        for k in k_values:
            obstacles = (build_obstacles(att_a, ca_a, att_b, ca_b, k,
                                        fixed_coords["ligand_heavy"])
                         if steric else None)
            err, params, chain, penalty = optimize_closure(
                att_a["cg"], att_a["cb"], ca_a, ca_b, att_b["cb"], k,
                restarts=restarts, seed=seed, obstacles=obstacles)
            built = np.array([xyz for _, xyz in chain])
            rows.append({
                "frag_a": att_a["name"], "frag_b": att_b["name"], "h_a": h_a, "h_b": h_b,
                "k": k, "d_cb": round(d_cb, 3), "d_ca": round(d_ca, 3),
                "cos_att": round(cos_att, 3), "theta_a": round(theta_a, 1),
                "theta_b": round(theta_b, 1), "closure": round(float(err), 4),
                "penalty": round(float(penalty), 4),
                "clash_frag": round(min_distance(built, probe), 3),
                "clash_lig": round(min_distance(built, fixed_coords["ligand_heavy"]), 3),
                "dihedrals": " ".join(f"{np.degrees(x) % 360:.0f}" for x in params),
            })
    return rows


def selected_attachments(outdir, frags, ligand, poses=None):
    """Attachment geometry for a set of placed poses.

    By default the poses flagged `selected` in energies.csv. `poses` overrides that with an
    explicit list of (fragment name, pose index), which is how a particular shell out of the
    starts x copies family gets analysed rather than just the default-order one.
    """
    new_molecules, _ = load_poses(outdir, frags, ligand)
    by_name = {f["name"]: i for i, f in enumerate(frags)}
    ie_of = {}
    with open(os.path.join(outdir, "energies.csv")) as f:
        rows = list(csv.DictReader(f))
    for row in rows:
        ie_of[(row["fragment"], int(row["pose"]))] = float(row["ie_kcal_mol"])
    if poses is None:
        chosen = [(r["fragment"], int(r["pose"]), float(r["ie_kcal_mol"]))
                  for r in rows if r.get("selected") == "1"]
    else:
        chosen = [(name, pose, ie_of[(name, pose)]) for name, pose in poses]
    out = []
    for name, pose, ie in chosen:
        frag = frags[by_name[name]]
        coords = np.asarray(new_molecules[by_name[name]][pose])[ligand["num_atoms"]:]
        att = pose_attachment(frag, coords)
        att["pose"] = pose
        att["ie"] = ie
        att["frag_coords"] = coords
        att["frag_symbols"] = list(frag["atoms"])
        att["charge"] = frag["charge"]
        out.append(att)
    return out


def parse_poses(spec):
    """frag:pose,frag:pose,... -> [(name, index), ...]; None stays None."""
    if not spec:
        return None
    out = []
    for item in spec.split(","):
        name, pose = item.strip().split(":")
        out.append((name, int(pose)))
    return out


def cmd_sweep(args):
    frags = define_fragments()
    ligand = load_ligand(args.ligand)
    atts = selected_attachments(args.outdir, frags, ligand, parse_poses(args.poses))
    print(f"{len(atts)} poses: " + ", ".join(f"{a['name']}{a['pose']}" for a in atts))
    lig_coords = np.asarray(ligand["coords"])
    lig_heavy = lig_coords[[i for i, sym in enumerate(ligand["atoms"]) if sym != "H"]]
    k_values = list(range(args.kmin, args.kmax + 1))

    default_name = "condense_pairs_steric.csv" if args.steric else "condense_pairs.csv"
    out_path = os.path.join(args.outdir, args.out or default_name)
    fields = ["frag_a", "pose_a", "frag_b", "pose_b", "h_a", "h_b", "k", "d_cb", "d_ca",
              "cos_att", "theta_a", "theta_b", "closure", "penalty", "clash_frag", "clash_lig",
              "dihedrals"]
    n_written = 0
    with open(out_path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields)
        writer.writeheader()
        for a, b in itertools.permutations(range(len(atts)), 2):
            att_a, att_b = atts[a], atts[b]
            others = np.vstack([att_a["frag_coords"], att_b["frag_coords"]])
            rows = pair_scan(att_a, att_b, k_values,
                             {"fragments": others, "ligand": lig_coords,
                              "ligand_heavy": lig_heavy},
                             restarts=args.restarts, seed=args.seed, steric=args.steric)
            for row in rows:
                row["pose_a"], row["pose_b"] = att_a["pose"], att_b["pose"]
                writer.writerow(row)
            fh.flush()          # so a partial sweep is readable while it runs
            n_written += len(rows)
            best = min(rows, key=lambda r: (r["closure"], r["penalty"]))
            print(f"{att_a['name']}{att_a['pose']} -> {att_b['name']}{att_b['pose']}  "
                  f"d_CA {best['d_ca']:.2f}  best closure {best['closure']:.3f} A "
                  f"(penalty {best['penalty']:.2f}) at k={best['k']}", flush=True)
    print(f"\n{n_written} rows -> {out_path}")
    return 0


def cmd_pair(args):
    frags = define_fragments()
    ligand = load_ligand(args.ligand)
    atts = {(a["name"], a["pose"]): a
            for a in selected_attachments(args.outdir, frags, ligand, parse_poses(args.poses))}
    att_a, att_b = atts[(args.frag_a, args.pose_a)], atts[(args.frag_b, args.pose_b)]
    others = np.vstack([att_a["frag_coords"], att_b["frag_coords"]])
    lig_coords = np.asarray(ligand["coords"])
    lig_heavy = lig_coords[[i for i, sym in enumerate(ligand["atoms"]) if sym != "H"]]
    rows = pair_scan(att_a, att_b, [args.k],
                     {"fragments": others, "ligand": lig_coords, "ligand_heavy": lig_heavy},
                     restarts=args.restarts, seed=args.seed, steric=args.steric)
    for row in sorted(rows, key=lambda r: r["closure"]):
        print(f"h_a {row['h_a']:>3} h_b {row['h_b']:>3}  d_CA {row['d_ca']:5.2f}  "
              f"theta_a {row['theta_a']:5.1f} theta_b {row['theta_b']:5.1f}  "
              f"closure {row['closure']:.3f}  penalty {row['penalty']:6.3f}  "
              f"clash frag {row['clash_frag']:.2f} lig {row['clash_lig']:.2f}")
    return 0


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="cmd", required=True)

    common = argparse.ArgumentParser(add_help=False)
    common.add_argument("outdir", help="run directory with poses/ and energies.csv")
    common.add_argument("--ligand", default="octinoxate")
    common.add_argument("--restarts", type=int, default=24)
    common.add_argument("--seed", type=int, default=0)
    common.add_argument("--poses", help="explicit shell as frag:pose,frag:pose,... "
                        "(default: the poses flagged selected in energies.csv)")
    common.add_argument("--steric", action="store_true",
                        help="make the backbone avoid both fragments, the ligand and itself")

    p_sweep = sub.add_parser("sweep", parents=[common], help="all ordered pairs of selected poses")
    p_sweep.add_argument("--kmax", type=int, default=4, help="largest spacer count (default 4)")
    p_sweep.add_argument("--kmin", type=int, default=0, help="smallest spacer count (default 0)")
    p_sweep.add_argument("--out", help="output CSV name inside outdir")
    p_sweep.set_defaults(func=cmd_sweep)

    p_pair = sub.add_parser("pair", parents=[common], help="one ordered pair in detail")
    p_pair.add_argument("frag_a")
    p_pair.add_argument("pose_a", type=int)
    p_pair.add_argument("frag_b")
    p_pair.add_argument("pose_b", type=int)
    p_pair.add_argument("--k", type=int, default=1)
    p_pair.set_defaults(func=cmd_pair)

    args = parser.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
