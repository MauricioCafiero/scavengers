"""Condense a whole path of placed fragments, not just a pair (step 4 of the plan).

A pair is easier than a chain, and the difference is not a detail. When two fragments are
condensed on their own, both ends are free: residue A's backbone can rotate about its CA-CB bond
(that is the parameter `t1`) and so can residue B's. In a chain, residue i's backbone frame is
already half-determined by the connection that arrived from residue i-1: its CA and CB are pinned
by the pose, and once the incoming N is placed, the outgoing C has nowhere left to go -- it is the
remaining tetrahedral slot, fixed by L chirality.

So the free parameters per segment drop from 2 + 2k to 1 + 2k, against the same four closure
conditions. A segment with k = 1 has 3 parameters against 4 conditions and is generically
unsolvable *in a chain* even where the same pair connects happily in isolation. **The pairwise
sweep is therefore an upper bound on what a chain can do, not a prediction of it**, and this
module is how that gap gets measured instead of assumed.

Usage:
    python code/condense_chain.py runs/octinoxate --path arginine:3 lysine:5 aspartic:21 --k 2 2
    python code/condense_chain.py runs/octinoxate --from-assign condense_pairs_geom.csv
"""
import argparse
import itertools
import os
import sys

import numpy as np
from scipy.optimize import minimize

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from condense import (A_C_N_CA, A_CA_C_N, A_CB_CA_C, A_CB_CA_N, A_N_CA_C, B_C_N,  # noqa: E402
                      B_C_O, B_CA_C, B_CA_HA, B_N_CA, B_N_H, OMEGA, build_obstacles,
                      closure_error, place_atom, selected_attachments, steric_penalty)
from condense_strain import (TETRA, A_CA_HA, EV_TO_KCAL, attached_atoms,  # noqa: E402
                             dihedral_angle, is_l_residue, kabsch_rmsd, planar_third, relax,
                             single_point, validate_structure, worst_contact)
from peptide_builder import ONE_LETTER, define_fragments, load_ligand, write_xyz  # noqa: E402


def outgoing_carbon(cg, cb, ca, n):
    """Residue i's carbonyl carbon, once its incoming nitrogen is known.

    CA and CB are pinned by the pose and N has just arrived, so C takes the remaining tetrahedral
    slot; L chirality decides which of the two it is. Nothing here is free, which is exactly why a
    chain is more constrained than a pair.
    """
    d0 = dihedral_angle(cg, cb, ca, n)
    c = place_atom(cg, cb, ca, B_CA_C, A_CB_CA_C, d0 + TETRA)
    if not is_l_residue(n, ca, c, cb):
        c = place_atom(cg, cb, ca, B_CA_C, A_CB_CA_C, d0 - TETRA)
    return c


def grow_from(ca, c, ref, k, params):
    """Grow from a residue whose CA and C are already placed, through k glycines, to the next N.

    `params` is [psi, (phi_i, psi_i) * k]; `ref` is the atom the outgoing psi is measured from,
    which is the residue's own CB. Returns the list of (label, xyz) starting at CA.
    """
    chain = [("CA", np.asarray(ca, dtype=float)), ("C", np.asarray(c, dtype=float))]
    n_next = place_atom(ref, ca, c, B_C_N, A_CA_C_N, params[0])
    chain.append(("N", n_next))
    ca_prev, c_prev = ca, c
    for i in range(k):
        phi, psi = params[1 + 2 * i], params[2 + 2 * i]
        ca_i = place_atom(ca_prev, c_prev, n_next, B_N_CA, A_C_N_CA, np.radians(OMEGA))
        c_i = place_atom(c_prev, n_next, ca_i, B_CA_C, A_N_CA_C, phi)
        n_i = place_atom(n_next, ca_i, c_i, B_C_N, A_CA_C_N, psi)
        chain += [("CA", ca_i), ("C", c_i), ("N", n_i)]
        ca_prev, c_prev, n_next = ca_i, c_i, n_i
    return chain


def n_params(ks):
    """a_0 for the first residue, then psi plus two dihedrals per glycine for each segment."""
    return 1 + sum(1 + 2 * k for k in ks)


def build_chain(sites, ks, params):
    """Walk the whole path. `sites` is one dict per fragment with cg, cb, ca.

    Returns (segments, closures): the grown atoms per segment, and the closure error at each
    junction. Each residue after the first inherits its carbonyl carbon from the incoming N, so
    only one dihedral per residue is free rather than two.
    """
    segments, closures = [], []
    cursor = 1
    c_i = place_atom(sites[0]["cg"], sites[0]["cb"], sites[0]["ca"],
                     B_CA_C, A_CB_CA_C, params[0])
    for i, k in enumerate(ks):
        take = 1 + 2 * k
        seg = grow_from(sites[i]["ca"], c_i, sites[i]["cb"], k, params[cursor:cursor + take])
        cursor += take
        segments.append(seg)
        nxt = sites[i + 1]
        closures.append(closure_error(seg, nxt["ca"], nxt["cb"]))
        if i + 1 < len(ks):
            c_i = outgoing_carbon(nxt["cg"], nxt["cb"], nxt["ca"], seg[-1][1])
    return segments, closures


def optimise_chain(sites, ks, restarts=200, seed=0, obstacles=None, w_steric=1.0,
                   w_closure=50.0, tol=0.02):
    """Minimise the total squared closure over every junction at once."""
    rng = np.random.default_rng(seed)
    npar = n_params(ks)

    def total(params, with_steric):
        segments, closures = build_chain(sites, ks, params)
        value = w_closure * float(np.sum(np.square(closures)))
        if with_steric and obstacles is not None:
            for seg, obs in zip(segments, obstacles):
                value += w_steric * steric_penalty(np.array([c for _, c in seg]), obs)
        return value

    best = (np.inf, None)
    extended = np.full(npar, np.pi)
    starts = [extended]
    for a0 in np.linspace(-np.pi, np.pi, 8, endpoint=False):
        start = extended.copy()
        start[0] = a0
        starts.append(start)
    starts += [rng.uniform(-np.pi, np.pi, npar) for _ in range(restarts)]

    for p0 in starts:
        res = minimize(lambda p: total(p, obstacles is not None), p0, method="L-BFGS-B")
        if res.fun < best[0]:
            best = (float(res.fun), res.x)
        if best[0] < w_closure * len(ks) * tol ** 2:
            break
    segments, closures = build_chain(sites, ks, best[1])
    return closures, best[1], segments


def build_chain_atoms(sites, ks, params, segments, drop_detached=False):
    """Assemble the whole condensed peptide for a path.

    Generalises the pair build: every fragment keeps its designed coordinates exactly, less the one
    hydrogen the backbone replaces, and the backbone runs through all of them. Residue 0 gets a
    neutral NH2 N-terminus and the last residue a primary amide, both arbitrary caps.

    Returns (symbols, coords, charge, fixed, groups) where `groups` locates each side chain and
    its attachment atoms, for the strain decomposition.
    """
    m = len(ks)
    symbols, coords, fixed, groups = [], [], [], []

    def add(sym, xyz, is_fixed=False):
        if is_fixed:
            fixed.append(len(symbols))
        symbols.append(sym)
        coords.append(np.asarray(xyz, dtype=float))

    charge = 0
    for i, site in enumerate(sites):
        ca = site["ca"]
        cb, cg = site["cb"], site["cg"]
        # nitrogen: the cap for residue 0, otherwise the one the previous segment delivered
        if i == 0:
            t1 = params[0]
            c_i = segments[0][1][1]
            n_i = place_atom(cg, cb, ca, B_N_CA, A_CB_CA_N, t1 + TETRA)
            ha_i = place_atom(cg, cb, ca, B_CA_HA, A_CA_HA, t1 - TETRA)
            if not is_l_residue(n_i, ca, c_i, cb):
                n_i = place_atom(cg, cb, ca, B_N_CA, A_CB_CA_N, t1 - TETRA)
                ha_i = place_atom(cg, cb, ca, B_CA_HA, A_CA_HA, t1 + TETRA)
        else:
            n_i = segments[i - 1][-1][1]
            d0 = dihedral_angle(cg, cb, ca, n_i)
            if i < m:
                # The carbonyl carbon was already placed by `outgoing_carbon`, which picked
                # whichever of the two tetrahedral slots L chirality required. HA has to go in the
                # *other* one, so find out which was taken rather than assuming: assuming it put
                # HA on top of the carbonyl carbon, 0.44 A away, whenever the minus slot was used.
                c_i = segments[i][1][1]
                plus = place_atom(cg, cb, ca, B_CA_C, A_CB_CA_C, d0 + TETRA)
                minus = place_atom(cg, cb, ca, B_CA_C, A_CB_CA_C, d0 - TETRA)
                took_plus = (np.linalg.norm(c_i - plus) < np.linalg.norm(c_i - minus))
                ha_i = place_atom(cg, cb, ca, B_CA_HA, A_CA_HA,
                                  d0 - TETRA if took_plus else d0 + TETRA)
            else:
                c_i = place_atom(cg, cb, ca, B_CA_C, A_CB_CA_C, d0 + TETRA)
                ha_i = place_atom(cg, cb, ca, B_CA_HA, A_CA_HA, d0 - TETRA)
                if not is_l_residue(n_i, ca, c_i, cb):
                    c_i = place_atom(cg, cb, ca, B_CA_C, A_CB_CA_C, d0 - TETRA)
                    ha_i = place_atom(cg, cb, ca, B_CA_HA, A_CA_HA, d0 + TETRA)

        add("N", n_i)
        if i == 0:                       # neutral NH2 N-terminus
            add("H", place_atom(cb, ca, n_i, B_N_H, 109.5, np.radians(60.0)))
            add("H", place_atom(cb, ca, n_i, B_N_H, 109.5, np.radians(180.0)))
        else:
            add("H", planar_third(n_i, segments[i - 1][-2][1], ca, B_N_H))
        i_ca = len(symbols)
        add("C", ca)
        add("H", ha_i)

        keep = attached_atoms(site) if drop_detached else set(range(len(site["frag_symbols"])))
        group, i_cb = [], None
        for j, (sym, xyz) in enumerate(zip(site["frag_symbols"], site["frag_coords"])):
            if j != site["h"] and j in keep:
                if j == site["site"]["c"]:
                    i_cb = len(symbols)
                group.append(len(symbols))
                add(sym, xyz, is_fixed=True)
        groups.append({"atoms": group, "cb": i_cb, "ca": i_ca, "charge": site["charge"]})
        charge += site["charge"]

        add("C", c_i)
        if i < m:
            add("O", planar_third(c_i, ca, segments[i][2][1], B_C_O))
            # the glycine linkers of segment i
            seg = segments[i]
            prev_c = c_i
            for j in range(ks[i]):
                n_g, ca_g, c_g = seg[2 + 3 * j][1], seg[3 + 3 * j][1], seg[4 + 3 * j][1]
                n_next = seg[5 + 3 * j][1]
                add("N", n_g)
                add("H", planar_third(n_g, prev_c, ca_g, B_N_H))
                add("C", ca_g)
                phi = params_of_segment(params, ks, i)[1 + 2 * j]
                add("H", place_atom(prev_c, n_g, ca_g, B_CA_HA, A_CA_HA, phi + TETRA))
                add("H", place_atom(prev_c, n_g, ca_g, B_CA_HA, A_CA_HA, phi - TETRA))
                add("C", c_g)
                add("O", planar_third(c_g, ca_g, n_next, B_C_O))
                prev_c = c_g
        else:                            # C-terminal primary amide
            n_cap = place_atom(n_i, ca, c_i, B_C_N, A_CA_C_N, np.pi)
            o_c = planar_third(c_i, ca, n_cap, B_C_O)
            add("O", o_c)
            add("N", n_cap)
            add("H", place_atom(o_c, c_i, n_cap, B_N_H, 120.0, 0.0))
            add("H", place_atom(o_c, c_i, n_cap, B_N_H, 120.0, np.pi))

    return symbols, np.array(coords), charge, fixed, groups


def params_of_segment(params, ks, i):
    """The slice of the parameter vector belonging to segment i."""
    cursor = 1
    for j in range(i):
        cursor += 1 + 2 * ks[j]
    return params[cursor:cursor + 1 + 2 * ks[i]]


def site_for(att, h):
    return {"cg": att["cg"], "cb": att["cb"], "ca": dict(att["ca"])[h], "name": att["name"],
            "pose": att["pose"], "site": att["site"], "frag_symbols": att["frag_symbols"],
            "frag_coords": att["frag_coords"], "h": h, "charge": att["charge"]}


def sidechains_only(symbols, xyz, groups):
    """Every side chain on its own, each recapped with the hydrogen the backbone displaced.

    Same idea as the pair version in `condense_strain`: it isolates the side-chain interactions the
    designed arrangement forbids, so they can be subtracted from the total and not be mistaken for
    backbone strain.
    """
    syms, pos = [], []
    for g in groups:
        for gi in g["atoms"]:
            syms.append(symbols[gi])
            pos.append(xyz[gi])
        cb, ca = np.asarray(xyz[g["cb"]]), np.asarray(xyz[g["ca"]])
        direction = ca - cb
        syms.append("H")
        pos.append(cb + 1.09 * direction / np.linalg.norm(direction))
    return syms, np.array(pos)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("outdir")
    parser.add_argument("--ligand", default="octinoxate")
    parser.add_argument("--path", nargs="+", required=True,
                        help="fragment:pose in order, e.g. arginine:3 lysine:5")
    parser.add_argument("--k", nargs="+", type=int, required=True,
                        help="linkers between consecutive path entries")
    parser.add_argument("--restarts", type=int, default=200)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--tol", type=float, default=0.2)
    parser.add_argument("--steric", action="store_true")
    parser.add_argument("--save", help="write the chain result to this JSON file")
    parser.add_argument("--max-combos", type=int, default=40,
                        help="cap on CA-candidate combinations tried (sampled beyond this)")
    parser.add_argument("--build", action="store_true", help="assemble the full-atom peptide")
    parser.add_argument("--relax", action="store_true",
                        help="also relax it with UMA and report the strain")
    parser.add_argument("--with-ligand", action="store_true",
                        help="hold the ligand fixed during the constrained relaxation")
    parser.add_argument("--drop-detached", action="store_true")
    parser.add_argument("--steps", type=int, default=400)
    parser.add_argument("--model", default="uma-s-1p2p1")
    args = parser.parse_args(argv)

    if len(args.k) != len(args.path) - 1:
        parser.error(f"need {len(args.path) - 1} linker counts for {len(args.path)} fragments")

    frags = define_fragments()
    ligand = load_ligand(args.ligand)
    atts = {(a["name"], a["pose"]): a for a in selected_attachments(args.outdir, frags, ligand)}
    chosen = []
    for entry in args.path:
        name, pose = entry.split(":")
        chosen.append(atts[(name, int(pose))])

    lig_xyz = np.asarray(ligand["coords"])
    lig_heavy = lig_xyz[[i for i, s in enumerate(ligand["atoms"]) if s != "H"]]

    sequence = ONE_LETTER[chosen[0]["name"]]
    for att, k in zip(chosen[1:], args.k):
        sequence += "G" * k + ONE_LETTER[att["name"]]
    print(f"path: {' -> '.join(e for e in args.path)}")
    print(f"linkers: {args.k}   sequence: {sequence}\n")

    # Every fragment offers ~3 places to graft its CA, and the design is free to choose, so the
    # combinations multiply as 3^n along the path. Exhaustive is right for short paths and
    # hopeless by six fragments (729), so past --max-combos this samples instead, always keeping
    # the all-first choice so the sampling can only improve on a fixed reference.
    candidate_sets = [[h for h, _ in att["ca"]] for att in chosen]
    all_combos = list(itertools.product(*candidate_sets))
    if len(all_combos) > args.max_combos:
        rng = np.random.default_rng(args.seed)
        first = all_combos[0]
        idx = rng.choice(len(all_combos), size=args.max_combos - 1, replace=False)
        combos = [first] + [all_combos[i] for i in sorted(idx)]
        print(f"{len(all_combos)} CA combinations; sampling {len(combos)}\n")
    else:
        combos = all_combos
    best = None
    for combo in combos:
        sites = [site_for(att, h) for att, h in zip(chosen, combo)]
        obstacles = None
        if args.steric:
            obstacles = [build_obstacles(sites[i], sites[i]["ca"], sites[i + 1],
                                         sites[i + 1]["ca"], k, lig_heavy)
                         for i, k in enumerate(args.k)]
        closures, params, _ = optimise_chain(sites, args.k, restarts=args.restarts,
                                             seed=args.seed, obstacles=obstacles)
        worst = float(np.max(closures))
        if best is None or worst < best[0]:
            best = (worst, combo, closures, params)
        print(f"  CA choice {combo}: worst junction {worst:.3f} A   "
              + " ".join(f"{c:.3f}" for c in closures), flush=True)

    worst, combo, closures, params = best
    print(f"\nbest CA choice {combo}: worst junction {worst:.3f} A")
    ok = worst < args.tol
    print(f"chain {'CLOSES' if ok else 'does NOT close'} at tolerance {args.tol} A")
    if not ok:
        bad = int(np.argmax(closures))
        print(f"  the junction that fails is {args.path[bad]} -> {args.path[bad + 1]} "
              f"with {args.k[bad]} linker(s), closure {closures[bad]:.3f} A")

    if args.save:
        import json
        with open(args.save, "w") as fh:
            json.dump({"path": args.path, "linkers": args.k, "sequence": sequence,
                       "ca_choice": list(combo), "closures": [round(float(c), 4) for c in closures],
                       "worst_junction": round(float(worst), 4),
                       "closes_at_tol": bool(worst < args.tol), "tol": args.tol}, fh, indent=2)
        print(f"saved chain result to {args.save}")

    if not (args.build or args.relax):
        return 0

    sites = [site_for(att, h) for att, h in zip(chosen, combo)]
    segments, _ = build_chain(sites, args.k, params)
    symbols, coords, charge, fixed, groups = build_chain_atoms(
        sites, args.k, params, segments, drop_detached=args.drop_detached)
    out_dir = os.path.join(args.outdir, "condense")
    os.makedirs(out_dir, exist_ok=True)
    tag = "chain_" + "_".join(f"{e.replace(':', '')}" for e in args.path) + \
          "_k" + "".join(str(x) for x in args.k)
    write_xyz(os.path.join(out_dir, tag + "_built.xyz"), symbols, coords,
              f"{sequence}, worst junction {worst:.3f} A, charge {charge}")
    print(f"\nbuilt {len(symbols)} atoms, charge {charge}, {len(fixed)} held "
          f"-> {tag}_built.xyz")
    complaints = validate_structure(symbols, coords)
    for complaint in complaints:
        print(f"  STRUCTURE: {complaint}")
    if not complaints:
        print("  structure audit: clean")
    if not args.relax:
        return 0

    import torch
    from fairchem.core import FAIRChemCalculator, pretrained_mlip
    device = "cuda" if torch.cuda.is_available() else "cpu"
    calculator = FAIRChemCalculator(pretrained_mlip.get_predict_unit(args.model, device=device),
                                    task_name="omol")

    if args.with_ligand:
        all_symbols = symbols + list(ligand["atoms"])
        all_coords = np.vstack([coords, lig_xyz])
        held_fixed = list(fixed) + list(range(len(symbols), len(all_symbols)))
        _, xyz_all, st_h, cv_h = relax(all_symbols, all_coords,
                                       charge + int(ligand.get("charge", 0)), calculator,
                                       fixed=held_fixed, steps=args.steps)
        xyz_held = xyz_all[:len(symbols)]
        e_held = single_point(symbols, xyz_held, charge, calculator)
    else:
        e_held, xyz_held, st_h, cv_h = relax(symbols, coords, charge, calculator,
                                             fixed=fixed, steps=args.steps)
    write_xyz(os.path.join(out_dir, tag + "_held.xyz"), symbols, xyz_held,
              f"side chains held, E = {e_held:.6f} eV")
    e_free, xyz_free, st_f, cv_f = relax(symbols, xyz_held, charge, calculator, steps=args.steps)
    write_xyz(os.path.join(out_dir, tag + "_free.xyz"), symbols, xyz_free,
              f"free relaxation, E = {e_free:.6f} eV")

    sc_charge = sum(g["charge"] for g in groups)
    sc_held = single_point(*sidechains_only(symbols, xyz_held, groups), sc_charge, calculator)
    sc_free = single_point(*sidechains_only(symbols, xyz_free, groups), sc_charge, calculator)
    strain = (e_held - e_free) * EV_TO_KCAL
    sidechain = (sc_held - sc_free) * EV_TO_KCAL
    backbone_atoms = [i for i in range(len(symbols))
                      if i not in set(fixed) and symbols[i] != "H"]
    lig_contact = float(np.linalg.norm(
        xyz_held[backbone_atoms][:, None, :] - lig_heavy[None, :, :], axis=-1).min())
    print(f"\nstrain {strain:8.2f} = side chains {sidechain:7.2f} + backbone "
          f"{strain - sidechain:7.2f} kcal/mol")
    print(f"  rmsd(held,free) {kabsch_rmsd(xyz_held, xyz_free):.2f} A   "
          f"worst contact {worst_contact(symbols, xyz_held):.2f} A   "
          f"ligand {lig_contact:.2f} A")
    print(f"  held {st_h} steps conv={cv_h}, free {st_f} steps conv={cv_f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
