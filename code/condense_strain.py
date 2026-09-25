"""Strain energy of a condensed fragment pair: the energetic half of method A.

`condense.py` measures whether a backbone can *reach* from one placed fragment to another. This
module measures what it costs. It builds the whole covalent molecule -- both side chains exactly
where the search put them, a real backbone between them, neutral caps at both ends -- and relaxes
it twice with UMA:

    E_held   fragment atoms fixed at their designed positions, backbone free
    E_free   nothing fixed

    strain = E_held - E_free      (kcal/mol, >= 0)

Low strain means the shell was realisable and the structure is now in hand. High strain means the
design was asking for something a peptide cannot do, and the pair that carries the strain is named.
The ligand is deliberately absent: it would add binding energy and confuse the readout, and a
backbone that runs through the ligand is already reported separately as `clash_lig`.

Residue chirality is L, checked against the Boltz-predicted peptides
(det[N-CA, C-CA, CB-CA] > 0 for every L residue in them).

Usage:
    python code/condense_strain.py runs/octinoxate --cases 12
    python code/condense_strain.py runs/octinoxate --pair arginine 3 glutamic 22 --k 1
"""
import argparse
import csv
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from condense import (A_CA_C_N, A_CB_CA_C, A_CB_CA_N, B_C_N, B_C_O, B_CA_C,  # noqa: E402
                      B_CA_HA, B_N_CA, B_N_H, angle_between, build_obstacles, closure_error,
                      grow_backbone, optimize_closure, place_atom, selected_attachments,
                      steric_penalty)
from attach_geom import connected_components, electron_count, perceive_bonds  # noqa: E402
from peptide_builder import define_fragments, load_ligand, write_xyz  # noqa: E402

EV_TO_KCAL = 23.06035
A_CA_HA = 108.0
TETRA = 2.0 * np.pi / 3.0


def dihedral_angle(p0, p1, p2, p3):
    """Dihedral p0-p1-p2-p3 in radians."""
    b0, b1, b2 = p0 - p1, p2 - p1, p3 - p2
    b1n = b1 / np.linalg.norm(b1)
    v = b0 - (b0 @ b1n) * b1n
    w = b2 - (b2 @ b1n) * b1n
    return float(np.arctan2(np.cross(b1n, v) @ w, v @ w))


def is_l_residue(n, ca, c, cb):
    """L amino acids have det[N-CA, C-CA, CB-CA] > 0."""
    return np.linalg.det(np.array([n - ca, c - ca, cb - ca])) > 0


def planar_third(centre, u1, u2, bond):
    """The third substituent of an sp2 centre: opposite the other two, in their plane."""
    d1 = (u1 - centre) / np.linalg.norm(u1 - centre)
    d2 = (u2 - centre) / np.linalg.norm(u2 - centre)
    d = -(d1 + d2)
    return centre + bond * d / np.linalg.norm(d)


def attached_atoms(att):
    """Indices of the fragment atoms actually bonded to the capping carbon.

    The tyrosine entry in the library is toluene plus a loose water molecule (see
    `attach_geom.py`). Left in, that water is frozen during the constrained relaxation and free to
    migrate during the unconstrained one, so it donates its solvation energy straight into the
    strain. This is what `--drop-detached` removes.
    """
    adj = perceive_bonds(att["frag_symbols"], att["frag_coords"])
    for comp in connected_components(adj):
        if att["site"]["c"] in comp:
            return set(comp)
    return set(range(len(att["frag_symbols"])))


def build_pair(att_a, h_a, att_b, h_b, k, params, drop_detached=False):
    """Assemble the condensed molecule for one pose pair.

    Both fragments keep their designed coordinates exactly, less the one hydrogen on each capping
    carbon that the backbone replaces. Returns (symbols, coords, charge, fixed, info).
    `fixed` indexes the fragment atoms, the ones held during the constrained relaxation.
    """
    ca_a = dict(att_a["ca"])[h_a]
    ca_b = dict(att_b["ca"])[h_b]
    cb_a, cg_a, cb_b, cg_b = att_a["cb"], att_a["cg"], att_b["cb"], att_b["cg"]

    chain, n_b = grow_backbone(cg_a, cb_a, ca_a, k, params)
    c_a = chain[1][1]
    t1 = params[0]

    # residue A: N and HA fill the two tetrahedral slots left on CA by CB and C, one dihedral
    # step either side of C. Which way round is the chirality.
    n_a = place_atom(cg_a, cb_a, ca_a, B_N_CA, A_CB_CA_N, t1 + TETRA)
    ha_a = place_atom(cg_a, cb_a, ca_a, B_CA_HA, A_CA_HA, t1 - TETRA)
    if not is_l_residue(n_a, ca_a, c_a, cb_a):
        n_a = place_atom(cg_a, cb_a, ca_a, B_N_CA, A_CB_CA_N, t1 - TETRA)
        ha_a = place_atom(cg_a, cb_a, ca_a, B_CA_HA, A_CA_HA, t1 + TETRA)

    # residue B: CA and CB are fixed by the pose, N comes from the chain, so C and HA follow
    d0 = dihedral_angle(cg_b, cb_b, ca_b, n_b)
    c_b = place_atom(cg_b, cb_b, ca_b, B_CA_C, A_CB_CA_C, d0 + TETRA)
    ha_b = place_atom(cg_b, cb_b, ca_b, B_CA_HA, A_CA_HA, d0 - TETRA)
    if not is_l_residue(n_b, ca_b, c_b, cb_b):
        c_b = place_atom(cg_b, cb_b, ca_b, B_CA_C, A_CB_CA_C, d0 - TETRA)
        ha_b = place_atom(cg_b, cb_b, ca_b, B_CA_HA, A_CA_HA, d0 + TETRA)

    symbols, coords, fixed = [], [], []

    def add(sym, xyz, is_fixed=False):
        if is_fixed:
            fixed.append(len(symbols))
        symbols.append(sym)
        coords.append(np.asarray(xyz, dtype=float))

    frag_span, frag_map = {}, {}

    def add_fragment(att, skip, tag):
        keep = attached_atoms(att) if drop_detached else set(range(len(att["frag_symbols"])))
        start = len(symbols)
        mapping = {}
        for i, (sym, xyz) in enumerate(zip(att["frag_symbols"], att["frag_coords"])):
            if i != skip and i in keep:
                mapping[i] = len(symbols)
                add(sym, xyz, is_fixed=True)
        frag_span[tag] = [i for i in range(start, len(symbols)) if symbols[i] != "H"]
        frag_map[tag] = {"all": sorted(mapping.values()), "cb": mapping[att["site"]["c"]],
                         "charge": att["charge"]}

    # ---- residue A. Both caps are neutral and arbitrary, so their orientation is chosen at the
    # end, once everything else is in place; they start as placeholders.
    i_n_a = len(symbols); add("N", n_a)
    cap_n = [len(symbols), len(symbols) + 1]
    add("H", n_a); add("H", n_a)
    i_ca_a = len(symbols); add("C", ca_a); add("H", ha_a)
    add_fragment(att_a, h_a, "a")
    add("C", c_a)
    # the carbonyl O of residue A points away from its CA and the next N
    add("O", planar_third(c_a, ca_a, chain[2][1], B_C_O))

    # ---- spacer glycines: chain is [CA_A, C_A, N_1, (CA_i, C_i, N_i+1) * k]
    prev_c = c_a
    for i in range(k):
        n_i = chain[2 + 3 * i][1]
        ca_i, c_i, n_next = (chain[3 + 3 * i][1], chain[4 + 3 * i][1], chain[5 + 3 * i][1])
        add("N", n_i); add("H", planar_third(n_i, prev_c, ca_i, B_N_H))
        add("C", ca_i)
        phi = params[2 + 2 * i]
        add("H", place_atom(prev_c, n_i, ca_i, B_CA_HA, A_CA_HA, phi + TETRA))
        add("H", place_atom(prev_c, n_i, ca_i, B_CA_HA, A_CA_HA, phi - TETRA))
        add("C", c_i); add("O", planar_third(c_i, ca_i, n_next, B_C_O))
        prev_c = c_i

    # ---- residue B
    add("N", n_b); add("H", planar_third(n_b, prev_c, ca_b, B_N_H))
    i_ca_b = len(symbols); add("C", ca_b); add("H", ha_b)
    add_fragment(att_b, h_b, "b")
    i_c_b = len(symbols); add("C", c_b)
    cap_c = [len(symbols) + j for j in range(4)]
    add("O", c_b); add("N", c_b); add("H", c_b); add("H", c_b)

    coords = np.array(coords)

    def place_nh2(angle):
        return [place_atom(cb_a, ca_a, n_a, B_N_H, 109.5, angle),
                place_atom(cb_a, ca_a, n_a, B_N_H, 109.5, angle + TETRA)]

    def place_amide(psi):
        n_cap = place_atom(n_b, ca_b, c_b, B_C_N, A_CA_C_N, psi)
        o = planar_third(c_b, ca_b, n_cap, B_C_O)
        return [o, n_cap,
                place_atom(o, c_b, n_cap, B_N_H, 120.0, 0.0),
                place_atom(o, c_b, n_cap, B_N_H, 120.0, np.pi)]

    def orient(place, own, exclude):
        """The cap orientation that keeps the furthest from every other atom."""
        keep = [i for i in range(len(coords)) if i not in set(own) | set(exclude)]
        others = coords[keep]
        best = (-np.inf, None)
        for angle in np.linspace(-np.pi, np.pi, 36, endpoint=False):
            xyz = np.array(place(angle))
            d = float(np.linalg.norm(xyz[:, None, :] - others[None, :, :], axis=-1).min())
            if d > best[0]:
                best = (d, xyz)
        return best[1]

    coords[cap_n] = orient(place_nh2, cap_n, [i_n_a, i_ca_a])
    coords[cap_c] = orient(place_amide, cap_c, [i_c_b, i_ca_b])

    frag_map["a"]["ca"] = i_ca_a
    frag_map["b"]["ca"] = i_ca_b
    info = {
        "sc": frag_map,
        "frag_a_idx": frag_span["a"], "frag_b_idx": frag_span["b"],
        "n_atoms": len(symbols),
        "closure": closure_error(chain, ca_b, cb_b),
        "n_b_ca_b": float(np.linalg.norm(n_b - ca_b)),
        "d_ca": float(np.linalg.norm(ca_a - ca_b)),
        "theta_a": angle_between(ca_b - ca_a, cb_a - ca_a),
        "theta_b": angle_between(ca_a - ca_b, cb_b - ca_b),
    }
    return symbols, coords, att_a["charge"] + att_b["charge"], fixed, info


def kabsch_rmsd(p, q):
    """RMSD after the best rigid superposition of p onto q."""
    p = np.asarray(p) - np.mean(p, axis=0)
    q = np.asarray(q) - np.mean(q, axis=0)
    v, _, w = np.linalg.svd(p.T @ q)
    d = np.sign(np.linalg.det(v @ w))
    rot = v @ np.diag([1.0, 1.0, d]) @ w
    return float(np.sqrt(np.mean(np.sum((p @ rot - q) ** 2, axis=1))))


# bonds each element may carry: sp2 carbon has 3 heavy neighbours, amide N 3, protonated
# lysine NZ 4, carbonyl O 1, hydroxyl O 2
ALLOWED_DEGREE = {"H": (1,), "C": (3, 4), "N": (2, 3, 4), "O": (1, 2)}


def validate_structure(symbols, coords, min_heavy=2.2, min_hydrogen=1.5):
    """Bond-count and contact audit of a built structure. Returns a list of complaints.

    This is what found the HA-on-top-of-the-carbonyl bug: closure was perfect, the steric penalty
    was clean and the chirality check passed, because none of them look at hydrogens. Counting
    bonds on every atom does.
    """
    adj = perceive_bonds(symbols, coords)
    problems = []
    pieces = connected_components(adj)
    if len(pieces) > 1:
        problems.append(f"{len(pieces)} disconnected pieces: {[len(p) for p in pieces]}")
    for i, sym in enumerate(symbols):
        if len(adj[i]) not in ALLOWED_DEGREE[sym]:
            neighbours = ", ".join(f"{symbols[j]}{j}" for j in adj[i])
            problems.append(f"{sym}{i} has {len(adj[i])} bonds ({neighbours})")
    xyz = np.asarray(coords)
    for i in range(len(symbols)):
        for j in range(i + 1, len(symbols)):
            if j in adj[i]:
                continue
            d = float(np.linalg.norm(xyz[i] - xyz[j]))
            limit = min_hydrogen if "H" in (symbols[i], symbols[j]) else min_heavy
            if d < limit:
                problems.append(f"{symbols[i]}{i}...{symbols[j]}{j} only {d:.2f} A apart")
    return problems


def worst_contact(symbols, coords, cutoff=1.9):
    """Closest heavy-atom pair separated by more than two bonds: a real clash, in angstroms.

    An earlier version of this just took the shortest non-bonded distance, which returned about
    2.25 A for every structure regardless of quality -- it was finding 1-3 neighbours (N...C across
    a bond angle sits near 2.4 A), not contacts. Bonded and 1-3 pairs have to be excluded through
    the connectivity for the number to mean anything.
    """
    adj = perceive_bonds(symbols, coords)
    near = [set([i]) | set(adj[i]) for i in range(len(symbols))]
    for i in range(len(symbols)):
        near[i] = near[i].union(*(set(adj[j]) for j in adj[i])) if adj[i] else near[i]
    heavy = [i for i, sym in enumerate(symbols) if sym != "H"]
    xyz = np.asarray(coords)
    best = np.inf
    for a_i, i in enumerate(heavy):
        for j in heavy[a_i + 1:]:
            if j in near[i]:
                continue
            d = float(np.linalg.norm(xyz[i] - xyz[j]))
            if d < best:
                best = d
    return best if np.isfinite(best) else float("nan")


def sidechain_system(symbols, xyz, info, tags):
    """One or both side chains on their own, each recapped with the hydrogen the backbone displaced.

    `tags` is ("a",), ("b",) or ("a", "b"). Returns (symbols, coords, charge).
    """
    syms, pos, charge = [], [], 0
    for tag in tags:
        sc = info["sc"][tag]
        for gi in sc["all"]:
            syms.append(symbols[gi])
            pos.append(xyz[gi])
        cb, ca = np.asarray(xyz[sc["cb"]]), np.asarray(xyz[sc["ca"]])
        direction = ca - cb
        syms.append("H")
        pos.append(cb + 1.09 * direction / np.linalg.norm(direction))
        charge += sc["charge"]
    return syms, np.array(pos), charge


def sidechain_terms(symbols, xyz, info, calculator):
    """(interaction, internal) energies of the side chains at one geometry, in eV.

    Two separate things change when the molecule is released, and lumping them together is what
    made an earlier version of this produce a *negative* "backbone" residual:

      interaction  E(A+B) - E(A) - E(B): how the two side chains see each other, which is where a
                   salt bridge forming on release shows up;
      internal     E(A) + E(B): each side chain's own relaxation, which has nothing to do with the
                   other one or with the backbone.
    """
    e_ab = single_point(*sidechain_system(symbols, xyz, info, ("a", "b")), calculator)
    e_a = single_point(*sidechain_system(symbols, xyz, info, ("a",)), calculator)
    e_b = single_point(*sidechain_system(symbols, xyz, info, ("b",)), calculator)
    return e_ab - e_a - e_b, e_a + e_b


def multiplicity(symbols, charge):
    """Spin multiplicity from the electron count: a doublet when odd, a singlet when even.

    Hardcoding 1 is wrong for the tyrosine fragment, which is a radical (59 electrons at charge 0)
    that the library declares as a singlet. It is the only entry affected.
    """
    return 2 if electron_count(symbols, charge) % 2 else 1


def single_point(symbols, coords, charge, calculator):
    from ase import Atoms
    atoms = Atoms(symbols, positions=coords)
    atoms.info.update({"charge": int(charge), "spin": multiplicity(symbols, charge)})
    atoms.calc = calculator
    return atoms.get_potential_energy()


def relax(symbols, coords, charge, calculator, fixed=None, fmax=0.05, steps=400):
    """BFGS relaxation; returns (energy in eV, coordinates, steps taken, converged)."""
    from ase import Atoms
    from ase.constraints import FixAtoms
    from ase.optimize import BFGS

    atoms = Atoms(symbols, positions=coords)
    atoms.info.update({"charge": int(charge), "spin": multiplicity(symbols, charge)})
    atoms.calc = calculator
    if fixed:
        atoms.set_constraint(FixAtoms(indices=list(fixed)))
    # logfile='-' not None: with None a run killed mid-relaxation leaves no record of how far it
    # got, and these relaxations are the expensive part of a sweep.
    opt = BFGS(atoms, logfile='-')
    converged = opt.run(fmax=fmax, steps=steps)
    return atoms.get_potential_energy(), atoms.get_positions(), opt.get_number_of_steps(), converged


def pick_cases(outdir, n_cases, kmax=4, pairs_csv=None, max_closure=None, tol=0.2):
    """Choose pair/k cases from the sweep, spread over the range of CA-CA distance.

    Per ordered pair, the case is the **smallest** k that connects, which is what a design would
    actually use; a pair that never connects contributes its closest attempt instead, because the
    strain of a connection that cannot be made is worth measuring too. Picking the best closure
    outright would instead favour large k, where the chain has slack to spare and the answer is
    uninteresting.
    """
    path = pairs_csv or os.path.join(outdir, "condense_pairs.csv")
    if not os.path.isabs(path) and not os.path.exists(path):
        path = os.path.join(outdir, path)
    with open(path) as fh:
        rows = [r for r in csv.DictReader(fh) if int(r["k"]) <= kmax]
    if max_closure is not None:
        rows = [r for r in rows if float(r["closure"]) <= max_closure]
    if not rows:
        raise SystemExit(f"no rows in {path} pass the filters")

    by_pair = {}
    for r in rows:
        key = (r["frag_a"], r["pose_a"], r["frag_b"], r["pose_b"])
        by_pair.setdefault(key, []).append(r)
    best = []
    for group in by_pair.values():
        connecting = [r for r in group if float(r["closure"]) < tol]
        if connecting:
            best.append(min(connecting, key=lambda r: (int(r["k"]), float(r["closure"]))))
        else:
            best.append(min(group, key=lambda r: float(r["closure"])))

    ordered = sorted(best, key=lambda r: float(r["d_ca"]))
    if n_cases >= len(ordered):
        return ordered
    idx = np.linspace(0, len(ordered) - 1, n_cases).round().astype(int)
    return [ordered[i] for i in sorted(set(idx))]


def run_cases(args, cases, atts, calculator, out_csv, xyz_dir, ligand_coords=None,
              ligand_full=None):
    fields = ["frag_a", "pose_a", "h_a", "frag_b", "pose_b", "h_b", "k", "n_atoms", "charge",
              "d_ca", "theta_a", "theta_b", "closure", "penalty", "e_held_ev", "e_free_ev",
              "strain_kcal", "sc_interaction_kcal", "sc_internal_kcal", "sidechain_kcal",
              "remainder_kcal", "rmsd_held_free",
              "frag_gap_held", "frag_gap_free", "contact_held", "contact_ligand",
              "steps_held", "steps_free", "conv_held", "conv_free"]
    done = set()
    if os.path.exists(out_csv):
        with open(out_csv) as fh:
            reader = csv.DictReader(fh)
            existing = reader.fieldnames or []
            for r in reader:
                done.add((r["frag_a"], r["pose_a"], r["h_a"], r["frag_b"], r["pose_b"],
                          r["h_b"], r["k"]))
        if existing != fields:
            # an older run wrote different columns; keep it and start a new file beside it
            base, ext = os.path.splitext(out_csv)
            n = 2
            while os.path.exists(f"{base}_v{n}{ext}"):
                n += 1
            out_csv = f"{base}_v{n}{ext}"
            print(f"existing CSV has different columns; writing {os.path.basename(out_csv)} "
                  f"and leaving the old one alone", flush=True)
            done = set()
            fh_out = open(out_csv, "w", newline="")
            writer = csv.DictWriter(fh_out, fieldnames=fields)
            writer.writeheader()
        else:
            fh_out = open(out_csv, "a", newline="")
            writer = csv.DictWriter(fh_out, fieldnames=fields)
    else:
        fh_out = open(out_csv, "w", newline="")
        writer = csv.DictWriter(fh_out, fieldnames=fields)
        writer.writeheader()

    with fh_out:
        for case in cases:
            key = (case["frag_a"], case["pose_a"], case["h_a"], case["frag_b"],
                   case["pose_b"], case["h_b"], case["k"])
            if key in done:
                print(f"skip {key} (already done)", flush=True)
                continue
            att_a = atts[(case["frag_a"], int(case["pose_a"]))]
            att_b = atts[(case["frag_b"], int(case["pose_b"]))]
            k = int(case["k"])
            ca_a = dict(att_a["ca"])[int(case["h_a"])]
            ca_b = dict(att_b["ca"])[int(case["h_b"])]
            obstacles = build_obstacles(att_a, ca_a, att_b, ca_b, k, ligand_coords)
            _, params, _, penalty = optimize_closure(att_a["cg"], att_a["cb"], ca_a, ca_b,
                                                    att_b["cb"], k, restarts=args.restarts,
                                                    seed=args.seed, obstacles=obstacles)
            symbols, coords, charge, fixed, info = build_pair(
                att_a, int(case["h_a"]), att_b, int(case["h_b"]), k, params,
                drop_detached=args.drop_detached)
            info["penalty"] = penalty
            tag = (f"{case['frag_a']}{case['pose_a']}h{case['h_a']}_"
                   f"{case['frag_b']}{case['pose_b']}h{case['h_b']}_k{k}")
            write_xyz(os.path.join(xyz_dir, tag + "_built.xyz"), symbols, coords,
                      f"condensed pair, closure {info['closure']:.3f} A, charge {charge}")
            for complaint in validate_structure(symbols, coords):
                print(f"    STRUCTURE: {complaint}", flush=True)
            print(f"{tag}: {info['n_atoms']} atoms, charge {charge}, closure "
                  f"{info['closure']:.3f} A, penalty {info['penalty']:.3f}, "
                  f"d_CA {info['d_ca']:.2f}", flush=True)

            if args.with_ligand and ligand_full is not None:
                # Relax the held structure with the ligand present and fixed, so the backbone
                # cannot relax into space the ligand occupies -- then take the peptide-only
                # energy at that geometry, so binding energy still stays out of the strain.
                lig_symbols, lig_xyz, lig_charge = ligand_full
                all_symbols = symbols + list(lig_symbols)
                all_coords = np.vstack([coords, lig_xyz])
                held_fixed = list(fixed) + list(range(len(symbols), len(all_symbols)))
                _, xyz_all, st_h, cv_h = relax(all_symbols, all_coords, charge + lig_charge,
                                               calculator, fixed=held_fixed, steps=args.steps)
                xyz_held = xyz_all[:len(symbols)]
                e_held = single_point(symbols, xyz_held, charge, calculator)
            else:
                e_held, xyz_held, st_h, cv_h = relax(symbols, coords, charge, calculator,
                                                     fixed=fixed, steps=args.steps)
            write_xyz(os.path.join(xyz_dir, tag + "_held.xyz"), symbols, xyz_held,
                      f"fragments fixed, E = {e_held:.6f} eV")
            e_free, xyz_free, st_f, cv_f = relax(symbols, xyz_held, charge, calculator,
                                                steps=args.steps)
            write_xyz(os.path.join(xyz_dir, tag + "_free.xyz"), symbols, xyz_free,
                      f"free relaxation, E = {e_free:.6f} eV")
            strain = (e_held - e_free) * EV_TO_KCAL
            rmsd = kabsch_rmsd(xyz_held, xyz_free)
            contact = worst_contact(symbols, xyz_held)
            # the ligand is kept out of the energy so that binding does not masquerade as strain,
            # but the relaxed backbone still has to leave room for it
            def frag_gap(xyz):
                """Closest approach between the two side chains: a salt bridge forming in the
                free relaxation shows up here as a collapse from the designed separation."""
                a = np.asarray(xyz)[info["frag_a_idx"]]
                b = np.asarray(xyz)[info["frag_b_idx"]]
                return float(np.linalg.norm(a[:, None, :] - b[None, :, :], axis=-1).min())

            gap_held, gap_free = frag_gap(xyz_held), frag_gap(xyz_free)

            # split the strain into what the side chains do and what is left over
            int_held, self_held = sidechain_terms(symbols, xyz_held, info, calculator)
            int_free, self_free = sidechain_terms(symbols, xyz_free, info, calculator)
            sc_interaction = (int_held - int_free) * EV_TO_KCAL
            sc_internal = (self_held - self_free) * EV_TO_KCAL
            sidechain = sc_interaction + sc_internal
            backbone_strain = strain - sidechain
            backbone_atoms = [i for i in range(len(symbols))
                              if i not in set(fixed) and symbols[i] != "H"]
            lig_contact = (float(np.linalg.norm(
                xyz_held[backbone_atoms][:, None, :] - ligand_coords[None, :, :], axis=-1).min())
                if ligand_coords is not None and backbone_atoms else float("nan"))
            print(f"    strain {strain:8.2f} = sc_int {sc_interaction:7.2f} + sc_self "
                  f"{sc_internal:7.2f} + rest {backbone_strain:7.2f} kcal/mol   "
                  f"gap {gap_held:.2f} -> {gap_free:.2f} A   "
                  f"rmsd {rmsd:.2f} A   contact {contact:.2f} A   ligand {lig_contact:.2f} A   "
                  f"(held {st_h} conv={cv_h}, free {st_f} conv={cv_f})", flush=True)

            writer.writerow({
                "frag_a": case["frag_a"], "pose_a": case["pose_a"], "h_a": case["h_a"],
                "frag_b": case["frag_b"], "pose_b": case["pose_b"], "h_b": case["h_b"],
                "k": k, "n_atoms": info["n_atoms"], "charge": charge,
                "d_ca": round(info["d_ca"], 3), "theta_a": round(info["theta_a"], 1),
                "theta_b": round(info["theta_b"], 1), "closure": round(info["closure"], 4),
                "penalty": round(float(info["penalty"]), 4),
                "e_held_ev": f"{e_held:.6f}", "e_free_ev": f"{e_free:.6f}",
                "strain_kcal": f"{strain:.3f}",
                "sc_interaction_kcal": f"{sc_interaction:.3f}",
                "sc_internal_kcal": f"{sc_internal:.3f}",
                "sidechain_kcal": f"{sidechain:.3f}",
                "remainder_kcal": f"{backbone_strain:.3f}", "rmsd_held_free": f"{rmsd:.3f}",
                "frag_gap_held": f"{gap_held:.3f}", "frag_gap_free": f"{gap_free:.3f}",
                "contact_held": f"{contact:.3f}", "contact_ligand": f"{lig_contact:.3f}",
                "steps_held": st_h, "steps_free": st_f,
                "conv_held": int(bool(cv_h)), "conv_free": int(bool(cv_f)),
            })
            fh_out.flush()


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("outdir")
    parser.add_argument("--ligand", default="octinoxate")
    parser.add_argument("--cases", type=int, default=12, help="how many pair/k cases to run")
    parser.add_argument("--kmax", type=int, default=4)
    parser.add_argument("--pairs-csv", help="sweep CSV to draw cases from "
                        "(default condense_pairs.csv inside outdir)")
    parser.add_argument("--tol", type=float, default=0.2,
                        help="closure that counts as connected when choosing cases")
    parser.add_argument("--max-closure", type=float, default=None,
                        help="only use cases whose closure is at most this, in A")
    parser.add_argument("--pair", nargs=4, metavar=("FRAG_A", "POSE_A", "FRAG_B", "POSE_B"),
                        help="one named pair instead of a spread of cases")
    parser.add_argument("--k", type=int, default=1, help="spacer count for --pair")
    parser.add_argument("--restarts", type=int, default=24)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--steps", type=int, default=400, help="max BFGS steps per relaxation")
    parser.add_argument("--model", default="uma-s-1p2p1")
    parser.add_argument("--out", help="CSV name inside outdir (default condense_strain.csv)")
    parser.add_argument("--with-ligand", action="store_true",
                        help="hold the ligand fixed during the constrained relaxation so the "
                             "backbone cannot relax into it; energies stay peptide-only")
    parser.add_argument("--drop-detached", action="store_true",
                        help="drop fragment atoms not bonded to the capping carbon "
                             "(the loose water in the tyrosine entry)")
    parser.add_argument("--build-only", action="store_true",
                        help="write the built structures and stop, no UMA")
    args = parser.parse_args(argv)

    frags = define_fragments()
    ligand = load_ligand(args.ligand)
    atts = {(a["name"], a["pose"]): a for a in selected_attachments(args.outdir, frags, ligand)}

    if args.pair:
        frag_a, pose_a, frag_b, pose_b = args.pair
        att_a, att_b = atts[(frag_a, int(pose_a))], atts[(frag_b, int(pose_b))]
        cases = []
        for h_a, ca_a in att_a["ca"]:
            for h_b, ca_b in att_b["ca"]:
                obstacles = build_obstacles(att_a, ca_a, att_b, ca_b, args.k,
                                            np.asarray(ligand["coords"]))
                err, _, _, _ = optimize_closure(att_a["cg"], att_a["cb"], ca_a, ca_b,
                                                att_b["cb"], args.k, restarts=args.restarts,
                                                seed=args.seed, obstacles=obstacles)
                cases.append({"frag_a": frag_a, "pose_a": pose_a, "h_a": str(h_a),
                              "frag_b": frag_b, "pose_b": pose_b, "h_b": str(h_b),
                              "k": str(args.k), "closure": err, "d_ca": 0.0})
        cases = [min(cases, key=lambda c: c["closure"])]
    else:
        cases = pick_cases(args.outdir, args.cases, kmax=args.kmax,
                           pairs_csv=args.pairs_csv, max_closure=args.max_closure,
                           tol=args.tol)
    print(f"{len(cases)} case(s) to run\n")

    xyz_dir = os.path.join(args.outdir, "condense")
    os.makedirs(xyz_dir, exist_ok=True)

    if args.build_only:
        for case in cases:
            att_a = atts[(case["frag_a"], int(case["pose_a"]))]
            att_b = atts[(case["frag_b"], int(case["pose_b"]))]
            ca_a = dict(att_a["ca"])[int(case["h_a"])]
            ca_b = dict(att_b["ca"])[int(case["h_b"])]
            k = int(case["k"])
            obstacles = build_obstacles(att_a, ca_a, att_b, ca_b, k,
                                        np.asarray(ligand["coords"]))
            _, params, _, penalty = optimize_closure(att_a["cg"], att_a["cb"], ca_a, ca_b,
                                                    att_b["cb"], k, restarts=args.restarts,
                                                    seed=args.seed, obstacles=obstacles)
            symbols, coords, charge, fixed, info = build_pair(
                att_a, int(case["h_a"]), att_b, int(case["h_b"]), k, params,
                drop_detached=args.drop_detached)
            info["penalty"] = penalty
            tag = (f"{case['frag_a']}{case['pose_a']}h{case['h_a']}_"
                   f"{case['frag_b']}{case['pose_b']}h{case['h_b']}_k{k}")
            write_xyz(os.path.join(xyz_dir, tag + "_built.xyz"), symbols, coords,
                      f"condensed pair, closure {info['closure']:.3f} A, charge {charge}")
            for complaint in validate_structure(symbols, coords):
                print(f"    STRUCTURE: {complaint}", flush=True)
            print(f"{tag}: {info['n_atoms']} atoms, charge {charge}, "
                  f"closure {info['closure']:.3f} A, penalty {info['penalty']:.3f}, "
                  f"fixed {len(fixed)}")
        return 0

    import torch
    from fairchem.core import FAIRChemCalculator, pretrained_mlip
    device = "cuda" if torch.cuda.is_available() else "cpu"
    predictor = pretrained_mlip.get_predict_unit(args.model, device=device)
    calculator = FAIRChemCalculator(predictor, task_name="omol")

    out_csv = os.path.join(args.outdir, args.out or "condense_strain.csv")
    lig_xyz = np.asarray(ligand["coords"])
    lig_heavy = lig_xyz[[i for i, sym in enumerate(ligand["atoms"]) if sym != "H"]]
    run_cases(args, cases, atts, calculator, out_csv, xyz_dir, ligand_coords=lig_heavy,
              ligand_full=(list(ligand["atoms"]), lig_xyz, int(ligand.get("charge", 0))))
    print(f"\nresults in {out_csv}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
