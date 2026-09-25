"""Attachment geometry for the fragment library: where the backbone would join each fragment.

Every fragment in `define_fragments` is a side-chain analogue capped with a methyl (or, for
serine, a methyl that is the whole side chain) where the backbone would attach. That capping
carbon marks roughly where CB sits, and the C-H bonds pointing away from the fragment mark the
directions CA could go. The rest of the pipeline throws this away -- `build_sequence` chains
fragment centroids -- which is why a design can ask for a backbone path no peptide can follow.

What this module provides, for a fragment or for a placed pose of one:

  * bond perception from distances (no RDKit needed; the fragments are small and unambiguous)
  * the capping carbon(s): a carbon with three hydrogens and exactly one heavy neighbour
  * a symmetry test over those candidates, so an ambiguous fragment (n-butane, isobutane) is
    reported as ambiguous-but-equivalent rather than as a choice that has to be made
  * per candidate, the CA positions obtained by replacing each of its hydrogens with a CA at
    1.53 A -- the discrete set of places the backbone can start from, with the fragment held
    exactly where the search put it

Usage:
    python code/attach_geom.py                    # table for the built-in library
    python code/attach_geom.py --pose runs/octinoxate/poses/octinoxate_w_leucine3.xyz \
                               --frag leucine --ligand-atoms 24
"""
import argparse
import itertools
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from peptide_builder import define_fragments, load_ligand  # noqa: E402

# enough for the library: C, N, O, H
COVALENT_RADIUS = {"H": 0.31, "C": 0.76, "N": 0.71, "O": 0.66, "S": 1.05}
ATOMIC_NUMBER = {"H": 1, "C": 6, "N": 7, "O": 8, "S": 16}
CB_CA_BOND = 1.53


def perceive_bonds(symbols, coords, tol=0.4):
    """Adjacency lists from interatomic distances: bonded if d < r1 + r2 + tol."""
    coords = np.asarray(coords, dtype=float)
    n = len(symbols)
    adj = [[] for _ in range(n)]
    for i, j in itertools.combinations(range(n), 2):
        cutoff = COVALENT_RADIUS[symbols[i]] + COVALENT_RADIUS[symbols[j]] + tol
        if np.linalg.norm(coords[i] - coords[j]) < cutoff:
            adj[i].append(j)
            adj[j].append(i)
    return adj


def connected_components(adj):
    """Atom index sets of the separate molecules in one coordinate block."""
    seen, comps = set(), []
    for start in range(len(adj)):
        if start in seen:
            continue
        stack, comp = [start], []
        seen.add(start)
        while stack:
            at = stack.pop()
            comp.append(at)
            for nb in adj[at]:
                if nb not in seen:
                    seen.add(nb)
                    stack.append(nb)
        comps.append(sorted(comp))
    return comps


def environment(symbols, coords, idx):
    """Rounded (symbol, distance) multiset seen from one atom: an order-independent
    fingerprint, so two atoms related by a symmetry of the fragment compare equal."""
    coords = np.asarray(coords, dtype=float)
    d = np.linalg.norm(coords - coords[idx], axis=1)
    return tuple(sorted((symbols[j], round(float(d[j]), 2)) for j in range(len(symbols)) if j != idx))


def attachment_sites(symbols, coords, adj=None):
    """Capping carbons: a carbon bonded to three hydrogens and one heavy atom.

    Returns one dict per candidate:
        c      index of the capping carbon (the approximate CB position)
        heavy  index of its single heavy neighbour (the approximate CG position)
        h      indices of its hydrogens, each of which is a place CA could go
        env    symmetry fingerprint of the capping carbon
    """
    if adj is None:
        adj = perceive_bonds(symbols, coords)
    sites = []
    for i, sym in enumerate(symbols):
        if sym != "C":
            continue
        hyd = [j for j in adj[i] if symbols[j] == "H"]
        heavy = [j for j in adj[i] if symbols[j] != "H"]
        if len(hyd) == 3 and len(heavy) == 1:
            sites.append({"c": i, "heavy": heavy[0], "h": hyd,
                          "env": environment(symbols, coords, i)})
    return sites


def sites_equivalent(sites):
    """True when every candidate has the same fingerprint, i.e. the choice does not matter."""
    return len(sites) > 1 and len({s["env"] for s in sites}) == 1


def ca_candidates(coords, site, bond=CB_CA_BOND):
    """Where CA can sit: one position per hydrogen on the capping carbon, that hydrogen
    replaced by a carbon at `bond` A along the same direction.

    Returns a list of (h_index, ca_xyz). Holding the fragment fixed leaves no rotation
    about the CG-CB bond, so this discrete set is the whole choice.
    """
    coords = np.asarray(coords, dtype=float)
    cb = coords[site["c"]]
    out = []
    for h in site["h"]:
        direction = coords[h] - cb
        direction /= np.linalg.norm(direction)
        out.append((h, cb + bond * direction))
    return out


def pose_attachment(frag, pose_coords, site_index=0):
    """Attachment geometry of one placed pose.

    Args:
        frag: fragment dictionary (for its atom symbols)
        pose_coords: that fragment's coordinates in the pose (fragment atoms only)
        site_index: which capping carbon, when a fragment has several equivalent ones
    Returns:
        dict with cb, cg, the CA candidates, and the attachment vector cb - cg
    """
    symbols = frag["atoms"]
    sites = attachment_sites(symbols, pose_coords)
    site = sites[site_index]
    coords = np.asarray(pose_coords, dtype=float)
    vec = coords[site["c"]] - coords[site["heavy"]]
    return {"name": frag["name"], "site": site, "cb": coords[site["c"]],
            "cg": coords[site["heavy"]], "attach_vector": vec / np.linalg.norm(vec),
            "ca": ca_candidates(coords, site), "n_sites": len(sites),
            "equivalent": sites_equivalent(sites)}


def electron_count(symbols, charge):
    """Electrons in a fragment. Odd means a radical, which must be a doublet, not a singlet."""
    return sum(ATOMIC_NUMBER[s] for s in symbols) - charge


def spin_is_consistent(symbols, charge, spin):
    """A closed-shell singlet needs an even electron count; a doublet an odd one."""
    odd = electron_count(symbols, charge) % 2 == 1
    return (spin == 2) if odd else (spin == 1)


def describe_library(frags):
    """One row per fragment: composition, disconnected pieces, capping carbons, ambiguity."""
    rows = []
    for frag in frags:
        symbols, coords = frag["atoms"], np.asarray(frag["coords"], dtype=float)
        adj = perceive_bonds(symbols, coords)
        comps = connected_components(adj)
        sites = attachment_sites(symbols, coords, adj)
        formula = "".join(f"{s}{symbols.count(s)}" for s in ("C", "H", "N", "O") if symbols.count(s))
        rows.append({
            "name": frag["name"], "n_atoms": len(symbols), "formula": formula,
            "pieces": [len(c) for c in comps],
            "sites": [s["c"] for s in sites],
            "equivalent": sites_equivalent(sites),
            "n_ca": [len(s["h"]) for s in sites],
            "electrons": electron_count(symbols, frag["charge"]),
            "spin_ok": spin_is_consistent(symbols, frag["charge"], frag["spin"]),
            "spin": frag["spin"],
        })
    return rows


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--pose", help="pose .xyz (ligand followed by one fragment)")
    parser.add_argument("--frag", help="fragment name in that pose")
    parser.add_argument("--ligand", default="octinoxate", help="ligand the pose was built around")
    args = parser.parse_args(argv)

    frags = define_fragments()
    by_name = {f["name"]: f for f in frags}

    if args.pose:
        import ase.io
        frag = by_name[args.frag]
        ligand = load_ligand(args.ligand)
        coords = ase.io.read(args.pose).get_positions()[ligand["num_atoms"]:]
        att = pose_attachment(frag, coords)
        print(f"{att['name']} pose {os.path.basename(args.pose)}")
        print(f"  capping carbons  {att['n_sites']}"
              + ("  (equivalent by symmetry)" if att["equivalent"] else ""))
        print(f"  CB  {np.round(att['cb'], 3)}")
        print(f"  CG  {np.round(att['cg'], 3)}")
        print(f"  CG->CB unit vector  {np.round(att['attach_vector'], 3)}")
        for h, ca in att["ca"]:
            print(f"  CA candidate (replacing H{h})  {np.round(ca, 3)}")
        return 0

    print(f"{'fragment':<14} {'atoms':>5} {'formula':<10} {'pieces':<10} "
          f"{'capping C':<12} {'CA options':<11} note")
    for row in describe_library(frags):
        note = ""
        if not row["spin_ok"]:
            note = (f"WRONG SPIN: {row['electrons']} electrons is odd, so this is a radical and "
                    f"needs spin 2, not {row['spin']}")
        elif len(row["pieces"]) > 1:
            note = f"NOT ONE MOLECULE: {len(row['pieces'])} pieces"
        elif row["equivalent"]:
            note = "candidates equivalent by symmetry"
        elif len(row["sites"]) > 1:
            note = "AMBIGUOUS: candidates are not equivalent"
        print(f"{row['name']:<14} {row['n_atoms']:>5} {row['formula']:<10} "
              f"{str(row['pieces']):<10} {str(row['sites']):<12} {str(row['n_ca']):<11} {note}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
