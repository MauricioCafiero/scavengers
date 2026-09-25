"""Wrap a molecule in explicit waters, so a desolvation energy can be computed.

Copied from the `solvation` class in `~/python_mac/UMADock/code/UMADock/UMADock.py` and modified
here; that repo is left untouched. The algorithm is the original one: take the molecule's bounding
box, pad it by a couple of water van der Waals radii, then throw waters at it at random
orientations, rejecting any that land too close to the solute or to a water already placed, and
stop early once placements keep failing.

Three changes from the original:

* it works on ASE `Atoms` in memory instead of round-tripping through xyz files, since nothing here
  needs the files;
* the rejection distance is a parameter. The original used the water vdW radius, 1.7 A, for both
  solute-water and water-water contacts, which permits real overlaps -- 1.7 A between two oxygens
  is not a contact, it is a clash -- and leaves the subsequent relaxation to sort it out. The
  default here is 2.4 A. Pass `min_contact=1.7` to reproduce UMADock's numbers exactly;
* the random seed is settable, so a desolvation energy is reproducible;
* the rejection test covers all three atoms of the candidate water. The original tested only the
  oxygen, so hydrogens could sit well inside the intended separation.
"""
import numpy as np
from ase import Atoms

WATER_VDW = 1.7
# rigid TIP3P-like geometry, the same coordinates the original used
H1 = np.array([0.580743, 0.0, 0.758810])
H2 = np.array([0.580743, 0.0, -0.758810])


def random_rotation(rng):
    """A uniformly random rotation matrix, as three successive axis rotations."""
    tx, ty, tz = rng.uniform(0, 2 * np.pi, 3)
    rx = np.array([[1, 0, 0], [0, np.cos(tx), -np.sin(tx)], [0, np.sin(tx), np.cos(tx)]])
    ry = np.array([[np.cos(ty), 0, -np.sin(ty)], [0, 1, 0], [np.sin(ty), 0, np.cos(ty)]])
    rz = np.array([[np.cos(tz), -np.sin(tz), 0], [np.sin(tz), np.cos(tz), 0], [0, 0, 1]])
    return rz @ ry @ rx


def water_at(origin, rng):
    """One water molecule, randomly oriented, its oxygen at `origin`."""
    rot = random_rotation(rng)
    return np.array([origin, origin + rot @ H1, origin + rot @ H2])


def add_waters(atoms, max_waters, water_radii=2.0, min_contact=2.4, stopping_criteria=10,
               seed=0, verbose=True, region=None, region_radius=2.8,
               attempts_per_water=60, stop_on_streak=True):
    """`atoms` plus up to `max_waters` waters placed around it.

    With `region` -- an array of points, normally the ligand's coordinates -- waters go only where
    the ligand would be, instead of all over the solute. That is what a cavity desolvation needs:
    water on the peptide's outer surface is present both before and after binding, so it cancels,
    and including it would mean extracting a small cavity term as the difference between two large
    numbers whose random water placements do not cancel at all. `region_radius` is how far from a
    region point a water may sit.

    Returns (solvated, n_added); the solute keeps atom indices 0..len(atoms)-1, so the waters can
    be split off again by slicing.
    """
    rng = np.random.default_rng(seed)
    solute = atoms.get_positions()
    region = None if region is None else np.asarray(region, dtype=float)
    if region is None:
        lo = solute.min(axis=0) - water_radii * WATER_VDW
        hi = solute.max(axis=0) + water_radii * WATER_VDW

    def sample_origin():
        """A candidate oxygen position: anywhere in the padded box, or -- with a region -- near
        one of its points. Sampling around the points beats rejecting box samples, because a
        molecular envelope fills little of its own bounding box."""
        if region is None:
            return rng.uniform(lo, hi)
        centre = region[rng.integers(len(region))]
        direction = rng.normal(size=3)
        direction /= np.linalg.norm(direction)
        return centre + direction * region_radius * rng.random() ** (1 / 3)

    placed = []                      # every water atom position so far
    recent_failures = []
    for _ in range(max_waters * attempts_per_water):
        if len(placed) // 3 >= max_waters:
            break
        if (stop_on_streak and len(recent_failures) >= stopping_criteria
                and all(recent_failures[-stopping_criteria:])):
            if verbose:
                print(f"  stopping early: {stopping_criteria} failed placements in a row")
            break
        origin = sample_origin()
        candidate = water_at(origin, rng)
        # test all three atoms of the candidate water, not just its oxygen. The original tested
        # the oxygen only, which lets the hydrogens sit inside the cutoff -- a 2.4 A oxygen test
        # still produced 2.24 A contacts.
        too_close = (len(solute) > 0 and np.min(np.linalg.norm(
            solute[None, :, :] - candidate[:, None, :], axis=-1)) < min_contact)
        if not too_close and placed:
            too_close = np.min(np.linalg.norm(
                np.array(placed)[None, :, :] - candidate[:, None, :], axis=-1)) < min_contact
        recent_failures.append(too_close)
        if not too_close:
            placed.extend(candidate)

    n_added = len(placed) // 3
    solvated = atoms.copy()
    if n_added:
        solvated += Atoms("OHH" * n_added, positions=np.array(placed))
    if verbose:
        print(f"  added {n_added}/{max_waters} waters")
    return solvated, n_added


def waters_for(atoms, per_heavy_atom=1.2, minimum=20):
    """A water count that scales with the solute instead of being fixed at 20.

    UMADock always asks for 20, which is about right for a drug-sized ligand and far too few to
    shell a 200-atom peptide. This scales with the heavy-atom count, which tracks surface area
    closely enough for the purpose.
    """
    heavy = sum(1 for s in atoms.get_chemical_symbols() if s != "H")
    return max(minimum, int(round(per_heavy_atom * heavy)))
