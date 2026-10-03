#!/usr/bin/env python3
"""Warm-restart (extend) a finished MD production leg from its OpenMM checkpoint.

The project's MD pipeline cannot resume: `omd run` writes `checkpoint.chk` via a
CheckpointReporter but `openmm_md.dynamics.run` never loads it and always
re-minimizes + re-equilibrates, so re-invoking it restarts from zero (see
HANDOFF/NEXT_STEPS "preempted leg restarts from step zero"). This script is the
missing piece: it rebuilds the *production* context exactly as dynamics.run left
it -- deserialize system.xml, add the NPT MonteCarloBarostat, add the k=0
protein-heavy positional restraint, LangevinMiddle(T, friction, dt) on the same
platform/precision -- then loadCheckpoint() and step N more ns, appending to
`<leg>/traj_ext.dcd` + `energy_ext.csv` with the step counter and clock preserved.

Because the context is reconstructed identically, the load is bit-exact: the
checkpoint's potential energy matches the last `energy.csv` row to the decimal, so
this is a true continuation (positions, velocities, box, RNG all restored), not a
re-thermalized restart. Verify the printed PE against the leg's last energy row
before trusting the extension.

It does NOT recompute MM/GBSA -- recombine the trajectories and rescore the whole
run with `md_mmgbsa_whole_run.sh` so the number stays comparable to every other
whole-run point.

    python code/md_warm_restart.py <structure> [--ns 5] [--leg prod_20ns] \
        [--platform OpenCL] [--outdir runs/octinoxate]

Example (what produced orig_f12_dock5's 20->25 ns extension, 2026-10-03):
    python code/md_warm_restart.py orig_f12_dock5 --ns 5
"""
import argparse
import os
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parent
OMD_SRC = os.environ.get("OMD_SRC", str(Path.home() / "python_mac/openmm/src"))
sys.path.insert(0, OMD_SRC)

from openmm import app, openmm as mm, unit, XmlSerializer  # noqa: E402
from openmm_md.config import Config  # noqa: E402
from openmm_md.dynamics import protein_heavy_indices, _new_sim  # noqa: E402


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("structure")
    ap.add_argument("--ns", type=float, default=5.0, help="nanoseconds to add (default 5)")
    ap.add_argument("--leg", default="prod_20ns", help="production dir holding checkpoint.chk")
    ap.add_argument("--platform", default="OpenCL")
    ap.add_argument("--outdir", default="runs/octinoxate")
    ap.add_argument("--dry", action="store_true", help="load + report energy, do not step")
    a = ap.parse_args(argv)

    cfg = Config()
    cfg.platform = a.platform
    M = REPO / a.outdir / "md" / a.structure
    P = M / a.leg
    steps = int(round(a.ns * 1e6 / cfg.timestep))  # ns -> 2 fs steps

    system = XmlSerializer.deserializeSystem(open(M / "system/system.xml").read())
    pdb = app.PDBFile(str(M / "system/complex.pdb"))
    topology, positions = pdb.topology, pdb.positions

    # Rebuild the production context exactly: xml forces -> barostat -> k=0 restraint.
    if cfg.pressure > 0:
        system.addForce(mm.MonteCarloBarostat(cfg.pressure * unit.atmosphere,
                                              cfg.temperature * unit.kelvin,
                                              cfg.barostat_interval))
    heavy = protein_heavy_indices(topology)
    if cfg.restrain_protein and heavy:
        r = mm.CustomExternalForce("k*((x-x0)^2+(y-y0)^2+(z-z0)^2)")
        r.addGlobalParameter("k", cfg.restraint_weight * unit.kilojoule_per_mole / unit.nanometer**2)
        for n in ("x0", "y0", "z0"):
            r.addPerParticleParameter(n)
        for i in heavy:
            p = positions[i].value_in_unit(unit.nanometer)
            r.addParticle(i, [float(p[0]), float(p[1]), float(p[2])])
        system.addForce(r)
        print(f"[warm] restraint force rebuilt for {len(heavy)} protein heavy atoms")

    sim = _new_sim(topology, system, cfg, a.platform)
    sim.loadCheckpoint(str(P / "checkpoint.chk"))
    if cfg.restrain_protein and heavy:
        sim.context.setParameter("k", 0.0)

    st = sim.context.getState(getEnergy=True)
    print(f"[warm] loaded: t = {st.getTime().value_in_unit(unit.picosecond):.1f} ps, "
          f"PE = {st.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole):.1f} kJ/mol "
          f"(check against the last energy.csv row before trusting this)")
    if a.dry:
        print("[warm] dry run OK, not stepping")
        return 0

    sim.reporters.append(app.StateDataReporter(
        str(P / "energy_ext.csv"), cfg.report_interval, step=True, time=True,
        potentialEnergy=True, kineticEnergy=True, totalEnergy=True,
        temperature=True, volume=True, density=True, separator=","))
    sim.reporters.append(app.DCDReporter(str(P / "traj_ext.dcd"), cfg.traj_interval))
    sim.reporters.append(app.CheckpointReporter(str(P / "checkpoint_ext.chk"), cfg.checkpoint_interval))
    print(f"[warm] extending {a.structure}/{a.leg} by {a.ns} ns ({steps} steps) on {a.platform} ...")
    sim.step(steps)
    print(f"[warm] DONE {a.ns} ns extension")
    return 0


if __name__ == "__main__":
    sys.exit(main())
