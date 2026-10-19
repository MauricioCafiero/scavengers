"""Write pocket-guided OpenFold3 input files from the Boltz forced-fold yamls.

A forced Boltz fold ties each designed residue to a named ligand atom and keeps
enforcing it (`force: true`). OpenFold3's nearest equivalent, `pocket_constraint`,
takes residues only -- so the hint's ligand atom is dropped, and the result biases
ligand placement rather than enforcing contact. These runs are therefore
**pocket-guided, not forced**; see NEXT_STEPS.md, "Forced vs unforced co-folds".

The residue list is read off the run's own yaml rather than typed, because a
mistyped position would not fail -- it would guide the ligand to the wrong face
and report a plausible-looking number. `max_distance` is the largest of that
yaml's contact distances, so no contact is tightened relative to Boltz.

    python code/pocket_inputs.py octinoxate s3_esm2_f4 orig_f12 \
        --out-dir runs/octinoxate/folds/inputs_pocket

Output is one `<struct>.txt` per structure in fold's plain-text format
(SEQUENCE/SMILES/POCKET/POCKET_MAX_DISTANCE), which `fold_batch.sh` feeds to the
app with `--input-file`.
"""
import argparse
import os
import re
import sys

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def read_forced_yaml(path):
    """Return (sequence, smiles, residues, max_distance) from a Boltz fold yaml.

    Parsed with pyyaml; the contact block is a list of `{contact: {token1: [chain,
    pos], ...}}` entries, and only `token1`'s position is usable here.
    """
    import yaml

    doc = yaml.safe_load(open(path))
    sequence = smiles = None
    for entry in doc["sequences"]:
        if "protein" in entry:
            sequence = entry["protein"]["sequence"]
        elif "ligand" in entry:
            smiles = entry["ligand"]["smiles"]
    contacts = [c["contact"] for c in doc.get("constraints", []) if "contact" in c]
    if not contacts:
        sys.exit(f"{path}: no contact constraints -- this is an unforced fold, so there "
                 f"is no hint set to carry over")
    residues = sorted({int(c["token1"][1]) for c in contacts})
    return sequence, smiles, residues, max(float(c["max_distance"]) for c in contacts)


def shell_pocket(system, name):
    """Return (residues, max_distance, shell_file) for `name` from its DESIGNED SHELL.

    The Boltz `contact` hints were never the design -- they are 4, 8 or 12 of the shell's slots,
    and every hint set is a strict subset of them (verified for all nine, 2026-10-19). That
    subsetting does not carry over to OpenFold3, because `pocket_constraint` is one residue SET
    scored by its own centroid, not N independent residue-atom restraints: naming 4 of 12 slots
    does not weaken the pocket, it MOVES it off the shell's centre. There is therefore no OF3
    analogue of the f4/f8/f12 dose, and the only faithful translation of "the design puts these
    residues on the ligand" is the whole slot set.

    Slots are the non-glycine positions of the shell's sequence -- the same rule
    pair_contacts.py:slots_of uses, and the shell is chosen by the same name map
    (pair_contacts.design_for), so a peptide's pocket and its pair-contact slots cannot drift
    apart. max_distance is the largest over EVERY dose of that peptide's yamls, so the pocket is
    no tighter than any Boltz run it is set beside.
    """
    import glob

    import yaml

    sys.path.insert(0, os.path.join(REPO, "code"))
    from pair_contacts import design_for, slots_of

    run_dir = os.path.join(REPO, "runs", system)
    shell = design_for(name)
    slots, seq = slots_of(run_dir, shell)
    base = re.sub(r"_(f\d+|control)$", "", name)
    dists = []
    for y in sorted(glob.glob(os.path.join(run_dir, "boltz", f"{base}_f*.yaml"))):
        doc = yaml.safe_load(open(y))
        dists += [float(c["contact"]["max_distance"])
                  for c in doc.get("constraints", []) or [] if "contact" in c]
    if not dists:
        sys.exit(f"{name}: no {base}_f*.yaml to take a max_distance from")
    return [i + 1 for i in slots], max(dists), shell, seq


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0],
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("system", help="molecule, e.g. octinoxate")
    p.add_argument("structures", nargs="+", help="forced fold names, e.g. s3_esm2_f4")
    p.add_argument("--out-dir", help="default runs/<system>/folds/inputs_pocket")
    p.add_argument("--from-shell", action="store_true",
                   help="take the pocket from the peptide's DESIGNED SHELL (all slots) instead of "
                        "from a yaml's contact subset, and read SEQUENCE/SMILES from the unguided "
                        "lane's own input file so the guided and unguided folds differ only by the "
                        "pocket. This is the campaign's setting -- see shell_pocket().")
    a = p.parse_args(argv)

    out = a.out_dir or os.path.join(REPO, "runs", a.system, "folds", "inputs_pocket")
    os.makedirs(out, exist_ok=True)

    for name in a.structures:
        if a.from_shell:
            # SEQUENCE/SMILES come from the unguided lane's own input, not from a yaml: that is
            # the exact pair the guided fold is compared against, so taking them from anywhere
            # else risks a silent one-residue difference between the two folds.
            src = os.path.join(REPO, "runs", a.system, "folds", "inputs", f"{name}.txt")
            if not os.path.exists(src):
                sys.exit(f"no unguided input at {src} -- the guided fold needs its twin to match")
            lines = dict(l.split(":", 1) for l in open(src) if ":" in l)
            sequence = lines.get("SEQUENCE", "").strip()
            smiles = lines.get("SMILES", "").strip()
            residues, max_distance, shell, shell_seq = shell_pocket(a.system, name)
            if len(shell_seq) != len(sequence):
                sys.exit(f"{name}: shell {shell} is {len(shell_seq)} residues but the fold is "
                         f"{len(sequence)} -- the slot positions would not line up")
            origin = f"{shell}, all {len(residues)} slots"
        else:
            yml = os.path.join(REPO, "runs", a.system, "boltz", f"{name}.yaml")
            if not os.path.exists(yml):
                sys.exit(f"no yaml at {yml}")
            sequence, smiles, residues, max_distance = read_forced_yaml(yml)
            origin = f"{name}.yaml contacts"
        if not smiles:
            sys.exit(f"{name}: no ligand -- a pocket constraint needs one")
        dest = os.path.join(out, f"{name}.txt")
        with open(dest, "w") as fh:
            fh.write(f"SEQUENCE: {sequence}\n")
            fh.write(f"SMILES: {smiles}\n")
            fh.write(f"POCKET: {','.join(str(r) for r in residues)}\n")
            fh.write(f"POCKET_MAX_DISTANCE: {max_distance}\n")
        print(f"{name}: {len(residues)} residues {residues}, max_distance {max_distance} "
              f"[{origin}] -> {os.path.relpath(dest, REPO)}")
    print("\nthese are pocket-guided inputs, not forced: OpenFold3 ranks ligand proposals "
          "against\nthese residues, it does not enforce contact, and the ligand atom of each "
          "hint is dropped.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
