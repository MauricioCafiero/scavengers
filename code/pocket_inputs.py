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


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0],
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("system", help="molecule, e.g. octinoxate")
    p.add_argument("structures", nargs="+", help="forced fold names, e.g. s3_esm2_f4")
    p.add_argument("--out-dir", help="default runs/<system>/folds/inputs_pocket")
    a = p.parse_args(argv)

    out = a.out_dir or os.path.join(REPO, "runs", a.system, "folds", "inputs_pocket")
    os.makedirs(out, exist_ok=True)

    for name in a.structures:
        yml = os.path.join(REPO, "runs", a.system, "boltz", f"{name}.yaml")
        if not os.path.exists(yml):
            sys.exit(f"no yaml at {yml}")
        sequence, smiles, residues, max_distance = read_forced_yaml(yml)
        if not smiles:
            sys.exit(f"{yml}: no ligand block -- a pocket constraint needs one")
        dest = os.path.join(out, f"{name}.txt")
        with open(dest, "w") as fh:
            fh.write(f"SEQUENCE: {sequence}\n")
            fh.write(f"SMILES: {smiles}\n")
            fh.write(f"POCKET: {','.join(str(r) for r in residues)}\n")
            fh.write(f"POCKET_MAX_DISTANCE: {max_distance}\n")
        print(f"{name}: {len(residues)} residues {residues}, max_distance {max_distance} "
              f"-> {os.path.relpath(dest, REPO)}")
    print("\nthese are pocket-guided inputs, not forced: OpenFold3 ranks ligand proposals "
          "against\nthese residues, it does not enforce contact, and the ligand atom of each "
          "hint is dropped.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
