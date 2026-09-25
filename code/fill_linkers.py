"""Replace the glycine linkers in a design with residues ESM2 thinks are probable.

`peptide_builder` uses glycine as a spacer wherever two fragments are too far apart to be
neighbours, which leaves runs of G that are floppy and low-complexity. This masks every
linker position and fills it with an ESM2 masked-language model, most-confident position
first, the way the GenMaskFill repo does it.

Glycine is not in the fragment library, so every G in a design sequence is a linker and the
mask set is unambiguous. The designed residues are never touched, so the fragment poses and
their total interaction energy still apply to the filled sequence.

Filled sequences are appended to the run's sequences.csv, so `boltz_check.py` picks them up
on its next pass and they can be compared with the poly-glycine originals.

Usage:
    python code/fill_linkers.py runs/octinoxate [--sequences S,S] [--limit N] [--dry-run]

ESM2 runs in process (pip install transformers); --genmask-venv falls back to another environment.
"""
import argparse
import csv
import json
import os
import subprocess
import sys
import tempfile

# ESM2 runs in process when transformers is installed; this is only the fallback, and only if set
GENMASK_VENV = os.environ.get("PEPTIDEBUILDER_GENMASK_VENV")
DEFAULT_MODEL = "facebook/esm2_t12_35M_UR50D"

# fills one mask at a time, choosing the position the model is most confident about but
# sampling the residue from the candidates above prob_cutoff, so a flat distribution does not
# collapse to its mode (GenMaskFill keeps every candidate over the cutoff rather than the argmax)
FILL_SCRIPT = r'''
import json, sys, torch
from transformers import AutoTokenizer, AutoModelForMaskedLM

in_path, out_path = sys.argv[1], sys.argv[2]
config = json.load(open(in_path))
model_name, cutoff, n_variants, seed = config["model"], config["cutoff"], config["variants"], config["seed"]

tokenizer = AutoTokenizer.from_pretrained(model_name)
model = AutoModelForMaskedLM.from_pretrained(model_name)
model.eval()

letters = "ACDEFGHIKLMNPQRSTVWY"
allowed = [tokenizer.convert_tokens_to_ids(a) for a in letters]
generator = torch.Generator().manual_seed(seed)

results = []
for sequence in config["sequences"]:
    for variant in range(n_variants):
        chars = list(sequence)
        todo = {i for i, c in enumerate(chars) if c == "G"}
        steps = []
        while todo:
            text = "".join("<mask>" if i in todo else c for i, c in enumerate(chars))
            enc = tokenizer(text, return_tensors="pt")
            with torch.no_grad():
                logits = model(**enc).logits[0]
            mask_rows = (enc["input_ids"][0] == tokenizer.mask_token_id).nonzero().flatten().tolist()
            order = sorted(todo)
            best = None
            for row, pos in zip(mask_rows, order):
                probs = torch.softmax(logits[row], dim=-1)[allowed]
                if best is None or float(probs.max()) > best[0]:
                    best = (float(probs.max()), pos, probs)
            confidence, pos, probs = best
            # sample among the candidates clearing the cutoff; fall back to the best one
            keep = probs >= cutoff
            if keep.sum() == 0:
                pick = int(probs.argmax())
            else:
                masked = torch.where(keep, probs, torch.zeros_like(probs))
                pick = int(torch.multinomial(masked / masked.sum(), 1, generator=generator))
            chars[pos] = letters[pick]
            todo.discard(pos)
            steps.append({"position": pos, "residue": letters[pick],
                          "probability": round(float(probs[pick]), 4),
                          "n_candidates": int(keep.sum())})
        results.append({"parent": sequence, "variant": variant, "filled": "".join(chars), "steps": steps})

json.dump(results, open(out_path, "w"))
'''


def fill(sequences, genmask_venv=None, model=DEFAULT_MODEL,
         cutoff=0.05, variants=1, seed=1):
    """Fill the glycine linkers of each sequence; returns a list of result dicts.

    transformers is an ordinary pip install, so this normally runs in process. `genmask_venv` is
    the original route -- the same script in another environment's interpreter -- kept as a
    fallback. Both run the identical script text, in process by exec rather than a second copy of
    the algorithm, so the two routes cannot drift apart and the same seed gives the same fills.
    """
    with tempfile.TemporaryDirectory() as tmp:
        in_path = os.path.join(tmp, "in.json")
        out_path = os.path.join(tmp, "out.json")
        script = os.path.join(tmp, "fill.py")
        json.dump({"sequences": list(sequences), "model": model, "cutoff": cutoff,
                   "variants": variants, "seed": seed}, open(in_path, "w"))
        open(script, "w").write(FILL_SCRIPT)
        try:
            import transformers  # noqa: F401
        except ImportError:
            if not genmask_venv:
                raise RuntimeError(
                    "transformers is not installed here and no --genmask-venv was given. Either\n"
                    "  pip install transformers    (or: uv pip install transformers)\n"
                    "or point --genmask-venv at an environment that has it.")
            proc = subprocess.run(
                [os.path.join(genmask_venv, "bin", "python"), script, in_path, out_path],
                capture_output=True, text=True)
            if not os.path.exists(out_path):
                raise RuntimeError(f"ESM2 fill failed:\n{proc.stdout[-2000:]}\n{proc.stderr[-2000:]}")
        else:
            saved = sys.argv
            sys.argv = [script, in_path, out_path]
            try:
                exec(compile(FILL_SCRIPT, "fill_linkers:FILL_SCRIPT", "exec"),
                     {"__name__": "__main__"})
            finally:
                sys.argv = saved
        return json.load(open(out_path))


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("outdir", help="a peptide_builder run directory, e.g. runs/octinoxate")
    parser.add_argument("--sequences", help="comma-separated sequences to fill (default: every design with a linker)")
    parser.add_argument("--limit", type=int, help="only the N best-scoring designs")
    parser.add_argument("--dry-run", action="store_true", help="show the filled sequences without writing them")
    parser.add_argument("--variants", type=int, default=1, help="filled sequences per design (default: 1)")
    parser.add_argument("--prob-cutoff", type=float, default=0.05, help="minimum candidate probability to sample from (default: 0.05)")
    parser.add_argument("--seed", type=int, default=1, help="sampling seed (default: 1)")
    parser.add_argument("--model", default=DEFAULT_MODEL, help=f"HF ESM2 checkpoint (default: {DEFAULT_MODEL})")
    parser.add_argument("--genmask-venv", default=GENMASK_VENV,
                        help="venv with transformers, only needed if it is not installed here "
                             "(or $PEPTIDEBUILDER_GENMASK_VENV)")
    args = parser.parse_args(argv)

    csv_path = os.path.join(args.outdir, "sequences.csv")
    with open(csv_path) as f:
        rows = list(csv.DictReader(f))
        fields = list(rows[0].keys())

    if args.sequences:
        wanted = set(args.sequences.split(","))
        chosen = [r for r in rows if r["sequence"] in wanted or r["sequence"].replace("cyclo-", "") in wanted]
    else:
        chosen = [r for r in rows if "G" in r["sequence"] and "esm2" not in r.get("start_residue", "")]
    chosen.sort(key=lambda r: float(r["total_ie_kcal_mol"]))
    if args.limit:
        chosen = chosen[:args.limit]
    if not chosen:
        sys.exit("no designs with glycine linkers to fill")

    print(f"filling linkers in {len(chosen)} designs with {args.model}\n")
    results = fill([r["sequence"].replace("cyclo-", "") for r in chosen],
                   args.genmask_venv, args.model, args.prob_cutoff, args.variants, args.seed)

    existing = {r["sequence"] for r in rows}
    new_rows = []
    for result in results:
        row = next(r for r in chosen if r["sequence"].replace("cyclo-", "") == result["parent"])
        cyclic = row["sequence"].startswith("cyclo-")
        filled = ("cyclo-" if cyclic else "") + result["filled"]
        n_linkers = len(result["steps"])
        print(f"{row['sequence']}\n  -> {filled}   ({n_linkers} linkers filled, "
              f"mean confidence {sum(s['probability'] for s in result['steps']) / n_linkers:.2f})")
        print(f"     {' '.join(s['residue'] + str(s['position']) for s in result['steps'])}")
        if filled in existing:
            print("     (already present, not added)")
            continue
        new = {k: row.get(k, "") for k in fields}
        new["sequence"] = filled
        new["start_residue"] = f"esm2_filled:{row['start_residue']}"
        new_rows.append(new)
        existing.add(filled)

    if args.dry_run:
        print(f"\ndry run: {len(new_rows)} sequences not written")
        return
    if new_rows:
        with open(csv_path, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=fields)
            writer.writeheader()
            writer.writerows(rows + new_rows)
        print(f"\nadded {len(new_rows)} filled sequences to {csv_path}")
        print("run boltz_check.py to co-fold them")


if __name__ == "__main__":
    main()
