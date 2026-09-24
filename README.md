# PeptideBuilder

A computational tool for designing and optimizing peptide structures through fragment-based molecular building and energy calculations using machine learning interatomic potentials.

## Overview

PeptideBuilder is a Python-based system that combines fragment-based molecular assembly with machine learning energy predictions to generate optimized peptide structures. The tool allows users to:

- Define amino acid fragments and ligand binding sites
- Randomly place and rotate fragments within a specified binding region
- Calculate interaction energies using the FAIRChem UMA-OMOL model
- Identify optimal fragment poses with spatial constraints
- Combine best poses to generate complete peptide sequences
- Estimate linker requirements between fragments

## Features

- **Fragment Library**: Pre-defined amino acid fragments including arginine, lysine, aspartic acid, glutamic acid, isoleucine, leucine, serine, tryptophan, tyrosine, and phenylalanine
- **Ligand Support**: Includes UV filter molecules (octinoxate, octocrylene, oxybenzone, avobenzone)
- **ML-Based Energy Calculation**: Uses FAIRChem's UMA-OMOL pre-trained model for accurate force field predictions
- **Pose Optimization**: BFGS geometry optimization for fragment placement
- **Spatial Analysis**: Calculates distances between fragment centers and suggests linker placement
- **3D Visualization**: Integration with py3Dmol for interactive molecular visualization
- **Flexible Configuration**: Supports customizable vdW distances, number of placement attempts, and charge/spin states

## Installation

Requires Python 3.10–3.13 (fairchem does not support 3.14 yet).

```bash
git clone https://github.com/lucia-71/peptidebuilder.git
cd peptidebuilder
uv venv --python 3.12 .venv
uv pip install --python .venv/bin/python -r requirements.txt
source .venv/bin/activate
```

The UMA weights are gated on Hugging Face: request access to `facebook/UMA` and log in once (`hf auth login`). The checkpoint is downloaded to `~/.cache/fairchem` on first use.

## Usage

```bash
python code/peptide_builder.py list                                   # built-in ligands and fragment library
python code/peptide_builder.py run octinoxate                         # built-in ligand, all fragments, 100 tries each
python code/peptide_builder.py run data/xyz_files/avobenzone.xyz -f R,K,W,Y -n 200 --seed 1
python code/peptide_builder.py run "CCCCC(CC)COC(=O)/C=C/c1ccc(OC)cc1" --name octinoxate_smi
python code/peptide_builder.py run octinoxate --all-starts --copies 2      # many sequences in one go
python code/peptide_builder.py combine runs/octinoxate --all-starts       # rebuild from saved poses, no sampling
```

The ligand can be a built-in name, an `.xyz` file, or a SMILES string. SMILES are embedded with RDKit (ETKDGv3, 5 conformers, MMFF minimisation with UFF fallback; lowest-energy conformer). XYZ and SMILES ligands are centred on their centroid, as the built-in ligands are.

| Option | Default | Meaning |
|---|---|---|
| `--name` | file stem / `ligand` | ligand name used in output file names |
| `--charge`, `--spin` | 0 / 1 (XYZ); formal charge / radicals + 1 (SMILES) | ligand charge and spin multiplicity |
| `-f, --fragments` | all | comma-separated names or one-letter codes |
| `-n, --tries` | 100 | placement attempts per fragment |
| `--box-pad` | 0.0 | distance (Å) added to the ligand bounding box |
| `--energy-cutoff` | 500 | unrelaxed IE (kcal/mol) below which a placement is relaxed and kept |
| `--clash` | 1.4 | minimum distance (Å) between atoms of selected fragments |
| `--copies` | 1 | how many times each residue may be placed |
| `--all-starts` | off | repeat the combination with each residue picking first, saving every unique sequence |
| `--cyclic` | off | close the chain into a ring, adding linkers for the last-to-first gap |
| `--linker-slack` | 0 | extra glycines per gap beyond the taut minimum, for backbone flexibility |
| `--model` | `uma-s-1p2p1` | fairchem model name (e.g. `uma-s-1p1`, `uma-m-1p1`) or checkpoint `.pt` path |
| `--device` | auto | `cpu`, `cuda`, or auto |
| `--seed` | none | random seed for reproducible sampling |
| `-o, --outdir` | `runs/<ligand>` | output directory |

### Pipeline

1. The search region is the ligand's bounding box (plus `--box-pad`).
2. For each fragment, `--tries` random rotations are placed at Gaussian-sampled origins around the box centre and scored with UMA-OMOL: IE = E(complex) − E(ligand) − E(fragment), in kcal/mol.
3. Placements under `--energy-cutoff` are relaxed with BFGS (fmax 0.05 eV/Å, ligand fixed), re-scored and kept.
4. Fragments are taken in library order. Each gets its lowest-IE pose that doesn't clash with fragments already selected. With `--copies N` the walk repeats N times, so a residue can be placed more than once (a different pose each time).
5. The selected fragments are chained from the first one, always stepping to the nearest unvisited centre. Every gap gets `round(d / 3.8) − 1` glycine linkers, which gives the one-letter sequence.
6. With `--all-starts` steps 4–5 are repeated with the order rotated so each residue in turn picks first. Because the first picker takes the best site, each start gives a different peptide. Every distinct sequence is saved; if two starts give the same sequence, the lower total IE is kept.

### Output Files

- `<outdir>/poses/<ligand>_w_<fragment><index>.xyz`: every kept pose, with IE in the comment line
- `<outdir>/combined.xyz`: ligand plus the selected pose of each placed fragment
- `<outdir>/energies.csv`: fragment, pose, IE (kcal/mol), and whether the pose was selected
- `<outdir>/sequence.txt`: the best peptide sequence and its fragment-to-fragment distances
- `<outdir>/sequences/<sequence>.xyz`: one structure per unique sequence
- `<outdir>/sequences.csv`: each unique sequence with its starting residue, fragment count and total IE
- `<outdir>/ligand.xyz`: the (centred) ligand geometry used

### Re-combining without re-sampling

Sampling is the expensive step and its poses are saved as they are found, so alternative peptides cost seconds:

```bash
python code/peptide_builder.py combine runs/octinoxate --all-starts --copies 2
```

`combine` reads `poses/*.xyz` and `energies.csv` from an earlier run and redoes only the selection, combination and sequence steps.

### Validating a design against Boltz-2

Three scripts check whether the designs survive contact with a folded structure.

**Requires the `boltz_local` repo**, a sibling of this one, which provides Boltz-2 co-folding on
Apple silicon through PyTorch's MPS backend (adapted from the
[CafChem](https://github.com/MauricioCafiero/CafChem) Boltz tools; Boltz-2 itself is
[jwohlwend/boltz](https://github.com/jwohlwend/boltz)). Set it up per its own README, then pass its
path with `--boltz-repo` (default `~/python_mac/boltz_local`). These scripts call its
`code/boltz_mps.py`, which runs the structure and affinity models in separate processes so the two
do not have to share 8 GB of unified memory. Nothing in that repo is modified.

`uma_binding.py` additionally needs pdbfixer, taken from any environment that has it via
`--fixer-venv` (default `~/python_mac/pocket_assist/venv`), so this repo needs no extra packages.

```bash
python code/boltz_check.py runs/octinoxate              # co-fold every sequence with the ligand
python code/overlay.py runs/octinoxate --save-overlays  # superpose the fold on the design
python code/uma_binding.py runs/octinoxate              # UMA interaction energy of each complex
```

- **`boltz_check.py`** writes one Boltz YAML per sequence and co-folds it with the ligand. The
  ligand SMILES is perceived from the run's own `ligand.xyz`, so it always matches the structure
  the fragments were sampled against. MSAs are skipped (`msa: empty`): these are designed peptides
  with no natural homologs, and it keeps them off the public MSA server. `cyclo-` sequences are
  folded as genuine rings (`cyclic: true`). Results go to `boltz/comparison.csv`.
- **`overlay.py`** superposes each folded complex onto its design using the ligand's rigid core
  (the flexible tail and the ring's symmetry are both handled), then measures how far each folded
  side chain sits from the fragment that designed it. Results go to `boltz/overlay.csv`.
- **`uma_binding.py`** recomputes each complex's interaction energy with UMA as
  E(complex) - E(peptide) - E(ligand) at a fixed geometry, so strain cancels. Boltz writes heavy
  atoms only, so hydrogens are added first: the peptide with pdbfixer at pH 7 (called in a separate
  environment via `--fixer-venv`, so nothing extra is needed here), the ligand with RDKit from its
  SMILES. Results go to `boltz/uma_binding.csv`.

#### Controls and linker substitution

```bash
python code/fill_linkers.py runs/octinoxate --variants 2      # replace glycine linkers with ESM2
python code/random_control.py runs/octinoxate --length 19     # 5 random peptides of the same length
python code/random_control.py runs/octinoxate --length 19 --shuffle RYGLSGIKWDKSFEGDGGE
```

- **`fill_linkers.py`** masks every glycine linker and fills it with an ESM2 masked-language
  model, most-confident position first, the way the
  [GenMaskFill](https://github.com/MauricioCafiero/CafChem) repo does it. Glycine is not in the
  fragment library, so every G in a design is a linker and the mask set is unambiguous; the designed
  residues are never touched, so the fragment poses and their total IE still apply. Residues are
  sampled from the candidates above `--prob-cutoff` rather than taken as the argmax: with this many
  masks the 35M checkpoint's per-position distribution is nearly flat (confidence 0.06-0.10), and
  the argmax collapses every linker to leucine. Filled sequences are appended to `sequences.csv`,
  so `boltz_check.py` picks them up. ESM2 runs in the GenMaskFill environment (`--genmask-venv`).
- **`random_control.py`** is the null model: random peptides of a given length, or with
  `--shuffle`, permutations of one design that hold length and composition fixed and vary only the
  order. Each goes through the same Boltz-2 and UMA path as the designs, and the results are
  compared with the designs of that length. Results go to `boltz/controls.csv`.

#### What the validation showed (octinoxate, 20 designs, uma-s-1p2p1)

A folded peptide does not reproduce the designed shell. The fragment search puts every side chain
in contact with the ligand; a real backbone can only reach a fraction of those positions at once.

- The ligand superposes well in every case (RMSD 0.07-0.86 A, mean 0.54 A), so the comparison frame
  is sound, but only **8 of 225 fragments (3.6%)** are reproduced within 3 A. Median displacement is
  12.0 A for the same residue type, 6.3 A for the nearest side chain of any type: the fold packs
  beside the ligand rather than around it.
- Predicted affinity does not track the fragment sum. Across 23 sequences the rank correlation
  between total fragment IE and Boltz pIC50 is **-0.42** — mildly inverted. Sequence length alone
  correlates better with predicted dG than the fragment energy does.
- The fragment sum is a sum of rigid gas-phase interaction energies with no desolvation or entropy,
  so its magnitude (-54 to -122 kcal/mol) is not comparable to a binding free energy (Boltz gives
  about -8 to -10 kcal/mol). Only its ranking could carry information.
- Boltz's affinity head reports confident numbers for complexes that are not in contact: the
  12-residue `RYGLGEFSDWKI` has its ligand 6.6 A away, yet scored pIC50 6.96 and binder probability
  0.75. Always read `overlay.csv` contact distances alongside an affinity.

Substituting the glycine linkers helps, consistently but modestly. For four 19-residue designs the
substitution improved Boltz's dG in 4 of 4 cases, by 0.01 to 0.38 kcal/mol (mean 0.22), and raised
the binder probability in 3 of 4. The design whose substitution changed fewest positions moved
least. Boltz's affinity head is nearly blind to the change in itself - one substituted sequence
differed from its parent by 0.003 log units of pIC50 despite 5 of 19 residues changing - so the
effect shows up through the fold rather than through sequence composition.

Against a null model there is no evidence the fragment search adds value. At length 19 the four
designs give UMA interaction energies of -16.2, -12.5, -10.1 and -7.2 kcal/mol (mean -11.5); five
random sequences of the same length give -40.8, -9.5, -7.7, -5.1 and +12.1 (mean -10.2). The means
are indistinguishable and the best random peptide is 2.5 times better than the best design. The
designs are more consistent (all four physical, in a narrow band) while random sampling has much
higher variance. n is far too small for significance, but there is no hint of an advantage either.

Two explanations for that were tested and refuted. The best random binder is more charged than the
designs (37% vs 23%), so a bias toward charged residues is not the cause. And the interaction
energy does not track contact area: Spearman +0.17 against contact pairs and +0.42 against buried
ligand atoms, the wrong sign, with two structures of near-identical burial differing by 35 kcal/mol.
At fixed length the energy is set by contact quality rather than quantity; what drives that is
still unexplained.

`--linker-slack` tests the flexibility limit. With the designed geometry and fragment IE held fixed
at -100.87 kcal/mol and only glycine added, one extra linker per gap clearly helps the fold reach
the ligand, but more does not:

| slack | residues | median displacement (any type) | pIC50 | binder prob |
|---|---|---|---|---|
| 0 | 12 | 11.6 A | 6.96 | 0.75 |
| 1 | 21 | **6.8 A** | 6.23 | 0.53 |
| 2 | 30 | 6.7 A | 6.02 | 0.17 |
| 3 | 39 | 8.5 A | 5.97 | 0.16 |

So slack has an optimum near 1 rather than being a monotonic lever, and predicted affinity falls
steadily as backbone is added. No setting reproduced the designed positions (0/10 within 3 A at
every level): flexibility is necessary but not sufficient, because nothing in the sequence
specifies the arrangement. Cyclisation is the constraint-based alternative (`--cyclic`), since a
ring restricts the backbone instead of merely loosening it.

### Python API and viewing

```python
import sys; sys.path.insert(0, "../code")   # from notebooks/
from peptide_builder import load_ligand, view_combined_poses, view_frag_pose   # viewers need py3Dmol

ligand = load_ligand("octinoxate")
view_combined_poses("../runs/octinoxate", ligand)
view_frag_pose("../runs/octinoxate", ligand, "serine", 0)
```

## Data Format

### Fragment Dictionary Structure

```python
fragment = {
    "num_atoms": int,           # Total number of atoms
    "name": str,                # Fragment identifier
    "atoms": [str, ...],        # Element symbols
    "coords": [[x, y, z], ...], # 3D coordinates (Angstroms)
    "charge": int,              # Net charge
    "spin": int,                # Spin multiplicity
    "size": float               # Spatial extent (calculated)
}
```

## Supported Molecules

### Amino Acids
Arginine, Lysine, Aspartic Acid, Glutamic Acid, Isoleucine, Leucine, Serine, Tryptophan, Tyrosine, Phenylalanine

### UV Filters (Example Binding Sites)
- Octinoxate (44 atoms)
- Octocrylene (56 atoms)
- Oxybenzone (29 atoms)
- Avobenzone (45 atoms)

## Computational Requirements

- **GPU**: Recommended (CUDA-compatible for faster calculations)
- **RAM**: Minimum 4GB (8GB+ recommended)
- **CPU**: Multi-core processor beneficial for parallel processing

## File Organization

```
peptidebuilder/
├── code/
│   └── peptide_builder.py      # ligand/fragment data, pipeline, notebook viewers and command line
├── notebooks/
│   └── frag_grow.ipynb         # original exploratory notebook
├── data/
│   ├── xyz_files/              # ligand and peptide–ligand structures
│   ├── md_videos/              # MD trajectories
│   └── log_files/              # MD logs
├── runs/                       # pipeline output (git-ignored)
├── requirements.txt
└── README.md
```

## Tips for Use

1. **Fragment Sampling**: Start with 100 placement attempts and increase if too few poses are kept
2. **Energy Thresholds**: Adjust `--energy-cutoff` to suit your binding site
3. **Steric Clashes**: Raise `--clash` to keep selected fragments further apart, or lower it if too few fragments get placed

## Limitations

- Designed for peptide/small molecule binding site decoration
- ML predictions valid within FAIRChem training domain
- Fragment orientations limited to random rotations (no directed sampling)
- Interaction energies exclude long-range electrostatics beyond ML cutoff

## Future Enhancements

- Support for custom fragment libraries
- Directed sampling algorithms (e.g., genetic algorithms)
- Integration with protein structure databases
- Machine learning model fine-tuning capability
- Molecular dynamics post-processing

## Citation

If you use PeptideBuilder in your research, please cite:

```
PeptideBuilder: Fragment-based peptide design using machine learning
[Your citation information here]
```

## License

[Specify your license here - e.g., MIT, Apache 2.0, etc.]

## Contact

For questions or issues, please visit the [GitHub repository](https://github.com/lucia-71/peptidebuilder) or open an issue.

## Acknowledgments

- FAIRChem team for the UMA-OMOL ML interatomic potential
- ASE (Atomic Simulation Environment) for geometry optimization
- py3Dmol for molecular visualization

---

**Last Updated**: June 2026
