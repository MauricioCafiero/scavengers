# Peptide strain against a Boltz apo fold — attempted 2026-09-26, does not work

Archived, not deleted: the code runs and the negative result is worth keeping. **Nothing here should be
quoted as a strain measurement.**

## The idea

Peptide strain was disabled in the main pipeline because its reference was meaningless: relaxing a free
peptide in vacuum collapses it into a compact hydrogen-bonded ball, giving 237.8 kcal/mol on a 33-mer
and not converging, so it measured collapse rather than strain.

The proposal was to replace that reference with the peptide's **apo fold** — Boltz's prediction for the
sequence with no ligand — hydrogens relaxed with heavy atoms fixed, exactly as the bound state is
treated. Then `peptide strain = E(peptide in the complex) − E(apo fold)`, one reference per sequence,
shared across that sequence's four constraint levels.

## Why it fails

Tested on shell 3's twelve folds. Two independent failure modes, both fatal.

**Vacuum still rewards compactness — the original disease, relocated.** `esm2`'s apo fold is an extended
100% helix at Rg 15.2 Å. Its forced folds collapse around the ligand to Rg 9.7–10.3 and come out
107–152 kcal/mol *lower* in energy, giving large negative "strain". Its control, which stays extended at
Rg 15.0, gives +15.9. The number tracks compaction, not conformational cost.

**The apo reference is only confident where it is useless.** Boltz's apo confidence tracks helicity
exactly:

| sequence | apo pTM | apo pLDDT | apo helix | bound helix |
|---|---|---|---|---|
| `orig` (glycine design) | 0.205 | 0.501 | 34% | 22/3/19/19% |
| `esm1` | 0.293 | 0.630 | 69% | 59/72/69/72% |
| `esm2` | 0.593 | 0.903 | 100% | 100/91/91/88% |

So the sequence with a trustworthy apo fold is the one whose fold barely changes on binding, where
strain is near zero and uninformative. The glycine designs, where a fold change would be interesting,
are intrinsically disordered — 64% glycine, no single apo structure to reference — and Boltz even gives
the free chain *more* helix (34%) than any of its bound folds, the signature of a guess. Their strains
came out +5.6, +82.2, +162.1 and **+478.6** kcal/mol at essentially constant Rg (7.8–8.1), i.e. noise.

Full results as measured: `runs/octinoxate/boltz/strain_peptide_shell3.csv`, range −152 to +479 kcal/mol,
8 of 12 negative.

## What would be needed

Implicit solvent. GB cancels the spurious vacuum preference for compact structures, which is the
dominant error term above, and a force field with minimisation is far less sensitive to the ångström-
scale coordinate error in a prediction than an MLIP single point is. That is MM/GBSA, and the machinery
exists in `~/python_mac/openmm` (`omd mmgbsa`, AmberTools present in the `openmm-md` conda env). Peptide
strain should be attempted there, as a minimised-with-GB comparison, not as a raw single-point
difference.

## What is here

| path | what |
|---|---|
| `code/peptide_apo.py` | `fold` (Boltz, no ligand, local) and `energy` (protonate, relax hydrogens only, single point) |
| `code/modal_peptide_ref.py` | runs the `energy` step on a Modal GPU |
| `code/strain_peptide.py` | `E(peptide in complex) − E(apo)`, from the peptide portion of `_complex_relaxed_h.xyz` |
| `runs/octinoxate/apo/` | the three shell-3 apo folds, with Boltz confidences |
| `runs/octinoxate/peptide_reference.json` | the three apo energies |
| `runs/octinoxate/boltz/strain_peptide_shell3.csv` | the twelve results — evidence, not measurements |
| `runs/octinoxate/apo_shell3_*/` | Modal staging, including one failed run |

The ligand-side equivalent — `code/ligand_reference.py` and `code/strain_global.py` — **does** work and
stays in the main pipeline. The difference is that the ligand is 44 atoms with a genuine relaxed global
minimum as its reference, and neither condition holds for a 34-residue peptide.
