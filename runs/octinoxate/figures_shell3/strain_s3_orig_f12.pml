# s3_orig_f12: the best-enclosed structure in the project, and net unfavourable.
#
#   enclosed 0.965, wrapped 1.00, 20/20 ligand atoms engaged, centroid separation 3.6 A
#   interaction -25.13, ligand strain 44.05, sum +18.92 kcal/mol
#
# The question this session answers: is that 44 kcal/mol real, or is it the potential misbehaving?
# It looks real. The bound ligand sits 1.35 A heavy-atom RMSD from its own lowest conformer, with four
# rotatable torsions displaced by more than 100 degrees -- one of them by 174, i.e. flipped. The
# reference conformer is loaded superposed so you can see which parts moved.
#
# No semicolons anywhere in this file. PyMOL splits commands on a semicolon even inside a comment, so
# text after one is executed.

cd /Users/cafierom/python_mac/peptidebuilder/runs/octinoxate/figures_shell3

load s3_orig_f12.cif, complex
load s3_orig_f12_ligand_bound.xyz, lig_bound
load s3_orig_f12_ligand_reference_aligned.xyz, lig_relaxed

hide everything
show sticks
set stick_radius, 0.15
bg_color black

# peptide green, the bound ligand red as everywhere else in this project
color green, complex and polymer
color red, complex and not polymer
color red, lig_bound

# the ligand's own lowest conformer, aligned onto the bound pose. Where yellow and red diverge is
# where binding has bent the molecule.
color yellow, lig_relaxed

# start on the comparison rather than the whole complex, since that is the point here
disable complex
zoom lig_bound, 3

# --- what to do next, to type at the prompt ---------------------------------------------------
# enable complex                      put the peptide back, to see how it encloses the ligand
# zoom complex, 5
# disable lig_relaxed                 just the bound pose
# show spheres, lig_bound             enclosure reads better as spheres inside the green cage
# set stick_radius, 0.08, complex     thin the peptide down so the ligand stays visible
#
# The four displaced torsions, by atom index in lig_bound:
#   20-0-1-2     bound  162  reference   60    diff 102
#   3-4-5-6      bound   53  reference  178    diff 125
#   3-4-7-8      bound   72  reference -103    diff 174
#   15-16-17-39  bound -168  reference   63    diff 130
# label them with:  label lig_bound and index 4+5+7+16+17, index
