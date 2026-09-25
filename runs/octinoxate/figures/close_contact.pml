# The 2.44 A oxygen contact, in orig_f12 -- the fold with interaction -45.86 kcal/mol.
# Gly14's backbone carbonyl O sits 2.44 A from ligand O25 and 2.57 A from C28. Both are O...O
# or O...C with no donor hydrogen available, so they are short without being hydrogen bonds --
# which is why the -45.86 needs checking before it is trusted.
load orig_f12.cif, f12
hide everything
# no cartoon: a ribbon drawn through the backbone sits near the atoms being judged and reads as
# extra bulk. Sticks only, so what is on screen is the atoms themselves.
show sticks, f12 and not polymer
color red, f12 and not polymer

select gly14_O, f12 and resi 14 and name O
select lig_O25, f12 and not polymer and name O25
select lig_C28, f12 and not polymer and name C28
# small markers only. At default scale these are full vdW radii, and any two atoms 2.44 A apart
# will appear to overlap regardless of whether the contact is real -- which prejudges the question.
set sphere_scale, 0.25
show spheres, gly14_O or lig_O25 or lig_C28
color magenta, gly14_O
color yellow, lig_O25 or lig_C28
distance d1, gly14_O, lig_O25
distance d2, gly14_O, lig_C28
show sticks, f12 and byres (polymer within 5 of (not polymer))
color green, f12 and polymer
set label_size, 18
# default axis-aligned orientation, not an auto-oriented one
reset
zoom (gly14_O or lig_O25 or lig_C28), 8
