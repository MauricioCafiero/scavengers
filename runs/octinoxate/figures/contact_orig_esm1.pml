# Leu6 backbone O to ligand O26, 2.55 A
# Sticks only, small markers, default orientation -- nothing implying volume.
load orig_esm1.cif, m
hide everything
show sticks, m and not polymer
color red, m and not polymer
select pep_atom, m and resi 6 and name O
select lig_atom, m and not polymer and name O26
set sphere_scale, 0.25
show spheres, pep_atom or lig_atom
color magenta, pep_atom
color yellow, lig_atom
distance dd, pep_atom, lig_atom
show sticks, m and byres (polymer within 5 of (not polymer))
color green, m and polymer
set label_size, 18
reset
zoom (pep_atom or lig_atom), 8
