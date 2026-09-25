# Arg1 CA to ligand O25, 2.56 A
# Sticks only, small markers, default orientation -- nothing implying volume.
load orig_f8.cif, m
hide everything
show sticks, m and not polymer
color red, m and not polymer
select pep_atom, m and resi 1 and name CA
select lig_atom, m and not polymer and name O25
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
