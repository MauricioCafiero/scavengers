# bg33_4: Boltz reference (green) vs Vina poses
# No alignment commands: every object is already in one frame.
load bg33_4_protein.pdb, receptor
load ligand_ref.sdf, boltz_ref
load poses.sdf, vina_poses
hide everything
show cartoon, receptor
color grey80, receptor
show sticks, boltz_ref
color green, boltz_ref
show sticks, vina_poses
color cyan, vina_poses
set all_states, on, vina_poses
set stick_radius, 0.12, vina_poses
set stick_radius, 0.20, boltz_ref
orient boltz_ref
set ray_opaque_background, 0
