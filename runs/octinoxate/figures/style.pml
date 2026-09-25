# Representation only. No loads, no alignment.
# Run it after loading whatever you want:  @style.pml
#
# Note for anyone editing this file: PyMOL splits commands on a semicolon even inside a comment,
# so text following one is executed as a command. Never put that character in this file.

set auto_show_cartoon, 0
set auto_show_lines, 0

hide everything
show sticks

bg_color black
color green, polymer
color red, not polymer
set stick_radius, 0.15

# yellow reads cleanly against black, green and red at once.
# A negative label_size is in Angstroms rather than screen points, so the text scales with the
# molecule instead of staying the same size as you zoom out.
set label_size, -0.7
set label_color, yellow

# every residue, on CA only, so one label per residue rather than one per atom
label polymer and name CA, resn+resi

# No zoom. Fitting all 33 residues on screen shrinks the molecule, and your current view is
# probably the one you want. Type "zoom polymer" yourself if you do want the whole thing framed.

# Adjustments, to type at the prompt:
#   set label_size, -0.5      smaller, still scaling with zoom
#   set label_size, -1.0      larger, still scaling with zoom
#   set label_color, wheat
#   set label_color, grey70
#   label none
#   label polymer and name CA, resi
