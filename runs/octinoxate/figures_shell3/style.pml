# Representation only. No loads, no alignment.
# Run it after loading whatever you want:  @style.pml

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

# No zoom: your current view is probably the one you want. Type "zoom polymer" if you want the
# whole thing framed.
#
# Adjustments, to type at the prompt:
#   set label_size, -0.5      smaller, still scaling with zoom
#   set label_size, -1.0      larger, still scaling with zoom
#   set label_color, wheat
#   label none
