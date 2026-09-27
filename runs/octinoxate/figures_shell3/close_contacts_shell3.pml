# Shell 3's two folds with positive interaction energy, and the contacts responsible.
#
# Everything this script loads is inside one of two groups, so the object panel shows two lines rather
# than five. Click the triangle on a group to expand it.
#
#   UNPHYSICAL        s3_orig_control  -- a real clash
#   LOOKS_BAD_ISNT    s3_esm1_f8       -- its shortest contact is a hydrogen bond
#
# No semicolons in this file. PyMOL splits commands on a semicolon even inside a comment, so text
# after one gets executed.

cd /Users/cafierom/python_mac/peptidebuilder/runs/octinoxate/figures_shell3

load s3_orig_control.cif, s3_orig_control_POSITIVE
load s3_esm1_f8.cif, s3_esm1_f8_POSITIVE

hide everything
show sticks
set stick_radius, 0.13
bg_color black
color green, polymer
color red, not polymer
set label_size, -0.8
set label_color, yellow
set dash_width, 2.5

# --- s3_orig_control: interaction +0.92, and the reason is unambiguous -------------------------
# TRP32 backbone carbonyl O to ligand O26 at 2.00 A. Both acceptors, no hydrogen between them.
# O...O van der Waals is about 3.0 A and a hydrogen bond needs 2.6-2.8 A plus a donor. Unphysical.
distance CLASH_2p00_backboneO_to_ligandO, s3_orig_control_POSITIVE and resi 32 and name O, s3_orig_control_POSITIVE and not polymer and name O26
color magenta, CLASH_2p00_backboneO_to_ligandO
label s3_orig_control_POSITIVE and resi 32 and name O, "TRP32 O (acceptor)"
label s3_orig_control_POSITIVE and not polymer and name O26, "O26 (acceptor)"

# --- s3_esm1_f8: interaction +2.00, but its shortest contact is not the problem ----------------
# TRP32 NE1 is the indole nitrogen and carries an N-H, so 2.57 A to O26 is a hydrogen bond.
distance HBOND_2p57_indoleNH_to_ligandO, s3_esm1_f8_POSITIVE and resi 32 and name NE1, s3_esm1_f8_POSITIVE and not polymer and name O26
color cyan, HBOND_2p57_indoleNH_to_ligandO
label s3_esm1_f8_POSITIVE and resi 32 and name NE1, "TRP32 NE1 (donor)"

# The actual suspect in that structure: a hydroxyl oxygen 2.60 A from a carbon, where van der Waals
# wants about 3.2 A.
distance CLASH_2p60_tyrOH_to_ligandC, s3_esm1_f8_POSITIVE and resi 1 and name OH, s3_esm1_f8_POSITIVE and not polymer and name C38
color magenta, CLASH_2p60_tyrOH_to_ligandC
label s3_esm1_f8_POSITIVE and resi 1 and name OH, "TYR1 OH"

# --- two groups, so the panel stays readable ---------------------------------------------------
group UNPHYSICAL_orig_control, s3_orig_control_POSITIVE CLASH_2p00_backboneO_to_ligandO
group LOOKS_BAD_ISNT_esm1_f8, s3_esm1_f8_POSITIVE HBOND_2p57_indoleNH_to_ligandO CLASH_2p60_tyrOH_to_ligandC

# start on the unambiguous one, with the other group collapsed and off
disable LOOKS_BAD_ISNT_esm1_f8
orient s3_orig_control_POSITIVE and resi 32
zoom s3_orig_control_POSITIVE and resi 32, 6

print ""
print "=============================================================="
print " UNPHYSICAL_orig_control   s3_orig_control, interaction +0.92"
print "   magenta 2.00 A : TRP32 backbone O to ligand O26"
print "   acceptor to acceptor, no donor. This one is a real clash."
print ""
print " LOOKS_BAD_ISNT_esm1_f8    s3_esm1_f8, interaction +2.00"
print "   cyan    2.57 A : TRP32 indole N-H to ligand O26  = H BOND"
print "   magenta 2.60 A : TYR1 OH to ligand C38           = the suspect"
print ""
print " showing UNPHYSICAL only. To switch:"
print "   disable UNPHYSICAL_orig_control"
print "   enable LOOKS_BAD_ISNT_esm1_f8"
print "   zoom s3_esm1_f8_POSITIVE and (resi 1 or resi 32), 6"
print "=============================================================="
print ""
