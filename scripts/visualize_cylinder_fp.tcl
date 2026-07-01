# visualize_cylinder_fp.tcl -- overlay a design, the real antibody, the cylinder, and the
# cylinder-flagged scaffold CAs, to see why the flags are not real steric clashes.
#
# Run from a cylinder_fp_probe.py output directory (which holds design.pdb, epitope_cas.pdb,
# antibody_aligned.pdb, flagged_cas.pdb, cylinder_frame.txt):
#     cd results/cylinder_fp/DP2_0804
#     vmd -e /path/to/episcaf/scripts/visualize_cylinder_fp.tcl
#
# What you should see: the antibody (blue) sits ABOVE the epitope in the outer part of the
# cylinder; the flagged scaffold CAs (orange) sit at the cylinder BASE, hugging the epitope,
# well clear of the antibody -- in the cylinder volume but not obstructing the paratope.

# --- read the cylinder frame written by the probe ---
set fp [open "cylinder_frame.txt" r]
while {[gets $fp line] >= 0} {
    set tok [split [string trim $line]]
    switch [lindex $tok 0] {
        base   { set cx [lindex $tok 1]; set cy [lindex $tok 2]; set cz [lindex $tok 3] }
        normal { set nx [lindex $tok 1]; set ny [lindex $tok 2]; set nz [lindex $tok 3] }
        R      { set R [lindex $tok 1] }
        H      { set H [lindex $tok 1] }
    }
}
close $fp
set nmag [expr {sqrt($nx*$nx+$ny*$ny+$nz*$nz)}]
set nx [expr {$nx/$nmag}]; set ny [expr {$ny/$nmag}]; set nz [expr {$nz/$nmag}]

# --- molecules ---
mol new design.pdb type pdb waitfor all
set d [molinfo top get id]; mol delrep 0 $d
mol addrep $d; mol modselect 0 $d {all}
mol modstyle 0 $d NewRibbons 0.3 10 1.5; mol modcolor 0 $d ColorID 6; mol modmaterial 0 $d AOChalky

mol new epitope_cas.pdb type pdb waitfor all
set e [molinfo top get id]; mol delrep 0 $e
mol addrep $e; mol modselect 0 $e {all}; mol modstyle 0 $e VDW 0.9 16; mol modcolor 0 $e ColorID 1

mol new flagged_cas.pdb type pdb waitfor all
set f [molinfo top get id]; mol delrep 0 $f
mol addrep $f; mol modselect 0 $f {all}; mol modstyle 0 $f VDW 0.9 16; mol modcolor 0 $f ColorID 3

mol new antibody_aligned.pdb type pdb waitfor all
set a [molinfo top get id]; mol delrep 0 $a
mol addrep $a; mol modselect 0 $a {all}
mol modstyle 0 $a QuickSurf 1.2 0.7 1.0 1.0; mol modcolor 0 $a ColorID 0
mol modmaterial 0 $a Transparent   ; # antibody as a translucent blue volume, so the
                                     # base gap below it is obvious

display projection Orthographic
display depthcue on
axes location off
color Display Background white
display resetview

# --- cylinder basis (two perpendicular vectors) ---
if {abs($nx) < 0.9} { set ax 1.0; set ay 0.0; set az 0.0 } else { set ax 0.0; set ay 1.0; set az 0.0 }
set p1x [expr {$ny*$az-$nz*$ay}]; set p1y [expr {$nz*$ax-$nx*$az}]; set p1z [expr {$nx*$ay-$ny*$ax}]
set p1mag [expr {sqrt($p1x*$p1x+$p1y*$p1y+$p1z*$p1z)}]
set p1x [expr {$p1x/$p1mag}]; set p1y [expr {$p1y/$p1mag}]; set p1z [expr {$p1z/$p1mag}]
set p2x [expr {$ny*$p1z-$nz*$p1y}]; set p2y [expr {$nz*$p1x-$nx*$p1z}]; set p2z [expr {$nx*$p1y-$ny*$p1x}]

# --- cyan cylinder mesh (base ring, top ring at H, verticals) ---
draw color cyan
set ntheta 24
for {set k 0} {$k < $ntheta} {incr k} {
    set t0 [expr {2*3.14159*$k/$ntheta}]; set t1 [expr {2*3.14159*($k+1)/$ntheta}]
    set bx0 [expr {$cx+$R*(cos($t0)*$p1x+sin($t0)*$p2x)}]
    set by0 [expr {$cy+$R*(cos($t0)*$p1y+sin($t0)*$p2y)}]
    set bz0 [expr {$cz+$R*(cos($t0)*$p1z+sin($t0)*$p2z)}]
    set bx1 [expr {$cx+$R*(cos($t1)*$p1x+sin($t1)*$p2x)}]
    set by1 [expr {$cy+$R*(cos($t1)*$p1y+sin($t1)*$p2y)}]
    set bz1 [expr {$cz+$R*(cos($t1)*$p1z+sin($t1)*$p2z)}]
    set tx0 [expr {$bx0+$H*$nx}]; set ty0 [expr {$by0+$H*$ny}]; set tz0 [expr {$bz0+$H*$nz}]
    set tx1 [expr {$bx1+$H*$nx}]; set ty1 [expr {$by1+$H*$ny}]; set tz1 [expr {$bz1+$H*$nz}]
    draw line "$bx0 $by0 $bz0" "$bx1 $by1 $bz1" width 2
    draw line "$tx0 $ty0 $tz0" "$tx1 $ty1 $tz1" width 2
    if {[expr {$k%6}]==0} { draw line "$bx0 $by0 $bz0" "$tx0 $ty0 $tz0" width 1 }
}

# --- blue approach-normal arrow ---
draw color blue
set arx [expr {$cx+$nx*($H+4)}]; set ary [expr {$cy+$ny*($H+4)}]; set arz [expr {$cz+$nz*($H+4)}]
set tipx [expr {$cx+$nx*($H+8)}]; set tipy [expr {$cy+$ny*($H+8)}]; set tipz [expr {$cz+$nz*($H+8)}]
draw cylinder "$cx $cy $cz" "$arx $ary $arz" radius 0.4
draw cone "$arx $ary $arz" "$tipx $tipy $tipz" radius 1.2

puts "silver ribbon = design scaffold | red = epitope | orange = cylinder-flagged scaffold CAs"
puts "blue points  = real antibody (from the AbDb complex, aligned) | cyan mesh = the cylinder"
puts "look: orange sits at the cylinder base near the epitope; blue antibody sits above it."
