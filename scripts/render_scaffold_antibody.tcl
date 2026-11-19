# render_scaffold_antibody.tcl -- headless VMD -> Tachyon scene for the homepage "antibody
# engaging the epitope" panel. The designed scaffold and the native antibody share one coordinate
# frame (RFdiffusion3 holds the epitope motif in native coordinates), so the antibody that binds
# the native epitope also sits on the scaffold's epitope -- a real co-registered render in the same
# ribbon art style as the antigen and scaffold panels.
#
# Auto-orients the pair: the epitope's long axis is laid flat along screen x, and the direction from
# the epitope to the antibody is stood up along screen +y, so the antibody reads as an upright Y
# above a mostly-flat epitope with the scaffold below. A final ROLL (degrees about the vertical
# axis) spins the antibody about that upright axis to face its Y toward the camera.
#
# Usage:
#   vmd -dispdev text -e render_scaffold_antibody.tcl -args DESIGN.pdb ANTIBODY.pdb OUT.dat ROLL [epitope resids...]
set argv $argv
set design [lindex $argv 0]
set abpdb  [lindex $argv 1]
set out    [lindex $argv 2]
set roll   [lindex $argv 3]
set resids [lrange $argv 4 end]
set rl [join $resids " "]

# ---- designed scaffold: grey ribbon + red epitope (ribbon + licorice) ----
mol new $design type pdb waitfor all
set md [molinfo top]
mol delrep 0 $md
mol representation NewCartoon 0.30 16 4.5
mol material AOChalky
mol color ColorID 2
mol selection "protein"
mol addrep $md
mol representation NewCartoon 0.36 16 4.5
mol material AOChalky
mol color ColorID 1
mol selection "protein and resid $rl"
mol addrep $md
mol representation Licorice 0.3 12 12
mol material AOShiny
mol color ColorID 1
mol selection "protein and resid $rl"
mol addrep $md

# ---- antibody: steel-blue ribbon (its own file holds only the antibody chain) ----
mol new $abpdb type pdb waitfor all
set ma [molinfo top]
mol delrep 0 $ma
color change rgb 20 0.34 0.44 0.60
mol representation NewCartoon 0.30 16 4.5
mol material AOChalky
mol color ColorID 20
mol selection "protein"
mol addrep $ma

catch {display projection Orthographic}
catch {axes location off}
catch {color Display Background white}
catch {display ambientocclusion on}
catch {display shadows on}
catch {display depthcue off}

# ---- build the orientation basis ----
# epitope CA centroid (E) on the design, antibody centroid (A)
set epi [atomselect $md "protein and resid $rl and name CA"]
set E [measure center $epi]
set abs [atomselect $ma "protein"]
set A [measure center $abs]

# epitope long axis = vector between the two most-distant epitope CAs (robust, no inertia API quirks)
set coords [$epi get {x y z}]
set n [llength $coords]
set best -1.0
set eaxis {1 0 0}
for {set i 0} {$i < $n} {incr i} {
  for {set j [expr {$i+1}]} {$j < $n} {incr j} {
    set v [vecsub [lindex $coords $j] [lindex $coords $i]]
    set d [veclength $v]
    if {$d > $best} {set best $d; set eaxis $v}
  }
}
set eaxis [vecnorm $eaxis]

# up = epitope -> antibody (stands the antibody vertically above the epitope)
set up [vecnorm [vecsub $A $E]]
# right = epitope long axis made perpendicular to up (so the epitope lies flat along x)
set d [vecdot $eaxis $up]
set right [vecnorm [vecsub $eaxis [vecscale $d $up]]]
# view (z, toward camera) and re-orthonormalized up
set view [vecnorm [veccross $right $up]]
set up   [vecnorm [veccross $view $right]]

# change-of-basis matrix (rows right, up, view): maps right->x, up->y, view->z
set m [list \
  [list [lindex $right 0] [lindex $right 1] [lindex $right 2] 0] \
  [list [lindex $up 0]    [lindex $up 1]    [lindex $up 2]    0] \
  [list [lindex $view 0]  [lindex $view 1]  [lindex $view 2]  0] \
  [list 0 0 0 1]]

# apply the same matrix to both molecules (keeps them registered), then roll about the vertical axis
foreach mm [list $md $ma] {
  set s [atomselect $mm "all"]
  $s move $m
  $s move [transaxis y $roll]
  $s delete
}
catch {display resetview}

render Tachyon $out
quit
