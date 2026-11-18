# render_scaffold_antibody.tcl -- headless VMD -> Tachyon scene for the homepage "tested for
# antibody binding" panel. The designed scaffold and the native antibody share one coordinate
# frame (RFdiffusion3 holds the epitope motif in native coordinates), so the antibody that binds
# the native epitope also sits on the scaffold's epitope -- a real co-registered render in the same
# ribbon art style as the antigen and scaffold panels.
#
# Usage:
#   vmd -dispdev text -e render_scaffold_antibody.tcl -args DESIGN.pdb ANTIBODY.pdb OUT.dat ROTX ROTY [epitope resids...]
set argv $argv
set design [lindex $argv 0]
set abpdb  [lindex $argv 1]
set out    [lindex $argv 2]
set rotx   [lindex $argv 3]
set roty   [lindex $argv 4]
set resids [lrange $argv 5 end]

# ---- designed scaffold: grey ribbon + red epitope (ribbon + licorice) ----
mol new $design type pdb waitfor all
set md [molinfo top]
mol delrep 0 $md
mol representation NewCartoon 0.30 16 4.5
mol material AOChalky
mol color ColorID 2
mol selection "protein"
mol addrep $md
if {[llength $resids] > 0} {
  set rl [join $resids " "]
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
}

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

# rotate BOTH molecules by the same matrix (keep them registered), then frame the pair
foreach m [list $md $ma] {
  set s [atomselect $m "all"]
  $s move [transaxis x $rotx]
  $s move [transaxis y $roty]
  $s delete
}
catch {display resetview}

render Tachyon $out
quit
