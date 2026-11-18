# render_islands.tcl -- headless VMD -> Tachyon: grey ribbon with TWO residue groups ("islands")
# highlighted in two colors, for the sequence-segments figure (two stretches of sequence that meet
# in the fold).
#
#   vmd -dispdev text -e render_islands.tcl -args IN.pdb OUT.dat CHAIN ROTX ROTY  r1 r2 ... -- s1 s2 ...
#   (residues before "--" = island 1 (red); after "--" = island 2 (orange))
set pdb   [lindex $argv 0]
set out   [lindex $argv 1]
set chain [lindex $argv 2]
set rotx  [lindex $argv 3]
set roty  [lindex $argv 4]
set sep   [lsearch $argv "--"]
set isl1  [lrange $argv 5 [expr {$sep - 1}]]
set isl2  [lrange $argv [expr {$sep + 1}] end]

mol new $pdb type pdb waitfor all
mol delrep 0 top

mol representation NewCartoon 0.30 16 4.5
mol material AOChalky
mol color ColorID 2
mol selection "chain $chain and protein"
mol addrep top

mol representation Licorice 0.3 12 12
mol material AOShiny
mol color ColorID 1
mol selection "chain $chain and protein and resid [join $isl1 { }]"
mol addrep top

mol representation Licorice 0.3 12 12
mol material AOShiny
mol color ColorID 3
mol selection "chain $chain and protein and resid [join $isl2 { }]"
mol addrep top

catch {display projection Orthographic}
catch {axes location off}
catch {color Display Background white}
catch {display ambientocclusion on}
catch {display shadows on}
catch {display depthcue off}

set sel [atomselect top "chain $chain and protein"]
$sel global
$sel move [transaxis x $rotx]
$sel move [transaxis y $roty]
catch {display resetview}

render Tachyon $out
quit
