# render_epitope_fig.tcl -- headless VMD -> Tachyon scene for episcaf web figures.
#
# Renders one chain of a structure as a grey NewCartoon ribbon, with a residue set
# highlighted as a red molecular surface region (the epitope or a candidate patch).
# Produces a Tachyon scene (.dat); ray-trace it with the external tachyon binary.
#
# Usage:
#   vmd -dispdev text -e render_epitope_fig.tcl -args IN.pdb OUT.dat CHAIN HLCOLOR ROTX ROTY [resids...]
#     HLCOLOR: 1=red 0=blue 7=green 11=purple ; use -1 for no highlight
# Consistent camera: same CHAIN/ROTX/ROTY across calls gives the same orientation,
# so candidate-patch panels differ only in which residues are highlighted.
set argv $argv
set pdb   [lindex $argv 0]
set out   [lindex $argv 1]
set chain [lindex $argv 2]
set hlcol [lindex $argv 3]
set rotx  [lindex $argv 4]
set roty  [lindex $argv 5]
set resids [lrange $argv 6 end]

mol new $pdb type pdb waitfor all
mol delrep 0 top

# grey ribbon for the whole chain
mol representation NewCartoon 0.30 16 4.5
mol material AOChalky
mol color ColorID 2
mol selection "chain $chain and protein"
mol addrep top

# highlighted residues: colored ribbon drawn over the grey, plus licorice sidechains
if {$hlcol >= 0 && [llength $resids] > 0} {
    set rl [join $resids " "]
    mol representation NewCartoon 0.36 16 4.5
    mol material AOChalky
    mol color ColorID $hlcol
    mol selection "chain $chain and protein and resid $rl"
    mol addrep top

    mol representation Licorice 0.3 12 12
    mol material AOShiny
    mol color ColorID $hlcol
    mol selection "chain $chain and protein and resid $rl"
    mol addrep top
}

catch {display projection Orthographic}
catch {axes location off}
catch {color Display Background white}
catch {display ambientocclusion on}
catch {display shadows on}
catch {display depthcue off}

# frame on the selected chain
set sel [atomselect top "chain $chain and protein"]
$sel global
catch {display resetview}
catch {rotate x by $rotx}
catch {rotate y by $roty}
catch {display resetview}
# resetview after rotating the view can re-fit; rotate the coordinates instead for stability
set m [transaxis x $rotx]
$sel move $m
set m2 [transaxis y $roty]
$sel move $m2
catch {display resetview}

render Tachyon $out
quit
