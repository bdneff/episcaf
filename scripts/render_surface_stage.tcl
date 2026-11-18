# render_surface_stage.tcl -- headless VMD -> Tachyon: a molecular SURFACE of one chain with a
# residue subset highlighted (for the seed-and-grow / directional-sampling doc figures).
#
#   vmd -dispdev text -e render_surface_stage.tcl -args IN.pdb OUT.dat CHAIN HLCOLOR ROTX ROTY [resids...]
#
# The whole chain is drawn as a grey QuickSurf (a smooth molecular surface); the given residues are
# drawn as colored VDW spheres that sit proud of the surface, so they read as residues on a curved
# surface. Same CHAIN/ROTX/ROTY across calls gives the same camera, so growth stages differ only in
# which residues are highlighted.
set pdb   [lindex $argv 0]
set out   [lindex $argv 1]
set chain [lindex $argv 2]
set hlcol [lindex $argv 3]
set rotx  [lindex $argv 4]
set roty  [lindex $argv 5]
set resids [lrange $argv 6 end]

mol new $pdb type pdb waitfor all
mol delrep 0 top

set hassel ""
if {$hlcol >= 0 && [llength $resids] > 0} {
    set rl [join $resids " "]
    set hassel "chain $chain and protein and resid $rl"
}

# grey molecular surface of the chain, EXCLUDING the highlighted residues (so the red patch shows)
mol representation QuickSurf 1.0 0.7 1.0 2.5
mol material AOChalky
mol color ColorID 2
if {$hassel ne ""} {
    mol selection "chain $chain and protein and not ($hassel)"
} else {
    mol selection "chain $chain and protein"
}
mol addrep top

# highlighted residues as a red surface patch (same surface style, carved out of the grey)
if {$hassel ne ""} {
    mol representation QuickSurf 1.0 0.7 1.0 2.5
    mol material AOShiny
    mol color ColorID $hlcol
    mol selection "$hassel"
    mol addrep top
}

catch {display projection Orthographic}
catch {axes location off}
catch {color Display Background white}
catch {display ambientocclusion on}
catch {display shadows on}
catch {display depthcue off}

set sel [atomselect top "chain $chain and protein"]
$sel global
set m [transaxis x $rotx]
$sel move $m
set m2 [transaxis y $roty]
$sel move $m2
catch {display resetview}

render Tachyon $out
quit
