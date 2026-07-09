# render_minibinder_story.tcl -- 3 click-through panels (same orientation) for a presentation:
#   fig1  8VDL: antigen (chain C) GREEN cartoon, C7 antibody (H/L) GRAY
#   fig2  same + the conserved hotspots F655/F656/E666 in RED licorice
#   fig3  the 8VDL CIDR (green) with our BindCraft minibinder (PURPLE) in place of the antibody
#         -- the minibinder's CIDR (bindcraft chain A) is RMSD-fit onto 8VDL chain C, carrying the
#            minibinder (chain B) into the binding site.
# Writes fig1/2/3 Tachyon scene (.dat) files; the driver ray-traces them at high res.
# Usage: vmd -dispdev text -e render_minibinder_story.tcl -args <outdir>

set outdir [lindex $argv 0]
set VD [mol new dp4_8vdl/data/8VDL.pdb waitfor all]
set MB [mol new dp4_8vdl/data/bindcraft-3.pdb waitfor all]

# --- RMSD-align the minibinder's CIDR (chain A) onto the 8VDL CIDR (chain C) ---
# both CIDRs are the same 76 resolved residues in the same order, so CA<->CA fit is 1:1.
set ref [atomselect $VD "chain C and name CA and protein"]
set tar [atomselect $MB "chain A and name CA and protein"]
set mov [atomselect $MB "all"]
$mov move [measure fit $tar $ref]

# --- global display / lighting (publication still) ---
catch {display projection Orthographic}
catch {display depthcue off}
catch {display ambientocclusion on}
catch {display shadows on}
catch {color Display Background white}
catch {axes location Off}

# Orient by rotating the COORDINATES (not the view), then `display resetview` fits EVERYTHING in
# frame. Rotating the *view* clips an elongated antigen+Fab off the fixed render canvas (that was the
# cut-off antibody). Both molecules are already aligned, so rotate both by the SAME matrix to keep
# them registered. One shared orientation -> the three panels superimpose for the click-through.
set R [transmult [transaxis z -20] [transmult [transaxis y 15] [transaxis x -90]]]
foreach m [list $VD $MB] { set s [atomselect $m all]; $s move $R; $s delete }
catch {display resetview}

proc clearreps {m} { while {[molinfo $m get numreps] > 0} { mol delrep 0 $m } }

# ColorIDs: 7 green, 2 gray, 1 red, 11 purple
set HOT {chain C and resid 655 656 666}

# ---- Figure 1: antigen green, antibody gray ----
clearreps $VD; clearreps $MB
mol representation NewCartoon 0.30 12.0 4.5
mol material AOChalky
mol color ColorID 7
mol selection {chain C and protein}
mol addrep $VD
mol color ColorID 2
mol selection {(chain H or chain L) and protein}
mol addrep $VD
render Tachyon $outdir/fig1_antibody.dat

# ---- Figure 2: conserved hotspots RED in both the ribbon and the sidechains ----
clearreps $VD; clearreps $MB
mol representation NewCartoon 0.30 12.0 4.5
mol material AOChalky
mol color ColorID 7
mol selection {chain C and protein and not (resid 655 656 666)}   ;# antigen, minus hotspots
mol addrep $VD
mol color ColorID 1
mol selection $HOT                                                ;# hotspot ribbon red
mol addrep $VD
mol color ColorID 2
mol selection {(chain H or chain L) and protein}                  ;# antibody gray
mol addrep $VD
mol representation Licorice 0.30 12.0 12.0
mol material AOShiny
mol color ColorID 1
mol selection $HOT                                                ;# hotspot sidechains red
mol addrep $VD
render Tachyon $outdir/fig2_hotspots.dat

# ---- Figure 3: antigen green + minibinder purple (antibody removed) ----
clearreps $VD; clearreps $MB
mol representation NewCartoon 0.30 12.0 4.5
mol material AOChalky
mol color ColorID 7
mol selection {chain C and protein}
mol addrep $VD
mol color ColorID 11
mol selection {chain B and protein}
mol addrep $MB
render Tachyon $outdir/fig3_minibinder.dat

quit
