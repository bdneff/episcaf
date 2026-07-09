# render_structure.tcl — clean cartoon of a structure -> a Tachyon SCENE file (headless VMD).
# Usage:
#   vmd -dispdev text -e render_structure.tcl -args <structure.pdb|.gro> <out.dat> [resid ...]
# Writes a Tachyon scene; analysis/render_structures.sh then ray-traces it at high
# resolution with the external `tachyon` binary and converts to PNG (macOS `sips`).
# (display resize is ignored in headless text mode, so we render via external Tachyon
#  with an explicit -res for print-quality output rather than render TachyonInternal.)
# Optional trailing resids are drawn as red licorice over the cartoon (e.g. epitope residues).

if {[llength $argv] < 2} {
    puts "usage: vmd -dispdev text -e render_structure.tcl -args <structure> <out.tga> \[resid...\]"
    quit
}
set struct [lindex $argv 0]
set outimg [lindex $argv 1]
set hilite [lrange $argv 2 end]

mol new $struct waitfor all
mol delrep 0 top

# base cartoon, colored by secondary structure, matte (ambient-occlusion) finish
mol representation NewCartoon 0.30 12.0 4.5
mol color Structure
mol selection {protein}
catch {mol material AOChalky}
mol addrep top

# optional: highlight alanine-scan / epitope residues as red licorice
if {[llength $hilite] > 0} {
    mol representation Licorice 0.30 12.0 12.0
    mol color ColorID 1
    mol selection "protein and resid $hilite"
    catch {mol material AOShiny}
    mol addrep top
}

# display / lighting for a publication-style still. Each wrapped in catch so an
# unsupported attribute (e.g. across VMD versions) can't silently abort the render.
catch {display projection Orthographic}
catch {display depthcue off}
catch {display ambientocclusion on}
catch {display shadows on}
catch {color Display Background white}
catch {axes location Off}
catch {display resetview}
catch {rotate x by -75}
catch {rotate y by 20}

render Tachyon $outimg
quit
