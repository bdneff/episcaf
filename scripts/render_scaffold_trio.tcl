# render_scaffold_trio.tcl -- headless VMD -> Tachyon scenes for the "generate scaffolds" figure.
#
# Loads several real episcaf designs that all present the SAME epitope motif, superposes
# them on the epitope CA atoms (so the red patch sits in an identical place and orientation),
# then writes one Tachyon scene per design with a SHARED camera. Each scene = grey NewCartoon
# fold + red epitope ribbon/licorice. Ray-trace each .dat with the external tachyon binary.
#
# Usage:
#   vmd -dispdev text -e render_scaffold_trio.tcl -args OUTPREFIX EPISEL ROTX ROTY PDB1 PDB2 PDB3 ...
#     EPISEL : VMD resid selection for the epitope motif, e.g. "35 69 70 71"
#     OUTPREFIX : scenes written as OUTPREFIX_0.dat, OUTPREFIX_1.dat, ...
# Same epitope CA set across inputs (identical contig) => clean superposition.
set argv $argv
set outpref [lindex $argv 0]
set episel  [lindex $argv 1]
set rotx    [lindex $argv 2]
set roty    [lindex $argv 3]
set pdbs    [lrange $argv 4 end]

set episelCA "protein and name CA and resid $episel"
set episelAll "protein and resid $episel"

# ---- load all, superpose on the first design's epitope CA ----
set mols {}
set refmol -1
foreach f $pdbs {
    mol new $f type pdb waitfor all
    set m [molinfo top]
    lappend mols $m
    if {$refmol < 0} {
        set refmol $m
        set ref [atomselect $refmol $episelCA]
    } else {
        set mob [atomselect $m $episelCA]
        if {[$mob num] == [$ref num] && [$ref num] > 0} {
            set T [measure fit $mob $ref]
            set all [atomselect $m "all"]
            $all move $T
        } else {
            puts "WARNING: epitope CA mismatch for $f ([$mob num] vs [$ref num]); not aligned"
        }
    }
}

# ---- representations for every molecule ----
foreach m $mols {
    mol delrep 0 $m
    # grey fold
    mol representation NewCartoon 0.30 16 4.5
    mol material AOChalky
    mol color ColorID 2
    mol selection "protein"
    mol addrep $m
    # epitope ribbon (red, slightly thicker) over the grey
    mol representation NewCartoon 0.36 16 4.5
    mol material AOChalky
    mol color ColorID 1
    mol selection $episelAll
    mol addrep $m
    # epitope sidechains as red licorice
    mol representation Licorice 0.3 12 12
    mol material AOShiny
    mol color ColorID 1
    mol selection $episelAll
    mol addrep $m
}

catch {display projection Orthographic}
catch {axes location off}
catch {color Display Background white}
catch {display ambientocclusion on}
catch {display shadows on}
catch {display depthcue off}

# ---- shared camera: rotate ALL coordinates by the same matrix, then frame once ----
foreach m $mols {
    set s [atomselect $m "all"]
    $s move [transaxis x $rotx]
    $s move [transaxis y $roty]
}
catch {display resetview}

# ---- render one scene per design (only that design visible) ----
set i 0
foreach m $mols {
    foreach mm $mols { mol off $mm }
    mol on $m
    render Tachyon ${outpref}_${i}.dat
    incr i
}
quit
