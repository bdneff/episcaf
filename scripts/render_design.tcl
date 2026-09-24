# silver scaffold ribbon + red EPITOPE ribbon (two NewCartoon reps)
set struct [lindex $argv 0]; set outimg [lindex $argv 1]; set hilite [lrange $argv 2 end]
mol new $struct waitfor all
mol delrep 0 top
mol representation NewCartoon 0.30 12.0 4.5
mol color ColorID 6
mol selection {protein}
catch {mol material AOChalky}
mol addrep top
if {[llength $hilite] > 0} {
  mol representation NewCartoon 0.34 12.0 4.5
  mol color ColorID 1
  mol selection "protein and resid $hilite"
  catch {mol material AOChalky}
  mol addrep top
}
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
