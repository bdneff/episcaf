set nat [lindex $argv 0]; set out [lindex $argv 1]; set epi [lindex $argv 2]
catch {material change opacity Transparent 0.17}
mol new $nat waitfor all
mol delrep 0 top
mol representation NewCartoon 0.30 12.0 4.5
mol color ColorID 6
mol selection {protein}
catch {mol material Transparent}
mol addrep top
mol new $epi waitfor all
mol delrep 0 top
mol representation NewCartoon 0.45 12.0 4.5
mol color ColorID 1
mol selection {protein}
catch {mol material Opaque}
mol addrep top
catch {display projection Orthographic}
catch {display shadows on}
catch {color Display Background white}
catch {axes location Off}
catch {display resetview}
catch {rotate x by -75}
render Tachyon $out
quit
