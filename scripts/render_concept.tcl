# Render a molecular scene from a pipe-delimited spec file (one NewCartoon rep per line).
# spec line:  pdbpath|colorID|material|scale|selection
# args: specfile out.dat rotx roty rotz
set spec [lindex $argv 0]; set out [lindex $argv 1]
set rx [lindex $argv 2]; set ry [lindex $argv 3]; set rz [lindex $argv 4]
catch {material change opacity Transparent 0.14}
catch {material change opacity Ghost 0.10}
set fp [open $spec r]
foreach line [split [read $fp] "\n"] {
  if {[string trim $line] eq ""} continue
  set p [split $line "|"]
  mol new [lindex $p 0] waitfor all
  mol delrep 0 top
  mol representation NewCartoon [lindex $p 3] 16.0 5.0
  mol color ColorID [lindex $p 1]
  mol selection [lindex $p 4]
  catch {mol material [lindex $p 2]}
  mol addrep top
}
close $fp
catch {display projection Orthographic}
catch {display depthcue off}
catch {display shadows on}
catch {display ambientocclusion on}
catch {color Display Background white}
catch {axes location Off}
catch {display resetview}
catch {rotate x by $rx}
catch {rotate y by $ry}
catch {rotate z by $rz}
render Tachyon $out
quit
