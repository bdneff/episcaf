#!/usr/bin/env bash
# Render the Epitope Scaffolder dropdown complexes with VMD: antigen grey, antibody blue,
# contact epitope red. Output -> docs/img/<ID>.png (served by GitHub Pages next to index.html).
set -euo pipefail
cd "$(dirname "$0")/.."
export VMDDIR=/Applications/VMD.app/Contents/vmd2/lib
VMD=$VMDDIR/vmd_MACOSXARM64; TACHYON=$VMDDIR/tachyon_MACOSXARM64
W=/tmp/webtool; mkdir -p "$W" docs/img docs/pdb

# ID  input  antigen_chain  antibody_chains(csv)  rotx roty rotz
targets=(
  "8PWW data/dp4_binding/structures/8PWW.cif A B   -75 20 0"
  "7OX3 data/dp4_binding/structures/7OX3.cif C A,B -75 20 0"
  "5FHX data/dp4_binding/structures/5FHX.cif A H,L -75 20 0"
  "8VDL dp4_8vdl/data/8VDL.pdb              C H,L  -75 20 0"
)
for t in "${targets[@]}"; do
  set -- $t; id=$1; inp=$2; ag=$3; ab=$4; rx=$5; ry=$6; rz=$7
  epi=$(/usr/bin/python3 scripts/webtool_prep.py "$inp" "$ag" "$ab" "docs/pdb/$id.pdb")
  {
    echo "docs/pdb/$id.pdb|6|AOChalky|0.30|protein and chain $ag"
    echo "docs/pdb/$id.pdb|0|AOChalky|0.30|protein and chain ${ab//,/ }"
    [ -n "$epi" ] && echo "docs/pdb/$id.pdb|1|AOChalky|0.50|protein and chain $ag and resid $epi"
  } > "$W/$id.spec"
  "$VMD" -dispdev text -e scripts/render_concept.tcl -args "$W/$id.spec" "$W/$id.dat" "$rx" "$ry" "$rz" >/dev/null 2>&1
  "$TACHYON" "$W/$id.dat" -res 1120 920 -aasamples 12 -format TARGA -o "$W/$id.tga" >/dev/null 2>&1
  sips -s format png "$W/$id.tga" --out "$W/$id.png" >/dev/null 2>&1
  magick "$W/$id.png" -background white -flatten -trim +repage -bordercolor white -border 26 -resize 960x "docs/img/$id.png"
  echo "  $id  antigen=$ag antibody=$ab  epitope=[$epi]  -> docs/img/$id.png ($(magick identify -format '%wx%h' docs/img/$id.png))"
done
echo "done."
