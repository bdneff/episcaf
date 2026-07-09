#!/bin/bash
# render_epitopes.sh -- 1x3 manuscript-ready render of the 8VDL CIDR antigen (chain C) with each
# epitope definition highlighted in red. Reuses the bcell-project pattern (render_structure.tcl ->
# external Tachyon at high res -> sips PNG -> magick trim), which \figorbox-style figures use.
#
# Runs LOCALLY (needs VMD + macOS `sips` + ImageMagick). VMD not at the default path? override:
#   VMD=/path/to/vmd_bin TACHYON=/path/to/tachyon_bin bash dp4_8vdl/scripts/render_epitopes.sh
# Output: dp4_8vdl/figures/epitope_definitions_1x3.png  (+ per-panel PNGs)
set -euo pipefail
cd "$(dirname "$0")/../.."          # repo root

VMD="${VMD:-/Applications/VMD.app/Contents/vmd2/lib/vmd_MACOSXARM64}"
TACHYON="${TACHYON:-/Applications/VMD.app/Contents/vmd2/lib/tachyon_MACOSXARM64}"
export VMDDIR="${VMDDIR:-/Applications/VMD.app/Contents/vmd2}"
RES="${RES:-1400}"
# a real font file (ImageMagick can't resolve 'Helvetica' by name on macOS); first that exists
FONT="${FONT:-}"
for f in /System/Library/Fonts/Supplemental/Arial.ttf /System/Library/Fonts/Helvetica.ttc; do
  [ -z "$FONT" ] && [ -f "$f" ] && FONT="$f"
done
TCL=dp4_8vdl/scripts/render_structure.tcl
FIG=dp4_8vdl/figures
PDB=dp4_8vdl/data/8VDL.pdb
CHAINC=$FIG/_8VDL_chainC.pdb
mkdir -p "$FIG"

# antigen chain C only, so `display resetview` frames the CIDR domain (not the whole C+H+L complex)
awk 'substr($0,1,6)=="ATOM  " && substr($0,22,1)=="C"' "$PDB" > "$CHAINC"

names=(epitope20 hotspots contact4A)
resids=("651 to 670" "655 656 666" "652 653 655 656 657 659 660 661 666 667 669 670 673")
labels=("Contiguous 651-670  (1 island)" "Hotspots F655/F656/E666  (3 res)" "Contact 4A  (13 res, 6 islands)")

panels=()
for i in 0 1 2; do
  n=${names[$i]}
  "$VMD" -dispdev text -e "$TCL" -args "$CHAINC" "/tmp/ep_${n}.dat" ${resids[$i]} >/dev/null 2>&1
  "$TACHYON" "/tmp/ep_${n}.dat" -res "$RES" "$RES" -aasamples 12 -format TARGA \
      -o "/tmp/ep_${n}.tga" >/dev/null 2>&1
  sips -s format png "/tmp/ep_${n}.tga" --out "$FIG/epitope_${n}.png" >/dev/null
  magick "$FIG/epitope_${n}.png" -trim +repage -bordercolor white -border 5% "$FIG/epitope_${n}.png"
  panels+=("$FIG/epitope_${n}.png")
  echo "  rendered $n -> $FIG/epitope_${n}.png"
done

# 1x3 montage with panel captions (fall back to no labels if no usable font was found)
if [ -n "$FONT" ]; then
  montage -font "$FONT" -pointsize 34 -background white -fill black \
      -label "${labels[0]}" "${panels[0]}" \
      -label "${labels[1]}" "${panels[1]}" \
      -label "${labels[2]}" "${panels[2]}" \
      -tile 3x1 -geometry +8+8 "$FIG/epitope_definitions_1x3.png"
else
  echo "  (no font found; montaging without labels)"
  montage "${panels[@]}" -tile 3x1 -geometry +8+8 -background white \
      "$FIG/epitope_definitions_1x3.png"
fi
rm -f "$CHAINC"
echo "wrote $FIG/epitope_definitions_1x3.png"
