#!/bin/bash
# render_minibinder_story.sh -- 3 SEPARATE presentation figures (click-through) of the 8VDL story:
#   fig1 antigen+antibody, fig2 + hotspots, fig3 antigen + our minibinder (RMSD-aligned).
# Same pipeline/skill as render_epitopes.sh (headless VMD -> Tachyon scene -> external tachyon PNG).
# Local only (VMD + macOS sips + ImageMagick). Override paths with VMD=/... TACHYON=/...
set -euo pipefail
cd "$(dirname "$0")/../.."          # repo root

VMD="${VMD:-/Applications/VMD.app/Contents/vmd2/lib/vmd_MACOSXARM64}"
TACHYON="${TACHYON:-/Applications/VMD.app/Contents/vmd2/lib/tachyon_MACOSXARM64}"
export VMDDIR="${VMDDIR:-/Applications/VMD.app/Contents/vmd2}"
RES="${RES:-1500}"
TCL=dp4_8vdl/scripts/render_minibinder_story.tcl
FIG=dp4_8vdl/figures
mkdir -p "$FIG"

"$VMD" -dispdev text -e "$TCL" -args /tmp >/dev/null 2>&1

for name in fig1_antibody fig2_hotspots fig3_minibinder; do
  "$TACHYON" "/tmp/${name}.dat" -res "$RES" "$RES" -aasamples 12 -format TARGA \
      -o "/tmp/${name}.tga" >/dev/null 2>&1
  sips -s format png "/tmp/${name}.tga" --out "$FIG/story_${name}.png" >/dev/null
  magick "$FIG/story_${name}.png" -trim +repage -bordercolor white -border 5% "$FIG/story_${name}.png"
  echo "  wrote $FIG/story_${name}.png"
done
echo "3 click-through figures ready in $FIG/ (story_fig1_antibody, story_fig2_hotspots, story_fig3_minibinder)"
