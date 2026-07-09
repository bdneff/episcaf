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

names=(fig1_antibody fig2_hotspots fig3_minibinder)
raws=()
for name in "${names[@]}"; do
  "$TACHYON" "/tmp/${name}.dat" -res "$RES" "$RES" -aasamples 12 -format TARGA \
      -o "/tmp/${name}.tga" >/dev/null 2>&1
  magick "/tmp/${name}.tga" "/tmp/${name}_raw.png"        # full canvas, NO per-image trim
  raws+=("/tmp/${name}_raw.png")
done

# ONE common crop box = TRUE union of all three panels' content boxes, so every figure keeps the same
# size AND pixel position (the VMD view is identical across panels). Superimposed in a slideshow,
# nothing moves. (Can't -flatten to get this: opaque white backgrounds hide lower layers, collapsing
# the "union" to the top image. So union the per-image trim boxes ourselves.)
X0=1000000; Y0=1000000; X1=0; Y1=0
for r in "${raws[@]}"; do
  bb=$(magick "$r" -format "%@" info:)          # WxH+X+Y content box
  w=${bb%%x*}; t=${bb#*x}; h=${t%%+*}; t=${t#*+}; x=${t%%+*}; y=${t##*+}
  (( x < X0 )) && X0=$x
  (( y < Y0 )) && Y0=$y
  (( x + w > X1 )) && X1=$((x + w))
  (( y + h > Y1 )) && Y1=$((y + h))
done
BOX="$((X1 - X0))x$((Y1 - Y0))+${X0}+${Y0}"
echo "  common crop box (true union of all 3): $BOX"
for name in "${names[@]}"; do
  magick "/tmp/${name}_raw.png" -crop "$BOX" +repage -bordercolor white -border 60 \
      "$FIG/story_${name}.png"
  echo "  wrote $FIG/story_${name}.png"
done
echo "3 click-through figures (identical size/position) ready in $FIG/"
identify -format "%f  %wx%h\n" "$FIG"/story_fig*.png
