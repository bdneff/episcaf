#!/usr/bin/env bash
# Real-structure concept renders for the "why it matters" slide (diagnostics + vaccine).
# Prep: scripts/dp4_concept_prep.py  (writes data/dp4_binding/concept/*.pdb, native frame).
# Output: data/dp4_binding/figs/{diagnostics,vaccine}.png
set -euo pipefail
cd "$(dirname "$0")/.."
export VMDDIR=/Applications/VMD.app/Contents/vmd2/lib
VMD=$VMDDIR/vmd_MACOSXARM64; TACHYON=$VMDDIR/tachyon_MACOSXARM64
FONT=/System/Library/Fonts/Supplemental/Arial.ttf
C=data/dp4_binding/concept; FIG=data/dp4_binding/figs; W=/tmp/concept; mkdir -p "$W"

/usr/bin/python3 scripts/dp4_concept_prep.py

# colorIDs: 6 silver, 1 red, 0 blue, 8 white-ish.  materials: AOChalky matte, Transparent faint.
epi_low="protein and resid 21 to 32 70 to 83"
epi_nat="protein and resid 205 to 216 279 to 292"

panel () {  # name specfile rotx roty rotz
  local n=$1 spec=$2 rx=$3 ry=$4 rz=$5
  "$VMD" -dispdev text -e scripts/render_concept.tcl -args "$spec" "$W/$n.dat" "$rx" "$ry" "$rz" >/dev/null 2>&1
  "$TACHYON" "$W/$n.dat" -res 1300 1300 -aasamples 12 -trans_vmd -format TARGA -o "$W/$n.tga" >/dev/null 2>&1
  sips -s format png "$W/$n.tga" --out "$W/$n.png" >/dev/null 2>&1
  magick "$W/$n.png" -background white -flatten -trim +repage -resize 1180x1180 \
         -gravity center -extent 1220x1220 -background white "$W/$n.png"
  echo "  panel $n"
}
panelfit () {  # name specfile rotx roty rotz  -> trimmed, natural aspect (portrait figures)
  local n=$1 spec=$2 rx=$3 ry=$4 rz=$5
  "$VMD" -dispdev text -e scripts/render_concept.tcl -args "$spec" "$W/$n.dat" "$rx" "$ry" "$rz" >/dev/null 2>&1
  "$TACHYON" "$W/$n.dat" -res 1300 1300 -aasamples 12 -trans_vmd -format TARGA -o "$W/$n.tga" >/dev/null 2>&1
  sips -s format png "$W/$n.tga" --out "$W/$n.png" >/dev/null 2>&1
  magick "$W/$n.png" -background white -flatten -trim +repage "$W/$n.png"
  echo "  panelfit $n"
}

# ---- EPITOPE DEFINITION: the real antibody gripping the epitope patch on the whole antigen ----
cat >"$W/epdef.spec" <<SPEC
$C/native_antigen.pdb|6|AOChalky|0.28|protein
$C/native_antigen.pdb|1|AOChalky|0.42|protein and resid 205 to 216 279 to 292
$C/scfv.pdb|0|AOChalky|0.28|protein
SPEC
panelfit epdef "$W/epdef.spec" -80 15 0
cp "$W/epdef.png" "$FIG/epitope_def.png"

# ---- DISCONTINUOUS EPITOPE: two sequence-distant segments held together in 3D (limitations slide) ----
# 205-216 (red) and 279-292 (orange) are ~63 residues apart in the chain but touch in space.
cat >"$W/disc.spec" <<SPEC
$C/native_antigen.pdb|6|Transparent|0.26|protein
$C/native_antigen.pdb|1|AOChalky|0.50|protein and resid 205 to 216
$C/native_antigen.pdb|3|AOChalky|0.50|protein and resid 279 to 292
SPEC
panelfit disc "$W/disc.spec" -80 15 0
cp "$W/disc.png" "$FIG/disc_epitope.png"

# ---- DIAGNOSTICS: the real scFv antibody docked on the scaffolded epitope ----
cat >"$W/diag.spec" <<SPEC
$C/aligned_low.pdb|6|AOChalky|0.30|protein
$C/aligned_low.pdb|1|AOChalky|0.45|$epi_low
$C/scfv.pdb|0|AOChalky|0.30|protein
SPEC
panel diag "$W/diag.spec" -80 15 0

# ---- VACCINE: conserved epitope as a patch on the whole antigen vs. isolated on the scaffold ----
cat >"$W/vac_a.spec" <<SPEC
$C/native_antigen.pdb|6|AOChalky|0.28|protein
$C/native_antigen.pdb|1|AOChalky|0.42|$epi_nat
SPEC
cat >"$W/vac_b.spec" <<SPEC
$C/aligned_low.pdb|6|AOChalky|0.30|protein
$C/aligned_low.pdb|1|AOChalky|0.45|$epi_low
SPEC
panel vac_a "$W/vac_a.spec" -80 15 0
panel vac_b "$W/vac_b.spec" -80 15 0

montage -font "$FONT" -pointsize 44 -background white -fill "#212121" \
  -label "the real antibody binds the scaffolded epitope" "$W/diag.png" \
  -tile 1x1 -geometry +0+10 "$FIG/diagnostics.png"
montage -font "$FONT" -pointsize 42 -background white -fill "#212121" \
  -label "conserved epitope: a patch on the antigen" "$W/vac_a.png" \
  -label "scaffolded: the dominant target" "$W/vac_b.png" \
  -tile 2x1 -geometry +16+12 "$FIG/vaccine.png"
echo "wrote $FIG/diagnostics.png and $FIG/vaccine.png"
