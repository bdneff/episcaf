#!/bin/bash
# render_webtool_figs.sh -- generate the episcaf web figures from real structures (headless VMD+Tachyon).
#
# Produces grey-ribbon renderings of real project structures with verified residue selections
# (epitope / candidate patches / design epitope), ray-traced and post-processed to trimmed,
# transparent-background PNGs under docs/img/fig/. Reproducible: no hand-drawn proteins.
#
#   bash scripts/render_webtool_figs.sh
#
# Structures: docs/pdb/8PWW.pdb (antigen+complex), docs/pdb/8pww_design.pdb (a real scaffold design).
# Selections: 8PWW:A contact epitope (scripts/webtool_prep.py); candidate patches from the growPatch
# demo data (docs/data/growpatch_8pww.js, patches #136/#8/#31); design epitope = design residues whose
# CA sits within 1.5 A of a native-epitope CA (RFdiffusion keeps the motif fixed in native coordinates).
set -e
cd "$(dirname "$0")/.."

export VMDDIR=/Applications/VMD.app/Contents/vmd2/lib
VMD=/Applications/VMD.app/Contents/vmd2/lib/vmd_MACOSXARM64
TACHYON=/Applications/VMD.app/Contents/vmd2/lib/tachyon_MACOSXARM64
TCL=scripts/render_epitope_fig.tcl
OUT=docs/img/fig
TMP="${SCRATCH:-/tmp}/episcaf_figs"
RES=1100
mkdir -p "$OUT" "$TMP"

# render NAME PDB CHAIN HLCOLOR ROTX ROTY RESIDS...
render () {
  name="$1"; pdb="$2"; chain="$3"; col="$4"; rx="$5"; ry="$6"; shift 6
  echo "  rendering $name ..."
  "$VMD" -dispdev text -e "$TCL" -args "$pdb" "$TMP/$name.dat" "$chain" "$col" "$rx" "$ry" "$@" > "$TMP/$name.vmd.log" 2>&1
  "$TACHYON" "$TMP/$name.dat" -res $RES $RES -aasamples 12 -format TARGA -o "$TMP/$name.tga" > "$TMP/$name.tac.log" 2>&1
  sips -s format png "$TMP/$name.tga" --out "$TMP/$name.png" > /dev/null 2>&1
  cp "$TMP/$name.png" "$OUT/$name.png"
}

EPI="205 209 212 213 216 278 279 281 282 284 285 286 288 289 292"
CAND_A="196 197 198 199 200 295 296 297 298 299 300 301 302 303 400 401 402 403"
CAND_B="329 330 331 332 333 334 335 336 337 338 339 340 341 342 343 344 345 346 347 348 349 350"
CAND_SEL="205 206 207 208 209 210 211 212 213 214 215 216 274 275 276 277 278 279 280 281 282 283 284 285 286 287 288 289 290 291 292"
DES_EPI="21 25 28 29 32 69 70 72 73 75 76 77 79 80 83"

# homepage: same antigen + camera, epitope and three candidate patches; plus the real design scaffold
render antigen_epitope   docs/pdb/8PWW.pdb        A 1 -75 20 $EPI
render cand_a            docs/pdb/8PWW.pdb        A 0 -75 20 $CAND_A
render cand_b            docs/pdb/8PWW.pdb        A 0 -75 20 $CAND_B
render cand_selected     docs/pdb/8PWW.pdb        A 1 -75 20 $CAND_SEL
render design_scaffold   docs/pdb/8pww_design.pdb A 1 -75 20 $DES_EPI

echo "post-processing (trim + transparent background) ..."
/Users/bneff/micromamba/envs/bcell-repair/bin/python scripts/trim_transparent.py "$OUT"/antigen_epitope.png "$OUT"/cand_a.png "$OUT"/cand_b.png "$OUT"/cand_selected.png "$OUT"/design_scaffold.png
echo "done -> $OUT"
