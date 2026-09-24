#!/usr/bin/env bash
# Rasterize the application-concept SVGs (diagnostics, vaccine) to PNG with headless Chrome @2x.
# Source: scripts/app_figs.py -> build/appfig_*.svg ; output: data/dp4_binding/figs/app_*.png
set -euo pipefail
cd "$(dirname "$0")/.."
CHROME="/Applications/Google Chrome.app/Contents/MacOS/Google Chrome"
FIG=data/dp4_binding/figs; B=build
/usr/bin/python3 scripts/app_figs.py

shot () {  # name  W  H
  local n=$1 W=$2 H=$3
  cat >"$B/$n.html" <<HTML
<!doctype html><html><head><meta charset="utf-8"><style>
html,body{margin:0;padding:0;background:#fff}
svg{display:block}
</style></head><body>$(cat "$B/appfig_$n.svg")</body></html>
HTML
  "$CHROME" --headless=new --disable-gpu --hide-scrollbars --force-device-scale-factor=2 \
    --default-background-color=FFFFFFFF --window-size="$W,$H" \
    --screenshot="$B/$n.png" "file://$PWD/$B/$n.html" >/dev/null 2>&1
  magick "$B/$n.png" -background white -flatten -trim +repage "$FIG/app_$n.png"
  echo "  app_$n.png  $(magick identify -format '%wx%h' "$FIG/app_$n.png")"
}
shot diagnostics 1040 860
shot vaccine 1460 860
echo "wrote $FIG/app_diagnostics.png and $FIG/app_vaccine.png"
