#!/usr/bin/env bash
set -euo pipefail

# Usage:
#   tools/run_rfd3_batch.sh <RUN_DIR> <N>
#
# Example:
#   tools/run_rfd3_batch.sh runs/dp2_rfd3_try1/run_20260218_105915 4

RUN_DIR="${1:?need RUN_DIR (e.g. runs/.../run_YYYYMMDD_HHMMSS)}"
N="${2:-4}"

MANIFEST="${RUN_DIR}/02_rfd3/inputs_manifest.csv"
INPUT_DIR="${RUN_DIR}/02_rfd3/inputs"
OUT_DIR="${RUN_DIR}/02_rfd3/outputs"
LOG_DIR="${RUN_DIR}/02_rfd3/logs"

mkdir -p "$OUT_DIR" "$LOG_DIR"

echo "RUN_DIR = $RUN_DIR"
echo "MANIFEST = $MANIFEST"
echo "N        = $N"
echo "OUT_DIR  = $OUT_DIR"
echo "LOG_DIR  = $LOG_DIR"

# Grab first N json paths (skip header)
mapfile -t JSONS < <(tail -n +2 "$MANIFEST" | head -n "$N" | awk -F, '{print $2}')

i=0
for js in "${JSONS[@]}"; do
  i=$((i+1))
  base=$(basename "$js" .json)
  log="${LOG_DIR}/${base}.log"

  echo
  echo "[$i/$N] rfd3 design --config $js"
  echo "log: $log"

  # Run from repo root so relative paths inside JSON resolve correctly
  ( cd "$(dirname "$RUN_DIR")/../.." && \
    rfd3 design --config "$js" --output_dir "$OUT_DIR" \
  ) 2>&1 | tee "$log"
done
