#!/usr/bin/env bash
set -euo pipefail

LABEL="${1:-dp2_rfd3}"
STAMP=$(date +"%Y%m%d_%H%M%S")
RUN_DIR="runs/${LABEL}/run_${STAMP}"

mkdir -p "runs/${LABEL}"
cp -r runs/.template "$RUN_DIR"

echo "$RUN_DIR"
