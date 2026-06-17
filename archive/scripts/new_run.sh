#!/usr/bin/env bash
set -euo pipefail

# Create a timestamped run directory and snapshot a dataset parquet into it.
#
# Usage:
#   scripts/new_run.sh <dataset.parquet> [run_name_prefix]
#
# Example:
#   scripts/new_run.sh datasets/dp2.parquet run_dp2

DATASET="${1:?need dataset parquet path}"
PREFIX="${2:-run}"
STAMP=$(date +"%Y%m%d_%H%M%S")
RUN_DIR="runs/${PREFIX}_${STAMP}"

python -m episcaf_pipeline init --dataset "$DATASET" --run_dir "$RUN_DIR"
echo "$RUN_DIR"
