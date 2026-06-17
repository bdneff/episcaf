#!/usr/bin/env bash
#SBATCH -J rfd3_design
#SBATCH -p gpu-a100
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH -t 00:10:00
#SBATCH -o %x_%A_%a.out
#SBATCH -e %x_%A_%a.err

set -euo pipefail

# Usage:
#   sbatch --array=1-N episcaf_pipeline/hpc/sbatch/rfd3_array.sbatch <RUN_DIR>
#
# RUN_DIR should be an absolute or relative path to your run folder:
#   runs/run_YYYYMMDD_HHMMSS

module load Mamba/25.3.0-3
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate ~/rfd3/env/rfd3_py312

RUN_DIR_IN="${1:?usage: sbatch --array=1-N rfd3_array.sbatch <RUN_DIR>}"
RUN_DIR="$(python - <<PY
import os
print(os.path.abspath("${RUN_DIR_IN}"))
PY
)"

MANIFEST="${RUN_DIR}/02_rfd3/inputs_manifest.csv"
OUT_DIR="${RUN_DIR}/02_rfd3/outputs"
LOG_DIR="${RUN_DIR}/02_rfd3/logs"
mkdir -p "$OUT_DIR" "$LOG_DIR"

TASK_ID="${SLURM_ARRAY_TASK_ID}"

read -r DESIGN_ID JSON_PATH < <(python - <<PY
import csv
manifest="${MANIFEST}"
task=int("${TASK_ID}")
with open(manifest, newline="") as f:
    r=csv.reader(f)
    header=next(r, None)
    rows=list(r)
if task < 1 or task > len(rows):
    raise SystemExit(f"No line for TASK_ID={task}. Manifest has {len(rows)} rows.")
design_id, json_path = rows[task-1][0].strip(), rows[task-1][1].strip()
print(design_id, json_path)
PY
)

LOG_FILE="${LOG_DIR}/${DESIGN_ID}.log"

echo "HOST      = $(hostname)"
echo "TASK_ID   = $TASK_ID"
echo "DESIGN_ID = $DESIGN_ID"
echo "JSON_PATH = $JSON_PATH"
echo "RUN_DIR   = $RUN_DIR"

set -x
rfd3 design +inputs="$JSON_PATH" +out_dir="$OUT_DIR" n_batches=1 2>&1 | tee "$LOG_FILE"
set +x

echo "Done."
