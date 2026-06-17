#!/usr/bin/env bash
#SBATCH -J rfd3_dp2
#SBATCH -p gpu-a100
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH -t 00:15:00
#SBATCH --array=1-10

set -euo pipefail

module load Mamba/25.3.0-3
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate ~/rfd3/env/rfd3_py312

RUN_DIR_IN="${1:?usage: sbatch [--array=1-N] sbatch/03_rfd3_array.sbatch <RUN_DIR>}"

# Absolute RUN_DIR (avoid SLURM spool surprises)
RUN_DIR="$(python - <<PY
import os
print(os.path.abspath("${RUN_DIR_IN}"))
PY
)"

MANIFEST="${RUN_DIR}/02_rfd3/inputs_manifest.csv"
LOG_DIR="${RUN_DIR}/02_rfd3/logs"
mkdir -p "$LOG_DIR"

TASK_ID="${SLURM_ARRAY_TASK_ID}"

# Robust CSV parsing (handles commas, whitespace, CRLF)
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
# Print space-separated for bash read
print(design_id, json_path)
PY
)

LOG_FILE="${LOG_DIR}/${DESIGN_ID}.log"

# Compute repo root from RUN_DIR (everything before /runs/)
REPO_ROOT="${RUN_DIR%%/runs/*}"

echo "HOST      = $(hostname)"
echo "TASK_ID   = $TASK_ID"
echo "DESIGN_ID = $DESIGN_ID"
echo "JSON_PATH (manifest) = $JSON_PATH"
echo "RUN_DIR   = $RUN_DIR"
echo "REPO_ROOT = $REPO_ROOT"
echo "LOG_FILE  = $LOG_FILE"

cd "$REPO_ROOT"

# Make JSON_PATH absolute relative to repo root (in case manifest has relative paths)
if [[ "$JSON_PATH" != /* ]]; then
  JSON_PATH="${REPO_ROOT}/${JSON_PATH}"
fi
echo "JSON_PATH (absolute) = $JSON_PATH"

if [[ ! -f "$JSON_PATH" ]]; then
  echo "ERROR: JSON file not found: $JSON_PATH"
  echo "PWD=$(pwd)"
  echo "Contents of RUN_DIR/02_rfd3/inputs:"
  ls -la "${RUN_DIR}/02_rfd3/inputs" | head -50 || true
  exit 3
fi

# Per-design output dir avoids collisions and makes debugging easy
OUT_DIR="${RUN_DIR}/02_rfd3/outputs/${DESIGN_ID}"
mkdir -p "$OUT_DIR"

echo "OUT_DIR   = $OUT_DIR"
echo "which rfd3 = $(which rfd3 || true)"
nvidia-smi || true

# Log the exact command line
set -x
rfd3 design \
  out_dir="$OUT_DIR" \
  inputs="$JSON_PATH" \
  skip_existing=False \
  dump_trajectories=False \
  prevalidate_inputs=True \
  2>&1 | tee "$LOG_FILE"
set +x

echo "Done."

