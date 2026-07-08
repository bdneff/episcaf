#!/usr/bin/env bash
#SBATCH -J rfd3_8vdl
#SBATCH -p gpu-v100
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH -t 2:00:00
# NOTE: no hardcoded --chdir (the minibinder original pointed at a dead ~/rfd3 path); submit from
# dp4_8vdl/ and pass the run dir as $1. Dry-run once on the current cluster before the full array.
#
# Usage (from dp4_8vdl/):
#   sbatch --array=1-N%200 scripts/03_rfd3_array.sbatch 02_rfd3/<run>
# where <run> is epitope20 or hotspots. Each task processes CHUNK_SIZE contigs (8 backbones each).

set -euo pipefail

module load Mamba/25.3.0-3
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate ~/rfd3/env/rfd3_py312

RUN_DIR="${1:?usage: sbatch --array=1-N scripts/03_rfd3_array.sbatch <run_dir e.g. 02_rfd3/epitope20>}"
RUN_DIR="$(python3 -c "import os; print(os.path.abspath('${RUN_DIR}'))")"

MANIFEST="${RUN_DIR}/inputs_manifest.csv"
OUT_ROOT="${RUN_DIR}/outputs"
LOG_DIR="${RUN_DIR}/logs"
mkdir -p "${OUT_ROOT}" "${LOG_DIR}"

CHUNK_SIZE="${CHUNK_SIZE:-1}"        # one contig (8 backbones) per task by default
TASK_ID="${SLURM_ARRAY_TASK_ID}"
START=$(( (TASK_ID - 1) * CHUNK_SIZE + 1 ))
END=$(( TASK_ID * CHUNK_SIZE ))

echo "HOST=$(hostname)  TASK=${TASK_ID}  CHUNK=${START}-${END}  RUN_DIR=${RUN_DIR}"

mapfile -t CHUNK < <(python3 - <<PY
import csv
with open("${MANIFEST}", newline="") as f:
    r = csv.reader(f); next(r); rows = list(r)
for row in rows[${START}-1:${END}]:
    print(row[0], row[1])
PY
)

if [[ ${#CHUNK[@]} -eq 0 ]]; then
    echo "No designs for task ${TASK_ID}, exiting"; exit 0
fi

nvidia-smi || true
for LINE in "${CHUNK[@]}"; do
    DESIGN_ID=$(echo "$LINE" | awk '{print $1}')
    JSON_PATH=$(echo "$LINE" | awk '{print $2}')
    OUT_DIR="${OUT_ROOT}/${DESIGN_ID}"
    DONE_FLAG="${OUT_DIR}/_DONE"
    if [[ -f "${DONE_FLAG}" ]]; then echo "SKIP (done): ${DESIGN_ID}"; continue; fi
    mkdir -p "${OUT_DIR}"
    echo "Running: ${DESIGN_ID}"
    rfd3 design \
        out_dir="${OUT_DIR}" \
        inputs="${JSON_PATH}" \
        skip_existing=False \
        dump_trajectories=False \
        prevalidate_inputs=True \
        && touch "${DONE_FLAG}" \
        || echo "FAILED: ${DESIGN_ID}"
done
echo "Task ${TASK_ID} complete."
