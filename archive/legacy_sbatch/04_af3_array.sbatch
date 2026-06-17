#!/usr/bin/env bash
#SBATCH -J af3_dp2
#SBATCH -p gpu-a100
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH -t 04:00:00
#SBATCH --array=1-208
#SBATCH -o %x_%A_%a.out
#SBATCH -e %x_%A_%a.err

set -euo pipefail

RUN_DIR_IN="${1:?usage: sbatch sbatch/04_af3_array.sbatch <RUN_DIR>}"

RUN_DIR="$(python - <<PY
import os
print(os.path.abspath("${RUN_DIR_IN}"))
PY
)"

IN_DIR="${RUN_DIR}/03_af3/inputs_all"
OUT_ROOT="${RUN_DIR}/03_af3/outputs"
LOG_DIR="${RUN_DIR}/03_af3/logs"
mkdir -p "$OUT_ROOT" "$LOG_DIR"

TASK_ID="${SLURM_ARRAY_TASK_ID}"

JSON_PATH="$(ls -1 "${IN_DIR}"/*.json | sed -n "${TASK_ID}p")"
if [[ -z "${JSON_PATH}" ]]; then
  echo "No JSON for TASK_ID=${TASK_ID} in ${IN_DIR}"
  exit 2
fi

PRED_ID="$(basename "$JSON_PATH" .json)"
OUT_DIR="${OUT_ROOT}/${PRED_ID}"
mkdir -p "$OUT_DIR"

LOG_FILE="${LOG_DIR}/${PRED_ID}.log"

echo "HOST    = $(hostname)"
echo "TASK_ID = ${TASK_ID}"
echo "JSON    = ${JSON_PATH}"
echo "OUT_DIR = ${OUT_DIR}"

set -x
singularity exec --nv \
  -B /home,/scratch,/tgen_labs,/ref_genomes \
  /tgen_labs/altin/alphafold3/containers/alphafold_3.0.1.sif \
  python /app/alphafold/run_alphafold.py \
    --json_path "${JSON_PATH}" \
    --model_dir /ref_genomes/alphafold/alphafold3/models \
    --db_dir /ref_genomes/alphafold/alphafold3/ \
    --output_dir "${OUT_DIR}" \
    --norun_data_pipeline \
    --num_diffusion_samples 1 \
  2>&1 | tee "${LOG_FILE}"
set +x

echo "Done."

