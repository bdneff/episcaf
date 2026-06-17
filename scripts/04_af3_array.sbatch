#!/usr/bin/env bash
#SBATCH -J af3_rfd3_mpnn
#SBATCH -p gpu-v100
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH -t 4:00:00
#SBATCH -o logs/rfd3_mpnn/03_af3/af3_%A_%a.out
#SBATCH -e logs/rfd3_mpnn/03_af3/af3_%A_%a.err
#SBATCH --chdir=/home/bneff/rfd3/repo_refactored

# Usage:
#   sbatch --array=1-N scripts/04_af3_array.sbatch <run_dir>
#
# Each array task processes CHUNK_SIZE JSONs sequentially on one GPU.
# V100 requires two extra flags vs A100:
#   - SINGULARITYENV_XLA_FLAGS to disable incompatible HLO pass
#   - --flash_attention_implementation xla (triton not supported on compute capability 7.x)

set -euo pipefail

module load singularity
export SINGULARITYENV_XLA_FLAGS="--xla_disable_hlo_passes=custom-kernel-fusion-rewriter"

RUN_DIR="${1:?usage: sbatch --array=1-N scripts/04_af3_array.sbatch <run_dir>}"
RUN_DIR="$(python3 -c "import os; print(os.path.abspath('${RUN_DIR}'))")"

IN_DIR="${RUN_DIR}/03_af3/inputs"
OUT_ROOT="${RUN_DIR}/03_af3/outputs"
LOG_DIR="logs/rfd3_mpnn/03_af3"
mkdir -p "${OUT_ROOT}" "${LOG_DIR}"

CHUNK_SIZE=100
TASK_ID="${SLURM_ARRAY_TASK_ID}"
START=$(( (TASK_ID - 1) * CHUNK_SIZE + 1 ))
END=$(( TASK_ID * CHUNK_SIZE ))

echo "HOST      = $(hostname)"
echo "TASK_ID   = ${TASK_ID}"
echo "CHUNK     = ${START}-${END}"

mapfile -t JSON_FILES < <(find "${IN_DIR}" -maxdepth 1 -name "*.json" | sort | sed -n "${START},${END}p")

if [[ ${#JSON_FILES[@]} -eq 0 ]]; then
    echo "No JSONs for task ${TASK_ID}, exiting"
    exit 0
fi

echo "Processing ${#JSON_FILES[@]} JSONs"

for JSON_PATH in "${JSON_FILES[@]}"; do
    PRED_ID="$(basename "${JSON_PATH}" .json)"
    OUT_DIR="${OUT_ROOT}/${PRED_ID}"
    DONE_FLAG="${OUT_DIR}/_DONE"

    if [[ -f "${DONE_FLAG}" ]]; then
        echo "SKIP (done): ${PRED_ID}"
        continue
    fi

    mkdir -p "${OUT_DIR}"
    echo "Running: ${PRED_ID}"

    singularity exec --nv \
        -B /home,/scratch,/tgen_labs,/ref_genomes \
        /tgen_labs/altin/alphafold3/containers/alphafold_3.0.1.sif \
        python /app/alphafold/run_alphafold.py \
            --json_path    "${JSON_PATH}" \
            --model_dir    /ref_genomes/alphafold/alphafold3/models \
            --db_dir       /ref_genomes/alphafold/alphafold3/ \
            --output_dir   "${OUT_ROOT}" \
            --norun_data_pipeline \
            --num_diffusion_samples 1 \
            --flash_attention_implementation xla \
        && touch "${DONE_FLAG}" \
        || echo "FAILED: ${PRED_ID}"

done

echo "Task ${TASK_ID} complete."
