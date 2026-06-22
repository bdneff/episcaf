#!/usr/bin/env bash
#SBATCH -J af3_rfd3_mpnn
#SBATCH -p gpu-v100
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH -t 4:00:00
#SBATCH -o af3_mpnn_%A_%a.out
#SBATCH -e af3_mpnn_%A_%a.err

# Usage:
#   sbatch --array=1-N scripts/stage04_af3_array.sbatch <run_dir>
#
# Each array task runs its CHUNK_SIZE JSONs through AlphaFold3 in ONE container launch via
# --input_dir, so the model loads ONCE per task (not once per JSON) and the GPU stays busy.
# (Earlier we launched the container per JSON, reloading the model ~100x/task and leaving the
# GPU idle during each load -- HPC flagged the idle GPUs.)
# V100 needs the XLA flag + --flash_attention_implementation xla (no triton on cc 7.x).

set -euo pipefail
module load singularity
export SINGULARITYENV_XLA_FLAGS="--xla_disable_hlo_passes=custom-kernel-fusion-rewriter"

RUN_DIR="${1:?usage: sbatch --array=1-N scripts/stage04_af3_array.sbatch <run_dir>}"
RUN_DIR="$(python3 -c "import os; print(os.path.abspath('${RUN_DIR}'))")"

IN_DIR="${RUN_DIR}/04_af3/inputs"
OUT_ROOT="${RUN_DIR}/04_af3/outputs"
STAGE_ROOT="${RUN_DIR}/04_af3/_task_inputs"
mkdir -p "${OUT_ROOT}" "${STAGE_ROOT}"

CHUNK_SIZE=100
TASK_ID="${SLURM_ARRAY_TASK_ID}"
START=$(( (TASK_ID - 1) * CHUNK_SIZE + 1 ))
END=$(( TASK_ID * CHUNK_SIZE ))

echo "HOST=$(hostname)  TASK=${TASK_ID}  CHUNK=${START}-${END}"

mapfile -t JSON_FILES < <(find "${IN_DIR}" -maxdepth 1 -name "*.json" | sort | sed -n "${START},${END}p")
[[ ${#JSON_FILES[@]} -eq 0 ]] && { echo "No JSONs for task ${TASK_ID}, exiting"; exit 0; }

# Stage only the not-yet-done JSONs into a per-task input dir. AF3 writes its output to
# OUT_ROOT/<name lowercased>/, so a JSON is "done" if that dir already holds a model CIF
# (also honor the old per-JSON _DONE flag from the earlier per-file version).
STAGE_DIR="${STAGE_ROOT}/task_${TASK_ID}"
rm -rf "${STAGE_DIR}"; mkdir -p "${STAGE_DIR}"
n_todo=0
for JSON_PATH in "${JSON_FILES[@]}"; do
    PRED_ID="$(basename "${JSON_PATH}" .json)"
    LC="${PRED_ID,,}"
    if compgen -G "${OUT_ROOT}/${LC}/"'*.cif' > /dev/null 2>&1 || [[ -f "${OUT_ROOT}/${PRED_ID}/_DONE" ]]; then
        continue
    fi
    ln -sf "${JSON_PATH}" "${STAGE_DIR}/${PRED_ID}.json"
    n_todo=$((n_todo+1))
done
echo "TASK ${TASK_ID}: ${n_todo}/${#JSON_FILES[@]} to run; model loads once for the batch"
[[ ${n_todo} -eq 0 ]] && { echo "all already done"; rm -rf "${STAGE_DIR}"; exit 0; }

# One container launch -> one model load -> all staged JSONs predicted, GPU busy throughout.
singularity exec --nv \
    -B /home,/scratch,/tgen_labs,/ref_genomes \
    /tgen_labs/altin/alphafold3/containers/alphafold_3.0.1.sif \
    python /app/alphafold/run_alphafold.py \
        --input_dir    "${STAGE_DIR}" \
        --model_dir    /ref_genomes/alphafold/alphafold3/models \
        --db_dir       /ref_genomes/alphafold/alphafold3/ \
        --output_dir   "${OUT_ROOT}" \
        --norun_data_pipeline \
        --num_diffusion_samples 1 \
        --flash_attention_implementation xla

rm -rf "${STAGE_DIR}"
echo "Task ${TASK_ID} complete."
