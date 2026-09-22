#!/bin/bash -l
#SBATCH --job-name=3hfm_alascan
#SBATCH --partition=gpu-v100
#SBATCH --gres=gpu:1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=12:00:00
#SBATCH --array=0-12%1
# User submits. Serial calculation; GPU allocation matches the existing working recipe.
set -euo pipefail
: "${ALASCAN_CODE:?Set absolute path to the staged alanine_scan directory}"
: "${ALASCAN_SOURCE:?Set absolute path to existing 3HFM MD out directory}"
: "${ALASCAN_OUT:?Set a fresh absolute output directory outside the MD directory}"
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate gmxMMPBSA
python "$ALASCAN_CODE/run_scan.py" --source "$ALASCAN_SOURCE" \
    --out "$ALASCAN_OUT" --task "$SLURM_ARRAY_TASK_ID"
