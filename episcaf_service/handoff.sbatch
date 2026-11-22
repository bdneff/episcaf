#!/usr/bin/env bash
#SBATCH --job-name=episcaf-pilot
#SBATCH --partition=gpu-a100
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=01:00:00
#SBATCH --no-requeue
set -euo pipefail
umask 077
module load Mamba/25.3.0-3
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$HOME/rfd3/env/rfd3_py312"
cd "$HOME/episcaf-advisor"
python -m episcaf_service.handoff_job
