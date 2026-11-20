#!/usr/bin/env bash
#SBATCH -J episcaf-lysozyme
#SBATCH -p gpu-a100
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --no-requeue
set -euo pipefail
umask 077
module load Mamba/25.3.0-3
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate ~/rfd3/env/rfd3_py312
python -m episcaf_service.lysozyme_job "$1"
