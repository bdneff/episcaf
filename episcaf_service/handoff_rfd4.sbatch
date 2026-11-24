#!/bin/bash
# The RFD4-scaffold cluster job: the data-mover worker sbatches this (one per request) with the
# worker inbox as SLURM_SUBMIT_DIR. Runs the RFD4 campaign engine on gpu-a100 in the RFD4 venv.
# Tools + weights live on tgen_labs (home is quota-bound); see episcaf_service/gcloud/rfd4_setup.md.
#
#SBATCH --job-name=episcaf-rfd4
#SBATCH --partition=gpu-a100
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=02:00:00         # <=2h -> axshort QOS; raise + set the right QOS for large campaigns
#SBATCH --no-requeue
#SBATCH --open-mode=append
set -euo pipefail
umask 077

TGEN=/tgen_labs/altin/alphafold3/workspace/episcaf_v2_bneff
source "$HOME/RFD4-Proteina/.venv/bin/activate"
export HF_HOME="$TGEN/weights/huggingface"
export TORCH_COMPILE_CACHE_BLOB="$TGEN/weights/rfd4_compile_cache.blob"
export EPISCAF_COLABFOLD="$TGEN/tools/localcolabfold/colabfold-conda/bin/colabfold_batch"
export EPISCAF_AF2_DATA="$TGEN/tools/localcolabfold/colabfold"
export EPISCAF_RFD4_CMD="rfd4-generate"

cd "$HOME/episcaf-advisor"
exec python -m episcaf_service.handoff_rfd4_job
