#!/bin/bash
# The RFD4-scaffold cluster job: the data-mover worker sbatches this (one per request) with the
# worker inbox as SLURM_SUBMIT_DIR. Runs the RFD4 campaign engine on gpu-a100 in the RFD4 venv.
# Tools + weights live on tgen_labs (home is quota-bound); see episcaf_service/gcloud/rfd4_setup.md.
#
#SBATCH --job-name=episcaf-rfd4
#SBATCH --partition=gpu-a100
#SBATCH --gres=gpu:4            # gpu-a100 nodes have 4 A100s; fold 4 designs at once (dial to gpu:1 if queueing)
#SBATCH --cpus-per-task=16
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
# AF2 is ~95% of the run and per-design independent. Fold N designs at once, one colabfold per GPU, so
# this must match --gres=gpu:N above. gpu-a100 nodes are 4-GPU -> 4 here. The data-mover worker still
# tracks this one job. (Cross-node scaling beyond one node's 4 GPUs is a future SLURM-array step.)
export EPISCAF_FOLD_WORKERS="${EPISCAF_FOLD_WORKERS:-4}"
# REVISIT (Brandon 2026-11-25): folds are short, so gpu-v100 is a lower-contention alternative to the
# in-demand gpu-a100 (slower per fold, but likely shorter queue). To try: --partition=gpu-v100 + the
# right QOS, keep --gres=gpu:N to this node type's GPU count.

cd "$HOME/episcaf-advisor"
exec python -m episcaf_service.handoff_rfd4_job
