#!/bin/bash
# Phase-0 smoke: prove ColabFold (AlphaFold2) folds a sequence on a Gemini GPU node and emits PAE.
# See episcaf_service/gcloud/rfd4_setup.md. Brandon submits this; it touches nothing in the pilot.
#
# EDIT the two marked lines (partition/account, colabfold_batch path), then:
#   sbatch episcaf_service/gcloud/af2_smoke.sbatch
#
#SBATCH --job-name=af2-smoke
#SBATCH --partition=gpu-a100        # <-- confirm Gemini GPU partition
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=00:40:00
#SBATCH --no-requeue
#SBATCH --output=af2_smoke_%j.log
set -euo pipefail
umask 077

nvidia-smi || true
COLABFOLD="$HOME/localcolabfold/colabfold-conda/bin/colabfold_batch"   # <-- set to your install

WORK="$SLURM_SUBMIT_DIR/runs/af2_smoke_$SLURM_JOB_ID"
mkdir -p "$WORK"
# A short ubiquitin-like test sequence (fast to fold; smoke only, not a real design):
cat > "$WORK/smoke.fasta" <<'FASTA'
>smoke
MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG
FASTA

# One model, one recycle for a quick proof of life (the real engine uses 5 models / 3 recycles, no templates):
"$COLABFOLD" --num-models 1 --num-recycle 1 "$WORK/smoke.fasta" "$WORK/out"

echo "=== af2 smoke output ==="
find "$WORK/out" -maxdepth 2 -type f | head -40
echo "AF2 smoke OK if a *_rank_001_*.pdb and a scores/PAE json are listed above."
echo "If it stalled at 'Running MSA', GPU-node outbound traffic to the MMseqs2 API is blocked (see setup §2)."
