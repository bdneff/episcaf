#!/bin/bash
# Phase-0 smoke: prove RFD4-Proteina generates one design on a Gemini GPU node.
# See episcaf_service/gcloud/rfd4_setup.md. Brandon submits this; it touches nothing in the pilot.
#
# EDIT the two marked lines below (partition/account, and the venv path OR the SIF path), then:
#   sbatch episcaf_service/gcloud/rfd4_smoke.sbatch
#
#SBATCH --job-name=rfd4-smoke
#SBATCH --partition=gpu-a100        # <-- confirm Gemini GPU partition
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=00:40:00
#SBATCH --no-requeue
#SBATCH --output=rfd4_smoke_%j.log
set -euo pipefail
umask 077

nvidia-smi || true
OUT="$SLURM_SUBMIT_DIR/runs/rfd4_smoke_$SLURM_JOB_ID"
mkdir -p "$OUT"

# inference_experiment=quickstart uses fast settings (one design, metrics off, no compile warm-up) but
# STILL needs an explicit inputs= (a YAML condition spec). We point it at the repo's bundled PD-L1
# example — this smoke just proves RFD4 conditions on a YAML spec (C_CRDGRP+C_SEQ), the same shape
# generate.py will build for our epitopes.
# ---- Route A (native venv): activate it, run from the RFD4-Proteina checkout ----
RFD4_REPO="$HOME/RFD4-Proteina"          # <-- set to your checkout
source "$RFD4_REPO/.venv/bin/activate"
cd "$RFD4_REPO"
cat > "$OUT/smoke_binder.yaml" <<'YAML'
binder:
  input: docs/examples/inputs/5o45_cropped.cif.zst
  contig: "50-120, /0, A17-131"
  conditions: {C_CRDGRP: {1: [{select: "sel('A/*/*/[N,CA,C,O]')"}]},
               C_SEQ: {true: [{select: "sel('A')"}]}}
YAML
rfd4-generate inference_experiment=quickstart inputs="$OUT/smoke_binder.yaml" \
    hf_ckpt_path=nvidia/RFD4-Proteina out_dir="$OUT"

# ---- Route B (Apptainer SIF): comment out Route A above and use this instead ----
# SIF="$HOME/.episcaf-pilot/rfd4-proteina-cu13.sif"   # <-- set to your SIF
# apptainer run --nv -v "${HF_HOME:-$HOME/.cache/huggingface}:/opt/rfd4/cache/huggingface" \
#     "$SIF" rfd4-generate inference_experiment=quickstart inputs=... hf_ckpt_path=nvidia/RFD4-Proteina out_dir="$OUT"

echo "=== smoke output ==="
find "$OUT" -maxdepth 3 -type f | head -40
echo "RFD4 smoke OK if a design PDB/CIF is listed above."
