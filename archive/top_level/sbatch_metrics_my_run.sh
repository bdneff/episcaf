#!/usr/bin/env bash
#SBATCH -J metrics_my_run
#SBATCH -p compute
#SBATCH -c 16
#SBATCH --mem=32G
#SBATCH -t 06:00:00
#SBATCH -o metrics_my_run.%j.out
#SBATCH -e metrics_my_run.%j.err

set -eo pipefail

# ---- paths ----
REPO="/home/bneff/rfd3/repo_refactored"
RUN_DIR="runs/run_test_rfd3_nompmn"
DP2="datasets/dp2.parquet"
TRUE_DIR="/tgen_labs/altin/alphafold3/workspace/episcaf-experiments/data/abdb/complex_pdbfiles/cleaned"
OUT_CSV="$RUN_DIR/04_filter/metrics_full.csv"

cd "$REPO"

source ~/.bashrc
conda activate /home/bneff/rfd3/env/rfd3_py312

export PYTHONPATH="$REPO:$PYTHONPATH"

python -u compute_metrics.py run \
  --run_dir      "$RUN_DIR" \
  --dp2_parquet  "$DP2" \
  --true_dir     "$TRUE_DIR" \
  --out_csv      "$OUT_CSV" \
  --clash_cutoff 4.0

echo "DONE: wrote $OUT_CSV"
