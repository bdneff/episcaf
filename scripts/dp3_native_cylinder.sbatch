#!/usr/bin/env bash
#SBATCH -J dp3_native_cyl
#SBATCH -p compute            # CPU-only (gemmi/numpy/scipy); set to your CPU partition
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH -t 8:00:00
#SBATCH -o dp3_native_cyl_%j.out
#SBATCH -e dp3_native_cyl_%j.err

# Recompute cylinder_native_aware for the FULL DP3 set (no 2.5A gate), so Figure 7 / the
# overlays can show plain AND native-aware side by side for both pools. Single process with the
# native complex cached per epitope; ~150k AF3 CIF loads, so ~1.5-3h. (Not parallelized -- if
# that is too slow we can shard later.)
#
# PRECONDITIONS (check first):
#   - the DP3 run's AF3 outputs still exist on disk at the af3_dir paths in metrics_cylinder_full.csv
#   - that CSV has an `af3_dir` column (the archived run did)
#
# Usage:
#   sbatch scripts/dp3_native_cylinder.sbatch <metrics_cylinder_full.csv> <dp2.parquet> <native_dir> <out.csv>
# Defaults use configs/paths.py locations:

set -euo pipefail
module load Mamba/25.3.0-3
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate ~/rfd3/env/rfd3_py312

WS=/tgen_labs/altin/alphafold3/workspace/episcaf_v2_bneff
METRICS="${1:-${WS}/runs/run_rfd3_mpnn/04_filter/metrics_cylinder_full.csv}"
DP2="${2:-${WS}/datasets/dp2.parquet}"
NATIVE_DIR="${3:-/tgen_labs/altin/alphafold3/workspace/episcaf-experiments/data/abdb/complex_pdbfiles/cleaned}"
OUT="${4:-${WS}/runs/run_rfd3_mpnn/04_filter/metrics_native_cyl_full.csv}"

echo "HOST=$(hostname)"
echo "METRICS=${METRICS}"
echo "OUT=${OUT}"

python3 scripts/dp3_native_cylinder.py \
    --metrics_csv  "${METRICS}" \
    --dp2_parquet  "${DP2}" \
    --native_dir   "${NATIVE_DIR}" \
    --out_csv      "${OUT}" \
    --exclude_dist 1.0
    # no --gate: process ALL designs (the old 44,862-design file used --gate 2.5)

echo "done -> ${OUT}"
