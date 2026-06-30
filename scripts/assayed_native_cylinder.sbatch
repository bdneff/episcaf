#!/usr/bin/env bash
#SBATCH -J assayed_native_cyl
#SBATCH -p compute            # CPU-only (gemmi/numpy/scipy); set to your CPU partition
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH -t 1:00:00
#SBATCH -o err_out_files/assayed_native_cyl_%j.out
#SBATCH -e err_out_files/assayed_native_cyl_%j.err
# NB: err_out_files/ must exist before submit (SLURM opens these before the script runs).

# Native-aware cylinder for the 403 ASSAYED scaffold designs (the DP3 binding set), so the
# cylinder surrogate can be correlated against experimental binding. Identical carve to
# dp3_native_cylinder.py; only ~400 designs, so minutes. AF3 dirs are resolved directly by the
# assay hash (no tree walk).
#
# PRECONDITIONS:
#   - results/assayed_cylinder_worklist.csv exists (python3 scripts/build_assayed_cylinder_worklist.py)
#   - the assayed AF3 outputs are at <af3_root>/<assay_scaffolded_epitope_id>/<hash>_model.cif.gz
#     (verified 403/403 under sourced_antibody_v1/no_antibody/af3_predictions).
#
# Usage (defaults already point at sourced_antibody_v1/no_antibody):
#   sbatch scripts/assayed_native_cylinder.sbatch [af3_root] [worklist] [dp2_or_403.parquet] [native_dir] [out.csv]

set -euo pipefail
module load Mamba/25.3.0-3
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate ~/rfd3/env/rfd3_py312

BASE=/tgen_labs/altin/alphafold3/workspace/episcaf-experiments/data/sourced_antibody_v1/no_antibody
AF3_ROOT="${1:-${BASE}/af3_predictions}"
WORKLIST="${2:-results/assayed_cylinder_worklist.csv}"
DP2="${3:-${BASE}/assay_scaffold_simple_metrics_403.parquet}"   # has all 403 hashes + epitope indices
NATIVE_DIR="${4:-/tgen_labs/altin/alphafold3/workspace/episcaf-experiments/data/abdb/complex_pdbfiles/cleaned}"
OUT="${5:-results/assayed_native_cyl.csv}"

echo "HOST=$(hostname)"
echo "AF3_ROOT=${AF3_ROOT}"
echo "OUT=${OUT}"

python3 scripts/assayed_native_cylinder.py \
    --worklist     "${WORKLIST}" \
    --dp2_parquet  "${DP2}" \
    --native_dir   "${NATIVE_DIR}" \
    --af3_root     "${AF3_ROOT}" \
    --out_csv      "${OUT}" \
    --exclude_dist 1.0

echo "done -> ${OUT}"
