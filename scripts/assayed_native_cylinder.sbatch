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
# dp3_native_cylinder.py; only ~400 designs, so minutes. The slow part is indexing the AF3
# output tree once -- if AF3_ROOT is huge the rglob dominates; point it as deep as possible.
#
# PRECONDITIONS:
#   - results/assayed_cylinder_worklist.csv exists (python3 scripts/build_assayed_cylinder_worklist.py)
#   - the assayed designs' AF3 outputs exist under AF3_ROOT, each design's *_model.cif in a dir
#     whose path contains its assay_scaffolded_epitope_id hash. If keyed differently, edit the
#     python call to add --key_col <col>.
#
# Usage:
#   sbatch scripts/assayed_native_cylinder.sbatch <af3_root> [worklist] [dp2.parquet] [native_dir] [out.csv]

set -euo pipefail
module load Mamba/25.3.0-3
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate ~/rfd3/env/rfd3_py312

WS=/tgen_labs/altin/alphafold3/workspace/episcaf_v2_bneff
AF3_ROOT="${1:?give the root dir of the assayed-design AF3 outputs}"
WORKLIST="${2:-results/assayed_cylinder_worklist.csv}"
DP2="${3:-${WS}/datasets/dp2.parquet}"
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
