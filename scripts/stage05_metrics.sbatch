#!/usr/bin/env bash
#SBATCH -J stage05_metrics
#SBATCH -p compute            # <-- SET to your CPU partition (run `sinfo -s` to list). This is
                              #     CPU-only (gemmi/numpy/scipy); do NOT use a GPU partition.
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH -t 6:00:00
#SBATCH -o err_out_files/stage05_metrics_%j.out
#SBATCH -e err_out_files/stage05_metrics_%j.err
# NB: err_out_files/ must exist before submit (SLURM opens these before the script runs).

# Single source-of-truth metrics extraction for the dual-island run. ONE job, many cores:
# stage05_extract_metrics.py shards the ~46 epitopes across --workers processes (each loads
# its native complex once), so an array adds nothing here -- just give it the cores.
#
# Walltime: the 2-epitope smoke (3,840 designs) finishes in ~1-2 min; the full 111k is ~30x
# that and I/O-bound on the networked FS, so ~30-90 min in practice. 6h is generous margin
# (sizing from throughput + margin, per the MPNN-timeout lesson in docs/WORKFLOW.md). Not
# checkpointed -- a kill restarts from scratch, which is fine at this runtime.
#
# Usage:
#   sbatch scripts/stage05_metrics.sbatch <run_dir> [out_parquet] [native_dir] [ledger_csv]
# e.g.
#   sbatch scripts/stage05_metrics.sbatch /scratch/bneff/episcaf/runs/dual_island_rfd3

set -euo pipefail
module load Mamba/25.3.0-3
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate ~/rfd3/env/rfd3_py312

RUN_DIR="${1:?usage: sbatch scripts/stage05_metrics.sbatch <run_dir> [out] [native_dir] [ledger]}"
RUN_DIR="$(python3 -c "import os,sys; print(os.path.abspath(sys.argv[1]))" "${RUN_DIR}")"
OUT="${2:-${RUN_DIR}/05_analysis/metrics_dual_island.parquet}"
NATIVE_DIR="${3:-/tgen_labs/altin/alphafold3/workspace/episcaf-experiments/data/abdb/complex_pdbfiles/cleaned}"
LEDGER="${4:-results/dual_island_designs.csv}"

echo "HOST=$(hostname)  CPUS=${SLURM_CPUS_PER_TASK}"
echo "RUN_DIR=${RUN_DIR}"
echo "OUT=${OUT}"

python3 scripts/stage05_extract_metrics.py \
    --af3_out_dir  "${RUN_DIR}/04_af3/outputs" \
    --mpnn_pdb_dir "${RUN_DIR}/03_mpnn/mpnn_pdbs" \
    --ledger       "${LEDGER}" \
    --native_dir   "${NATIVE_DIR}" \
    --out          "${OUT}" \
    --workers      "${SLURM_CPUS_PER_TASK}"

echo "stage05 metrics complete -> ${OUT}"
