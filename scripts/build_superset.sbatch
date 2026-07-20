#!/usr/bin/env bash
#SBATCH -J build_superset
#SBATCH -p compute            # CPU-only (pandas + gemmi reads of AF3 cifs); no GPU needed
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G             # C1 is 140,716 rows; C3 82,712; C2 111,322 -- all held in pandas
#SBATCH -t 8:00:00
#SBATCH -o err_out_files/build_superset_%j.out
#SBATCH -e err_out_files/build_superset_%j.err
# NB: err_out_files/ must exist before submit (SLURM opens these before the script runs).

# Build the ALL-DESIGNS superset of the DP4 scaffolded arms in one pass, on the cluster
# (John, 2026-07-16: "a super-set version of dp4_library.csv ... include all designs").
#
# WHY THE CLUSTER: C2's metrics only exist here, and --sequences reads each passing design's AF3
# chain A from runs/*/04_af3/outputs. Those runs were copied to durable $WS on 2026-07-17, so the
# job reads them there via --af3-remap (the metrics still record the old /scratch paths).
#
# Output: ONE dp4_superset.csv (~335k rows) written straight to durable $WS.
#
# Usage:
#   mkdir -p err_out_files
#   sbatch scripts/build_superset.sbatch [M1] [M2] [M3] [OUT]
# where M1/M2/M3 are the C1/C2/C3 metrics files. The defaults read all three from the durable
# workspace $WS -- nothing needs staging. Pass positional args to point at another copy.

set -euo pipefail

WS=${WS:-/tgen_labs/altin/alphafold3/workspace/episcaf_v2_bneff}
REPO=${REPO:-/scratch/bneff/episcaf}

# All three verified 2026-07-16 to be byte-identical (md5) to the local copies the selection ran on:
#   C1  metrics_whole_epitope.csv     md5 dea529a2..  140,716 designs  (local name adds a _103 suffix)
#   C3  metrics_12mer.csv             md5 e199faf6..   82,712 designs
# C1 and C2 were read from /scratch when this job first ran; both runs have since been copied to $WS
# (2026-07-17), so the defaults now point there -- /scratch is swept on ~30 days. Same bytes either way.
M1=${1:-${WS}/runs/whole_epitope_rfd3/05_analysis/metrics_whole_epitope.csv}
M2=${2:-${WS}/runs/dual_island_rfd3/05_analysis/metrics_dual_island.parquet}
M3=${3:-${WS}/run_12mer_scaffolding/06_score/metrics_12mer.csv}
OUT=${4:-${WS}/dp4_superset.csv}
LIB=${LIB:-${REPO}/data/libraries/dp4_library.csv}

set +u
source "$(conda info --base)/etc/profile.d/conda.sh" 2>/dev/null || true
conda activate "${CONDA_ENV:-${HOME}/rfd3/env/rfd3_py312}"   # has pandas + gemmi
set -u

echo "HOST=$(hostname)"
for f in "${M1}" "${M2}" "${M3}" "${LIB}"; do
    [[ -e "${f}" ]] || { echo "ERROR: missing input ${f}"; exit 1; }
done
mkdir -p "$(dirname "${OUT}")"
cd "${REPO}"

# The C1/C2 metrics bake in ABSOLUTE af3_dir paths pointing at the /scratch run dirs they were built
# from. Those runs now live under $WS, so remap the prefix; without it --sequences silently resolves
# nothing once /scratch is swept (the script now hard-fails on that rather than writing blanks).
REMAP=${REMAP:-/scratch/bneff/episcaf/runs:${WS}/runs}

# C1 writes (overwrites) the file; C2 and C3 append. --sequences reads AF3 chain A for every
# global-passing design that did not already ship (selected designs get their sequence copied
# verbatim from dp4_library.csv, so the superset agrees with the library by construction).
echo "=== C1 $(date) ==="
python scripts/stage06_superset.py --component C1 --metrics-csv "${M1}" \
    --library "${LIB}" --sequences --af3-remap "${REMAP}" --out "${OUT}"

echo "=== C2 $(date) ==="
python scripts/stage06_superset.py --component C2 --metrics-csv "${M2}" \
    --library "${LIB}" --sequences --af3-remap "${REMAP}" --append --out "${OUT}"

echo "=== C3 $(date) ==="
python scripts/stage06_superset.py --component C3 --metrics-csv "${M3}" \
    --library "${LIB}" --append --out "${OUT}"

echo "=== done $(date) ==="
wc -l "${OUT}"
echo "superset -> ${OUT}"
