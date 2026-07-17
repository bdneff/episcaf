#!/usr/bin/env bash
#SBATCH -J copy_run_to_tgen
#SBATCH -p data-mover        # gemini data-transfer partition (14-day limit, node gemini-data2)
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH -t 12:00:00
#SBATCH -o err_out_files/copy_run_to_tgen_%j.out
#SBATCH -e err_out_files/copy_run_to_tgen_%j.err
# NB: err_out_files/ must exist before submit (SLURM opens these before the script runs).

# Copy a scratch run directory to durable /tgen_labs storage before scratch autopurges.
#
# SAFE BY CONSTRUCTION:
#   - COPY ONLY. No --delete anywhere, so nothing on /tgen_labs is ever removed (project rule).
#   - Does NOT touch the scratch source. Freeing scratch is a SEPARATE manual step you run
#     yourself AFTER verifying this copy (see the printed instructions at the end).
#   - Resumable + idempotent: rsync skips already-copied files, so if this job hits the
#     walltime, just resubmit the same command and it continues where it left off.
#   - Excludes 04_af3/_task_inputs (transient per-task staging symlinks; not worth copying).
#
# Usage:
#   sbatch scripts/copy_run_to_tgen.sbatch [SRC_RUN_DIR] [DEST_PARENT]
# Defaults copy the dual-island run into the episcaf_v2_bneff workspace runs/ dir:
#   sbatch scripts/copy_run_to_tgen.sbatch

set -euo pipefail

SRC="${1:-/scratch/bneff/episcaf/runs/dual_island_rfd3}"
DEST_PARENT="${2:-/tgen_labs/altin/alphafold3/workspace/episcaf_v2_bneff/runs}"

SRC="${SRC%/}"
NAME="$(basename "${SRC}")"
DEST="${DEST_PARENT%/}/${NAME}"

echo "HOST=$(hostname)"
echo "SRC : ${SRC}"
echo "DEST: ${DEST}"
[[ -d "${SRC}" ]] || { echo "ERROR: source ${SRC} does not exist"; exit 1; }
mkdir -p "${DEST}"

# -a archive (perms/times/symlinks), -W whole-file (no delta math; faster local->networked),
# -h human sizes, --info progress+stats. Trailing slash on SRC copies its CONTENTS into DEST.
echo "=== rsync starting $(date) ==="
rsync -aWh --info=progress2,stats2 --exclude '_task_inputs/' "${SRC}/" "${DEST}/"
echo "=== rsync finished $(date) ==="

# Integrity check: a dry-run second pass should find nothing left to transfer.
# NB: grep exits 1 when it matches nothing -- which is the SUCCESS case here (0 files to transfer).
# Under `set -o pipefail` + `set -e` that non-zero would kill the script right here, before the verdict
# prints, so a perfect copy would look like a silent death. `|| true` swallows the no-match exit; the
# `wc -l` still yields 0. (This bit the first several copy jobs: they all ended at the "verifying" line.)
echo "=== verifying (dry-run; files still differing should be 0) ==="
REMAIN="$( { rsync -aWn --itemize-changes --exclude '_task_inputs/' "${SRC}/" "${DEST}/" \
             | grep -E '^>f' || true; } | wc -l | tr -d ' ')"
echo "files still needing transfer: ${REMAIN}"

echo "=== sizes (will differ by the excluded _task_inputs) ==="
du -sh "${SRC}" "${DEST}" 2>/dev/null || true

if [[ "${REMAIN}" -eq 0 ]]; then
    echo "COPY VERIFIED COMPLETE."
    echo "To free scratch AFTER you have eyeballed ${DEST}, run (manually, not in this job):"
    echo "    rm -rf ${SRC}"
else
    echo "WARNING: ${REMAIN} files still differ -- resubmit this job to finish before deleting anything."
fi
