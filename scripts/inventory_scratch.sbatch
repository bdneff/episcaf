#!/usr/bin/env bash
#SBATCH -J inventory_scratch
#SBATCH -p compute
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH -t 4:00:00
#SBATCH -o err_out_files/inventory_scratch_%j.out
#SBATCH -e err_out_files/inventory_scratch_%j.err
# NB: err_out_files/ must exist before submit (SLURM opens these before the script runs).

# Characterize what is on scratch vs what is already durable on /tgen_labs, so the
# migration (scripts/copy_run_to_tgen.sbatch) can be pointed only at what is actually
# missing. Walking these trees takes many minutes -- run it as a job, not in a terminal.
#
# READ-ONLY. Copies nothing, deletes nothing. It only measures and writes one report.
#
# Usage:
#   mkdir -p err_out_files
#   sbatch scripts/inventory_scratch.sbatch [SCRATCH_RUNS] [WS]
# Report lands at ${WS}/migration/inventory_<jobid>.txt (durable) and is echoed to the
# job's .out file.

set -euo pipefail

SCRATCH_RUNS="${1:-/scratch/bneff/episcaf/runs}"
WS="${2:-/tgen_labs/altin/alphafold3/workspace/episcaf_v2_bneff}"

REPORT_DIR="${WS}/migration"
mkdir -p "${REPORT_DIR}"
REPORT="${REPORT_DIR}/inventory_${SLURM_JOB_ID:-manual}.txt"

# Everything below goes to both the report and the job log.
exec > >(tee "${REPORT}") 2>&1

echo "=== episcaf scratch inventory ==="
echo "host    : $(hostname)"
echo "date    : $(date)"
echo "scratch : ${SCRATCH_RUNS}"
echo "ws      : ${WS}"
echo

echo "=== scratch runs (size, largest last) ==="
du -sh "${SCRATCH_RUNS}"/* 2>/dev/null | sort -h || echo "(none)"
echo
echo "scratch runs TOTAL:"
du -sh "${SCRATCH_RUNS}" 2>/dev/null || true
echo

echo "=== already durable in ${WS}/runs ==="
du -sh "${WS}"/runs/* 2>/dev/null | sort -h || echo "(none)"
echo

# The decision this report exists to answer: for each scratch run, is there a
# same-named dir in $WS/runs, and does rsync think anything still needs to move?
echo "=== per-run status: MISSING / needs N files / IN SYNC ==="
for src in "${SCRATCH_RUNS}"/*; do
    [[ -d "${src}" ]] || continue
    name="$(basename "${src}")"
    dest="${WS}/runs/${name}"
    if [[ ! -d "${dest}" ]]; then
        echo "  MISSING       ${name}  -> not in \$WS/runs at all; must copy"
        continue
    fi
    remain="$(rsync -aWn --itemize-changes --exclude '_task_inputs/' "${src}/" "${dest}/" 2>/dev/null \
              | grep -cE '^>f' || true)"
    if [[ "${remain}" -eq 0 ]]; then
        echo "  IN SYNC       ${name}  -> safe to rm from scratch"
    else
        echo "  NEEDS ${remain} files  ${name}  -> resume copy"
    fi
done
echo

echo "=== stray scratch junk (loose logs at /scratch/bneff top level) ==="
ls -lah /scratch/bneff/*.out /scratch/bneff/*.err 2>/dev/null || echo "(none)"
echo

echo "=== report written to ${REPORT} ==="
