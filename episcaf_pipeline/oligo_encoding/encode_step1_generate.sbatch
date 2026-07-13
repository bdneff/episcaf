#!/bin/bash
#SBATCH --time=4:00:00
#SBATCH --mem=50G
#SBATCH -c 4
#SBATCH --job-name=dp4_encode1
#SBATCH --partition=compute
##SBATCH --mail-type=END,FAIL          # uncomment + set your address to get notified
##SBATCH --mail-user=you@example.org
#
# ============================================================================
# DP4 oligo encoding -- STEP 1 of 2: generate candidate DNA encodings
# ============================================================================
# For every peptide in the input, sample many codon-level DNA encodings, weighted
# by the codon table and targeting a GC ratio. Writes candidate sequences (OUT_SEQS)
# and their codon-usage ratios (OUT_RATIO); step 2 then scores these with the neural
# net and keeps the best per peptide. This mirrors the DP3 recipe exactly -- only the
# input file changes for DP4.
#
# RUN FROM a working directory that contains the codon weights and the named-peptides
# input file, e.g.:
#     cd <rundir>                       # holds codon_weights_updated.csv + DP4_named_peptides.csv
#     sbatch .../encode_step1_generate.sbatch
# Override any default by exporting the variable first, e.g.  INPUT=DP4_named_peptides.csv
#
# Encoder tool: ekelley's cluster install (referenced, not vendored). If it moves,
# set TOOL_DIR. See README.md in this directory.
# ============================================================================
set -euo pipefail

# conda activate pepseq_encoding   # the env ekelley used for the encoder

TOOL_DIR=${TOOL_DIR:-/home/ekelley/bin/Library-Design/oligo_encoding}
LIB=${LIB:-DP4}
INPUT=${INPUT:-${LIB}_named_peptides.csv}      # name,seq  (NO header) -- from the 06_library export
WEIGHTS=${WEIGHTS:-codon_weights_updated.csv}  # aa,codon,weight,index  (same table as DP3)
OUT_SEQS=${OUT_SEQS:-out_seqs}                 # -> step 2 -s
OUT_RATIO=${OUT_RATIO:-output_ratio}           # -> step 2 -r
TRIALS=${TRIALS:-10000}    # candidate encodings generated per peptide   (DP3: 10000)
N_OUT=${N_OUT:-300}        # encodings written per peptide                (DP3: 300)
GC=${GC:-0.55}             # target GC ratio                             (DP3: 0.55)
CORES=${CORES:-2}          # encoder threads                             (DP3: 2)

# Encoder binary name differs by install: ekelley's compiled it as `main`; the GitHub master
# Makefile (`make optimized`) produces `oligo_encoding`. Auto-detect, or set BIN to override.
BIN=${BIN:-}
if [ -z "${BIN}" ]; then
    for c in "${TOOL_DIR}/main" "${TOOL_DIR}/oligo_encoding"; do
        [ -x "${c}" ] && BIN="${c}" && break
    done
fi
[ -n "${BIN}" ] || { echo "ERROR: no encoder binary (main / oligo_encoding) in ${TOOL_DIR}; set TOOL_DIR or BIN" >&2; exit 1; }

echo "HOST=$(hostname)  BIN=${BIN}  INPUT=${INPUT}"
"${BIN}" \
    -i "${INPUT}" \
    -p "${WEIGHTS}" \
    -s "${OUT_SEQS}" \
    -r "${OUT_RATIO}" \
    -t "${TRIALS}" \
    -n "${N_OUT}" \
    -g "${GC}" \
    -c "${CORES}"
echo "done -> ${OUT_SEQS}, ${OUT_RATIO}  (feed these to encode_step2_select.sbatch)"
