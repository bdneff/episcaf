#!/bin/bash
#SBATCH --time=4:00:00
#SBATCH --mem=50G
#SBATCH -c 4
#SBATCH --job-name=dp4_encode2
#SBATCH --partition=compute
##SBATCH --mail-type=END,FAIL
##SBATCH --mail-user=you@example.org
#
# ============================================================================
# DP4 oligo encoding -- STEP 2 of 2: neural-net selection of the best encoding
# ============================================================================
# Score the step-1 candidate encodings with the pretrained H2O deep-learning model
# (it learned which encodings synthesize/perform best) and keep the single best
# encoding per peptide (smallest NN deviation). Mirrors the DP3 recipe.
#
# Run in the SAME working directory as step 1 (it reads out_seqs / output_ratio).
#     sbatch .../encode_step2_select.sbatch
# ============================================================================
set -euo pipefail

# Erin's env (pepseq_encoding.yml) carries the H2O runtime needed to load the model.
source "$(conda info --base)/etc/profile.d/conda.sh" 2>/dev/null || true
conda activate "${CONDA_ENV:-pepseq_encoding}"

TOOL_DIR=${TOOL_DIR:-/home/ekelley/bin/Library-Design/oligo_encoding}
MODEL=${MODEL:-${TOOL_DIR}/DeepLearning_model_R_1539970074840_1_20181019}   # the DP3 model
LIB=${LIB:-DP4}
OUT_SEQS=${OUT_SEQS:-out_seqs}        # from step 1
OUT_RATIO=${OUT_RATIO:-output_ratio}  # from step 1
OUT=${OUT:-${LIB}_best_encodings}     # selected best encoding(s) per peptide
SUBSAMPLE=${SUBSAMPLE:-300}           # candidates subsampled before scoring   (DP3: 300)
NN=${NN:-1}                           # encodings kept per peptide             (DP3: 1 = single best)

# NN selector script name differs by install: ekelley's is `encoding_with_nn.py`; the GitHub master
# ships `oligo_encoding.py`. Auto-detect, or set SEL to override. (If master's flags differ from the
# DP3 recipe below, the smoke test will surface it and we adjust -- the two are the same lineage.)
SEL=${SEL:-}
if [ -z "${SEL}" ]; then
    for c in "${TOOL_DIR}/encoding_with_nn.py" "${TOOL_DIR}/oligo_encoding.py"; do
        [ -f "${c}" ] && SEL="${c}" && break
    done
fi
[ -n "${SEL}" ] || { echo "ERROR: no NN selector (encoding_with_nn.py / oligo_encoding.py) in ${TOOL_DIR}; set TOOL_DIR or SEL" >&2; exit 1; }

echo "HOST=$(hostname)  SEL=${SEL}  MODEL=${MODEL}  OUT=${OUT}"
python "${SEL}" \
    -m "${MODEL}" \
    -s "${OUT_SEQS}" \
    -r "${OUT_RATIO}" \
    -o "${OUT}" \
    --subsample "${SUBSAMPLE}" \
    --read_per_loop 10 \
    -n "${NN}"
echo "done -> ${OUT}"

# -------------------------------------------------------------------------
# OPEN ITEM (see README "The order-file gap"): this writes ${OUT} (best encodings).
# The DP3 *order file* that goes to synthesis (DP3_order_file.csv, columns
# "Seq ID,nucleotide_encoding_with_twist_adapters") is a FURTHER reformat + Twist-
# adapter step that is NOT in ekelley's two scripts. Confirm with Erin how DP3 made
# it (likely encoding_with_nn.py's --adapter default already flanks the oligo, and
# the order file is a column rename of ${OUT}) before treating the run as complete.
# -------------------------------------------------------------------------
