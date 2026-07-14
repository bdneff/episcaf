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

# Conda's activate.d scripts reference unbound vars (MKL_INTERFACE_LAYER, ...), so relax `set -u`
# around activation. The env (pepseq_encoding.yml) carries the H2O runtime for the model.
set +u
source "$(conda info --base)/etc/profile.d/conda.sh" 2>/dev/null || true
conda activate "${CONDA_ENV:-pepseq_encoding}"
set -u

# REQUIRED (no default): the dir with the built `main` (or `oligo_encoding`), oligo_encoding.py, the
# model. Pass it INLINE on the sbatch line (same value as step 1).
TOOL_DIR=${TOOL_DIR:?set TOOL_DIR inline, e.g. TOOL_DIR=/home/bneff/Library-Design/oligo_encoding sbatch ...}
MODEL=${MODEL:-${TOOL_DIR}/DeepLearning_model_R_1539970074840_1_20181019}   # the DP3 model
LIB=${LIB:-DP4}
OUT_SEQS=${OUT_SEQS:-out_seqs}        # from step 1
OUT_RATIO=${OUT_RATIO:-output_ratio}  # from step 1
OUT=${OUT:-${LIB}_best_encodings}     # selected best encoding(s) per peptide
SUBSAMPLE=${SUBSAMPLE:-300}           # candidates subsampled before scoring   (DP3: 300)
NN=${NN:-1}                           # encodings kept per peptide             (DP3: 1 = single best)

# Twist adapters, PINNED EXPLICITLY -- do not rely on the tool's default. The two encoder installs
# ship DIFFERENT defaults, and the difference is silent: master's oligo_encoding.py defaults to the
# 19-mers CCTATACTTCCAAGGCGCA|GGTGACTCTCTGTCTTGGC, while the DP3 order file that was actually
# synthesized carries 20-mers (a leading A and a trailing T more). Left to the default, a master-built
# encoder emits 347-nt oligos where DP3 shipped 349. Verified against examples/DP3_order_file.csv:
# 6000/6000 rows are 349 nt = 20 + 309 (103 codons) + 20, and the 309-nt core translates exactly to
# the peptide. These are those 20-mers.
ADAPTER=${ADAPTER:-ACCTATACTTCCAAGGCGCA|GGTGACTCTCTGTCTTGGCT}

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
echo "ADAPTER=${ADAPTER}"
python "${SEL}" \
    -m "${MODEL}" \
    -s "${OUT_SEQS}" \
    -r "${OUT_RATIO}" \
    -o "${OUT}" \
    --subsample "${SUBSAMPLE}" \
    --read_per_loop 10 \
    --adapter "${ADAPTER}" \
    -n "${NN}"
echo "done -> ${OUT}"

# -------------------------------------------------------------------------
# The order file (RESOLVED -- there is no missing script). ${OUT} already carries the
# adapter-flanked oligo in its "Nucleotide Encoding w/ Adapters" column, so the synthesis
# order file is a two-column slice + rename of ${OUT}, not a further encoding step:
#     "Seq ID" -> Seq ID
#     "Nucleotide Encoding w/ Adapters" -> nucleotide_encoding_with_twist_adapters
# Emit it with scripts/stage07_order_file.py (which re-verifies length + translation).
# -------------------------------------------------------------------------
