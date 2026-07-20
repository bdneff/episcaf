#!/bin/bash
#SBATCH --job-name=case_encode
#SBATCH --time=02:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=2
#SBATCH --output=logs/case_encode_%j.out
#SBATCH --error=logs/case_encode_%j.err
# Rebuild the case-encoded scaffoldEPITOPE for the selected whole-epitope designs (C1 + C5),
# healing the lost token->design mapping by matching each design's sequence to its dp2 row.
# Reads design PDBs -> needs the cluster. Submit from the repo checkout after `git pull`:
#     mkdir -p logs && sbatch scripts/case_encode_selected.sbatch
# Walltime sized ~generously (reading ~4.2k small PDBs + one dp2 load ~ 15 min); see memory
# cluster-job-walltime. Output feeds build_c6_mutants.py and library assembly.
set -eo pipefail        # NOT -u: the cluster's /etc/bashrc uses an unbound var (BASHRCSOURCED)
source ~/.bashrc
conda activate ~/rfd3/env/rfd3_py312

WS=${WS:-/tgen_labs/altin/alphafold3/workspace/episcaf_v2_bneff}
cd "${REPO:-/scratch/bneff/episcaf}"        # the git checkout: disposable, re-clone anywhere
mkdir -p logs results

# dp2 that carries OUR rerun tokens (assay_scaffolded_epitope_id). Try the $WS copy first; if the log
# reports 0 token coverage, override: DP2=/home/bneff/rfd3/repo_refactored/datasets/dp2.parquet (run's dp2).
DP2="${DP2:-$WS/datasets/dp2.parquet}"
echo "DP2: $DP2"

# The selections' mpnn_pdb paths point to the OLD repo_refactored home dir (gone); the files live in
# the persistent workspace. Remap to $WS by default; override with REMAP="" to disable.
if [ -z "${REMAP:-}" ]; then
  REMAP="--pdb-remap /home/bneff/rfd3/repo_refactored/runs/run_rfd3_mpnn:$WS/runs/run_rfd3_mpnn"
fi
echo "REMAP: $REMAP"

for pair in \
  "results/dp4_C1_whole_epitope_ranked.top20.csv:results/dp4_C1_scaffoldEPITOPE.csv" \
  "results/dp4_C5_titration.csv:results/dp4_C5_scaffoldEPITOPE.csv" ; do
  sel="${pair%%:*}"; out="${pair##*:}"
  echo "=== case-encoding $sel -> $out ==="
  python3 scripts/case_encode_selected.py \
      --selection "$sel" \
      --dp2 "$DP2" \
      --token-col token --pdb-col mpnn_pdb --target-col id \
      $REMAP \
      --out "$out"
done
echo "DONE. Commit results/dp4_C{1,5}_scaffoldEPITOPE.csv and push, then Claude pulls them back."
