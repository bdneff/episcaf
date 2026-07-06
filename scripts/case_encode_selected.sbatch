#!/bin/bash
#SBATCH --job-name=case_encode
#SBATCH --time=02:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=2
#SBATCH --output=logs/case_encode_%j.out
#SBATCH --error=logs/case_encode_%j.err
# Rebuild the case-encoded scaffoldEPITOPE for the selected whole-epitope designs (C1 + C5),
# healing the lost token->design mapping by matching each design's sequence to its dp2 row.
# Reads design PDBs -> needs the cluster. Submit from /scratch/bneff/episcaf after `git pull`:
#     mkdir -p logs && sbatch scripts/case_encode_selected.sbatch
# Walltime sized ~generously (reading ~4.2k small PDBs + one dp2 load ~ 15 min); see memory
# cluster-job-walltime. Output feeds build_c6_mutants.py and library assembly.
set -euo pipefail
source ~/.bashrc
conda activate ~/rfd3/env/rfd3_py312

WS=/tgen_labs/altin/alphafold3/workspace/episcaf_v2_bneff
cd /scratch/bneff/episcaf
mkdir -p logs results

# If the mpnn_pdb paths in the selection are stale (older repo_refactored home path), enable a
# prefix remap to the persistent workspace copy, e.g.:
#   REMAP="--pdb-remap /home/bneff/rfd3/repo_refactored/runs/run_rfd3_mpnn:$WS/runs/run_rfd3_mpnn"
REMAP="${REMAP:-}"

for pair in \
  "results/dp4_C1_whole_epitope_ranked.top20.csv:results/dp4_C1_scaffoldEPITOPE.csv" \
  "results/dp4_C5_titration.csv:results/dp4_C5_scaffoldEPITOPE.csv" ; do
  sel="${pair%%:*}"; out="${pair##*:}"
  echo "=== case-encoding $sel -> $out ==="
  python3 scripts/case_encode_selected.py \
      --selection "$sel" \
      --dp2 "$WS/datasets/dp2.parquet" \
      --id-col id --pdb-col mpnn_pdb --target-col id \
      $REMAP \
      --out "$out"
done
echo "DONE. Commit results/dp4_C{1,5}_scaffoldEPITOPE.csv and push, then Claude pulls them back."
