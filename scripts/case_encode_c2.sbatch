#!/bin/bash
#SBATCH --job-name=case_encode_c2
#SBATCH --time=02:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=2
#SBATCH --output=logs/case_encode_c2_%j.out
#SBATCH --error=logs/case_encode_c2_%j.err
# Case-encode C2 (single-island): reads the dual-island run's af3 outputs; epitope positions come from
# the local n_flank/island_size in the selection. Submit from /scratch/bneff/episcaf after `git pull`:
#     mkdir -p logs && sbatch scripts/case_encode_c2.sbatch
# The af3_dir paths point to /scratch/.../dual_island_rfd3/04_af3/outputs -- if that run was swept off
# /scratch, first `ls` one af3_dir; if gone, point AF3REMAP at the archived/relocated copy, e.g.:
#     AF3REMAP="--af3-remap /scratch/bneff/episcaf/runs/dual_island_rfd3:$WS/runs/dual_island_rfd3"
#     sbatch --export=ALL,AF3REMAP="$AF3REMAP" scripts/case_encode_c2.sbatch   # (no spaces issue: single arg)
set -eo pipefail
source ~/.bashrc
conda activate ~/rfd3/env/rfd3_py312
cd /scratch/bneff/episcaf
mkdir -p logs results
AF3REMAP="${AF3REMAP:-}"
echo "AF3REMAP: ${AF3REMAP:-<none>}"
python3 scripts/case_encode_c2.py \
    --selection results/dp4_C2_single_island_ranked.top20.csv \
    $AF3REMAP \
    --out results/dp4_C2_scaffoldEPITOPE.csv
echo "DONE. Then commit results/dp4_C2_scaffoldEPITOPE.csv + push."
echo "NOTE: also re-run 'sbatch scripts/case_encode_selected.sbatch' to refresh C5's scaffoldEPITOPE"
echo "      for the 56-set (the committed dp4_C5_scaffoldEPITOPE.csv predates the C5 re-run)."
