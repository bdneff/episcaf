#!/bin/bash
#SBATCH --job-name=case_encode_c2
#SBATCH --time=02:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=2
#SBATCH --output=logs/case_encode_c2_%j.out
#SBATCH --error=logs/case_encode_c2_%j.err
# Case-encode C2 (single-island): reads the dual-island run's af3 outputs; epitope positions come from
# the local n_flank/island_size in the selection. Submit from the repo checkout after `git pull`:
#     mkdir -p logs && sbatch scripts/case_encode_c2.sbatch
set -eo pipefail
source ~/.bashrc
conda activate ~/rfd3/env/rfd3_py312

WS=${WS:-/tgen_labs/altin/alphafold3/workspace/episcaf_v2_bneff}
cd "${REPO:-/scratch/bneff/episcaf}"        # the git checkout: disposable, re-clone anywhere
mkdir -p logs results

# The selection's af3_dir paths record the OLD /scratch home of the dual-island run; the run itself was
# copied to durable $WS on 2026-07-17 (/scratch is swept on ~30 days). Remap by default, same as
# case_encode_selected.sbatch does for --pdb-remap. Override with AF3REMAP="" to read the raw paths.
if [ -z "${AF3REMAP+set}" ]; then
  AF3REMAP="--af3-remap /scratch/bneff/episcaf/runs/dual_island_rfd3:$WS/runs/dual_island_rfd3"
fi
echo "AF3REMAP: ${AF3REMAP:-<none>}"
python3 scripts/case_encode_c2.py \
    --selection results/dp4_C2_single_island_ranked.top20.csv \
    $AF3REMAP \
    --out results/dp4_C2_scaffoldEPITOPE.csv
echo "DONE. Then commit results/dp4_C2_scaffoldEPITOPE.csv + push."
echo "NOTE: also re-run 'sbatch scripts/case_encode_selected.sbatch' to refresh C5's scaffoldEPITOPE"
echo "      for the 56-set (the committed dp4_C5_scaffoldEPITOPE.csv predates the C5 re-run)."
