#!/usr/bin/env bash
set -euo pipefail

ROOT="/tgen_labs/altin/alphafold3/workspace/episcaf-experiments/data/sourced_antibody_v1/no_antibody"
OUTCSV="../runs/dp2_lawson_validation_full.csv"

sbatch --partition=compute \
       --nodes=1 \
       --ntasks=1 \
       --cpus-per-task=8 \
       --mem=32G \
       --time=04:00:00 \
       --job-name=dp2_lawson_val \
       --output=dp2_lawson_val.%j.out \
       --wrap="source ~/.bashrc && \
               conda activate ~/rfd3/env/rfd3_py312 && \
               cd $(pwd)/scripts && \
               python validate_dp2_rmsd_full.py \
                 --dp2 ../datasets/dp2.parquet \
                 --root ${ROOT} \
                 --out_csv ${OUTCSV} \
                 --progress_every 500"
echo "Submitted. Output will be in dp2_lawson_val.<jobid>.out"
