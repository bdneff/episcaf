#!/bin/bash
#SBATCH --job-name=rfd3_cif_to_pdb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=2:00:00
#SBATCH --output=logs/rfd3_mpnn/01_cif_to_pdb_%j.log

mkdir -p logs/rfd3_mpnn

/home/bneff/rfd3/env/rfd3_py312/bin/python scripts/01_rfd3_cif_to_fixed_pdb.py \
    --metrics_csv runs/run_test_rfd3_nompmn/04_filter/metrics_full.csv \
    --dp2_parquet datasets/dp2.parquet \
    --outdir      runs/run_rfd3_mpnn/01_fixed_pdbs \
    --n_workers   16

