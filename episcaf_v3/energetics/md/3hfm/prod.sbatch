#!/bin/bash
#SBATCH --job-name=3hfm_prod
#SBATCH --partition=gpu-v100
#SBATCH --gres=gpu:1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=24G
#SBATCH --time=12:00:00
# 3HFM holo production (20 ns). Run AFTER equil.sbatch has written out/md.tpr + EQUIL_OK:
#   sbatch --dependency=afterok:<equil_jobid> prod.sbatch
set -e
module load Gromacs
export GMX_MAXBACKUP=-1
cd "$SLURM_SUBMIT_DIR/out"
echo "=== 3HFM HOLO PROD start $(date) ==="
gmx mdrun -deffnm md -nb gpu -ntmpi 1 -ntomp 8
echo "=== 3HFM HOLO PROD done $(date) ==="
