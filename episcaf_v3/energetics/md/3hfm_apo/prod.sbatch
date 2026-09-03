#!/bin/bash
#SBATCH --job-name=3hfm_apo_prod
#SBATCH --partition=gpu-v100
#SBATCH --gres=gpu:1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=24G
#SBATCH --time=8:00:00
# 3HFM apo production (20 ns, matches holo). After equil writes out/md.tpr + EQUIL_OK:
#   sbatch --dependency=afterok:<equil_jobid> prod.sbatch
set -e
module load Gromacs
export GMX_MAXBACKUP=-1
cd "$SLURM_SUBMIT_DIR/out"
echo "=== 3HFM APO PROD start $(date) ==="
gmx mdrun -deffnm md -nb gpu -ntmpi 1 -ntomp 8
echo "=== 3HFM APO PROD done $(date) ==="
