#!/bin/bash -l
#SBATCH --job-name=3hfm_gbsa
#SBATCH --partition=gpu-v100
#SBATCH --gres=gpu:1
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=3:00:00
# Per-residue MM-GBSA on the 3HFM holo trajectory (Decision D2), MPI-parallel over frames (~8x).
# CPU-only work; the GPU idles -- `--gres=gpu:1` just lets the job land on the gpu-v100 partition
# (swap to a CPU partition if Gemini has one). Uses the conda gmx_MMPBSA env + its own gmx/sander.
#
# Submit FROM the holo run dir (so $SLURM_SUBMIT_DIR is md/3hfm/out, which has md.tpr/md.xtc/
# topol.top/gbsa_index.ndx):
#     cd /scratch/bneff/episcaf_run/episcaf_v3/energetics/md/3hfm/out
#     sbatch /scratch/bneff/episcaf_run/episcaf_v3/energetics/mmgbsa/run_mmgbsa.sbatch

# activate the env (sets AMBERHOME etc. that sander needs) -- the likely first-run snag
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate gmxMMPBSA

set -e
cd "$SLURM_SUBMIT_DIR"
rm -rf _GMXMMPBSA_*                 # clear any leftover temp from an interrupted run
echo "=== 3HFM MM-GBSA start $(date) on $SLURM_NTASKS ranks ==="
mpirun -np "$SLURM_NTASKS" gmx_MMPBSA -O \
    -i "$SLURM_SUBMIT_DIR/../../../mmgbsa/mmpbsa.in" \
    -cs md.tpr -ci gbsa_index.ndx -cg 17 18 \
    -ct md.xtc -cp topol.top \
    -o FINAL_RESULTS_MMGBSA.dat -eo FINAL_RESULTS_MMGBSA.csv \
    -do FINAL_DECOMP_MMGBSA.dat -deo FINAL_DECOMP_MMGBSA.csv -nogui
echo "=== 3HFM MM-GBSA done $(date) ==="
