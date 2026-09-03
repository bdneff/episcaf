#!/bin/bash -l
#SBATCH --job-name=3hfm_gbsa
#SBATCH --partition=gpu-v100
#SBATCH --gres=gpu:1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=4:00:00
# Per-residue MM-GBSA on the 3HFM holo trajectory (Decision D2), SERIAL.
# We tried MPI (mpirun -np 8) and the ranks ran independently (mpi4py not sharing a communicator
# under SLURM -> COM.prmtop race). Serial is proven to work; mmpbsa.in uses interval=100 (~20 frames)
# for a fast first pass. For the final figure, sort out MPI or lower interval and give more walltime.
# CPU-only; GPU idles (--gres just lands it on gpu-v100). Uses the conda gmx_MMPBSA env.
#
# Submit FROM the holo run dir:
#   cd /scratch/bneff/episcaf_run/episcaf_v3/energetics/md/3hfm/out
#   sbatch /scratch/bneff/episcaf_run/episcaf_v3/energetics/mmgbsa/run_mmgbsa.sbatch

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate gmxMMPBSA

set -e
cd "$SLURM_SUBMIT_DIR"
rm -rf _GMXMMPBSA_*                 # clear leftover temp from the interrupted/failed runs
echo "=== 3HFM MM-GBSA (serial) start $(date) ==="
gmx_MMPBSA -O \
    -i "$SLURM_SUBMIT_DIR/../../../mmgbsa/mmpbsa.in" \
    -cs md.tpr -ci gbsa_index.ndx -cg 17 18 \
    -ct md.xtc -cp topol.top \
    -o FINAL_RESULTS_MMGBSA.dat -eo FINAL_RESULTS_MMGBSA.csv \
    -do FINAL_DECOMP_MMGBSA.dat -deo FINAL_DECOMP_MMGBSA.csv -nogui
echo "=== 3HFM MM-GBSA done $(date) ==="
