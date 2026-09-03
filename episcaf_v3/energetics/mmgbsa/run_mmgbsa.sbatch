#!/bin/bash -l
#SBATCH --job-name=3hfm_gbsa
#SBATCH --partition=gpu-v100
#SBATCH --gres=gpu:1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=12:00:00
# Per-residue MM-GBSA on the 3HFM holo trajectory (Decision D2), SERIAL, 100 frames (interval=20).
# Denser resample of the SAME 20 ns trajectory (was 20 frames): tests whether Lys96 recovers and the
# false-positive zone empties. ~5x the 20-frame cost -> 12 h walltime. If this times out, raise
# interval back toward 100 (fewer frames) or, if the v100 partition allows, bump --time.
# Serial because MPI ranks ran independently under SLURM (mpi4py communicator issue).
# PBC FIX: gmx_MMPBSA does NOT correct PBC; a frame with the complex split across the box makes
# sander overflow (energy prints '****' -> parse crash). -merge all made the protein one molecule,
# so `gmx trjconv -pbc whole` reassembles it; we feed that whole trajectory as -ct.
# CPU-only; GPU idles. Uses the conda gmx_MMPBSA env.
#
# Submit FROM the holo run dir:
#   cd /scratch/bneff/episcaf_run/episcaf_v3/energetics/md/3hfm/out
#   sbatch /scratch/bneff/episcaf_run/episcaf_v3/energetics/mmgbsa/run_mmgbsa.sbatch

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate gmxMMPBSA

set -e
cd "$SLURM_SUBMIT_DIR"
# clear intermediates from any prior/failed run so this starts clean
rm -rf _GMXMMPBSA_* COM_traj_*.xtc COM.prmtop REC.prmtop LIG.prmtop reference.frc

# PBC-correct: make the (merged) protein whole in every frame
if [ ! -f md_whole.xtc ]; then
    printf '0\n' | gmx trjconv -s md.tpr -f md.xtc -o md_whole.xtc -pbc whole
fi

echo "=== 3HFM MM-GBSA (serial, PBC-whole) start $(date) ==="
gmx_MMPBSA -O \
    -i "$SLURM_SUBMIT_DIR/../../../mmgbsa/mmpbsa.in" \
    -cs md.tpr -ci gbsa_index.ndx -cg 17 18 \
    -ct md_whole.xtc -cp topol.top -cr complex.pdb \
    -o FINAL_RESULTS_MMGBSA.dat -eo FINAL_RESULTS_MMGBSA.csv \
    -do FINAL_DECOMP_MMGBSA.dat -deo FINAL_DECOMP_MMGBSA.csv -nogui
echo "=== 3HFM MM-GBSA done $(date) ==="
