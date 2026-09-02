#!/bin/bash
#SBATCH --job-name=3hfm_eq
#SBATCH --partition=gpu-v100
#SBATCH --gres=gpu:1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=24G
#SBATCH --time=0:40:00
# 3HFM holo equilibration: build + EM/NVT/NPT the full Fab+antigen complex, then build md.tpr.
# Submit from this directory (episcaf_v3/energetics/md/3hfm) on Gemini:  sbatch equil.sbatch
# The MD protocol (configs/*.mdp) is copied from bcell_epitope (Decision D3); see README.md.
set -e
module load Gromacs
export GMX_MAXBACKUP=-1
cd "$SLURM_SUBMIT_DIR"
mkdir -p out && cd out
echo "=== 3HFM HOLO EQUIL start $(date) ==="

# Keep the whole complex: antigen (chain Y) + antibody (H, L). Drop crystal water / heteroatoms.
awk '($0 ~ /^ATOM/ && (substr($0,22,1)=="Y" || substr($0,22,1)=="H" || substr($0,22,1)=="L")) || $0 ~ /^TER/ {print} END{print "END"}' \
    ../structures/3hfm.pdb > complex.pdb

# NOTE (see README, "known first hiccup"): the Fab has an H-L INTERCHAIN disulfide. If pdb2gmx does
# not form it across chains, add `-merge all` (merges chains into one moleculetype so inter-chain SS
# bonds form). Left off here to stay closest to the bcell template; check pdb2gmx.log for SS bonds.
gmx pdb2gmx -f complex.pdb -o complex.gro -p topol.top -i posre.itp \
    -ff amber99sb-ildn -water tip3p -ignh

gmx editconf -f complex.gro -o box.gro -c -d 1.2 -bt cubic
gmx solvate -cp box.gro -cs spc216.gro -o solv.gro -p topol.top
gmx grompp -f ../configs/ions.mdp -c solv.gro -p topol.top -o ions.tpr -maxwarn 2
echo SOL | gmx genion -s ions.tpr -o ions.gro -p topol.top -pname NA -nname CL -neutral -conc 0.15

gmx grompp -f ../configs/em.mdp  -c ions.gro -p topol.top -o em.tpr  -maxwarn 2
gmx mdrun  -deffnm em  -nb gpu -ntmpi 1 -ntomp 8
gmx grompp -f ../configs/nvt.mdp -c em.gro  -r em.gro  -p topol.top -o nvt.tpr -maxwarn 2
gmx mdrun  -deffnm nvt -nb gpu -ntmpi 1 -ntomp 8
gmx grompp -f ../configs/npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -o npt.tpr -maxwarn 2
gmx mdrun  -deffnm npt -nb gpu -ntmpi 1 -ntomp 8
gmx grompp -f ../configs/md.mdp  -c npt.gro -t npt.cpt -p topol.top -o md.tpr -maxwarn 2

touch EQUIL_OK
echo "=== 3HFM HOLO EQUIL done $(date) ==="
