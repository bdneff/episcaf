#!/bin/bash
#SBATCH --job-name=3hfm_apo_eq
#SBATCH --partition=gpu-v100
#SBATCH --gres=gpu:1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=24G
#SBATCH --time=0:40:00
# 3HFM APO lysozyme equilibration -- the antigen (chain Y) alone, antibody stripped.
# SAME protocol as the holo run (configs/*.mdp are copied from 3hfm/): the apo/holo pair must be
# identical except for the antibody, or the bound-vs-apo RMSF/entropy comparison is confounded.
# Submit from this dir; best from a /scratch clone (see 3hfm/README on the /tgen_labs NFS issue).
set -e
module load Gromacs
export GMX_MAXBACKUP=-1
cd "$SLURM_SUBMIT_DIR"
mkdir -p out
cp configs/*.mdp out/
# keep ONLY the antigen (chain Y = lysozyme); drop antibody chains H,L and water/heteroatoms.
awk '($0 ~ /^ATOM/ && substr($0,22,1)=="Y") || $0 ~ /^TER/ {print} END{print "END"}' \
    structures/3hfm.pdb > out/antigen.pdb
cd out

echo "=== 3HFM APO EQUIL start $(date) ==="
# single chain -> no -merge. Lysozyme's 4 intra-chain disulfides form automatically; confirm in the
# slurm stdout (Linking CYS-6/127, 30/115, 64/80, 76/94).
gmx pdb2gmx -f antigen.pdb -o apo.gro -p topol.top -i posre.itp \
    -ff amber99sb-ildn -water tip3p -ignh

gmx editconf -f apo.gro -o box.gro -c -d 1.2 -bt cubic
gmx solvate -cp box.gro -cs spc216.gro -o solv.gro -p topol.top
gmx grompp -f ions.mdp -c solv.gro -p topol.top -o ions.tpr -maxwarn 2
echo SOL | gmx genion -s ions.tpr -o ions.gro -p topol.top -pname NA -nname CL -neutral -conc 0.15

gmx grompp -f em.mdp  -c ions.gro -p topol.top -o em.tpr  -maxwarn 2
gmx mdrun  -deffnm em  -nb gpu -ntmpi 1 -ntomp 8
gmx grompp -f nvt.mdp -c em.gro  -r em.gro  -p topol.top -o nvt.tpr -maxwarn 2
gmx mdrun  -deffnm nvt -nb gpu -ntmpi 1 -ntomp 8
gmx grompp -f npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -o npt.tpr -maxwarn 2
gmx mdrun  -deffnm npt -nb gpu -ntmpi 1 -ntomp 8
gmx grompp -f md.mdp  -c npt.gro -t npt.cpt -p topol.top -o md.tpr -maxwarn 2

touch EQUIL_OK
echo "=== 3HFM APO EQUIL done $(date) ==="
