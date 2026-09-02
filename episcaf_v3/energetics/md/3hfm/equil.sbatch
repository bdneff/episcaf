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
# Submit from this directory (episcaf_v3/energetics/md/3hfm):  sbatch equil.sbatch
# Best run from a /scratch clone: /tgen_labs is a network FS and gave cross-node "file not found"
# on freshly written files; /scratch is coherent. The MD protocol (configs/*.mdp) is copied from
# bcell_epitope (Decision D3); see README.md.
set -e
module load Gromacs
export GMX_MAXBACKUP=-1
cd "$SLURM_SUBMIT_DIR"
mkdir -p out

# Copy the frozen inputs INTO the run dir first, then work entirely with local paths -- no
# cross-directory lookups mid-run (this is what broke on /tgen_labs), and a fast, clear failure
# here if an input is missing rather than deep inside solvate.
cp configs/*.mdp out/
# Keep the whole complex: antigen (chain Y) + antibody (H, L). Drop crystal water / heteroatoms.
awk '($0 ~ /^ATOM/ && (substr($0,22,1)=="Y" || substr($0,22,1)=="H" || substr($0,22,1)=="L")) || $0 ~ /^TER/ {print} END{print "END"}' \
    structures/3hfm.pdb > out/complex.pdb
cd out

echo "=== 3HFM HOLO EQUIL start $(date) ==="

# -merge all: merge the chains into one moleculetype so pdb2gmx forms the Fab H-L INTERCHAIN
# disulfide (L:CYS214 <-> H:CYS215). Without it, chains are processed separately and that bond is
# missed (confirmed on the first run). Check pdb2gmx.log links CYS-214 (L) to CYS-215 (H).
gmx pdb2gmx -f complex.pdb -o complex.gro -p topol.top -i posre.itp \
    -ff amber99sb-ildn -water tip3p -ignh -merge all

gmx editconf -f complex.gro -o box.gro -c -d 1.2 -bt cubic
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
echo "=== 3HFM HOLO EQUIL done $(date) ==="
