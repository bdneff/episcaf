# 3HFM — the D1 validation pilot (holo MD)

3HFM is the HyHEL-10 Fab bound to hen egg white lysozyme. It is the v3 validation pilot (Decision
D1): a textbook hot-spot system with experimental alanine-scan ΔΔG in SKEMPI (18 measurements, up to
+6.5 kcal/mol), and lysozyme is the system `bcell_epitope` already runs MD on. We run the holo
complex, then (later) decompose per-residue energies and check that our ranking recovers the known
hot spots before the definition drives any design.

## The structure
`structures/3hfm.pdb` — the raw RCSB deposition (`https://files.rcsb.org/download/3HFM.pdb`).
Chains: **Y = lysozyme (antigen)**, **H / L = HyHEL-10 Fab (antibody)**. No missing residues
(REMARK 465 empty); the only heteroatom is one crystal water, dropped at build time. So no PDBFixer
repair is needed for this structure (unlike some AbDb complexes in bcell).

## The MD protocol (Decision D3 — adopted from bcell as the starting template)
`configs/*.mdp` are copied verbatim from `bcell_epitope`'s frozen stack
(`campaign1_apo_features/md/1AKI/apo/configs/`), retitled for this complex, with production set to
**20 ns** (bcell's holo length). In one line: AMBER99SB-ILDN + TIP3P, cubic box 1.2 nm clearance,
0.15 M NaCl, PME (`rcoulomb=rvdw=1.0`), V-rescale thermostat at 300 K, C-rescale barostat for NPT
equilibration and Parrinello-Rahman for production, 2 fs step with LINCS on h-bonds. `em` (steep) →
`nvt` (100 ps, position-restrained) → `npt` (100 ps, restrained) → `md` (20 ns). This is a
**decision, logged in `docs/DECISIONS.md` (D3)**: it is bcell's protocol taken as a starting point,
to revisit if the energetics call for it — not an unexamined default.

## How to run it (on Gemini — Brandon submits; Claude does not drive SLURM)
1. `git pull` on Gemini, `cd episcaf_v3/energetics/md/3hfm`.
2. `sbatch equil.sbatch` → builds + equilibrates the solvated complex, writes `out/md.tpr` and
   `out/EQUIL_OK`. (`module load Gromacs` = GROMACS 2023.2-dev; `gpu-v100`.)
3. `sbatch --dependency=afterok:<equil_jobid> prod.sbatch` → 20 ns production, `out/md.xtc`.
4. Move `out/` off `/scratch` promptly — scratch is purged.

## Known first hiccup to check
The Fab has an **H–L interchain disulfide**. If `pdb2gmx` does not form disulfides across chains,
add `-merge all` to the `pdb2gmx` call in `equil.sbatch` (merges chains into one moleculetype so
inter-chain SS bonds form). Check `out/pdb2gmx.log` for the expected disulfides (lysozyme has 4;
the Fab has intra- and inter-chain SS bonds). This is the one place the holo build departs from
bcell's single-chain apo template, so it is the most likely thing to need a second pass.

## What comes next (not built yet)
From the production trajectory: Decision D2 (how we compute a defensible per-residue energy, including
solvent) → rank residues → check the ranking recovers 3HFM's SKEMPI hot spots. Only then does the
energetic definition earn the right to drive contig design.
