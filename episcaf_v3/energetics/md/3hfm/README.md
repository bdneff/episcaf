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

## How to run it (Brandon submits; Claude does not drive SLURM)
**Run from a `/scratch` clone, not `/tgen_labs`.** The network FS gave cross-node "file not found"
on freshly written inputs (`../configs/*.mdp` invisible to the compute node while `../structures/`
was fine); `/scratch` is coherent. So: `git clone` the repo under `/scratch/$USER`, then:
1. `cd episcaf_v3/energetics/md/3hfm`.
2. `sbatch equil.sbatch` → copies the frozen inputs into `out/`, builds + equilibrates the solvated
   complex, writes `out/md.tpr` and `out/EQUIL_OK`. (`module load Gromacs`; `gpu-v100`.)
3. `sbatch --dependency=afterok:<equil_jobid> prod.sbatch` → 20 ns production, `out/md.xtc`.
4. Copy the keepers (`out/md.xtc md.tpr md.log topol.top`) to `$WS` when done — `/scratch` is purged.

## Disulfides — verify (built in, but confirm)
The Fab has an **H–L interchain disulfide** (L:CYS214 ↔ H:CYS215). `equil.sbatch` now passes
`-merge all` to `pdb2gmx` so the chains are one moleculetype and that bond forms. Confirm in
`out/pdb2gmx.log` that it links CYS-214 (L) to CYS-215 (H), on top of the intra-chain bonds
(lysozyme's 4: 6–127, 30–115, 64–80, 76–94; two each in H and L). The first run *without* `-merge`
formed every intra-chain bond but missed the interchain one — hence this fix.

## What comes next (not built yet)
From the production trajectory: Decision D2 (how we compute a defensible per-residue energy, including
solvent) → rank residues → check the ranking recovers 3HFM's SKEMPI hot spots. Only then does the
energetic definition earn the right to drive contig design.
