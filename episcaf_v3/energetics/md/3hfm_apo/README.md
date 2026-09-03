# Apo lysozyme run (3HFM antigen alone) — for the entropy / RMSF comparison

The per-residue entropic signal for Category-2 residues comes from how a residue's fluctuations
*change on binding* — so it needs the unbound (apo) state to compare against the holo run. This is
that apo half: lysozyme (chain Y) alone, antibody stripped. It is only meaningful **paired with the
holo run at an identical protocol**, which is why `configs/*.mdp` are copied verbatim from `3hfm/`.

## What it is
- Antigen only: chain Y (lysozyme), antibody chains H/L stripped by the `awk` in `equil.sbatch`.
- Same force field / water / box / ions / equilibration / **20 ns** production as the holo run.
- Single chain, so **no `-merge`**; lysozyme's 4 intra-chain disulfides (6–127, 30–115, 64–80,
  76–94) form automatically — confirm them in the equil slurm stdout.

## Run (Gemini; Brandon submits, from a /scratch clone)
```bash
cd /scratch/bneff/episcaf_run && git pull
cd episcaf_v3/energetics/md/3hfm_apo
sbatch equil.sbatch
sbatch --dependency=afterok:<equil_jobid> prod.sbatch
```

## Residue mapping (for the bound-vs-apo comparison)
- **apo** lysozyme = residues **1–129** (single chain).
- **holo** lysozyme = residues **430–558** (it was last in the `-merge all` complex; see
  `holo_ie_mean.csv`'s `resid` column).
- Both correspond to lysozyme 1–129; the RMSF/entropy analysis subtracts the 429 offset to align them.

## Length — 20 ns now, extend later if needed
20 ns matches the existing holo run, so we get a matched pair for **RMSF** immediately (backbone RMSF
is reasonably converged at this length). **Configurational entropy** (quasi-harmonic) needs longer —
~100 ns and ideally replicates. If the 20 ns RMSF looks worth pursuing, extend **both** bound and apo
to ~100 ns together (`gmx convert-tpr -s md.tpr -extend <ps> -o md.tpr; gmx mdrun -deffnm md -cpi
md.cpt`) so they stay matched. Decide after seeing the 20 ns signal.

## Downstream
Per-residue RMSF (bound vs apo) and ΔRMSF; optionally per-residue quasi-harmonic entropy. The
formalism and the open questions (including the sign-ambiguity caveat and dynamic coupling) are in the
manuscript section "Open questions: entropy and dynamics" (`sections/entropy_questions.tex`).
