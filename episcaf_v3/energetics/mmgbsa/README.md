# Per-residue MM-GBSA on the 3HFM interface (Decision D2 — the desolvation step)

Raw interaction energy over-called a charged residue (Arg73) because it ignores desolvation;
LJ-only fixed that but over-called well-packed residues. MM-GBSA is the principled fix: it keeps
Coulomb + LJ but adds the GB polar-desolvation term, so a charged residue's raw contact is netted
against the cost of shedding its water. Single-trajectory, so it runs on the **existing 20 ns
trajectory** (`md/3hfm/out/`) — no new simulation, no apo run needed.

**Untested first draft** (I can't run gmx_MMPBSA locally). Expect to iterate on two things: the env
install and the index groups. Brandon runs it; Claude stages.

## 1. Environment (one-time install)
`gmx_MMPBSA` wraps AmberTools' MMPBSA.py. Install it in its own conda env:
```bash
mamba create -n gmxMMPBSA -c conda-forge gmx_MMPBSA    # or: conda create ...
conda activate gmxMMPBSA
gmx_MMPBSA --help                                       # sanity
```
It needs a `gmx` on PATH too (its own or the module). See how bcell handled Amber tooling if this
env is already around.

## 2. Index groups (the one manual bit)
Define the two binding partners. After `-merge all` the residues are numbered continuously, antigen
last (this matches the `resid` column in `holo_ie_mean.csv`, which ran 430–558 for the antigen):
- **antibody** = residues 1–429 (chains L 1–214, H 215–429)
- **antigen (lysozyme)** = residues 430–558
```bash
cd md/3hfm/out
gmx make_ndx -f md.tpr -o gbsa_index.ndx
# at the prompt:
#   ri 1-429        then   name <N> antibody     (use the number it just created)
#   ri 430-558      then   name <M> antigen
#   q
```

## 3. Run — use the committed sbatch (the working recipe)
Everything below (env activate, PBC-whole fix, groups, flags) is baked into
`energetics/mmgbsa/run_mmgbsa.sbatch`. Submit it from the holo run dir on Gemini:
```bash
cd /scratch/bneff/episcaf_run/episcaf_v3/energetics/md/3hfm/out
sbatch /scratch/bneff/episcaf_run/episcaf_v3/energetics/mmgbsa/run_mmgbsa.sbatch
```
It runs **serial** (MPI ranks ran independently under SLURM here — mpi4py communicator issue), does a
`gmx trjconv -pbc whole` first (gmx_MMPBSA does not correct PBC; a split complex overflows sander),
and passes `-cs md.tpr -ci gbsa_index.ndx -cg 17 18 -ct md_whole.xtc -cp topol.top -cr complex.pdb`.
`mmpbsa.in` sets igb=5, 0.15 M salt, per-residue decomposition (`idecomp=1`); `interval` is the
sampling dial — currently `20` → **100 frames** over the 20 ns trajectory (the first pass used
`interval=100` → 20 frames, which was too noisy; see `docs/DECISIONS.md` D2).

## 4. Downstream — extract the per-residue ΔG and run the same false-positive test
gmx_MMPBSA's own decomp CSV writer **crashes** on our `-merge all` topology (H/L chains share residue
numbers), so `FINAL_DECOMP_MMGBSA.csv` never gets written — but the per-residue energies it already
computed sit in the sander `TDC` lines of the decomp mdouts it leaves behind (the crash happens at the
*writing* stage, so the intermediates are not cleaned up). Grep those two files (the same step that
made the committed `complex_tdc.txt` / `ligand_tdc.txt`), then run the parser and the plot:
```bash
# from the run dir, after the job -- confirm the mdout names first:  ls _GMXMMPBSA_*mdout
grep '^TDC' _GMXMMPBSA_complex_gb.mdout > ../complex_tdc.txt      # complex decomp mdout
grep '^TDC' _GMXMMPBSA_ligand_gb.mdout  > ../ligand_tdc.txt       # ligand  decomp mdout

python ../../mmgbsa_decomp_to_csv.py \
    --complex ../complex_tdc.txt --ligand ../ligand_tdc.txt \
    --resnames ../holo_ie_mean.csv --out ../mmgbsa_perres.csv

python ../../plot_ie_vs_ddg.py --ie ../mmgbsa_perres.csv --ddg ../../skempi_3hfm_ddg.csv \
    --channel mmgbsa_dg --out ../../../manuscript/figures/ie_vs_ddg_3hfm_mmgbsa.png
```
The parser cross-checks that the antigen residues match between the two files, so a wrong mdout is
caught. The test is the same: **does the desolvation-inclusive ΔG empty the false-positive zone** —
demote Arg73 and Leu75/Arg21 while keeping K96/K97/Y20/D101 — and, at 100 frames, does Lys96 recover?

## Gotchas (status)
- **Topology split — resolved.** `-merge all` gives one moleculetype; the decomp writer chokes on the
  duplicate H/L residue numbers, hence the direct `TDC` parse above rather than gmx_MMPBSA's CSV.
- **PBC — resolved.** gmx_MMPBSA does not correct PBC; a split complex overflows sander (`****`). The
  sbatch runs `gmx trjconv -pbc whole` first and feeds `md_whole.xtc`.
- **Serial — required here.** MPI ranks ran independently under SLURM (mpi4py communicator issue), so
  the sbatch is `ntasks=1`. Per-residue GB over 100 frames on this system is the reason for the 12 h
  walltime; if it times out, raise `interval` (fewer frames) in `mmpbsa.in`.
- **Frames.** `interval` is the sampling dial: `100`→20 frames (first pass, too noisy), `20`→100
  frames (current), `10`→200 frames.
