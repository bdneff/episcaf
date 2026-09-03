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

## 3. Run (CPU / MPI, ~hours; config is ../../mmgbsa/mmpbsa.in)
```bash
mpirun -np 8 gmx_MMPBSA -O \
    -i ../../mmgbsa/mmpbsa.in \
    -cs md.tpr -ci gbsa_index.ndx -cg antibody antigen \
    -ct md.xtc -cp topol.top \
    -o FINAL_RESULTS_MMGBSA.dat  -eo FINAL_RESULTS_MMGBSA.csv \
    -do FINAL_DECOMP_MMGBSA.dat  -deo FINAL_DECOMP_MMGBSA.csv \
    -nogui
```
`-cg` takes the receptor and ligand group names (or numbers). `mmpbsa.in` sets igb=5, 0.15 M salt,
per-residue decomposition (`idecomp=1`), 100 frames (`interval=20`).

## 4. Downstream — same false-positive test
`FINAL_DECOMP_MMGBSA.csv` has per-residue TOTAL ΔG contributions. The antigen residues are numbered
430–558; subtract 429 to get lysozyme numbering (= `ag_res_idx`). Reshape those into a small csv with
columns `ag_res_idx,resname,mmgbsa_total` and feed it to the existing test:
```bash
python ../plot_ie_vs_ddg.py --ie <that csv> --ddg ../skempi_3hfm_ddg.csv \
    --channel mmgbsa_total --out ie_vs_ddg_mmgbsa.png
```
(Claude will write the small parser from the gmx_MMPBSA output once we see its exact format.) The
test is the same: **does the desolvation-inclusive ΔG empty the false-positive zone** — demote Arg73
and Leu75/Arg21 while keeping K96/K97/Y20/D101?

## First-run gotchas to watch
- **Topology split.** `-merge all` gives one moleculetype; gmx_MMPBSA rebuilds sub-topologies from the
  index groups, which is standard, but if it complains about the ligand not being a separate molecule,
  that's the thing to sort (paste me the error).
- **Frames.** `interval=20` → ~100 frames for a first pass; lower it for more sampling once it works.
- **Runtime/mem.** Per-residue GB over ~100 frames on this size system is CPU-heavy — give it several
  hours and enough MPI ranks.
