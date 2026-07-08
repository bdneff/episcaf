# Scaffolded-epitope controls (DP4 component C6)

Mutation controls that measure, per scaffolded design, **the contribution of each epitope island**
(alanine-substitute one island at a time) and **the contribution of the scaffold** in holding the
islands in the right shape and orientation (a structure-disrupting scaffold mutation). These are
**not new designs** — they take an existing design's sequence and substitute residues in a
prespecified way, carried over from the DP3 controls.

## The key idea: case-encoded sequence (`scaffoldEPITOPE`)
Each design is a string where **UPPERCASE = epitope residue, lowercase = scaffold**, and a
**contiguous run of uppercase letters = one epitope island**. That casing *is* the annotation — no
residue-index map or graft join is needed; island membership is read straight off the string. See
`docs/CASE_ENCODING.md` for how we recover this casing for our designs.

## Reference implementation (`reference_dp3/` — the original DP3 code, do not edit)
- `scaffoldedEpitopeControls.R` — driver: reads `scaffolded403.csv`, emits 5 flavors per design.
- `functions/mutateToA.R` — `mutateToA(seq, n)`: take the n-th uppercase run (island n), substitute it
  to alanine. DP3 mutated **every other** residue; returns `""` if island n doesn't exist (single-island
  designs → the island-2 flavor is empty and filtered out).
- `functions/mutate6.R` — `mutate6(seq, "PPDDGG", n)`: replace `n` lowercase 6-mers (each ≥4 residues
  from any uppercase/epitope position) with the structure-disrupting hexamer `PPDDGG`.
- `scaffolded403.csv` — the 403 assayed DP3 designs (input); `DP3_annot.csv` — the full DP3 library
  annotation (6,000 rows) with the `scaffoldEPITOPE` column and `numCIDRs` (island count). `DP3_annot.csv`
  is also our reference instance of the 8-column annotated format used across DP4.

DP3 flavors per design: `original`, `_scaffoldMutX1`, `_scaffoldMutX4`, `_epitopeIsland1Mut>A`,
`_epitopeIsland2Mut>A`.

## DP4 modifications
1. **Alanine-substitute EVERY epitope residue of the focal island** (not every other) — cleaner readout.
2. **Keep only `_scaffoldMutX4`** (drop `_scaffoldMutX1`); X4 was the stronger perturbation and saves space.
3. **Graceful fallback 4 → 3 → 2 → 1** (`--scaffold-min 1`, the default): when a scaffold can't fit 4
   disjoint hexamers ≥4 from the epitope, degrade the count (`_scaffoldMutX3`/`X2`/`X1`) rather than
   dropping the control; a design loses its scaffold control only if it can't place even one hexamer.
   Set `--scaffold-min 3` for the "≥3" policy, or `4` for the old drop-if-not-4 behavior.

DP4 flavors per design: `_scaffoldMutX4`, `_epitopeIsland1Mut>A`, `_epitopeIsland2Mut>A` (island 2 only
for dual-island). The unmutated design is not emitted here — it is already shipped as C1, and C6 is the
controlled comparison against it.

## The port: `build_c6_mutants.py`
Ports the scheme to Python with the two DP4 modifications, **seeded** (so the library is exactly
regenerable) and **fail-soft** (a design that can't fit four disjoint scaffold hexamers is logged and
skipped rather than crashing the run). Input is the case-encoded `scaffoldEPITOPE` per selected design,
produced by `scripts/case_encode_selected.py` (C1/C5) / `scripts/case_encode_c2.py` (C2).

```bash
python episcaf_pipeline/scaffolded_epitope_controls/build_c6_mutants.py \
    --input results/dp4_C1_scaffoldEPITOPE.csv \
    --id-col token --target-col target --seq-col scaffoldEPITOPE \
    --drop-targets 2h32,4xwo,7a3t \
    --out results/dp4_C6_controls.csv
```

Current build (C1 top-20 over the 56 mAbs): **3,100 constructs** — island1→A 1,120, island2→A 860,
scaffold disruption 1,120 (`scaffoldMutX4` 1,034 + `scaffoldMutX3` 53 + `scaffoldMutX2` 33; **every design
gets a scaffold control** — 86 fell back below 4, none dropped). This build is off the **current 104-mer
C1 pool** and will be **regenerated off the native-103 C1 rerun** (see `docs/DP4_LIBRARY.md`, the
C1-at-103 redo) once that pool is scored; at that point no 104→103 trim is needed. Full component
context: `docs/DP4_LIBRARY.md`.
