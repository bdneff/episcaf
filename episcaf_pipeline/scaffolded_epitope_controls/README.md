# Scaffolded-epitope controls (DP4 component C6)

Mutation controls that measure, for each scaffolded design, **the contribution of each epitope
island** (alanine-substitute one island at a time) and **the contribution of the scaffold** in
imposing the correct island shape/orientation (structure-disrupting scaffold mutation). These are
**not new designs** — they take an *existing* design's sequence and substitute residues in a
prespecified way (John Altin, Slack 2026-07-06: *"this would basically just be taking existing
designs and substituting alanines in a prespecified way — yeah exactly"*).

## The key idea: case-encoded sequence (`scaffoldEPITOPE`)
Each design is stored as a string where **UPPERCASE = epitope residues, lowercase = scaffold**.
A **contiguous run of uppercase letters = one epitope island**. That casing *is* the annotation —
no residue-index map or graft join is needed; island membership is read straight off the string.

## Reference implementation (`reference_dp3/`, John's original DP3 code — do not edit)
- `scaffoldedEpitopeControls.R` — driver: reads `scaffolded403.csv`, emits 5 flavors per design.
- `functions/mutateToA.R` — `mutateToA(seq, n)`: take the n-th uppercase run (island n), substitute
  it to Alanine. DP3 mutated **every other** residue in the run; returns `""` if island n doesn't
  exist (single-island designs → island-2 flavor is empty and filtered out).
- `functions/mutate6.R` — `mutate6(seq, "PPDDGG", n)`: replace `n` lowercase 6-mers (each ≥4 residues
  away from any uppercase/epitope position) with the structure-disrupting hexamer `PPDDGG`. `n=1` →
  `_scaffoldMutX1`, `n=4` → `_scaffoldMutX4`.
- `scaffolded403.csv` — the 403 assayed DP3 designs (input); `DP3_annot.csv` — the full DP3 library
  annotation (6,000 rows) incl. the `scaffoldEPITOPE` column and `numCIDRs` (island count).

DP3 flavors per design: `original`, `_scaffoldMutX1`, `_scaffoldMutX4`,
`_epitopeIsland1Mut>A`, `_epitopeIsland2Mut>A`.

## DP4 modifications (John, 2026-07-06)
1. **Alanine-substitute EVERY epitope residue of the focal island** (not every other) — cleaner.
2. **Keep only `_scaffoldMutX4`** (drop `_scaffoldMutX1`); X4 had the bigger/cleaner effect and saves space.

So DP4 flavors per design: `original` (= C1, already shipped), `_scaffoldMutX4`,
`_epitopeIsland1Mut>A`, `_epitopeIsland2Mut>A` (island 2 only for dual-island).

## What the DP4 port needs
A Python `build_c6_mutants.py` (TODO) that: takes our selected designs' **case-encoded
`scaffoldEPITOPE`** strings, applies all-residue island→A (ports `mutateToA` with mod #1) and
`_scaffoldMutX4` (ports `mutate6`), and writes the DP2-format rows. The one input to produce is the
case-encoded string per selected design (full sequence with epitope positions uppercased, from the
epitope mask) — the same per-design sequence step that library assembly needs anyway.
