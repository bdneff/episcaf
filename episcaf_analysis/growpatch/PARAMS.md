# growPatch parameters and provenance

The fitted dials of the no-antibody patch enumerator, with where each number comes from, so they
can be changed deliberately. All are `core.py` module constants and every public function takes
them as keyword arguments — nothing is hard-coded past the defaults.

These, plus the eyeview `footprint_radius` and the exhaustiveness knob `nms_jaccard` / `top_k`, are
the user-facing dials of the combined production entry `read.enumerate_candidates`. The
seed-and-grow-vs-eyeview benchmark and the decision to ship both (unioned) are in
`results/growpatch/ENUMERATOR_DECISION.md`.

Source: John Altin's `segmentation/growPatch.py` (AltinLab/epitope-universe @ `upstream/main`,
commit `668de1d`, 2026-11-11). Our `core.py` is a byte-parity port of the algorithm (verified:
identical candidate set to his code on 8pww chain A), with the structure reading and the freesasa
surface step lifted into `read.py`.

## The dials

| constant | default | role |
|---|---|---|
| `SASA_CUT` | **0.2** | a residue is "surface" if its relative SASA ≥ this (applied in `read.py`) |
| `PATCH_SIZE` | **16** | raw residues a seed grows to before gap-fill |
| `GAP` | 10 | sequence gaps ≤ this are bridged in gap-fill |
| `MIN_RES` | 12 | min residues in a kept, gap-filled candidate |
| `MAX_ISLANDS` | 3 | max sequence islands (summed over chains) in a candidate |
| `MIN_ISLAND_LEN` | 3 | gap-filled islands shorter than this are dropped unless protected |
| `CONTACT` | 4.5 Å | heavy-atom distance for the residue adjacency graph |
| `GROW_LINK` | 5.0 Å | heavy-atom distance a new residue must be within to join a patch |
| `NMS_JACCARD` | 0.8 | two candidates with Jaccard ≥ this are near-duplicates (one dropped) |

## Provenance of the two fitted dials (SASA_CUT, PATCH_SIZE)

Fit on a stratified random sample of 30 reference epitopes from
`reference/sabdb_epitopes_4A.csv.gz` (single-chain, gap-10 island count in {1,2,3}, 10 per island
count, seed 42), scored against the design goal directly — not a proxy:

1. maximize the mean fraction of each epitope's (un-gap-filled) residues covered by its single best
   candidate,
2. minimize that candidate's size,
3. minimize the total number of candidates the tool outputs (post-NMS, across the chain).

Swept `SASA_CUT ∈ [0.05..0.25] × PATCH_SIZE ∈ [10..28]`:

- **PATCH_SIZE is the dominant lever.** Mean best-single coverage rises 65% (size 10) → 85% (12) →
  90% (14–16), then flattens (92–96% through 28) for +15–20 residues/design and no consistent gain.
- **SASA_CUT barely moves coverage/size** (86–90% / 23–25 residues across the whole range) but
  strongly affects candidate count: 0.05→0.2 cuts mean candidates 92.2 → 54.4 (−41%) while coverage
  drops only 89.5% → 88.8%. `SASA_CUT=0.2` is the point past which further count reduction starts
  costing real coverage (0.25 → 36.5 candidates but coverage falls to 86.3%).

**PATCH_SIZE revision 14 → 16:** re-tested on a disjoint sample of 71 single-chain epitopes drawn
from the SASA-filtered table; mean best-single coverage rises smoothly with size (82.7% @14, 86.9%
@16, 91.1% @20) with no knee — 16 beat 14 in 32/71 epitopes and lost in none, for +2.9 residues in
the best design. 14 was a point on a smooth curve, not a plateau.

## Code-vs-report reconciliation (resolved 2026-11-11)

Earlier our notes flagged a split — the report used `SASA_CUT=0.2` while the committed reference
code defaulted `0.10`. Resolved: `upstream/main` was updated 2026-11-11 and now defaults
**`SASA_CUT=0.2, PATCH_SIZE=16`** (the `sasa-epitope-filter` work merged). Report and code agree;
no divergence remains. We adopt `0.2 / 16`.

## relSASA engine and the surface cut

relSASA is computed in `read.py` with **freesasa**, normalizing each residue's absolute SASA by the
**Tien et al. 2013 theoretical max-ASA** table (not freesasa's own built-in relative field, which uses
a different reference). This matches John's `growPatch.py` exactly — verified: `read.enumerate_patches`
reproduces his candidate set on 8pww:A (108 candidates over the 159-residue surface set at cut 0.2).

- **`SASA_CUT=0.2` is the fitted value**, from the coverage/count sweep (above). The b-cell project
  uses **0.25** as a generic "exposed" cutoff; that is the textbook convention, not a fit, so we keep
  the fitted 0.2. Worth re-confirming on our own structures (the CLI's sweep does this).
- **relSASA-engine alternative (dependency reduction, not yet adopted):** the same Tien-2013 relSASA
  is computable with **Biopython `Bio.PDB.SASA.ShrakeRupley`** (as in
  `bcell_epitope/.../recompute_apo_baselines.py`), which would drop the freesasa dependency entirely
  (biopython is in every env; freesasa is an extra). The two engines' absolute areas are **not**
  byte-identical (freesasa Lee-Richards vs Biopython Shrake-Rupley, probe 1.4 Å, n_points 100), so the
  surface set at a fixed cut could shift and the whole enumeration would need re-validation. Deferred
  as a measured swap: compare the two surface sets on a panel; if they agree ≥99% at cut 0.2, drop
  freesasa. Until then, freesasa is the parity engine.

## To re-sweep

John's driver sweeps with `--sasa-cuts` / `--patch-sizes` and writes a `growpatch_grid_<tag>.csv`.
Our enumerator (added in `read.py` / the CLI, next step) will expose the same sweep so the fit is
reproducible on our structures before we treat `0.2 / 16` as final for the episcaf set.
