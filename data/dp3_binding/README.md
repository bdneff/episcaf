# DP3 assay binding data

Ground-truth PepSeq binding intensities for the DP3 (known-mAb) designs, from John Altin
(2026-06). These are the experimental binding labels the scorer has been waiting for: they
let us ask, on real data, which AlphaFold3 design metrics actually track antibody binding.

## Files (raw, as received)
- `20260114_IM226_DP2_mAbs.csv` — assay run **IM226**, 6 antibodies: `6o9i 8cz8 8jnk 6xxv
  5fhx 7ox3`.
- `20260202_IM229_DP2.csv` — assay run **IM229**, 2 antibodies: `8db4 8pww`.

Both files contain the **same 1000 DP2 library members** with identical design columns;
they differ only in their `NoAb` baseline column and the per-antibody intensity columns.
8 antibodies total. Two were dropped from the original 10: **4xwo** (synthesized antibody
too low yield) and **7a3t** (epitope too small, RMSD ineffective) — their designs are still
in the library (`Target` 4xwo_5P, 7a3t_0P) but have no usable Ab column.

## Columns
- `library_member`, `sequence` (full 103-mer construct), `category`, `designedSequence`
  (the bare scaffold/peptide), `designedSequenceLength`, `Model`, `Design_ID`, `Target`, `Kd`.
- `NoAb_…` — the no-antibody baseline (scatterplot x-axis). One per run.
- `X<pdb>_…` — intensity for antibody `<pdb>`. John's readout is
  `log10(1+Ab)` vs `log10(1+NoAb)`, designs targeted to that Ab highlighted.

`category=scaffoldedAbEpitope` (403 rows) are our designs; the rest (`pMHCbinder`,
`published_binder`, `published_pMHCbinder`) are other library content / controls.

## Join to design metrics
`designedSequence` matches `scaffolded_epitope_seq` in `dp2.parquet` **1-to-1** (403/403
verified). `scripts/build_dp3_binding_join.py` merges the two runs, attaches the AF3 metrics
(`overall_rmsd`, `epitope_chunk_rmsd_vs_mpnn`, `mean_pae`, `af3_n_clash_res`, `is_pass`), and
resolves each design's **cognate** antibody (Target `7ox3_0P` -> Ab `7ox3`), writing
`results/dp3_binding_metrics.csv`.

## Caveat for fitting
Every assayed scaffold design already **passed the 4-filter** (`is_pass=True` for all 377
with a cognate Ab), so the metrics span only a narrow passing range (e.g. overall_rmsd
0.35–1.95, mean_pae 3.0–5.0, clash all 0) and the pass/fail filter has no variance here.
Pooled correlations are confounded by between-antibody offsets; within-antibody the signal
is weak. DP3 gives a starting prior, not a fit — the DP4 metric-space sampling is what
provides the variance to learn what predicts binding.
