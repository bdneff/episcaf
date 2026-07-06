# Recovering epitope positions for our reran designs (token → dp2 → case encoding)

This documents how to get, for any of our DP4-selected designs, **which residues are the epitope** —
needed to case-encode sequences for the C6 controls (`episcaf_pipeline/scaffolded_epitope_controls/`)
and to attach sequences at library assembly. Written because reconstructing this cost a long detour;
the recipe below is the authoritative source of truth.

## The problem
Our selections (C1 whole-epitope `dp4_C1_*`, C5 `dp4_C5_titration`) identify designs by a hash
`token` (col `token` in `metrics_native_cyl_full.csv`). The epitope-position annotation is **not** in
those tables — it has to be recovered.

## The answer (from `compute_metrics.py::run_metrics`, the driver that made our metrics)
We **reran RFdiffusion3 using Lawson's contigs** (the same DP2 library members; manuscript
`sec:methods`). An island sits at residues `N+1 … N+ℓ` of an `N`-flank contig, so each design's
epitope span is a property of its **contig**, not its sequence. The metrics `run` mode therefore:

1. joins each design to `dp2` by **`token == dp2.assay_scaffolded_epitope_id`**
   (`compute_metrics.py` lines ~700, 708), and
2. reads the epitope positions as **`dp2.scaffolded_epitope_chunk_resindices`** — 0-based indices into
   the design's chain A (line ~712: `rfd3_epi_ris = parse_index_list(r["scaffolded_epitope_chunk_resindices"])`).

Because the contigs are identical to Lawson's, `scaffolded_epitope_chunk_resindices` (a contig
property) is valid for **our** reran designs **even though the MPNN scaffold sequences differ** from
Lawson's. The chunk spans are contiguous → each island is one contiguous block of positions.

## Which `dp2` to use (this is the trap)
- Join against the dp2 the run used, which carries **our** tokens: `datasets/dp2.parquet` in the run's
  repo (`/home/bneff/rfd3/repo_refactored/datasets/dp2.parquet`, and/or `$WS/datasets/dp2.parquet`).
- **Do NOT** use the local `known_antigen/analysis/full_run/dp2.parquet` for the token join — it's a
  different/Lawson-token version; only ~149 of our tokens match it. (Verify token coverage before using
  any dp2: `set(sel.token) & set(dp2.assay_scaffolded_epitope_id)` should cover ~all selected designs.)

## Why the sequence-match attempt failed (don't repeat it)
`dp2.scaffolded_epitope_seq` is **Lawson's** design sequence, not ours (different MPNN run). So matching
our PDB-read sequence to dp2 by sequence fails for all but coincidences. Positions match by **contig**;
sequences do **not**. Join by token; take positions from dp2; read the sequence from our own PDB.

## Recipe (`scripts/case_encode_selected.py`)
Per selected design: `token → dp2 row` (by `assay_scaffolded_epitope_id`) → `scaffolded_epitope_chunk_resindices`;
read the design chain-A sequence from its `mpnn_pdb`; **uppercase those positions, lowercase the rest**
→ the case-encoded `scaffoldEPITOPE` (UPPERCASE = epitope, lowercase = scaffold; each contiguous
uppercase run = one island). PDB paths point to the stale `repo_refactored` home dir → remap to
`$WS/runs/run_rfd3_mpnn` (baked into `case_encode_selected.sbatch`).

## Provenance chain (multiple dp2 / run structures — do not confuse)
- `dp2` = Lawson's DP2 designs (RFD1). We reran RFD3 → different scaffolds, same contigs.
- The run that produced `metrics_native_cyl_full.csv`: `compute_metrics.py run --run_dir <run> --dp2_parquet datasets/dp2.parquet`
  (see `$WS/sbatch_metrics_my_run.sh`). It joins run designs to dp2 by token and pulls chunk positions.
- Local `known_antigen/analysis/data/metrics_native_cyl_full.csv` is our C1/C5 selection source (token-keyed).
