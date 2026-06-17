"""Schema contract for the episcaf RFD3→AF3 pipeline.

This project treats a Parquet file as the *design ledger* (source-of-truth).
Downstream stages append/emit artifacts and (optionally) results.

Contract:
- Stages 02–04 may only depend on SAFE_COLS.
- Stage 05+ may read/produce RESULT_COLS (nullable).

If you add columns, update these lists and keep the above invariants true.
"""

SAFE_COLS = ['antigen_seq', 'heavy_seq', 'light_seq', 'epitope_seq', 'id', 'resolution', 'light_original_segid', 'heavy_original_segid', 'antigen_original_segid', 'r_work', 'r_free', 'antigen_ncbi_taxonomy_id', 'light_ncbi_taxonomy_id', 'heavy_ncbi_taxonomy_id', 'epitope_resindices', 'epitope_boolmask', 'contig_string', 'epitope_chunks', 'contig_gaps', 'contig_min', 'contig_id', 'contig_length', 'epitope_chunk_resindices', 'scaffolded_epitope_chunk_resindices', 'scaffolded_epitope_resindices', 'rfd_id', 'mpnn_id', 'scaffolded_epitope_seq', 'assay_scaffolded_epitope_seq', 'assay_scaffolded_epitope_id', 'assay_scaffolded_epitope_chunk_resindices', 'assay_scaffolded_epitope_resindices']

RESULT_COLS = ['overall_rmsd', 'epitope_chunk_rmsd_vs_mpnn', 'mean_pae', 'mpnn_clash_resindices', 'af3_clash_resindices']

# Minimal keys that uniquely identify a design through generation.
DESIGN_KEYS = ["id", "contig_id", "rfd_id", "mpnn_id"]

