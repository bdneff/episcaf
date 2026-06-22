"""Absolute data locations on the cluster — the ONLY place paths live.

The repo holds no data. Edit WORKSPACE if the lab workspace moves, or point these
at a local copy when developing off-cluster. Confirm the exact metrics filenames
against the cluster (06_score / 04_filter contents) before relying on them.
"""
from pathlib import Path

WORKSPACE = Path("/tgen_labs/altin/alphafold3/workspace/episcaf_v2_bneff")

DATA = {
    "antibody_runs":  WORKSPACE / "runs",                   # DP3 / mAb set   (~312 G)
    "twelvemer_runs": WORKSPACE / "run_12mer_scaffolding",  # tiled 12mer set (~114 G)
    "datasets":       WORKSPACE / "datasets",               # dp2.parquet etc.
}

# Lawson's DP3 design ledger (verified on cluster 2026-06-22). The DP3 set is the 59-epitope
# subset of the 1,134 AbDb complexes in ABDB_CLEANED_PDB_DIR (<=2 islands, antigen >103 aa).
DP2_PARQUET = DATA["datasets"] / "dp2.parquet"

# canonical scoring inputs (VERIFY exact filenames on the cluster)
METRICS_12MER    = DATA["twelvemer_runs"] / "06_score" / "metrics_12mer.csv"
METRICS_ANTIBODY = DATA["antibody_runs"] / "run_rfd3_mpnn" / "04_filter" / "metrics_cylinder_full.csv"

# AbDb cleaned antibody-antigen complex PDBs, one per ledger id, named "<id>.pdb"
# (e.g. 4xwo_5P.pdb). These are the per-epitope antigen coordinate source for the design
# run: pass this as stage03 --pdb_dir so each row resolves <id>.pdb. The antigen is chain A
# (the contig A-segments and epitope_resindices index chain A); RFD3 fixes only the chain-A
# epitope atoms and ignores the unreferenced antibody chains. Lives under a different
# workspace root than WORKSPACE above. [verify on cluster: chain A == antigen]
ABDB_CLEANED_PDB_DIR = Path(
    "/tgen_labs/altin/alphafold3/workspace/episcaf-experiments/data/abdb/complex_pdbfiles/cleaned"
)
