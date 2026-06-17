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

# canonical scoring inputs (VERIFY exact filenames on the cluster)
METRICS_12MER    = DATA["twelvemer_runs"] / "06_score" / "metrics_12mer.csv"
METRICS_ANTIBODY = DATA["antibody_runs"] / "run_rfd3_mpnn" / "04_filter" / "metrics_cylinder_full.csv"
