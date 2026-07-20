"""Absolute data locations on the cluster — the ONLY place paths live.

The repo holds no data. Edit WORKSPACE if the lab workspace moves. Cluster metrics
filenames are confirmed below; local-dev copies (used by the figure scripts) live in
sibling dirs outside the repo and are centralized in LOCAL_METRICS at the bottom.
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

# Off-cluster (local-dev) copy of the same dp2.parquet, recorded 2026-06-25. 24 MB; lives
# outside the v2 repo under the sibling known_antigen analysis dir (data stays out of git;
# see .gitignore). Verified == the cluster ledger: 59 unique patch ids == data/sequences/
# dp3_mab_antigens.fasta exactly (Lawson scaffolded all 59; "57" is the unique-PROTEIN count).
# Use locally for analysis/plots; on the cluster use DP2_PARQUET above.
DP2_PARQUET_LOCAL = Path(
    "/Users/bneff/Desktop/projects/episcaf/known_antigen/analysis/full_run/dp2.parquet"
)

# canonical scoring inputs (filenames confirmed on cluster 2026-07: metrics_12mer.csv,
# metrics_cylinder_full.csv). METRICS_ANTIBODY points at the original 104-mer whole-epitope run;
# the shipped C1 is the native-103 redo (METRICS_C1 below).
METRICS_12MER    = DATA["twelvemer_runs"] / "06_score" / "metrics_12mer.csv"
METRICS_ANTIBODY = DATA["antibody_runs"] / "run_rfd3_mpnn" / "04_filter" / "metrics_cylinder_full.csv"

# The three metrics tables the SHIPPED DP4 library was selected from, in their durable homes. C1 and
# C2 were built on /scratch and copied here 2026-07-17 (/scratch is swept on ~30 days); md5-verified
# byte-identical to the local copies the selection ran on. NOTE these tables' own `af3_dir` column
# still records the old /scratch paths -- scripts that follow it take an OLD:NEW --af3-remap.
METRICS_C1 = DATA["antibody_runs"] / "whole_epitope_rfd3" / "05_analysis" / "metrics_whole_epitope.csv"
METRICS_C2 = DATA["antibody_runs"] / "dual_island_rfd3" / "05_analysis" / "metrics_dual_island.parquet"
METRICS_C3 = METRICS_12MER

# Where the runs those metrics came from now live, and the prefix swap that redirects a stale
# af3_dir/mpnn_pdb path onto them. Keep in sync with scripts/build_superset.sbatch REMAP.
SCRATCH_RUNS_OLD = "/scratch/bneff/episcaf/runs"
AF3_REMAP = f"{SCRATCH_RUNS_OLD}:{DATA['antibody_runs']}"

# ---------------------------------------------------------------------------------------------
# Local-dev copies (OFF-cluster). The manuscript figure scripts read metrics from sibling analysis
# dirs that live OUTSIDE the repo (data stays out of git; see data/README.md). Centralized here so
# the figure scripts stop hardcoding laptop paths. On the cluster these same tables live under
# WORKSPACE; a newcomer without the sibling dirs cannot regenerate the DP3/12-mer figures (noted in
# manuscript/figures/FIGURES.md).
LOCAL_ROOT = Path("/Users/bneff/Desktop/projects/episcaf")   # sibling of the v2 repo
LOCAL_METRICS = {
    "rfd1_mpnn_lawson":     LOCAL_ROOT / "known_antigen/analysis/full_run/metrics_full_rfd1_mpnn_LAWSON.csv",
    "rfd3_cylinder_full":   LOCAL_ROOT / "known_antigen/analysis/data/metrics_cylinder_full.csv",
    "native_cyl_full":      LOCAL_ROOT / "known_antigen/analysis/data/metrics_native_cyl_full.csv",
    "whole_epitope_103":    LOCAL_ROOT / "known_antigen/analysis/data/metrics_whole_epitope_103.csv",
    "metrics_12mer":        LOCAL_ROOT / "12mer_tiling/analysis/data/metrics_12mer.csv",
    "composite_12mer_top5": LOCAL_ROOT / "12mer_tiling/analysis/data/composite_12mer_top5_allscored.csv",
}

# AbDb cleaned antibody-antigen complex PDBs, one per ledger id, named "<id>.pdb"
# (e.g. 4xwo_5P.pdb). These are the per-epitope antigen coordinate source for the design
# run: pass this as stage03 --pdb_dir so each row resolves <id>.pdb. The antigen is chain A
# (the contig A-segments and epitope_resindices index chain A); RFD3 fixes only the chain-A
# epitope atoms and ignores the unreferenced antibody chains. Lives under a different
# workspace root than WORKSPACE above. [verify on cluster: chain A == antigen]
ABDB_CLEANED_PDB_DIR = Path(
    "/tgen_labs/altin/alphafold3/workspace/episcaf-experiments/data/abdb/complex_pdbfiles/cleaned"
)
