"""Build the DP4 binder modeling table (designs x features + binding label + epitope group).

One place that reconstructs the same table the DP4 analyses used (scripts/dp4_binder_ml*.py), so the
shippable model and the feature evaluation read identical data:

  universe : scaffolded designs in data/libraries/dp4_library.csv whose epitope was assayed
             (present in epitope_struct_features.csv), with a composite score.
  label    : bound = library_member appears in any epitope's hitIDs
             (data/dp4_binding/john/scaffoldedEpitopeSummary.csv).
  group    : epitope = target.split("_")[0].lower()  (grouped-by-epitope CV uses this).

Feature blocks mirror the DP4 ablation (findings: design 0.74 -> +SS 0.78 -> +native/dyn 0.81):
  DESIGN  per-design AF3/metric columns
  SS      secondary-structure counts (helix, strand, n_islands)
  NATIVE  native-context composition/size/exposure from epitope_struct_features.csv
  DYN     the dynamics subset (B-factor z, GNM dfit, dSASA) -- weakest, confounded by bound state
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
LIB = ROOT / "data/libraries/dp4_library.csv"
SUMM = ROOT / "data/dp4_binding/john/scaffoldedEpitopeSummary.csv"
STRUCTF = ROOT / "data/dp4_binding/epitope_struct_features.csv"
NATIVE_TIER1 = ROOT / "results/feature_model/native_tier1_features.csv"

DESIGN = ["cylinder_clashes", "epitope_rmsd", "overall_rmsd", "epitope_pae", "scaffold_pae",
          "mean_pae", "ptm", "af3_clashes", "island_index"]
SS = ["helix", "strand", "n_islands"]
# NATIVE composition is defined for the ablation but EXCLUDED from the baseline: adding it under
# grouped-by-epitope CV drops ROC-AUC 0.812 -> 0.766 (too many epitope-level terms for 56 groups).
NATIVE = ["epi_size", "epi_rel_sasa", "epi_rg", "epi_frac_hydrophobic", "epi_frac_charged",
          "epi_frac_aromatic"]
DYN = ["epi_bfac_z", "epi_gnm_dfit", "epi_dsasa"]

# Shipped baseline: reproduces the published DP4 number (grouped-CV ROC-AUC 0.812, findings ~0.81).
BASELINE_FEATURES = DESIGN + SS + DYN

# Recommended clean-static alternative (point BindingModel here): the three static/apo Tier-1 features
# that MATCH-or-BEAT the confounded dynamics block (0.817 vs 0.812) with better causal footing --
# see results/feature_model/model_v1_findings.md.
CLEAN_STATIC = ["epi_coulomb", "frame_ratio", "ag_near15"]
CLEAN_STATIC_FEATURES = DESIGN + SS + CLEAN_STATIC

ID_COL = "library_member"
LABEL_COL = "bound"
GROUP_COL = "epitope"


def _hit_ids(summ):
    hits = set()
    for ids in summ["hitIDs"].dropna():
        for t in str(ids).split(","):
            t = t.strip()
            if t.startswith("DP4_"):
                hits.add("DP4_" + t.split("_")[1])
    return hits


def load_table(with_tier1=True):
    """Return the modeling DataFrame (one row per assayed scaffolded design)."""
    lib = pd.read_csv(LIB, low_memory=False)
    summ = pd.read_csv(SUMM, encoding="latin-1")
    summ.columns = [c.strip() for c in summ.columns]
    summ = summ.rename(columns={"#islands": "n_islands"})
    summ["epitope"] = summ["epitope"].astype(str).str.lower()
    st = pd.read_csv(STRUCTF)
    st["epitope"] = st["epitope"].astype(str).str.lower()

    epi = (summ[["epitope", "helix", "strand", "n_islands"]]
           .merge(st[["epitope"] + NATIVE + DYN], on="epitope", how="inner"))
    hits = _hit_ids(summ)
    lib["epitope"] = lib["target"].astype(str).str.split("_").str[0].str.lower()
    d = lib[lib["composite"].notna() & lib["epitope"].isin(set(epi["epitope"]))].copy()
    d[LABEL_COL] = d[ID_COL].isin(hits).astype(int)
    d = d.merge(epi, on="epitope", how="left")

    if with_tier1 and NATIVE_TIER1.exists():
        t1 = pd.read_csv(NATIVE_TIER1)
        t1["epitope"] = t1["epitope"].astype(str).str.lower()
        d = d.merge(t1.drop(columns=[c for c in ("epi_size_check",) if c in t1.columns]),
                    on="epitope", how="left")
    return d


if __name__ == "__main__":
    d = load_table()
    y = d[LABEL_COL].values
    print(f"designs={len(d):,} bound={int(y.sum())} epitopes={d[GROUP_COL].nunique()} "
          f"base rate={y.mean():.3f}")
    missing = [c for c in BASELINE_FEATURES if c not in d.columns]
    print("missing baseline features:", missing or "none")
