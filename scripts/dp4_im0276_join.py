#!/usr/bin/env python3
"""Join the IM0276 DP4 assay counts to the design metrics, and run the readouts that are
attributable given the pooled-mAb design.

The assay row ids are DP4_<n>_<m>, where <n> is the library_member number (contiguous 1..36000)
and <m> is a secondary tag. We parse <n> and merge on library_member. That is the whole join;
every metric is already a column in data/libraries/dp4_library.csv.

What is attributable:
  - 8pww is the ONE individually assayed antibody -> a cognate readout is possible.
  - CIDRa1.* is the minibinder / PfEMP1 readout -> handled separately (its own targets).
  - Baselines: DP4_NoAbNoSer_* (no antibody, no serum) is the enrichment baseline.
What is blocked:
  - DP4_mAbs_Pool1..8 pool the antibodies, and the pool->antibody composition is NOT in the
    demux output, so cognate log-enrichment cannot be attributed to a specific antibody. Raise
    with John before attempting the DP3-style per-antibody analysis.

Column semantics (dilutions, the _JA suffix) are read off the header, not documented; this script
prints the full classification so they can be confirmed.

Run:  /usr/bin/python3 scripts/dp4_im0276_join.py
"""
import re, sys
from pathlib import Path
import numpy as np, pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
ASSAY = ROOT / "data/dp4_binding/table_raw_full_IM0276_DP4.tsv"
LIB = ROOT / "data/libraries/dp4_library.csv"
OUT = ROOT / "data/dp4_binding/dp4_im0276_joined.csv"

def classify(cols):
    """assay condition column -> group."""
    groups = {}
    pat = [
        ("baseline_noabnoser", r"NoAbNoSer"),
        ("baseline_noabjaser", r"NoAbJASer"),
        ("baseline_noprot",    r"NoProt"),
        ("ab_8pww",            r"^DP4_8pww_"),
        ("mab_pool",           r"^DP4_mAbs_Pool"),
        ("minibinder_cidr",    r"^DP4_CIDRa1"),
        ("serum_cvax",         r"^DP4_CVAX_"),
        ("serum_abrc",         r"^DP4_ABRC_"),
        ("serum_vf",           r"^DP4_VF_"),
    ]
    for c in cols:
        for name, rx in pat:
            if re.search(rx, c, re.I):
                groups.setdefault(name, []).append(c); break
        else:
            groups.setdefault("other", []).append(c)
    return groups

def main():
    lib = pd.read_csv(LIB, low_memory=False)
    print(f"library: {lib.shape[0]:,} rows, {lib.shape[1]} cols")
    for k in ["library_member", "target", "category", "composite"]:
        print(f"  has {k}: {k in lib.columns}")
    print(f"  target examples: {lib['target'].dropna().astype(str).unique()[:6].tolist()}")
    cog8 = lib["target"].astype(str).str.lower().str.startswith("8pww")
    print(f"  8pww-cognate designs in library: {int(cog8.sum())}")

    assay = pd.read_csv(ASSAY, sep="\t", low_memory=False)
    idcol = assay.columns[0]
    print(f"\nassay: {assay.shape[0]:,} rows, {assay.shape[1]} cols; id column = {idcol!r}")
    assay["library_member"] = "DP4_" + assay[idcol].astype(str).str.split("_").str[1]

    cond = [c for c in assay.columns if c not in (idcol, "library_member")]
    groups = classify(cond)
    print("\ncondition columns by group (confirm these inferences):")
    for g in sorted(groups):
        print(f"  {g:20s} {len(groups[g]):3d}   e.g. {groups[g][:3]}")

    m = assay.merge(lib, on="library_member", how="inner")
    print(f"\nmerge: {len(m):,}/{len(assay):,} assay rows matched a library design "
          f"({100*len(m)/len(assay):.1f}%)")

    # ---- baseline + 8pww cognate readout (the one attributable antibody) ----
    base_cols = groups.get("baseline_noabnoser", [])
    ab_cols   = sorted(groups.get("ab_8pww", []))
    print(f"\nbaseline (NoAbNoSer) columns: {base_cols}")
    print(f"8pww dilution columns: {ab_cols}")
    if base_cols and ab_cols:
        l1p = lambda s: np.log10(1.0 + s)
        base = l1p(m[base_cols].astype(float)).mean(axis=1)   # mean log10(1+NoAbNoSer)
        for c in ab_cols:
            m[f"le_{c}"] = l1p(m[c].astype(float)) - base       # cognate log-enrichment
        # pick 1X as the representative dilution if present, else the middle one
        rep = next((c for c in ab_cols if re.search(r"_1X(_|$)", c)), ab_cols[len(ab_cols)//2])
        le = f"le_{rep}"
        cog = m[m["target"].astype(str).str.lower().str.startswith("8pww")]
        print(f"\n8pww readout (dilution {rep}):")
        print(f"  cognate designs assayed: {len(cog)}")
        if len(cog):
            print(f"  median cognate log-enrichment: {cog[le].median():+.3f}")
            print(f"  fraction above baseline (>0):  {(cog[le] > 0).mean():.2f}")
            for metric in ["composite", "cylinder_clashes", "epitope_rmsd"]:
                if metric in cog and cog[metric].notna().sum() > 3:
                    r, p = stats.spearmanr(cog[metric], cog[le], nan_policy="omit")
                    print(f"  Spearman({metric}, enrichment) = {r:+.2f} (p={p:.2g}, n={cog[metric].notna().sum()})")

    keep = ["library_member", idcol, "target", "category", "composite", "rank_in_group",
            "cylinder_clashes", "epitope_rmsd", "mean_pae", "is_global_pass"]
    keep = [c for c in keep if c in m.columns] + cond + [c for c in m.columns if c.startswith("le_")]
    m[keep].to_csv(OUT, index=False)
    print(f"\nwrote {OUT.relative_to(ROOT)}  ({len(m):,} rows, {len(keep)} cols)")
    print("\nBLOCKED until pool composition arrives: cognate attribution for DP4_mAbs_Pool1..8 "
          f"({len(groups.get('mab_pool', []))} columns).")

if __name__ == "__main__":
    main()
