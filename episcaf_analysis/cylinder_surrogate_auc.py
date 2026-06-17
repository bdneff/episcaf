#!/usr/bin/env python3
"""
cylinder_surrogate_auc.py

How well does the cylinder (the no-antibody accessibility surrogate) predict the REAL
antibody clash on DP3, where the ground truth af3_n_clash_res exists? Reports the AUC for
predicting clash-free (af3_n_clash_res==0) from cylinder_ca_clashes and cylinder_native_aware,
plus Spearman correlation with the true clash count. This is the data-driven justification
for how much weight the cylinder term should carry.

Input: metrics_native_cyl.csv (cols af3_n_clash_res, cylinder_ca_clashes, cylinder_native_aware),
produced on Gemini by scripts/add_native_cylinder.py --exclude_dist 1.0.
"""
import numpy as np, pandas as pd
from pathlib import Path
from scipy.stats import rankdata

CSV = Path("/Users/bneff/Desktop/projects/episcaf/known_antigen/analysis/data/metrics_native_cyl.csv")

def auc_low_is_pos(score, y):           # low cylinder -> clash-free; rank on -score
    r = rankdata(-score); npos = y.sum(); nneg = len(y) - npos
    return (r[y == 1].sum() - npos*(npos+1)/2) / (npos*nneg)

def main():
    d = pd.read_csv(CSV, low_memory=False)
    tc = pd.to_numeric(d.af3_n_clash_res, errors="coerce")
    ca = pd.to_numeric(d.cylinder_ca_clashes, errors="coerce")
    na = pd.to_numeric(d.cylinder_native_aware, errors="coerce")
    m = tc.notna() & ca.notna() & na.notna()
    tc, ca, na = tc[m].values, ca[m].values, na[m].values
    y = (tc == 0).astype(int)
    print(f"DP3 valid n={len(y):,}; truly clash-free {int(y.sum()):,} ({100*y.mean():.1f}%)")
    print(f"AUC clash-free  | cylinder_ca_clashes  = {auc_low_is_pos(ca,y):.3f}")
    print(f"AUC clash-free  | cylinder_native_aware= {auc_low_is_pos(na,y):.3f}")
    print(f"Spearman vs true clash count | ca={pd.Series(ca).corr(pd.Series(tc),method='spearman'):.3f}"
          f"  native_aware={pd.Series(na).corr(pd.Series(tc),method='spearman'):.3f}")
    print(f"median cylinder (clash-free / clashing): "
          f"ca {np.median(ca[y==1]):.0f}/{np.median(ca[y==0]):.0f}, "
          f"native_aware {np.median(na[y==1]):.0f}/{np.median(na[y==0]):.0f}")

if __name__ == "__main__":
    main()
