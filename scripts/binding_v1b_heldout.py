#!/usr/bin/env python3
"""Held-out confirmation: v1 (bound-state dynamics block) vs v1b (clean-static features).

The v1b switch is for causal cleanliness, not a score gain, so the bar is NON-INFERIORITY: v1b must
not lose AUC out-of-sample. This runs a repeated, PAIRED grouped-by-epitope held-out comparison -- on
each identical train/test split both feature sets are fit and scored on the same held-out epitopes --
and reports the paired difference with a CI and a verdict.

Run: /Users/bneff/micromamba/envs/bcell-repair/bin/python scripts/binding_v1b_heldout.py
"""
import os
import sys
import warnings

import numpy as np
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedGroupKFold

warnings.simplefilter("ignore")
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from episcaf_analysis.binding_model import data as D, model as M   # noqa: E402

REPEATS, SPLITS = 40, 5


def main():
    d = D.load_table()
    y = d[D.LABEL_COL].to_numpy()
    g = d[D.GROUP_COL].to_numpy()
    Xv1 = d[D.BASELINE_FEATURES].to_numpy(dtype=float)
    Xv1b = d[D.CLEAN_STATIC_FEATURES].to_numpy(dtype=float)
    print(f"designs={len(d):,} bound={int(y.sum())} epitopes={d[D.GROUP_COL].nunique()} "
          f"base_rate={y.mean():.3f}")
    print(f"v1  = design+SS+dynamics ({len(D.BASELINE_FEATURES)} feat): {D.DYN}")
    print(f"v1b = design+SS+clean-static ({len(D.CLEAN_STATIC_FEATURES)} feat): {D.CLEAN_STATIC}")

    a1, a1b, diffs = [], [], []
    for s in range(REPEATS):
        for tr, te in StratifiedGroupKFold(SPLITS, shuffle=True, random_state=s).split(Xv1, y, g):
            if len(np.unique(y[te])) < 2:
                continue
            p1 = M._pipeline().fit(Xv1[tr], y[tr]).predict_proba(Xv1[te])[:, 1]
            p1b = M._pipeline().fit(Xv1b[tr], y[tr]).predict_proba(Xv1b[te])[:, 1]
            u1, u1b = roc_auc_score(y[te], p1), roc_auc_score(y[te], p1b)
            a1.append(u1); a1b.append(u1b); diffs.append(u1b - u1)

    a1, a1b, diffs = map(np.array, (a1, a1b, diffs))
    md, sem = diffs.mean(), diffs.std(ddof=1) / np.sqrt(len(diffs))
    lo, hi = md - 1.96 * sem, md + 1.96 * sem
    print(f"\nfolds evaluated: {len(diffs)}  ({REPEATS} repeats x {SPLITS} folds)")
    print(f"  v1   held-out ROC-AUC : {a1.mean():.3f} +/- {a1.std():.3f}")
    print(f"  v1b  held-out ROC-AUC : {a1b.mean():.3f} +/- {a1b.std():.3f}")
    print(f"  paired diff (v1b-v1)  : {md:+.4f}  95% CI [{lo:+.4f}, {hi:+.4f}]")
    print(f"  folds where v1b >= v1 : {(diffs >= 0).mean()*100:.0f}%")

    # non-inferiority margin: v1b must not be worse by more than 0.01 AUC
    MARGIN = -0.01
    if lo > MARGIN and md >= 0:
        verdict = "ADOPT v1b -- non-inferior (and causally cleaner: drops the confounded dynamics block)"
    elif lo > MARGIN:
        verdict = "ADOPT v1b -- within the non-inferiority margin; switch justified by causal cleanliness"
    else:
        verdict = "KEEP v1 -- v1b is meaningfully worse out-of-sample (lower CI below the -0.01 margin)"
    print(f"\nVERDICT: {verdict}")


if __name__ == "__main__":
    main()
