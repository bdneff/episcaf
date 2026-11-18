#!/usr/bin/env python
"""Reproduce the binding-enrichment figures quoted in the episcaf web docs.

The site states a base binder rate, the fraction of designs passing the structural
filters, the binder rate among those passers, and a top-ranked-band enrichment. This
script regenerates all of them from the DP4 assay table the deployed binding model is
fitted to, so the numbers on Methods and the scaffold tool are reproducible rather than
asserted.

    python scripts/webtool_enrichment.py

Needs an env with pandas + scikit-learn (e.g. bcell-repair). Reads only; prints a table.
Source of truth: episcaf_analysis.binding_model.data.load_table (the same universe used
for model fitting) and BindingModel.cross_val_auc (grouped-by-epitope out-of-fold preds).
"""
import sys
import warnings
from pathlib import Path

warnings.filterwarnings("ignore")
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))  # repo root on path
import numpy as np

from episcaf_analysis.binding_model import data as D, model as M

# Four structural filters (assay-table column names); matches the documented pass gate:
# epitope RMSD <= 1 A, overall RMSD <= 2 A, mean PAE < 5 A, AF3 clashing residues == 0.
FILTER_COLS = ["epitope_rmsd", "overall_rmsd", "mean_pae", "af3_clashes"]


def main():
    df = D.load_table()
    y = df["bound"].values.astype(int)
    n, n_bind = len(df), int(y.sum())
    print(f"universe            : {n} designs, {n_bind} binders, "
          f"{df['epitope'].nunique()} epitopes")
    print(f"base binder rate    : {n_bind / n:.1%}")

    have = df[FILTER_COLS].notna().all(axis=1)
    passf = (have & (df.epitope_rmsd <= 1) & (df.overall_rmsd <= 2)
             & (df.mean_pae < 5) & (df.af3_clashes == 0)).values
    npass = int(passf.sum())
    print(f"four-filter passers : {npass}/{n} = {npass / n:.1%} of designs")
    print(f"  binder rate among passers: {y[passf].mean():.1%} "
          f"({int(y[passf].sum())}/{npass})")

    cv = M.BindingModel().cross_val_auc(df, repeats=5)
    oof = cv["oof"]  # grouped-by-epitope out-of-fold predicted probabilities
    print(f"grouped-CV ROC-AUC  : {cv['roc_mean']:.3f} +/- {cv['roc_sd']:.3f}")
    order = np.argsort(-oof)
    for frac in (0.05, 0.10, 0.25):
        k = max(1, int(round(n * frac)))
        print(f"  top {int(frac * 100):>2d}% by binding model (n={k}): "
              f"binder rate {y[order[:k]].mean():.1%}")


if __name__ == "__main__":
    main()
