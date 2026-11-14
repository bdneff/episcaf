"""Evaluate the binding-model feature set: baseline ablation, then the Tier-1 NEW-feature keep test.

Two parts:
 (1) ablation of the existing blocks (design -> +SS -> +dynamics -> +native composition) to fix the
     baseline empirically under grouped-by-epitope CV;
 (2) for each NEW Tier-1 native feature (results/feature_model/native_tier1_features.csv): its max
     |correlation| with the baseline features (orthogonality) and its MARGINAL ROC-AUC when added to
     the baseline, same repeated grouped CV. Keep rule (FEATURE_CATALOG.md): keep only features that
     add meaningful marginal AUC (a real gain, not ~0.001) and are not near-redundant.

Run:  <env>/bin/python scripts/binding_feature_eval.py
"""
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from episcaf_analysis.binding_model import data as D          # noqa: E402
from episcaf_analysis.binding_model.model import BindingModel  # noqa: E402
from episcaf_analysis.binding_model.native_features import NEW_FEATURES  # noqa: E402

REPEATS = 5


def cv(df, feats):
    m = BindingModel(features=feats)
    r = m.cross_val_auc(df, repeats=REPEATS)
    return r["roc_mean"], r["roc_sd"], r["pr_mean"]


def main():
    df = D.load_table()
    y = df[D.LABEL_COL].values
    print(f"designs={len(df):,} bound={int(y.sum())} epitopes={df[D.GROUP_COL].nunique()} "
          f"base rate={y.mean():.3f}  (repeated grouped-by-epitope CV, {REPEATS} shuffles)\n")

    # ---- part 1: ablation of existing blocks ----
    blocks = [("design", D.DESIGN),
              ("+ secondary structure", D.DESIGN + D.SS),
              ("+ dynamics", D.DESIGN + D.SS + D.DYN),
              ("+ native composition", D.DESIGN + D.SS + D.DYN + D.NATIVE)]
    print(f"{'feature block':26s} {'nfeat':>5s} {'ROC-AUC':>14s} {'PR-AUC':>7s}")
    best = None
    for name, feats in blocks:
        rm, rs, pm = cv(df, feats)
        print(f"{name:26s} {len(feats):5d}   {rm:.3f}+/-{rs:.3f}   {pm:.3f}")
        if best is None or (name == "+ dynamics"):
            pass
    # baseline = design + SS + dynamics (the set that reproduces the published ~0.81)
    BASELINE = D.DESIGN + D.SS + D.DYN
    base_rm, base_rs, base_pm = cv(df, BASELINE)
    print(f"\nBASELINE for keep test = design+SS+dynamics ({len(BASELINE)} feat): "
          f"ROC-AUC {base_rm:.3f}+/-{base_rs:.3f}\n")

    # ---- part 2: NEW-feature orthogonality + marginal AUC ----
    Xb = df[BASELINE].to_numpy(float)
    print(f"{'new feature':16s} {'max|r| vs base':>14s} {'+feat ROC-AUC':>14s} {'marginal':>9s}")
    results = []
    for f in NEW_FEATURES:
        v = df[f].to_numpy(float)
        ok = np.isfinite(v)
        # orthogonality: max |Pearson r| with any baseline feature (finite rows)
        cors = []
        for j in range(Xb.shape[1]):
            cj = Xb[ok, j]
            m = np.isfinite(cj)
            if m.sum() > 2 and np.std(v[ok][m]) > 0 and np.std(cj[m]) > 0:
                cors.append(abs(np.corrcoef(v[ok][m], cj[m])[0, 1]))
        maxr = max(cors) if cors else float("nan")
        rm, rs, pm = cv(df, BASELINE + [f])
        marg = rm - base_rm
        results.append((f, maxr, rm, rs, marg))
        print(f"{f:16s} {maxr:14.3f} {rm:8.3f}+/-{rs:.3f} {marg:+9.4f}")

    # also: all kept-candidate features together
    KEEP_THRESH = 0.005
    kept = [f for f, maxr, rm, rs, marg in results if marg >= KEEP_THRESH and maxr < 0.9]
    print(f"\nkeep rule: marginal ROC-AUC >= {KEEP_THRESH} and max|r| < 0.9")
    print("kept:", kept or "(none)")
    if kept:
        rm, rs, pm = cv(df, BASELINE + kept)
        print(f"baseline + all kept ({len(BASELINE)+len(kept)} feat): "
              f"ROC-AUC {rm:.3f}+/-{rs:.3f}  (baseline {base_rm:.3f}, gain {rm-base_rm:+.4f})")

    # full-set reference (everything), for the record
    allf = BASELINE + NEW_FEATURES
    rm, rs, pm = cv(df, allf)
    print(f"baseline + ALL new ({len(allf)} feat): ROC-AUC {rm:.3f}+/-{rs:.3f}")

    # ---- part 3: marginal on the LEANER design+SS baseline (does the signal exist before the
    #      existing native/dynamics block absorbs it?) + the catalog's grouped candidates ----
    LEAN = D.DESIGN + D.SS
    lean_rm, _, _ = cv(df, LEAN)
    print(f"\n--- marginal on LEAN baseline design+SS ({len(LEAN)} feat, ROC-AUC {lean_rm:.3f}) ---")
    print(f"{'new feature/group':20s} {'+feat ROC-AUC':>14s} {'marginal':>9s}")
    groups = [(f, [f]) for f in NEW_FEATURES]
    groups += [("crowding_trio", ["ag_out25", "ag_out25_hi", "ag_near15"]),
               ("charge_pair", ["epi_net_charge", "epi_abs_charge"])]
    for name, feats in groups:
        rm, rs, pm = cv(df, LEAN + feats)
        print(f"{name:20s} {rm:8.3f}+/-{rs:.3f} {rm-lean_rm:+9.4f}")


if __name__ == "__main__":
    main()
