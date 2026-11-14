"""Shippable interpretable binding-model scorer.

An L2-regularized logistic regression (the DP4 result: binding is additive and mildly non-linear, so
the model stays interpretable, not a black box) over design + secondary-structure + native-context
features. Honest evaluation is grouped by epitope (test antibodies never seen in training).

Public API:
  BindingModel(features=...)          choose the feature list (defaults to data.BASELINE_FEATURES)
  .fit(df)                            fit the imputer+scaler+logistic pipeline
  .predict(df) -> DataFrame           design_id + binding_pred probability, keyed by library_member
  .cross_val_auc(df, repeats=5)       repeated StratifiedGroupKFold ROC/PR-AUC (mean, sd)
  .save(path) / BindingModel.load(path)
  .coefficients()                     standardized logistic weights (what the model leans on)

Fit + save the tracked artifact (results-tab ranker):
  <env>/bin/python -m episcaf_analysis.binding_model.model
"""
from __future__ import annotations

import pickle
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, roc_auc_score
from sklearn.model_selection import StratifiedGroupKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

from . import data as D

ROOT = Path(__file__).resolve().parents[2]
ARTIFACT = Path(__file__).resolve().parent / "binding_model_v1.pkl"


def _pipeline(C=1.0):
    """Median-impute -> standardize -> L2 logistic (class_weight balanced for the 8% base rate)."""
    return Pipeline([
        ("impute", SimpleImputer(strategy="median")),
        ("scale", StandardScaler()),
        # LogisticRegression defaults to L2; C is the inverse L2 strength (the interpretable dial).
        ("logit", LogisticRegression(C=C, max_iter=5000, class_weight="balanced")),
    ])


@dataclass
class BindingModel:
    # v1b (adopted 2026-11-13): the clean-static feature set. A held-out, paired, grouped-by-epitope
    # comparison (scripts/binding_v1b_heldout.py) had v1b beat the bound-state dynamics block by
    # +0.013 ROC-AUC, 95% CI [+0.008,+0.018] -- reliably better out-of-sample AND causally cleaner
    # (nothing downstream of binding). Use D.BASELINE_FEATURES to get the old dynamics-block v1.
    features: list = field(default_factory=lambda: list(D.CLEAN_STATIC_FEATURES))
    C: float = 1.0
    pipeline: Pipeline = None
    id_col: str = D.ID_COL

    def fit(self, df):
        X = df[self.features].to_numpy(dtype=float)
        y = df[D.LABEL_COL].to_numpy()
        self.pipeline = _pipeline(self.C).fit(X, y)
        return self

    def predict(self, df):
        """Return a DataFrame with the design id and a `binding_pred` probability (P hit)."""
        if self.pipeline is None:
            raise RuntimeError("model is not fitted")
        X = df[self.features].to_numpy(dtype=float)
        p = self.pipeline.predict_proba(X)[:, 1]
        out = pd.DataFrame({"design_id": df[self.id_col].to_numpy(), "binding_pred": p})
        return out

    def cross_val_auc(self, df, repeats=5, n_splits=5):
        """Repeated grouped-by-epitope CV. Returns dict with roc/pr mean+sd and the OOF predictions."""
        X = df[self.features].to_numpy(dtype=float)
        y = df[D.LABEL_COL].to_numpy()
        g = df[D.GROUP_COL].to_numpy()
        rocs, prs = [], []
        oof_last = np.zeros(len(y))
        for s in range(repeats):
            oof = np.zeros(len(y))
            for tr, te in StratifiedGroupKFold(n_splits, shuffle=True, random_state=s).split(X, y, g):
                m = _pipeline(self.C).fit(X[tr], y[tr])
                oof[te] = m.predict_proba(X[te])[:, 1]
            rocs.append(roc_auc_score(y, oof))
            prs.append(average_precision_score(y, oof))
            oof_last = oof
        return dict(roc_mean=float(np.mean(rocs)), roc_sd=float(np.std(rocs)),
                    pr_mean=float(np.mean(prs)), pr_sd=float(np.std(prs)), oof=oof_last)

    def coefficients(self):
        logit = self.pipeline.named_steps["logit"]
        return dict(sorted(zip(self.features, logit.coef_[0]), key=lambda t: t[1]))

    def save(self, path=ARTIFACT):
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        with open(path, "wb") as fh:
            pickle.dump({"features": self.features, "C": self.C, "pipeline": self.pipeline,
                         "id_col": self.id_col}, fh)
        return path

    @classmethod
    def load(cls, path=ARTIFACT):
        with open(path, "rb") as fh:
            st = pickle.load(fh)
        m = cls(features=st["features"], C=st["C"], id_col=st.get("id_col", D.ID_COL))
        m.pipeline = st["pipeline"]
        return m


def main():
    df = D.load_table()
    model = BindingModel()
    cv = model.cross_val_auc(df)
    print(f"designs={len(df):,} bound={int(df[D.LABEL_COL].sum())} "
          f"epitopes={df[D.GROUP_COL].nunique()} base rate={df[D.LABEL_COL].mean():.3f}")
    print(f"features ({len(model.features)}): {model.features}")
    print(f"grouped-CV ROC-AUC {cv['roc_mean']:.3f} +/- {cv['roc_sd']:.3f}   "
          f"PR-AUC {cv['pr_mean']:.3f} +/- {cv['pr_sd']:.3f}")
    model.fit(df)
    path = model.save()
    print(f"\nsaved artifact: {path.relative_to(ROOT)}")
    pred = model.predict(df.head(5))
    print("\nexample predictions:")
    print(pred.to_string(index=False))


if __name__ == "__main__":
    main()
