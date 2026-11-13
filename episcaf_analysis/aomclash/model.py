"""AOMclash scorer: score antibody-clash propensity from a prepared feature frame.

Port of John Altin's `aomclash/aom_model.py` (AltinLab/epitope-universe). The map, the epitope frame
and the per-design feature extraction live in `features.py`; this module is only the scorer. It never
reads structures and never trains a map, so it runs anywhere the feature matrix exists. Building that
matrix needs the AF3 design structures, which are cluster-only ($WS); this module is staged to run
there, not here.

The shipped scorer: one feature
-------------------------------
We ship a SINGLE interpretable metric: `s_p1_free` -- the antibody-occupancy map read at every
scaffold atom, summed over only the atoms NOT inside native-antigen space (`free`, dag > 4 A). That
is "how much of the design sticks into antibody territory, beyond what mimicking the native antigen
already explains." `fit_prior` / `score_prior` below are the shipped path: an isotonic calibration of
`s_p1_free` onto the AF3 clash count. Because it is one monotone feature, the calibration only sets a
clash-count SCALE -- for ranking designs within a target (the deployment use) the raw sum ranks
identically. There is no learned multi-feature regressor to overfit, and the whole scorer reproduces
from (occupancy map + native antigen) alone.

Why not the 111 (or 39) features? John's ablation (`aom_ablation.py`, `results/ablation.csv`) and the
8-Nov report put the single carved-occupancy feature at ~0.87-0.89 and the full model only ~0.02
higher. Adding ~38 highly-correlated features for +0.02, on a set this size, buys overfitting risk,
not signal. So the 39-feature gradient-boosted regressor below (`FEATURES_39`, `gbm`, `fit_final`,
`leave_one_group_out`'s `aomclash` head) is KEPT as the tested, rejected alternative and the
head-to-head baseline -- not the shipped path. [Margins are from John's AF2 campaign; confirm the
single-feature result holds on our AF3 (DP3/DP4) labels when the matrix exists, but shipping does not
gate on the 39-feature model.] See `results/aomclash/feature_set_decision.md`.

Fold definition  [DECISION / partly UNVERIFIED against John's latest code]
--------------------------------------------------------------------------
We default to **leave-one-EPITOPE-out**: pass a `groups` array keyed by epitope identity so every fold
is tested on an epitope (and same-sequence epitopes) it never trained on. This follows the 8-Nov
AOMclash report and our project decision.

NOTE: John's newest *committed* code (`upstream/feat/clash-method`, tip 2026-11-04 -- `aom_model.py`
there is byte-identical to `feat/aomclash`) groups by `biological_antigen_key`, i.e. leave-one-
ANTIGEN-out, and keeps the 70 histograms (109 model inputs). The report's epitope-out / 111 is not in
any committed branch. We follow the report; switching back to antigen folds is a one-arg change (pass
an antigen-keyed `groups`), and restoring the full bank is `feats=None`. Flagged so the choice is not
silent.

Label
-----
AF3 clashing-residue count. `af3_n_clash_res` in stage05 metrics; renamed `af3_clashes` in the shipped
libraries (`data/libraries/dp4_*.csv`). `resolve_label` accepts either.
"""
from __future__ import annotations

import numpy as np

# --- the 39-feature subset, by group (names must match features.design_feats / target_features) ---
OCCUPANCY = [
    "s_p1", "s_p2", "s_p1_free", "s_p1_near",
    "n_p05", "n_p15", "n_p30", "n_p50", "n_p15_free", "n_p30_free",
    "r_sum", "r_sum_free", "r_n30", "r_n50", "r_n70",
]                                                                    # 15
OUTWARD_COUNTS = ["b_near20", "b_out25", "b_out15", "b_out35_hi"]    # 4
ANTIGEN_CONTACT = ["ag3", "ag4", "ag6", "ag_free_out25"]            # 4
SCAFFOLD_SHAPE = ["z_mean", "z_p90", "z_max", "r_mean", "ca_rg",
                  "n_scaf", "n_seg", "seg_max", "seg_sum", "apr"]    # 10
EPITOPE_DESCRIPTORS = ["frame_ratio", "n_epi", "n_isl",
                       "ag_out25", "ag_out25_hi", "ag_near15"]       # 6

FEATURES_39 = (OCCUPANCY + OUTWARD_COUNTS + ANTIGEN_CONTACT
               + SCAFFOLD_SHAPE + EPITOPE_DESCRIPTORS)

# what we deliberately leave out of the model (computed by the extractor, not fed to the regressor)
DROPPED_HISTOGRAMS = [f"h{j}" for j in range(35)] + [f"hf{j}" for j in range(35)]   # 70
DROPPED_CYLINDER = ["cyl_aware", "cyl_plain"]                                        # 2

assert len(FEATURES_39) == 39, len(FEATURES_39)
assert len(set(FEATURES_39)) == 39, "duplicate feature name"

# John's exact HistGradientBoostingRegressor hyperparameters (aom_model.gbm)
GBM_KW = dict(max_iter=300, learning_rate=0.05, max_leaf_nodes=15,
              min_samples_leaf=40, l2_regularization=1.0, random_state=0)

# candidate AF3 clash-label column names, most-specific first
LABEL_CANDIDATES = ("af3_n_clash_res", "af3_clashes")

# the single feature we ship as the scorer: the carved-occupancy sum (occupancy over free atoms)
SHIPPED_FEATURE = "s_p1_free"


def gbm():
    """A fresh regressor with John's hyperparameters."""
    from sklearn.ensemble import HistGradientBoostingRegressor
    return HistGradientBoostingRegressor(**GBM_KW)


def resolve_label(df, label=None):
    """Return the AF3 clash-count column name present in `df` (or the one asked for)."""
    if label is not None:
        if label not in df.columns:
            raise KeyError(f"label column {label!r} not in frame")
        return label
    for c in LABEL_CANDIDATES:
        if c in df.columns:
            return c
    raise KeyError(f"no AF3 clash label; looked for {LABEL_CANDIDATES}")


def stage_weights(stage):
    """Weight rows so each stage (e.g. post-RFD3 vs post-AF) contributes equally in total.

    Port of aom_model's `w = where(stage==0, .5/n0, .5/n1); w/=w.mean()`, generalized to k stages.
    For a single stage this is all ones.
    """
    stage = np.asarray(stage)
    uniq = np.unique(stage)
    w = np.zeros(len(stage), dtype=float)
    for s in uniq:
        m = stage == s
        w[m] = (1.0 / len(uniq)) / m.sum()
    return w / w.mean()


def _weights(df, stage_col):
    if stage_col is None or stage_col not in df.columns:
        return np.ones(len(df), dtype=float)
    return stage_weights(df[stage_col].to_numpy())


# ------------------------------------------------------------------ the shipped scorer (one feature)

def fit_prior(df, label=None, prior_feature=SHIPPED_FEATURE, stage_col="stage"):
    """Fit the SHIPPED scorer: isotonic calibration of the carved-occupancy sum onto AF3 clash count.

    Returns a fitted sklearn `IsotonicRegression` mapping `prior_feature` -> expected clashing-residue
    count. This *is* the whole shipped model -- one monotone feature, no multi-feature regressor;
    pickle it as the artifact. For ranking designs within a target the calibration is cosmetic (it is
    monotone), so the raw `prior_feature` ranks identically; the fit only puts the score on a
    clash-count scale comparable to `af3_n_clash_res` on the antibody arm.
    """
    from sklearn.isotonic import IsotonicRegression
    label = resolve_label(df, label)
    if prior_feature not in df.columns:
        raise KeyError(f"{prior_feature!r} not in frame")
    x = df[prior_feature].to_numpy(dtype=float)
    y = df[label].to_numpy(dtype=float)
    w = _weights(df, stage_col)
    return IsotonicRegression(out_of_bounds="clip").fit(x, y, sample_weight=w)


def score_prior(iso, values):
    """Predicted clashing-residue count for each design's carved-occupancy sum (the shipped scorer)."""
    return iso.predict(np.asarray(values, dtype=float))


# ------------------------------------------------------------------ the rejected 39-feature baseline

def leave_one_group_out(df, groups, feats=FEATURES_39, label=None, stage_col="stage",
                        prior_feature=SHIPPED_FEATURE):
    """Leave-one-group-out predictions (port of aom_model.leave_one_antigen_out, group-agnostic).

    Every row is predicted by a model that never trained on its group. `groups` is one value per row
    (epitope identity for leave-one-epitope-out; pass an antigen key for antigen folds).

    Returns {"aomclash": preds, "aom_prior_sum": preds}:
      aomclash      -- the 39-feature gradient-boosted regressor (expm1 of its log1p prediction).
      aom_prior_sum -- isotonic calibration of the single carved-occupancy sum `prior_feature`
                       (the ~0.87 map-only baseline). Omitted (NaN) if that column is absent.
    """
    feats = list(feats)
    label = resolve_label(df, label)
    missing = [c for c in feats if c not in df.columns]
    if missing:
        raise KeyError(f"feature frame is missing {len(missing)} columns, e.g. {missing[:6]}")

    y = df[label].to_numpy(dtype=float)
    ly = np.log1p(y)
    g = np.asarray(groups)
    if len(g) != len(df):
        raise ValueError("groups length != number of rows")
    w = _weights(df, stage_col)

    pred = np.full(len(df), np.nan)
    prior = np.full(len(df), np.nan)
    do_prior = prior_feature in df.columns
    if do_prior:
        from sklearn.isotonic import IsotonicRegression

    X = df[feats]
    for grp in np.unique(g):
        te = g == grp
        tr = ~te
        if not tr.any():
            continue
        model = gbm().fit(X[tr], ly[tr], sample_weight=w[tr])
        pred[te] = np.expm1(model.predict(X[te]))
        if do_prior:
            iso = IsotonicRegression(out_of_bounds="clip").fit(
                df[prior_feature].to_numpy()[tr], y[tr], sample_weight=w[tr])
            prior[te] = iso.predict(df[prior_feature].to_numpy()[te])
    return {"aomclash": pred, "aom_prior_sum": prior}


def fit_final(df, feats=FEATURES_39, label=None, stage_col="stage"):
    """Train the shipping regressor on ALL rows (no held-out fold). Returns the fitted estimator.

    Target is log1p(clash); use `predict()` to get a residue count back.
    """
    feats = list(feats)
    label = resolve_label(df, label)
    ly = np.log1p(df[label].to_numpy(dtype=float))
    w = _weights(df, stage_col)
    return gbm().fit(df[feats], ly, sample_weight=w)


def predict(model, df, feats=FEATURES_39):
    """Predicted clashing-residue count for each row (expm1 of the model's log1p output)."""
    return np.expm1(model.predict(df[list(feats)]))


def per_group_spearman(groups, y_true, y_pred):
    """Median within-group Spearman rho (the deployment metric: ranking designs within a target)."""
    from scipy.stats import spearmanr
    g = np.asarray(groups)
    y_true = np.asarray(y_true, dtype=float)
    y_pred = np.asarray(y_pred, dtype=float)
    rhos = []
    for grp in np.unique(g):
        m = (g == grp) & np.isfinite(y_pred) & np.isfinite(y_true)
        if m.sum() >= 3 and np.ptp(y_true[m]) > 0 and np.ptp(y_pred[m]) > 0:
            rhos.append(spearmanr(y_pred[m], y_true[m]).statistic)
    return float(np.nanmedian(rhos)) if rhos else float("nan")


def epitope_groups(df, target_col="target"):
    """Derive an epitope key per row from a `target` column (matches dp4_binder_ml: prefix before '_').

    This is the default leave-one-EPITOPE-out grouping. Adjust if the feature frame carries an explicit
    epitope/patch id column instead.
    """
    if target_col not in df.columns:
        raise KeyError(f"no {target_col!r} column to derive epitope groups from; "
                       "pass `groups` explicitly")
    return df[target_col].astype(str).str.split("_").str[0].str.lower().to_numpy()


def _main(argv=None):
    """CLI: evaluate the 39-feature model leave-one-group-out on a prepared feature frame.

    The frame (CSV/Parquet/pickle) must carry the 39 feature columns, an AF3 clash label, a grouping
    column (or a `target` column to derive epitope groups), and optionally a `stage` column.
    """
    import argparse
    from pathlib import Path
    import pandas as pd

    ap = argparse.ArgumentParser(description=_main.__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--frame", required=True, help="feature matrix: .csv / .parquet / .pkl")
    ap.add_argument("--label", default=None, help=f"clash label column (default: first of {LABEL_CANDIDATES})")
    ap.add_argument("--group-col", default=None,
                    help="grouping column for folds; omit to derive epitope groups from --target-col")
    ap.add_argument("--target-col", default="target", help="column to derive epitope groups from")
    ap.add_argument("--stage-col", default="stage", help="stage column for equal-weighting (optional)")
    ap.add_argument("--out", default=None, help="write per-row predictions CSV here")
    a = ap.parse_args(argv)

    p = Path(a.frame)
    if p.suffix == ".parquet":
        df = pd.read_parquet(p)
    elif p.suffix in (".pkl", ".pickle"):
        df = pd.read_pickle(p)
    else:
        df = pd.read_csv(p, low_memory=False)

    groups = df[a.group_col].to_numpy() if a.group_col else epitope_groups(df, a.target_col)
    label = resolve_label(df, a.label)
    preds = leave_one_group_out(df, groups, label=label, stage_col=a.stage_col)

    from scipy.stats import spearmanr
    y = df[label].to_numpy(dtype=float)
    for name, pr in preds.items():
        m = np.isfinite(pr)
        if not m.any():
            print(f"{name:16s}  (not computed)")
            continue
        pooled = spearmanr(pr[m], y[m]).statistic
        within = per_group_spearman(groups[m], y[m], pr[m])
        mae = np.abs(pr[m] - y[m]).mean()
        print(f"{name:16s}  within-group rho={within:.3f}  pooled rho={pooled:.3f}  MAE={mae:.2f}  (n={m.sum()})")

    if a.out:
        out = df.copy()
        for name, pr in preds.items():
            out[name] = pr
        out.to_csv(a.out, index=False)
        print("wrote", a.out)


if __name__ == "__main__":
    _main()
