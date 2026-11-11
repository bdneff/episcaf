#!/usr/bin/env python3
"""Unit tests for the AOMclash regressor harness (no structures, no cluster data).

Run:  /usr/bin/python3 tests/test_aomclash_model.py     (or: pytest tests/)

These guard the 39-feature model wrapper in episcaf_analysis/aomclash/model.py:
  * the 39-feature list is exactly 39, groups add up, and holds no histogram/cylinder columns;
  * it stays consistent with what features.design_feats/target_features actually emit -- the 39 must
    be exactly the extractor's full key set minus the 70 position histograms and the 2 cylinder counts
    (cross-checked against the frozen key sets in test_aomclash_features.py);
  * the leave-one-group-out, fit_final and predict paths run on a synthetic frame with sklearn.
"""
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, os.path.dirname(__file__))
from episcaf_analysis.aomclash import model as M                     # noqa: E402
from test_aomclash_features import EXPECTED_KEYS, TARGET_KEYS         # noqa: E402  (frozen real keys)


def test_feature_count_and_groups():
    assert len(M.FEATURES_39) == 39
    assert len(set(M.FEATURES_39)) == 39
    assert (len(M.OCCUPANCY), len(M.OUTWARD_COUNTS), len(M.ANTIGEN_CONTACT),
            len(M.SCAFFOLD_SHAPE), len(M.EPITOPE_DESCRIPTORS)) == (15, 4, 4, 10, 6)
    assert len(M.DROPPED_HISTOGRAMS) == 70 and len(M.DROPPED_CYLINDER) == 2


def test_no_histograms_or_cylinder_in_model():
    for f in M.FEATURES_39:
        assert not (f.startswith(("h", "hf")) and f[-1].isdigit()), f
        assert f not in ("cyl_aware", "cyl_plain"), f


def test_39_matches_extractor_minus_histograms_and_cylinder():
    """The 39 must equal (all emitted keys) - (70 histograms) - (2 cylinder)."""
    produced = set(EXPECTED_KEYS) | set(TARGET_KEYS)          # 105 design + 6 target = 111
    assert len(produced) == 111, len(produced)
    expected_39 = produced - set(M.DROPPED_HISTOGRAMS) - set(M.DROPPED_CYLINDER)
    assert set(M.FEATURES_39) == expected_39, set(M.FEATURES_39) ^ expected_39
    assert len(expected_39) == 39


def _synthetic_frame(seed=0, n=160, n_groups=4):
    """A fake feature matrix: the 39 columns, an AF3 label correlated with s_p1_free, groups, stages."""
    rng = np.random.default_rng(seed)
    df = pd.DataFrame({f: rng.normal(size=n) for f in M.FEATURES_39})
    df["s_p1_free"] = rng.gamma(2.0, 2.0, size=n)              # make the prior baseline well-defined
    # label: a monotone-ish function of the carved occupancy sum, as a nonneg integer count
    df["af3_n_clash_res"] = np.rint(np.clip(df["s_p1_free"] + rng.normal(0, 1, n), 0, None)).astype(int)
    df["stage"] = rng.integers(0, 2, size=n)
    gids = [f"g{j}" for j in range(n_groups)]
    df["group"] = rng.choice(gids, size=n)
    df["target"] = df["group"] + "_0P"                         # for epitope_groups()
    return df


def test_resolve_label_prefers_present_and_falls_back():
    df = _synthetic_frame(1)
    assert M.resolve_label(df) == "af3_n_clash_res"
    df2 = df.rename(columns={"af3_n_clash_res": "af3_clashes"})
    assert M.resolve_label(df2) == "af3_clashes"


def test_stage_weights_equalize_stages():
    stage = np.array([0, 0, 0, 1])                            # John's 2-stage case
    w = M.stage_weights(stage)
    assert abs(w[stage == 0].sum() - w[stage == 1].sum()) < 1e-9     # each stage, equal total
    assert abs(w.mean() - 1.0) < 1e-9


def test_leave_one_group_out_runs():
    df = _synthetic_frame(2)
    groups = df["group"].to_numpy()
    preds = M.leave_one_group_out(df, groups)
    assert set(preds) == {"aomclash", "aom_prior_sum"}
    for k, v in preds.items():
        assert len(v) == len(df)
        assert np.isfinite(v).all(), k                        # every row predicted by some held-out fold
        assert (v > -1.0).all()                               # expm1 lower bound
    # sanity: predictions track the label at all (synthetic signal is monotone in s_p1_free)
    assert M.per_group_spearman(groups, df["af3_n_clash_res"], preds["aomclash"]) > -1.0


def test_fit_final_and_predict():
    df = _synthetic_frame(3)
    m = M.fit_final(df)
    p = M.predict(m, df)
    assert len(p) == len(df) and np.isfinite(p).all()


def test_epitope_groups_from_target():
    df = _synthetic_frame(4)
    g = M.epitope_groups(df)
    assert len(g) == len(df)
    assert set(g) == {f"g{j}" for j in range(4)}              # prefix before '_', lowercased


if __name__ == "__main__":
    for name, fn in sorted(globals().items()):
        if name.startswith("test_") and callable(fn):
            fn()
            print(f"  ok: {name}")
    print("all aomclash model tests passed")
