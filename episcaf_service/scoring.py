"""Score a completed run's metrics table (the /results endpoint's work).

Pure CPU: reads the AF3 metrics CSV the cluster wrote (already carrying `aomclash_pred` from
`aomclash_score_designs.py`), attaches the interpretable binding-model probability, and applies the
config-driven composite scorer. Reuses the shipped artifacts -- no reimplementation:
  - `episcaf_analysis.binding_model.BindingModel` -> `binding_pred`
  - `episcaf_analysis.score.score(df, preset)` + `episcaf_analysis.presets.PRESETS`  -> `composite`
"""
from __future__ import annotations

import os

import pandas as pd

from episcaf_analysis.presets import PRESETS
from episcaf_analysis.score import score as _composite


def _read_metrics(path):
    return pd.read_parquet(path) if str(path).endswith(".parquet") else pd.read_csv(path)


def score_run(run_dir, preset="twelvemer", topk=0):
    """Return (ranked_df, info). Expects $RUN_DIR/05_analysis/metrics.csv to exist.

    `info` reports which steps ran so the caller/UI can be honest about partial scoring.
    """
    metrics = os.path.join(run_dir, "05_analysis", "metrics.csv")
    if not os.path.exists(metrics):
        raise FileNotFoundError(f"no metrics yet at {metrics} -- the cluster run is not scored")
    df = _read_metrics(metrics)
    info = {"n_designs": int(len(df)), "binding_pred": False, "composite": False, "preset": preset}

    # 1. interpretable binding probability (skip cleanly if its features are not all present).
    try:
        from episcaf_analysis.binding_model.model import BindingModel
        bm = BindingModel.load()
        if set(bm.features).issubset(df.columns):
            pred = bm.predict(df)
            df = df.merge(pred, left_on=bm.id_col, right_on="design_id", how="left")
            info["binding_pred"] = True
    except Exception as e:                                   # never fail the whole request on scoring
        info["binding_pred_error"] = f"{type(e).__name__}: {e}"

    # 2. config-driven composite (renormalizes over whichever metric columns are present).
    if preset not in PRESETS:
        raise ValueError(f"unknown preset {preset!r}; choose from {sorted(PRESETS)}")
    try:
        df = _composite(df, PRESETS[preset])
        info["composite"] = True
        sort_key = "composite" if "composite" in df.columns else None
    except Exception as e:
        info["composite_error"] = f"{type(e).__name__}: {e}"
        sort_key = None

    if sort_key is None:
        sort_key = "binding_pred" if "binding_pred" in df.columns else df.columns[0]
    df = df.sort_values(sort_key, ascending=False, kind="mergesort").reset_index(drop=True)
    if topk and topk > 0:
        df = df.head(topk)
    info["ranked_by"] = sort_key
    return df, info
