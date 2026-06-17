#!/usr/bin/env python3
"""Config-driven composite scorer for episcaf designs.

Pipeline:  load CSV -> gate -> per-metric transform (within scope) -> weighted sum
           -> top-k per group -> write.

Each metric is transformed by its own activation (percentile / minmax / zscore /
sigmoid / identity), oriented so higher = better, then weighted-summed. The dials
live in presets.py. See docs/REORG.md for the design rationale (and the v2 hook for
fitting weights against DP3 `is_pass`).

Usage:
  python -m episcaf_analysis.score --preset twelvemer \
      --metrics-csv /path/metrics_12mer.csv --out /path/composite_12mer_top5.csv
  # or, run as a plain script from inside episcaf_analysis/:
  python score.py --preset antibody --metrics-csv ... --out ...
"""
from __future__ import annotations

import argparse
import sys

import numpy as np
import pandas as pd

try:                                   # works as a package module …
    from episcaf_analysis.presets import PRESETS
except ImportError:                    # … or as a plain script in this dir
    from presets import PRESETS


def _transform_series(x: pd.Series, spec: dict) -> pd.Series:
    """Return a 'higher = better' score for one metric over one population."""
    better = spec.get("better", "low")
    kind = spec.get("transform", "percentile")
    x = pd.to_numeric(x, errors="coerce")

    if kind == "percentile":
        r = x.rank(pct=True)
        return r if better == "high" else (1.0 - r)
    if kind == "minmax":
        lo, hi = x.min(), x.max()
        z = (x - lo) / (hi - lo) if hi > lo else x * 0.0
        return z if better == "high" else (1.0 - z)
    if kind == "zscore":
        mu, sd = x.mean(), x.std(ddof=0)
        z = (x - mu) / sd if sd > 0 else x * 0.0
        return z if better == "high" else -z
    if kind == "sigmoid":
        midpoint = float(spec.get("midpoint", x.median()))
        k = float(spec.get("k", 1.0))
        s = 1.0 if better == "low" else -1.0
        return 1.0 / (1.0 + np.exp(s * k * (x - midpoint)))
    if kind == "identity":
        return x if better == "high" else -x
    raise ValueError(f"unknown transform: {kind!r}")


def _scoped(df: pd.DataFrame, col: str, spec: dict, scope: str, antigen_col: str) -> pd.Series:
    if scope == "per_antigen" and antigen_col in df.columns:
        return df.groupby(antigen_col)[col].transform(lambda s: _transform_series(s, spec))
    return _transform_series(df[col], spec)


def score(df: pd.DataFrame, preset: dict) -> pd.DataFrame:
    df = df.copy()

    # 1. gate (lower-better) ----------------------------------------------------
    gcol, gthr = preset["gate"]
    if gcol not in df.columns:
        sys.exit(f"[score] gate column {gcol!r} not in CSV. columns: {list(df.columns)}")
    n0 = len(df)
    df = df[pd.to_numeric(df[gcol], errors="coerce") <= gthr].copy()
    print(f"[score] gate {gcol} <= {gthr}: {len(df)}/{n0} rows kept")

    # 2. per-metric transforms within scope ------------------------------------
    scope = preset.get("scope", "pooled")
    antigen_col = preset.get("antigen_col", "antigen")
    metrics = preset["metrics"]
    present = {c: s for c, s in metrics.items() if c in df.columns}
    missing = [c for c in metrics if c not in df.columns]
    if missing:
        print(f"[score] WARNING: metrics absent from CSV, dropped: {missing}")
    if not present:
        sys.exit("[score] no scoring metrics present in CSV.")
    wsum = sum(s["weight"] for s in present.values())  # renormalize over present metrics

    composite = pd.Series(0.0, index=df.index)
    for col, spec in present.items():
        sc = _scoped(df, col, spec, scope, antigen_col)
        df[f"score_{col}"] = sc
        composite = composite + (spec["weight"] / wsum) * sc.fillna(0.0)
    df["composite"] = composite

    # 3. select top-k per group ------------------------------------------------
    sel = preset.get("select") or {}
    group, topk = sel.get("group"), sel.get("topk")
    if group and topk:
        if group not in df.columns:
            print(f"[score] WARNING: select group {group!r} absent; taking global top-{topk}")
            out = df.nlargest(topk, "composite")
        else:
            out = (df.sort_values("composite", ascending=False)
                     .groupby(group, sort=False).head(topk))
        print(f"[score] selected top-{topk} per {group!r}: {len(out)} rows")
    else:
        out = df
    return out.sort_values("composite", ascending=False).reset_index(drop=True)


def main() -> None:
    ap = argparse.ArgumentParser(description="config-driven composite scorer (dials in presets.py)")
    ap.add_argument("--preset", required=True, choices=sorted(PRESETS))
    ap.add_argument("--metrics-csv", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--topk", type=int, default=None, help="override preset select topk")
    ap.add_argument("--gate", type=float, default=None, help="override preset gate threshold")
    args = ap.parse_args()

    preset = {k: (v.copy() if isinstance(v, dict) else v) for k, v in PRESETS[args.preset].items()}
    if args.gate is not None:
        preset["gate"] = (preset["gate"][0], args.gate)
    if args.topk is not None:
        sel = dict(preset.get("select") or {})
        sel["topk"] = args.topk
        preset["select"] = sel

    df = pd.read_csv(args.metrics_csv, low_memory=False)
    out = score(df, preset)
    out.to_csv(args.out, index=False)
    print(f"[score] wrote {len(out)} rows -> {args.out}")


if __name__ == "__main__":
    main()
