#!/usr/bin/env python3
"""Unit tests for the scorer transforms + selection (no data files required).

Run:  python tests/test_scoring.py     (or: pytest tests/)
"""
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "episcaf_analysis"))
from score import _transform_series, score  # noqa: E402


def test_percentile_orientation():
    x = pd.Series([1.0, 2.0, 3.0, 4.0])
    lo = _transform_series(x, dict(better="low", transform="percentile"))
    hi = _transform_series(x, dict(better="high", transform="percentile"))
    assert lo.iloc[0] > lo.iloc[-1]     # smallest is best when better=low
    assert hi.iloc[0] < hi.iloc[-1]     # largest is best when better=high


def test_minmax_bounds():
    x = pd.Series([0.0, 5.0, 10.0])
    z = _transform_series(x, dict(better="high", transform="minmax"))
    assert abs(z.min()) < 1e-9 and abs(z.max() - 1.0) < 1e-9


def test_sigmoid_monotonic_and_bounded():
    x = pd.Series(np.linspace(0, 10, 11))
    s = _transform_series(x, dict(better="low", transform="sigmoid", midpoint=5, k=1))
    assert (s.diff().dropna() < 0).all()   # better=low => score decreases as x grows
    assert 0.0 < s.min() and s.max() < 1.0


def test_end_to_end_gate_and_select():
    df = pd.DataFrame({
        "id":                    ["a", "a", "a", "b", "b"],
        "antigen":               ["x"] * 5,
        "epitope_chunk_rmsd":    [0.5, 1.0, 3.0, 0.8, 1.2],   # the 3.0 row is gated out
        "epitope_pae":           [2, 3, 9, 4, 5],
        "overall_rmsd":          [1, 2, 8, 2, 3],
        "cylinder_native_aware": [5, 8, 40, 6, 9],
    })
    preset = dict(
        gate=("epitope_chunk_rmsd", 2.5), scope="pooled", antigen_col="antigen",
        select=dict(group="id", topk=1),
        metrics={
            "epitope_chunk_rmsd": dict(weight=0.5, better="low", transform="percentile"),
            "epitope_pae":        dict(weight=0.5, better="low", transform="percentile"),
        },
    )
    out = score(df, preset)
    assert len(out) == 2                  # one winner per id (a and b)
    assert set(out["id"]) == {"a", "b"}
    assert "composite" in out.columns


def test_missing_metric_is_dropped_and_renormalized():
    df = pd.DataFrame({
        "id": ["a", "b"], "antigen": ["x", "x"],
        "epitope_chunk_rmsd": [1.0, 2.0], "epitope_pae": [3.0, 4.0],
    })
    preset = dict(
        gate=("epitope_chunk_rmsd", 2.5), scope="pooled", antigen_col="antigen",
        select=dict(group="id", topk=5),
        metrics={
            "epitope_chunk_rmsd": dict(weight=0.5, better="low", transform="percentile"),
            "epitope_pae":        dict(weight=0.25, better="low", transform="percentile"),
            "cylinder_native_aware": dict(weight=0.25, better="low", transform="percentile"),  # absent
        },
    )
    out = score(df, preset)               # should not raise; cylinder dropped, weights renormalized
    assert len(out) == 2 and "composite" in out.columns


if __name__ == "__main__":
    for name, fn in sorted(globals().items()):
        if name.startswith("test_") and callable(fn):
            fn()
            print(f"  ok: {name}")
    print("all scoring tests passed")
