#!/usr/bin/env python3
"""
stage06_sample_c5.py -- DP4 component C5: metric-space titration set.

Instead of shipping the TOP designs per mAb (that's C1), C5 deliberately samples designs SPREAD
across the metric space -- including ones the filters would reject -- so binding read off this
spread lets us learn which metrics predict enrichment and fit the composite weights (manuscript
sec:dp4c5, Q2). Method: farthest-point sampling (FPS, max-min) over the standardized 4 scoring
axes, PER mAb. Accessibility axis is the native-aware CYLINDER (not the real clash) on purpose --
so what we learn transfers to the no-antibody setting.

Standardization: each axis -> its percentile rank over the kept pool (scale-free, bounded [0,1],
same transform the scorer uses). FPS is deterministic (seed = the design nearest each mAb's
centroid, then greedily add the point maximizing min-distance to the selected set) -> reproducible.

Pool = the known-Ab whole-epitope designs (same table as C1), over the 57 assayable mAbs
(drop 4xwo_5P low-yield, 7a3t_0P epitope-too-small). Per-mAb quotas sum to --total (even split).

Usage:
  python scripts/stage06_sample_c5.py \
      --metrics-csv known_antigen/analysis/full_run/metrics_native_cyl_full.csv \
      --total 3000 --out results/dp4_C5_titration.csv
"""
from __future__ import annotations
import argparse
from pathlib import Path
import numpy as np
import pandas as pd

AXES = ["epitope_chunk_rmsd", "epitope_pae", "overall_rmsd", "cylinder_native_aware"]
DROP_IDS = {"4xwo_5p", "7a3t_0p"}   # assay-dropped: low yield / epitope too small


def fps(X: np.ndarray, k: int) -> list[int]:
    """Farthest-point (max-min) sampling of k rows of X. Deterministic: seed = nearest-centroid."""
    n = len(X)
    if k >= n:
        return list(range(n))
    c = X.mean(0)
    start = int(np.argmin(((X - c) ** 2).sum(1)))
    sel = [start]
    d = ((X - X[start]) ** 2).sum(1)          # min sq-dist to the selected set
    while len(sel) < k:
        i = int(np.argmax(d))
        sel.append(i)
        d = np.minimum(d, ((X - X[i]) ** 2).sum(1))
    return sel


def quotas(ids: list, total: int) -> dict:
    """Even per-group allocation summing to exactly `total`."""
    n = len(ids)
    base, rem = divmod(total, n)
    return {g: base + (1 if i < rem else 0) for i, g in enumerate(ids)}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--metrics-csv", required=True)
    ap.add_argument("--group", default="id")
    ap.add_argument("--total", type=int, default=3000)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    df = pd.read_csv(args.metrics_csv, low_memory=False)
    df = df[~df[args.group].astype(str).str.lower().isin(DROP_IDS)].copy()
    for a in AXES:
        df[a] = pd.to_numeric(df[a], errors="coerce")
    n0 = len(df)
    df = df.dropna(subset=AXES).reset_index(drop=True)
    print(f"[c5] pool {n0} -> {len(df)} after dropping NaN axes; mAbs = {df[args.group].nunique()}")

    # percentile-standardize each axis over the whole kept pool
    P = np.column_stack([df[a].rank(pct=True).to_numpy() for a in AXES])

    groups = sorted(df[args.group].unique())
    q = quotas(groups, args.total)
    picks = []
    for g in groups:
        idx = df.index[df[args.group] == g].to_numpy()
        k = min(q[g], len(idx))
        sel_local = fps(P[idx], k)
        for order, li in enumerate(sel_local):
            picks.append((idx[li], g, order))
    sel = pd.DataFrame(picks, columns=["row", args.group, "fps_order"]).set_index("row")
    out = df.loc[sel.index].copy()
    out["fps_order"] = sel["fps_order"].to_numpy()
    out["category"] = "metricSpaceTitration"

    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.out, index=False)
    print(f"[c5] sampled {len(out)} designs over {len(groups)} mAbs "
          f"(per-mAb {min(q.values())}-{max(q.values())})")
    # spread check: does the sample cover the pool's range on each axis?
    print("[c5] axis coverage (sampled range / full-pool range):")
    for a in AXES:
        fp = df[a]; sp = out[a]
        cov = (sp.max() - sp.min()) / (fp.max() - fp.min() + 1e-9)
        print(f"       {a:24s} pool[{fp.min():.2f},{fp.max():.2f}] "
              f"sample[{sp.min():.2f},{sp.max():.2f}]  cover {cov:.0%}")
    print(f"[c5] wrote {args.out}")


if __name__ == "__main__":
    main()
