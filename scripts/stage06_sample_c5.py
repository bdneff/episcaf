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
      --metrics-csv known_antigen/analysis/data/metrics_native_cyl_full.csv \
      --total 3000 --out results/dp4_C5_titration.csv
"""
from __future__ import annotations
import argparse
from pathlib import Path
import numpy as np
import pandas as pd

BASE_AXES = ["epitope_chunk_rmsd", "epitope_pae", "overall_rmsd"]
# Accessibility is SPLIT: half the sample spans the REAL AF3 clash (af3_n_clash_res, available because
# the antibody is known), half spans the native-aware cylinder surrogate, so the titration calibrates
# BOTH accessibility measures against binding (Brandon, 2026-07-13). Each half runs its own FPS over
# BASE_AXES + its accessibility axis; the two selections are merged and de-duplicated.
ACCESS_AXES = ["af3_n_clash_res", "cylinder_native_aware"]
# Canonical known-Ab exclusion set -> 56 mAbs: 4xwo (low yield), 7a3t (4-residue epitope),
# 2h32 (not a standard antibody case). See docs/DP4_LIBRARY.md.
DROP_IDS = {"4xwo_5p", "7a3t_0p", "2h32_0p"}


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


def sample_half(df, access, total, group):
    """Per-mAb FPS over BASE_AXES + one accessibility axis. Returns {original df index -> fps_order}."""
    axes = BASE_AXES + [access]
    d = df.dropna(subset=axes)
    oi = d.index.to_numpy()
    gv = d[group].to_numpy()
    P = np.column_stack([d[a].rank(pct=True).to_numpy() for a in axes])   # percentile-standardized
    q = quotas(sorted(pd.unique(gv)), total)
    out = {}
    for g in q:
        loc = np.where(gv == g)[0]
        for order, li in enumerate(fps(P[loc], min(q[g], len(loc)))):
            out[int(oi[loc[li]])] = order
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--metrics-csv", required=True)
    ap.add_argument("--group", default="id")
    ap.add_argument("--total", type=int, default=3000)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    df = pd.read_csv(args.metrics_csv, low_memory=False)
    df = df[~df[args.group].astype(str).str.lower().isin(DROP_IDS)].reset_index(drop=True)
    for a in BASE_AXES + ACCESS_AXES:
        df[a] = pd.to_numeric(df[a], errors="coerce")
    print(f"[c5] pool {len(df)} designs; mAbs = {df[args.group].nunique()}")

    # split the quota: half spanning the real AF3 clash, half the native-aware cylinder. The cylinder
    # half is drawn from the pool with the af3-half removed, so the two halves are DISJOINT and sum to
    # exactly --total (a clean "half to each" with no design counted twice).
    h1 = args.total // 2
    sel_af3 = sample_half(df, "af3_n_clash_res", h1, args.group)
    sel_cyl = sample_half(df.drop(index=list(sel_af3)), "cylinder_native_aware", args.total - h1, args.group)
    tag = {}
    for i in sel_af3:
        tag.setdefault(i, []).append("af3_clash")
    for i in sel_cyl:
        tag.setdefault(i, []).append("cylinder")
    rows = sorted(tag)
    out = df.loc[rows].copy()
    out["access_sampled"] = ["+".join(tag[i]) for i in rows]
    out["category"] = "metricSpaceTitration"

    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.out, index=False)
    n_both = sum(1 for i in rows if len(tag[i]) == 2)
    print(f"[c5] sampled {len(sel_af3)} (af3-clash) + {len(sel_cyl)} (cylinder) "
          f"-> {len(out)} unique ({n_both} picked by both halves)")
    print("[c5] axis coverage (sampled range / full-pool range):")
    for a in BASE_AXES + ACCESS_AXES:
        fp, sp = df[a].dropna(), out[a].dropna()
        cov = (sp.max() - sp.min()) / (fp.max() - fp.min() + 1e-9)
        print(f"       {a:24s} pool[{fp.min():.2f},{fp.max():.2f}] "
              f"sample[{sp.min():.2f},{sp.max():.2f}]  cover {cov:.0%}")
    print(f"[c5] wrote {args.out}")


if __name__ == "__main__":
    main()
