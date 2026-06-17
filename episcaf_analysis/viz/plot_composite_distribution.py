#!/usr/bin/env python3
"""
plot_composite_distribution.py

Diagnostics for the composite score, from the --scored_csv that
apply_composite_filter.py writes (needs columns: composite, is_pass, and the
epitope group column, default 'id').

Panels:
  (1) histogram of composite over all designs, four-filter passers overlaid;
  (2) scatter of composite vs a chosen metric, colored by pass/fail;
  (3) per-epitope capture curve: as you raise the per-epitope top-k, how
      coverage (targets keeping >=1 passer) and capacity-adjusted per-epitope
      recall scale -- the curve that tells you whether k=15 is the right cut.

Example
-------
    python scripts/plot_composite_distribution.py \
        --scored_csv runs/run_rfd3_mpnn/04_filter/scored.csv \
        --out_png    runs/run_rfd3_mpnn/04_filter/composite_distribution.png
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

Y_METRIC_ALIASES = ["mean_pae", "epitope_pae", "epitope_chunk_rmsd"]
KS = [1, 3, 5, 10, 15, 20, 30, 50, 100]


def first_present(cols, candidates):
    for c in candidates:
        if c in cols:
            return c
    return None


def per_epitope_curve(df, group_col, comp_col, passer, ks):
    """For each k: coverage (passer-targets keeping >=1 passer in their top-k)
    and capacity-adjusted per-epitope recall (mean over passer-targets of
    kept / min(n_passers, k))."""
    d = df.copy()
    d["_rank"] = d.groupby(group_col)[comp_col].rank(ascending=False, method="first")
    d["_pass"] = passer.values
    pg = d[d["_pass"]]
    passer_groups = pg[group_col].astype(str).unique()
    n_pg = len(passer_groups)
    out = []
    grouped = list(pg.groupby(pg[group_col].astype(str)))
    for k in ks:
        covered = 0
        per_ep = []
        for _, grp in grouped:
            kept = int((grp["_rank"] <= k).sum())
            if kept > 0:
                covered += 1
            per_ep.append(kept / min(len(grp), k))
        out.append((k,
                    covered / n_pg if n_pg else np.nan,
                    float(np.mean(per_ep)) if per_ep else np.nan))
    return out, n_pg


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--scored_csv", required=True)
    ap.add_argument("--out_png", required=True)
    ap.add_argument("--group_col", default="id")
    ap.add_argument("--y_metric", default=None)
    ap.add_argument("--top_k", type=int, default=15)
    args = ap.parse_args()

    df = pd.read_csv(args.scored_csv, low_memory=False)
    if "composite" not in df.columns:
        raise SystemExit("scored_csv has no 'composite' column")
    comp = pd.to_numeric(df["composite"], errors="coerce")
    has_pass = "is_pass" in df.columns and df["is_pass"].sum() > 0
    passer = (df["is_pass"].astype(bool) if "is_pass" in df.columns
              else pd.Series(False, index=df.index))
    y_col = args.y_metric or first_present(df.columns, Y_METRIC_ALIASES)

    fig, axes = plt.subplots(1, 3, figsize=(17, 4.6))

    # panel 1: distribution
    ax = axes[0]
    ax.hist(comp.dropna(), bins=60, color="#9bb8d3", alpha=0.85,
            label=f"all designs (n={comp.notna().sum():,})")
    if has_pass:
        pc = comp[passer].dropna()
        ax.hist(pc, bins=60, color="#c0392b", alpha=0.65,
                label=f"passers (n={len(pc):,})")
        ax.axvline(pc.median(), color="#c0392b", ls="--", lw=1.5,
                   label=f"passer median = {pc.median():.3f}")
    ax.set_xlabel("composite score"); ax.set_ylabel("count")
    ax.set_title("Composite score distribution"); ax.legend(fontsize=8)

    # panel 2: scatter
    ax = axes[1]
    if y_col and y_col in df.columns:
        y = pd.to_numeric(df[y_col], errors="coerce")
        ax.scatter(comp[~passer], y[~passer], s=4, alpha=0.22,
                   color="#9bb8d3", label="fail", rasterized=True)
        if has_pass:
            ax.scatter(comp[passer], y[passer], s=14, alpha=0.9, color="#c0392b",
                       edgecolor="k", linewidth=0.2, label="pass")
        ax.set_xlabel("composite score"); ax.set_ylabel(y_col)
        ax.set_title(f"composite vs {y_col}"); ax.legend(fontsize=8)
    else:
        ax.text(0.5, 0.5, "no y-metric column", ha="center")

    # panel 3: per-epitope capture curve
    ax = axes[2]
    if has_pass and args.group_col in df.columns:
        curve, n_pg = per_epitope_curve(df, args.group_col, "composite", passer, KS)
        ks = [c[0] for c in curve]
        cov = [c[1] for c in curve]
        rec = [c[2] for c in curve]
        ax.plot(ks, cov, "-o", color="#2c7fb8", label="coverage (targets w/ a passer)")
        ax.plot(ks, rec, "-s", color="#c0392b",
                label="per-epitope recall (capacity-adj.)")
        ax.axvline(args.top_k, color="gray", ls="--", lw=1,
                   label=f"current top_k = {args.top_k}")
        ax.set_xscale("log"); ax.set_xlabel("per-epitope top-k")
        ax.set_ylabel("fraction"); ax.set_ylim(0, 1.02)
        ax.set_title(f"Per-epitope capture ({n_pg} passer-targets)")
        ax.legend(fontsize=8)
        print(f"\nPer-epitope capture ({n_pg} passer-targets):")
        print("  k      coverage   per-ep recall")
        for k, c, r in curve:
            print(f"  {k:4d}    {c:.3f}      {r:.3f}")
    else:
        ax.text(0.5, 0.5, "need is_pass + group_col", ha="center")

    fig.tight_layout()
    out = Path(args.out_png); out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=150)
    print(f"\nWrote {out}")

    if has_pass:
        order = comp.rank(pct=True)
        n_pass = int(passer.sum())
        print(f"\nGlobal capture (top X% by composite), {n_pass} passers:")
        for frac in (0.01, 0.02, 0.05, 0.10, 0.20):
            captured = int(((order >= 1.0 - frac) & passer).sum())
            print(f"  top {frac*100:4.0f}%  ->  {captured}/{n_pass} ({captured/n_pass:.0%})")


if __name__ == "__main__":
    main()
