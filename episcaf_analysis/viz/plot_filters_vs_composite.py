#!/usr/bin/env python3
"""
plot_filters_vs_composite.py

2x2 scatter grid: each of the four validation filters vs the composite score,
from the --scored_csv apply_composite_filter.py writes. Four-filter passers are
highlighted and each filter's threshold is drawn as a dashed line, so you can see
how the composite tracks each underlying metric and where the passers sit.

A good composite shows passers packed into the high-composite / passing-region
corner of every panel; a filter whose panel shows passers spread across the
composite axis is one the composite is NOT capturing well.

Example
-------
    python scripts/plot_filters_vs_composite.py \
        --scored_csv runs/run_rfd3_mpnn/04_filter/scored.csv \
        --out_png    runs/run_rfd3_mpnn/04_filter/filters_vs_composite.png
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# (column, threshold, passing side, label)
FILTERS = [
    ("overall_rmsd",        2.0, "below", "overall_rmsd \u2264 2"),
    ("epitope_chunk_rmsd",  1.0, "below", "epitope_chunk_rmsd \u2264 1"),
    ("mean_pae",            5.0, "below", "mean_pae < 5"),
    ("af3_n_clash_res",     0.0, "below", "af3_n_clash_res = 0"),
]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--scored_csv", required=True)
    ap.add_argument("--out_png", required=True)
    ap.add_argument("--dedupe_col", default="token",
                    help="collapse to best-composite row per design (default: token)")
    ap.add_argument("--clip_pct", type=float, default=99.0,
                    help="clip each y-axis at this percentile so outliers don't "
                         "squash the bulk (default: 99)")
    args = ap.parse_args()

    df = pd.read_csv(args.scored_csv, low_memory=False)
    if "composite" not in df.columns:
        raise SystemExit("scored_csv has no 'composite' column")
    df["composite"] = pd.to_numeric(df["composite"], errors="coerce")
    if args.dedupe_col in df.columns:
        before = len(df)
        df = (df.sort_values("composite", ascending=False)
                .drop_duplicates(args.dedupe_col, keep="first").reset_index(drop=True))
        print(f"  deduped {before:,} -> {len(df):,} designs by {args.dedupe_col!r}")

    passer = (df["is_pass"].astype(bool) if "is_pass" in df.columns
              else pd.Series(False, index=df.index))
    comp = df["composite"]

    fig, axes = plt.subplots(2, 2, figsize=(12, 9))
    for ax, (col, thr, side, label) in zip(axes.ravel(), FILTERS):
        if col not in df.columns:
            ax.text(0.5, 0.5, f"no column {col!r}", ha="center")
            continue
        y = pd.to_numeric(df[col], errors="coerce")
        ax.scatter(comp[~passer], y[~passer], s=5, alpha=0.2,
                   color="#9bb8d3", label="fail", rasterized=True)
        if passer.any():
            ax.scatter(comp[passer], y[passer], s=18, alpha=0.9, color="#c0392b",
                       edgecolor="k", linewidth=0.2, zorder=3, label="pass")
        ax.axhline(thr, color="k", ls="--", lw=1.2, alpha=0.7)
        # clip y so a long tail (esp. clash count) doesn't flatten the plot
        hi = np.nanpercentile(y, args.clip_pct)
        ax.set_ylim(min(0, np.nanmin(y)) - 0.02 * abs(hi),
                    max(hi, thr * 1.5 if thr else hi) * 1.05 + 1e-9)
        ax.set_xlabel("composite score")
        ax.set_ylabel(col)
        ax.set_title(label)
        ax.legend(fontsize=8, loc="upper right")

    fig.suptitle("Filters vs composite score (passers in red, threshold dashed)",
                 y=1.00, fontsize=12)
    fig.tight_layout()
    out = Path(args.out_png)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=150)
    print(f"Wrote {out}")

    # quick correlation of composite vs each filter, on passers and overall
    print("\nSpearman corr of composite vs filter (all / passers):")
    for col, *_ in FILTERS:
        if col not in df.columns:
            continue
        y = pd.to_numeric(df[col], errors="coerce")
        m = comp.notna() & y.notna()
        r_all = comp[m].corr(y[m], method="spearman")
        mp = m & passer
        r_p = (comp[mp].corr(y[mp], method="spearman")
               if mp.sum() > 5 else float("nan"))
        print(f"  {col:22s}  all={r_all:+.3f}   passers={r_p:+.3f}")


if __name__ == "__main__":
    main()
