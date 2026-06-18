#!/usr/bin/env python3
"""
plot_island_split.py  --  DP3 metric distributions split by epitope island count (1 vs 2).

Answers John's question: do single-island ("category #1") epitopes scaffold differently from
two-island ones? Splits the DP3 design set by epitope_chunks (1 vs 2 islands) and overlays the
per-metric distributions for the four-filter metrics.

Source: dp2.parquet (Lawson's RFD1+MPNN deposited set). Designs without an AF3 result
(mean_pae NaN) are excluded. N is very uneven (549 single-island vs ~124k two-island with AF3
results), so distributions are density-normalized.

Regenerate:
  python episcaf_analysis/viz/plot_island_split.py \
      --dp2 ../known_antigen/analysis/full_run/dp2.parquet \
      --out manuscript/figures/island_split_metrics.png
"""
from __future__ import annotations
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams.update({
    "font.size": 14,
    "axes.titlesize": 17,
    "axes.labelsize": 16,
    "xtick.labelsize": 14,
    "ytick.labelsize": 14,
    "legend.fontsize": 14,
    "figure.titlesize": 18,
    "axes.linewidth": 1.1,
})

# (column, short title, x-axis unit, threshold, filter-direction, x-clip)
METRICS = [
    ("epitope_chunk_rmsd_vs_mpnn", "Epitope RMSD", "Å",        1.0, "≤", (0, 8)),
    ("overall_rmsd",               "Overall RMSD", "Å",        2.0, "≤", (0, 15)),
    ("mean_pae",                   "Global PAE",   "PAE",      5.0, "<", (0, 30)),
    ("n_clash",                    "AF3 clashes",  "residues", 0.0, "=", (0, 40)),
]
GROUPS = [(1, "1 island", "#1f77b4"), (2, "2 islands", "#d62728")]


def clashlen(x):
    try:
        return len(x)
    except TypeError:
        return 0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dp2", type=Path,
                    default=Path("../known_antigen/analysis/full_run/dp2.parquet"))
    ap.add_argument("--out", type=Path,
                    default=Path("manuscript/figures/island_split_metrics.png"))
    args = ap.parse_args()

    df = pd.read_parquet(args.dp2)
    df["n_clash"] = df["af3_clash_resindices"].apply(clashlen)
    df = df[df["mean_pae"].notna()].copy()          # AF3 result required

    # four-filter pass per design, for the per-group pass rate annotation
    df["_pass"] = ((df.epitope_chunk_rmsd_vs_mpnn <= 1) & (df.overall_rmsd <= 2)
                   & (df.mean_pae < 5) & (df.n_clash == 0))

    fig, axes = plt.subplots(1, len(METRICS), figsize=(5 * len(METRICS), 4.6))
    for ax, (col, title, unit, thr, direction, clip) in zip(axes, METRICS):
        bins = np.linspace(clip[0], clip[1], 41)
        for k, name, color in GROUPS:
            v = df.loc[df.epitope_chunks == k, col].dropna().clip(*clip)
            ax.hist(v, bins=bins, density=True, histtype="step", lw=2.5, color=color,
                    label=name)
            ax.axvline(v.median(), color=color, ls=":", lw=2.0, alpha=0.9)
        ax.axvline(thr, color="0.35", ls="--", lw=2.0)
        # short threshold tag near the filter line, off the canvas detail goes in the caption
        ax.annotate(f"filter {direction}{thr:g}", xy=(thr, 0.97), xycoords=("data", "axes fraction"),
                    fontsize=12, color="0.35", ha="left", va="top", rotation=90,
                    xytext=(3, 0), textcoords="offset points")
        ax.set_title(title)
        ax.set_xlabel(unit)
        ax.set_yticks([])
    axes[0].set_ylabel("density")
    axes[0].legend(frameon=False, loc="center right")
    fig.suptitle("DP3 metric distributions by epitope island count", y=0.99)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    args.out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.out, dpi=150)
    print(f"wrote {args.out}")
    for k, name, _ in GROUPS:
        sub = df[df.epitope_chunks == k]
        print(f"  {name}: {len(sub):,} designs, pass {100*sub._pass.mean():.2f}%, "
              f"median epi-RMSD {sub.epitope_chunk_rmsd_vs_mpnn.median():.2f}, "
              f"median PAE {sub.mean_pae.median():.1f}")


if __name__ == "__main__":
    main()
