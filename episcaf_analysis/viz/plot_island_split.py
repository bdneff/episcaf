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

# (column, label, threshold, filter-direction, x-clip)
METRICS = [
    ("epitope_chunk_rmsd_vs_mpnn", "Epitope chunk RMSD (Å)", 1.0, "<=", (0, 8)),
    ("overall_rmsd",               "Overall RMSD (Å)",       2.0, "<=", (0, 15)),
    ("mean_pae",                   "Mean PAE, global",       5.0, "<",  (0, 30)),
    ("n_clash",                    "AF3 clashing residues",  0.0, "==", (0, 40)),
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

    fig, axes = plt.subplots(1, len(METRICS), figsize=(4 * len(METRICS), 3.6))
    for ax, (col, label, thr, direction, clip) in zip(axes, METRICS):
        bins = np.linspace(clip[0], clip[1], 41)
        for k, name, color in GROUPS:
            v = df.loc[df.epitope_chunks == k, col].dropna().clip(*clip)
            n = len(v)
            ax.hist(v, bins=bins, density=True, histtype="step", lw=2, color=color,
                    label=f"{name} (n={n:,})")
            ax.axvline(v.median(), color=color, ls=":", lw=1.2, alpha=0.8)
        ax.axvline(thr, color="0.3", ls="--", lw=1.2)
        ax.set_xlabel(label)
        ax.set_yticks([])
        ax.set_title(f"filter: {direction} {thr:g}", fontsize=9, color="0.3")
    axes[0].set_ylabel("density")
    axes[0].legend(fontsize=8, frameon=False)

    pr = {k: 100 * df.loc[df.epitope_chunks == k, "_pass"].mean() for k, _, _ in GROUPS}
    fig.suptitle("DP3 (RFD1+MPNN, dp2) metric distributions by epitope island count  "
                 f"— four-filter pass: 1-island {pr[1]:.1f}%, 2-island {pr[2]:.1f}%  "
                 "(density-normalized; dotted = median, dashed = filter)",
                 fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
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
