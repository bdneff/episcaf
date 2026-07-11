#!/usr/bin/env python3
"""
plot_whole_epitope_103.py -- characterize the native-103 whole-epitope (C1) design set, in the same
style as plot_island_split.py. Produces two figures on our RFD3-103 run:

  metrics_103_island_split.png -- four-filter metric distributions split by island count (1 vs 2)
  metrics_104_vs_103.png       -- the same metrics, RFD3-104 (old C1 reproduction) vs RFD3-103 (ours),
                                  i.e. what the PepSeq 104->103 length constraint costs us

Inputs (local sibling data dir $D = ../.. from the repo, see filesystem-map):
  103 metrics : known_antigen/analysis/data/metrics_whole_epitope_103.csv  (our run; scp'd from Gemini)
  104 metrics : known_antigen/analysis/data/metrics_native_cyl_full.csv    (old C1 reproduction)
  island count: episcaf_v2/results/whole_epitope_designs.csv (epitope_chunks per id, joined in)

Run:
  python episcaf_analysis/viz/plot_whole_epitope_103.py \
      --m103 ../known_antigen/analysis/data/metrics_whole_epitope_103.csv \
      --m104 ../known_antigen/analysis/data/metrics_native_cyl_full.csv \
      --ledger results/whole_epitope_designs.csv --outdir manuscript/figures
"""
from __future__ import annotations
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams.update({"font.size": 14, "axes.titlesize": 17, "axes.labelsize": 16,
                     "xtick.labelsize": 14, "ytick.labelsize": 14, "legend.fontsize": 13,
                     "figure.titlesize": 18, "axes.linewidth": 1.1})

# (column, short title, x-unit, threshold, filter-direction, x-clip)
METRICS = [
    ("epitope_chunk_rmsd", "Epitope RMSD", "Å",        1.0, "≤", (0, 8)),
    ("overall_rmsd",       "Overall RMSD", "Å",        2.0, "≤", (0, 15)),
    ("mean_pae",           "Global PAE",   "PAE",      5.0, "<", (0, 30)),
    ("af3_n_clash_res",    "AF3 clashes",  "residues", 0.0, "=", (0, 40)),
]
FOUR = ["epitope_chunk_rmsd", "overall_rmsd", "mean_pae", "af3_n_clash_res"]


def four_filter(df):
    return ((df.epitope_chunk_rmsd <= 1) & (df.overall_rmsd <= 2)
            & (df.mean_pae < 5) & (df.af3_n_clash_res == 0))


def overlay(dplot, gcol, groups, title, out):
    fig, axes = plt.subplots(1, len(METRICS), figsize=(5 * len(METRICS), 4.6))
    for ax, (col, t, unit, thr, direction, clip) in zip(axes, METRICS):
        bins = np.linspace(clip[0], clip[1], 41)
        for val, name, color in groups:
            v = dplot.loc[dplot[gcol] == val, col].dropna().clip(*clip)
            ax.hist(v, bins=bins, density=True, histtype="step", lw=2.5, color=color, label=name)
            ax.axvline(v.median(), color=color, ls=":", lw=2.0, alpha=0.9)
        ax.axvline(thr, color="0.35", ls="--", lw=2.0)
        ax.annotate(f"filter {direction}{thr:g}", xy=(thr, 0.97), xycoords=("data", "axes fraction"),
                    fontsize=16, color="0.35", ha="left", va="top", rotation=90,
                    xytext=(3, 0), textcoords="offset points")
        ax.set_title(t); ax.set_xlabel(unit); ax.set_yticks([])
    axes[0].set_ylabel("density")
    axes[0].legend(frameon=False, loc="center right")
    fig.suptitle(title, y=0.99)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=150)
    print(f"wrote {out}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--m103", type=Path, required=True)
    ap.add_argument("--m104", type=Path, required=True)
    ap.add_argument("--ledger", type=Path, default=Path("results/whole_epitope_designs.csv"))
    ap.add_argument("--outdir", type=Path, default=Path("manuscript/figures"))
    args = ap.parse_args()

    d103 = pd.read_csv(args.m103, low_memory=False)
    chunks = pd.read_csv(args.ledger).drop_duplicates("id").set_index("id")["epitope_chunks"]
    d103["epitope_chunks"] = d103.id.map(chunks)
    valid103 = d103[d103.mean_pae.notna()].copy()          # distributions need an AF3 result

    # --- Figure 1: island split on the 103 set ---
    overlay(valid103, "epitope_chunks",
            [(1, "1 island", "#1f77b4"), (2, "2 islands", "#d62728")],
            "Native-103 metric distributions by epitope island count",
            args.outdir / "metrics_103_island_split.png")

    # --- Figure 2: 104 vs 103 (the PepSeq length constraint) ---
    d104 = pd.read_csv(args.m104, low_memory=False)
    ids56 = set(d103.id.unique())
    d104 = d104[d104.id.isin(ids56)].copy()                 # same 56 epitopes, fair comparison
    valid103["dataset"] = "103"
    d104["dataset"] = "104"
    both = pd.concat([d104[d104.mean_pae.notna()], valid103], ignore_index=True)
    overlay(both, "dataset",
            [("104", "104-mer (RFD3, Lawson contigs)", "#7f7f7f"),
             ("103", "103-mer (RFD3, our run)", "#2ca02c")],
            "RFD3 metric distributions: 104-mer vs 103-mer (PepSeq length constraint)",
            args.outdir / "metrics_104_vs_103.png")

    # --- printed numbers for the manuscript ---
    print("\n[103] four-filter pass rates:")
    print(f"  overall: {four_filter(d103).sum()}/{len(d103)} = {100*four_filter(d103).mean():.2f}%")
    for k in (1, 2):
        s = d103[d103.epitope_chunks == k]
        print(f"  {k}-island: {len(s):>6} designs  pass {100*four_filter(s).mean():.2f}%  "
              f"med epiRMSD {s.epitope_chunk_rmsd.median():.2f}  med PAE {s.mean_pae.median():.1f}")
    print(f"[104] four-filter pass (same 56 epitopes): "
          f"{four_filter(d104).sum()}/{len(d104)} = {100*four_filter(d104).mean():.2f}%")


if __name__ == "__main__":
    main()
