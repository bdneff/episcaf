#!/usr/bin/env python3
"""
plot_c5_titration.py -- how well the C5 metric-space titration spans each sampled dimension.

C5 is a farthest-point (max-min) sample over four standardized scoring axes; its purpose is to SPREAD
designs across the metric space (including the tails the top-designs pile avoids) so binding read off
the spread can calibrate the scorer. This makes one panel per axis: the full whole-epitope pool (grey,
peaked at its mode) vs the C5 sample (blue) -- a good titration is flatter and reaches further into the
tails than the pool. Coverage % (sampled range / pool range) is annotated per panel.

Inputs (local sibling data dir):
  pool   : known_antigen/analysis/data/metrics_whole_epitope_103.csv  (the C1-103 pool C5 samples from)
  sample : results/dp4_C5_titration.csv                                (the 3,000 sampled designs)

Run:
  python episcaf_analysis/viz/plot_c5_titration.py \
      --pool ../known_antigen/analysis/data/metrics_whole_epitope_103.csv \
      --sample results/dp4_C5_titration.csv --out manuscript/figures/c5_titration_coverage.png
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
                     "xtick.labelsize": 13, "ytick.labelsize": 13, "legend.fontsize": 13,
                     "figure.titlesize": 18, "axes.linewidth": 1.1})

import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "scripts"))
from stage06_sample_c5 import fps, quotas, DROP_IDS   # noqa: E402  (reuse the exact FPS)

AXES = [
    ("epitope_chunk_rmsd",    "Epitope RMSD (Å)"),
    ("epitope_pae",           "Epitope PAE"),
    ("overall_rmsd",          "Overall RMSD (Å)"),
    ("af3_n_clash_res",       "AF3 clash (real)"),
    ("cylinder_native_aware", "Native-aware cylinder"),
]


def span_1d(pool, col, total, lo, hi, nbins=25, seed=0):
    """The 'uniform span of THIS axis alone' reference: a stratified-uniform sample of ~`total` designs,
    an equal quota drawn from each value-bin of [lo,hi], so it is FLAT by construction wherever the pool
    has designs (and only falls off where a bin runs dry -- the sparse extremes). This is what spanning
    one dimension uniformly looks like, against which the joint C5 sample (which tracks pool density on
    any single axis) is compared."""
    v = pool[col].dropna()
    edges = np.linspace(lo, hi, nbins + 1)
    per = max(1, total // nbins)
    out = []
    for i in range(nbins):
        m = v[(v >= edges[i]) & (v < edges[i + 1])]
        if len(m):
            out.append(m.sample(min(per, len(m)), random_state=seed))
    return pd.concat(out) if out else v.iloc[[]]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pool", type=Path, required=True)
    ap.add_argument("--sample", type=Path, required=True)
    ap.add_argument("--out", type=Path, default=Path("manuscript/figures/c5_titration_coverage.png"))
    args = ap.parse_args()

    pool = pd.read_csv(args.pool, low_memory=False)
    samp = pd.read_csv(args.sample, low_memory=False)

    fig, axes = plt.subplots(1, len(AXES), figsize=(4.4 * len(AXES), 4.6))
    for ax, (col, label) in zip(axes, AXES):
        p = pool[col].dropna()
        s = samp[col].dropna()
        hi = np.nanpercentile(p, 99)               # clip the far outlier tail so the spread is visible
        lo = float(min(p.min(), s.min()))
        r = span_1d(pool, col, len(samp), lo, hi)  # uniform-span-of-this-axis reference (~flat)
        bins = np.linspace(lo, hi, 41)
        ax.hist(p.clip(lo, hi), bins=bins, density=True, color="0.7", alpha=0.55, label="full pool")
        ax.hist(s.clip(lo, hi), bins=bins, density=True, histtype="step", lw=2.6, color="#1f77b4",
                label="C5 sample (joint)")
        ax.hist(r.clip(lo, hi), bins=bins, density=True, histtype="step", lw=2.4, color="#d62728",
                label="span this axis alone")
        cover = 100 * (s.max() - s.min()) / (p.max() - p.min())
        ax.set_title(label)
        ax.set_yticks([])
        ax.annotate(f"spans {cover:.0f}%", xy=(0.96, 0.96), xycoords="axes fraction",
                    ha="right", va="top", fontsize=13, color="#1f77b4")
    axes[0].set_ylabel("density")
    axes[3].legend(frameon=False, loc="center right", fontsize=11)
    fig.suptitle("C5 metric-space titration: sample coverage of each sampled axis "
                 f"(n={len(samp):,} vs pool {len(pool):,})", y=0.99)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    args.out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.out, dpi=150)
    print(f"wrote {args.out}")
    for col, label in AXES:
        p, s = pool[col].dropna(), samp[col].dropna()
        print(f"  {label:28s} pool[{p.min():.2f},{p.max():.2f}] sample[{s.min():.2f},{s.max():.2f}]  "
              f"spans {100*(s.max()-s.min())/(p.max()-p.min()):.0f}%")


if __name__ == "__main__":
    main()
