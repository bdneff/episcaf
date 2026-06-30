#!/usr/bin/env python3
"""Do the filter metrics correlate with binding? DP3 cognate scaffolded designs.

Binding readout = cognate_log_enrichment (log10(1+Ab) - log10(1+NoAb)) from
results/dp3_binding_metrics.csv. We look at the three metrics that have variance in the
assayed (all-passing) set: overall_rmsd, epitope_chunk_rmsd_vs_mpnn, mean_pae
(af3_n_clash_res is all-zero here, so it is dropped).

Three views per metric, because pooling across antibodies is confounded by each mAb's own
baseline:
  (row 1) POOLED / marginalized over all designs   -- one cloud, one trend line
  (row 2) colored by antibody, per-Ab trend lines  -- shows the between-Ab structure
  + a WITHIN-Ab (fixed-effects) correlation: both variables group-mean-centered per
    antibody, which removes the baseline offsets and asks the marginalized question without
    the confound.

Writes the figure to the scratchpad (exploratory) and prints all correlations.
Run: python scripts/plot_dp3_metric_binding.py [out.png]
"""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
CSV = ROOT / "results" / "dp3_binding_metrics.csv"
CYL = ROOT / "results" / "assayed_native_cyl.csv"   # post-hoc cylinder on the assayed structures
OUT = Path(sys.argv[1]) if len(sys.argv) > 1 else ROOT / "manuscript" / "figures" / "dp3_metric_binding.png"

METRICS = [("overall_rmsd", "overall RMSD (Å)"),
           ("epitope_chunk_rmsd_vs_mpnn", "epitope RMSD (Å)"),
           ("mean_pae", "mean PAE"),
           ("cylinder_ca_clashes", "cylinder (plain)"),
           ("cylinder_native_aware", "cylinder (native-aware)")]
Y = "cognate_log_enrichment"


def corr_line(ax, x, y, color, label=None, lw=1.6):
    m = x.notna() & y.notna()
    x, y = x[m], y[m]
    if len(x) < 3:
        return None
    r, p = stats.pearsonr(x, y)
    rho, _ = stats.spearmanr(x, y)
    xs = np.linspace(x.min(), x.max(), 50)
    b1, b0 = np.polyfit(x, y, 1)
    ax.plot(xs, b0 + b1 * xs, color=color, lw=lw, label=label)
    return r, rho, p, len(x)


def main() -> None:
    d = pd.read_csv(CSV)
    cyl = pd.read_csv(CYL)[["designedSequence", "cylinder_ca_clashes", "cylinder_native_aware"]]
    d = d.merge(cyl, on="designedSequence", how="left")
    s = d[(d.category == "scaffoldedAbEpitope") & d.cognate_ab.notna()].copy()
    print(f"cognate scaffolded designs: {len(s)}  (antibodies: {sorted(s.cognate_ab.unique())})")
    print(f"cylinder merged for {s.cylinder_native_aware.notna().sum()}/{len(s)} of them")

    # within-Ab (fixed-effect) centering: subtract per-Ab mean from metric and from y
    g = s.groupby("cognate_ab")
    s["y_c"] = s[Y] - g[Y].transform("mean")
    for col, _ in METRICS:
        s[col + "_c"] = s[col] - g[col].transform("mean")

    abs_sorted = s.cognate_ab.value_counts().index.tolist()
    cmap = plt.get_cmap("tab10")
    abcol = {ab: cmap(i % 10) for i, ab in enumerate(abs_sorted)}

    ncol = len(METRICS)
    fig, axes = plt.subplots(2, ncol, figsize=(5 * ncol, 9))
    print("\n=== correlations of metric vs cognate_log_enrichment ===")
    for j, (col, lab) in enumerate(METRICS):
        # row 1: pooled
        ax = axes[0, j]
        ax.scatter(s[col], s[Y], s=14, c="0.55", linewidths=0)
        res = corr_line(ax, s[col], s[Y], "black", lw=2.0)
        r, rho, p, n = res
        # within-Ab fixed-effect correlation
        rw, pw = stats.pearsonr(s[col + "_c"], s["y_c"])
        rhow, _ = stats.spearmanr(s[col + "_c"], s["y_c"])
        ax.set_title(f"POOLED  Pearson r={r:+.2f} (p={p:.1g})\nSpearman rho={rho:+.2f}",
                     fontsize=10)
        ax.set_xlabel(lab); ax.set_ylabel("cognate log-enrichment")
        print(f"\n{col}:")
        print(f"  POOLED      n={n:3d}  Pearson r={r:+.3f} (p={p:.2g})  Spearman rho={rho:+.3f}")
        print(f"  WITHIN-Ab   (fixed-effect, centered)  Pearson r={rw:+.3f} (p={pw:.2g})  Spearman rho={rhow:+.3f}")

        # row 2: colored by antibody + per-Ab lines for n>=12
        ax2 = axes[1, j]
        for ab in abs_sorted:
            gg = s[s.cognate_ab == ab]
            ax2.scatter(gg[col], gg[Y], s=16, color=abcol[ab], linewidths=0,
                        label=f"{ab} (n={len(gg)})")
            if len(gg) >= 12:
                rr = corr_line(ax2, gg[col], gg[Y], abcol[ab], lw=1.6)
                if rr:
                    print(f"  within {ab:5s} n={rr[3]:3d}  Pearson r={rr[0]:+.3f}  Spearman rho={rr[1]:+.3f}")
        ax2.set_title(f"by antibody  (within-Ab fixed-effect r={rw:+.2f})", fontsize=10)
        ax2.set_xlabel(lab); ax2.set_ylabel("cognate log-enrichment")
        if j == ncol - 1:
            ax2.legend(fontsize=7, loc="best", framealpha=0.9)

    fig.suptitle("DP3: do the filter metrics predict binding? (cognate scaffolded designs)",
                 fontsize=13)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, dpi=140)
    print(f"\nwrote {OUT}")


if __name__ == "__main__":
    main()
