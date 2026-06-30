#!/usr/bin/env python3
"""Descriptive analysis of the DP3 assay binding data + the canonical scatter figure.

Reads results/dp3_binding_metrics.csv (built by build_dp3_binding_join.py) and prints every
number the manuscript "Experimental binding data" section cites, then renders
manuscript/figures/dp3_binding_scatter.png: per-antibody log10(1+NoAb) vs log10(1+Ab), all
library members in grey, the cognate scaffolded designs highlighted (John's red-highlight view).

Run:
    python scripts/analyze_dp3_binding.py
"""
from __future__ import annotations
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams.update({"font.size": 14, "axes.titlesize": 16, "axes.labelsize": 15, "xtick.labelsize": 12, "ytick.labelsize": 12, "legend.fontsize": 12, "figure.titlesize": 18})  # paper-legible fonts

ROOT = Path(__file__).resolve().parents[1]
CSV = ROOT / "results" / "dp3_binding_metrics.csv"
FIG = ROOT / "manuscript" / "figures" / "dp3_binding_scatter.png"
METRICS = ["overall_rmsd", "epitope_chunk_rmsd_vs_mpnn", "mean_pae", "af3_n_clash_res"]

# antibody -> (Ab intensity col, that run's NoAb col, run label). Mirrors the assay CSVs.
RUNS = {
    "6o9i": ("X6o9i_1X_0.1pmol_beadsFirst", "NoAb_0.1pmol_beadsFirst", "IM226"),
    "8cz8": ("X8cz8_1X_0.1pmol_beadsFirst", "NoAb_0.1pmol_beadsFirst", "IM226"),
    "8jnk": ("X8jnk_1X_0.1pmol_beadsFirst", "NoAb_0.1pmol_beadsFirst", "IM226"),
    "6xxv": ("X6xxv_1X_0.1pmol_beadsFirst", "NoAb_0.1pmol_beadsFirst", "IM226"),
    "5fhx": ("X5fhx_1X_0.1pmol_beadsFirst", "NoAb_0.1pmol_beadsFirst", "IM226"),
    "7ox3": ("X7ox3_1X_0.1pmol_beadsFirst", "NoAb_0.1pmol_beadsFirst", "IM226"),
    "8db4": ("X8db4_10X_1pmol_beadsFirst", "NoAb_0.1pmol_BeadsLast_B", "IM229"),
    "8pww": ("X8pww_1X_1pmol_beadsFirst",  "NoAb_0.1pmol_BeadsLast_B", "IM229"),
}
log = lambda x: np.log10(1.0 + x)


def main() -> None:
    d = pd.read_csv(CSV)
    print(f"library members: {len(d)}   unique designedSequence: {d.designedSequence.nunique()}")
    print("\ncategory counts:")
    print(d.category.value_counts().to_string())

    scaf = d[d.category == "scaffoldedAbEpitope"]
    cog = scaf[scaf.cognate_ab.notna()]
    print(f"\nscaffolded designs: {len(scaf)};  with a usable cognate Ab: {len(cog)}")
    dropped = scaf[scaf.cognate_ab.isna()].Target.value_counts()
    print(f"no-Ab (dropped 4xwo/7a3t) designs: {len(scaf)-len(cog)}  -> {dropped.to_dict()}")
    print(f"all assayed scaffolds passed the 4-filter: is_pass all True = {cog.is_pass.all()}")

    print("\nmetric range within the assayed (all-passing) set:")
    for c in METRICS:
        print(f"  {c:28s} min={cog[c].min():.3f} median={cog[c].median():.3f} max={cog[c].max():.3f}")

    print("\nper-antibody (cognate scaffolded designs):")
    rows = []
    for ab, g in cog.groupby("cognate_ab"):
        med = g.cognate_log_enrichment.median()
        frac_up = (g.cognate_log_enrichment > 0).mean()
        rows.append((ab, RUNS[ab][2], len(g), med, frac_up))
    rows.sort(key=lambda r: -r[2])
    print(f"  {'Ab':6s} {'run':6s} {'n':>4s} {'med_enrich':>11s} {'frac>baseline':>13s}")
    for ab, run, n, med, fr in rows:
        print(f"  {ab:6s} {run:6s} {n:4d} {med:11.3f} {fr:13.2f}")

    print("\nSpearman rho vs cognate_log_enrichment:")
    print(f"  POOLED (n={len(cog)}, confounded by between-Ab offsets):")
    for c in METRICS:
        v = cog[[c, "cognate_log_enrichment"]].corr("spearman").iloc[0, 1]
        print(f"    {c:28s} {v:+.3f}")
    print("  WITHIN-Ab (the honest view; only n>=12 Abs shown):")
    for ab, g in cog.groupby("cognate_ab"):
        if len(g) >= 12:
            parts = " ".join(
                f"{c.split('_')[0]}={g[[c,'cognate_log_enrichment']].corr('spearman').iloc[0,1]:+.2f}"
                for c in ["overall_rmsd", "epitope_chunk_rmsd_vs_mpnn", "mean_pae"])
            print(f"    {ab} (n={len(g)}): {parts}")

    print(f"\nKd present only for published controls: {d.Kd.notna().sum()} rows "
          f"(none scaffolded: {scaf.Kd.notna().sum()})")

    _figure(d, cog)


def _figure(d: pd.DataFrame, cog: pd.DataFrame) -> None:
    from matplotlib.lines import Line2D
    order = sorted(RUNS, key=lambda ab: -(cog.cognate_ab == ab).sum())
    fig, axes = plt.subplots(2, 4, figsize=(12.5, 6.8))
    for ax, ab in zip(axes.flat, order):
        ab_col, noab_col, run = RUNS[ab]
        x, y = log(d[noab_col]), log(d[ab_col])
        ax.scatter(x, y, s=6, c="0.78", linewidths=0)
        m = d.cognate_ab == ab
        ax.scatter(x[m], y[m], s=18, c="crimson", edgecolors="darkred", linewidths=0.3)
        lim = [0, max(x.max(), y.max()) * 1.02]
        ax.plot(lim, lim, "--", c="0.4", lw=0.8)  # NoAb diagonal
        ax.set_xlim(lim); ax.set_ylim(lim)
        ax.set_title(f"{ab}  ({run}, n={int(m.sum())})", fontsize=15)  # cognate count in title
        ax.set_xlabel(r"$\log_{10}(1+\mathrm{NoAb})$")
        ax.set_ylabel(r"$\log_{10}(1+\mathrm{Ab})$")
    fig.suptitle("DP3 binding: cognate scaffolded designs vs library background",
                 fontsize=18)
    handles = [
        Line2D([0], [0], marker="o", color="none", markerfacecolor="0.78",
               markersize=10, label="all library members"),
        Line2D([0], [0], marker="o", color="none", markerfacecolor="crimson",
               markeredgecolor="darkred", markersize=10, label="cognate scaffolded designs"),
    ]
    fig.legend(handles=handles, loc="lower center", ncol=2, fontsize=15,
               frameon=False, bbox_to_anchor=(0.5, -0.01))
    fig.tight_layout(rect=[0, 0.05, 1, 0.97])
    FIG.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG, dpi=150)
    print(f"\nwrote {FIG}")


if __name__ == "__main__":
    main()
