#!/usr/bin/env python3
"""Analyze DP4 assay binding: does our scoring predict binding? (Q1-Q6 of the plan.)

Reads the join output `results/dp4_binding_metrics.csv` (from build_dp4_binding_join.py) and
answers the questions in docs/DP4_RESULTS_ANALYSIS.md, then writes two figures. It is the DP4
analogue of DP3's `analyze_dp3_binding.py` + `plot_dp3_metric_binding.py`, merged into one command.

WHY THIS IS ASSAY-LAYOUT-INDEPENDENT (and could be written + verified before the data landed):
the join step already resolved the raw, unknown assay columns into a FIXED schema --
`cognate_ab`, `cognate_log_enrichment` (= log10(1+Ab) - log10(1+NoAb)), `cognate_log_ab/noab`,
plus every design metric and score (`composite`, `rank_in_group`, `epitope_rmsd`, `overall_rmsd`,
`mean_pae`, `cylinder_clashes`, `af3_clashes`). This script only ever touches those stable names,
so it does NOT depend on how the assay CSV was laid out or which antibodies were on the panel --
it groups by whatever `cognate_ab` values appear. Smoke-tested on a synthetic join output
2026-08-05 (fixed schema); run it for real once build_dp4_binding_join.py has produced the CSV.

The readout and the deconfounding rule (carry over from DP3, do NOT skip):
  binding = cognate_log_enrichment; NEVER trust a pooled correlation -- subtract each antibody's
  mean from BOTH the metric and the enrichment (within-antibody fixed-effect centering) before
  correlating. Pooling mixes in each antibody's baseline offset and manufactures signal.

Arm <-> category strings (the plan's C-labels are these `category` values):
  C1 scaffoldedAbEpitope · C2 scaffoldedSingleIsland · C3 scaffoldedPolyclonal · C4 tiled30mer
  C5 metricSpaceTitration · C6 scaffoldedEpitopeControl · 8VDL scaffolded8VDL · minibinder.
  composite/metrics are present on C1/C2/C5/8VDL; BLANK on C3/C4/C6/minibinder (never folded / no
  antibody), so metric regressions run only where the axis is non-null.

Run:
    python scripts/analyze_dp4_binding.py                       # after the join has been built
    python scripts/analyze_dp4_binding.py --metrics-csv <path>  # point at a specific join output
"""
from __future__ import annotations
import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams.update({"font.size": 13, "axes.titlesize": 15, "axes.labelsize": 13,
                     "xtick.labelsize": 11, "ytick.labelsize": 11, "legend.fontsize": 10})

ROOT = Path(__file__).resolve().parents[1]
CSV = ROOT / "results" / "dp4_binding_metrics.csv"
FIG_SCATTER = ROOT / "manuscript" / "figures" / "dp4_binding_scatter.png"
FIG_METRIC = ROOT / "manuscript" / "figures" / "dp4_metric_binding.png"

Y = "cognate_log_enrichment"
# Axes to test against binding. Higher composite = better -> expect +; lower rank = better -> expect -;
# the RMSD/PAE/clash axes are costs -> expect -. Only rows where the axis is non-null are used.
SCORE_AXES = [("composite", "composite score"), ("rank_in_group", "rank in group")]
METRIC_AXES = [("epitope_rmsd", "epitope RMSD (Å)"), ("overall_rmsd", "overall RMSD (Å)"),
               ("mean_pae", "mean PAE"), ("cylinder_clashes", "cylinder clashes"),
               ("af3_clashes", "AF3 clashes")]
ARM = {"scaffoldedAbEpitope": "C1", "scaffoldedSingleIsland": "C2", "scaffoldedPolyclonal": "C3",
       "tiled30mer": "C4", "metricSpaceTitration": "C5", "scaffoldedEpitopeControl": "C6",
       "scaffolded8VDL": "8VDL", "minibinder": "mini"}


def within_center(df: pd.DataFrame, cols: list[str], by: str = "cognate_ab") -> pd.DataFrame:
    """Group-mean-center each col (and Y) within antibody -- the fixed-effect deconfounding."""
    g = df.groupby(by)
    out = df.copy()
    out["_yc"] = df[Y] - g[Y].transform("mean")
    for c in cols:
        out[c + "_c"] = df[c] - g[c].transform("mean")
    return out


def corr_pair(x: pd.Series, y: pd.Series):
    m = x.notna() & y.notna()
    x, y = x[m], y[m]
    if len(x) < 3 or x.nunique() < 2 or y.nunique() < 2:
        return None
    r, p = stats.pearsonr(x, y)
    rho, _ = stats.spearmanr(x, y)
    return r, rho, p, len(x)


def report(d: pd.DataFrame) -> pd.DataFrame:
    print(f"rows: {len(d):,}   columns: {len(d.columns)}")
    if "category" in d:
        print("\nby arm:")
        vc = d["category"].value_counts()
        for cat, n in vc.items():
            print(f"  {ARM.get(cat, '?'):5s} {cat:26s} {n:,}")

    if "cognate_ab" not in d or Y not in d:
        sys.exit("\nERROR: join output lacks cognate_ab / cognate_log_enrichment -- rebuild with "
                 "build_dp4_binding_join.py.")

    cog = d[d["cognate_ab"].notna() & d[Y].notna()].copy()
    print(f"\ncognate rows (an antibody resolved AND an enrichment value): {len(cog):,}")
    if len(cog) == 0:
        print("!! No cognate rows. Most likely the join's antibody-column detection didn't match, "
              "or the 8VDL antibody is named C7 (target says 8VDL) -- re-run the join with --ab-map. "
              "See docs/DP4_RESULTS_ANALYSIS.md. Nothing to correlate; stopping the cognate analysis.")
        return cog
    print("antibodies:", ", ".join(f"{ab}(n={n})" for ab, n
                                    in cog["cognate_ab"].value_counts().items()))

    # per-antibody enrichment
    print("\nper-antibody cognate enrichment:")
    print(f"  {'Ab':8s} {'n':>5s} {'median':>8s} {'frac>0':>7s}")
    for ab, g in cog.groupby("cognate_ab"):
        print(f"  {ab:8s} {len(g):5d} {g[Y].median():8.3f} {(g[Y] > 0).mean():7.2f}")

    # Q1-Q3: does score / cylinder / metric predict binding? pooled vs within-antibody
    axes = SCORE_AXES + METRIC_AXES
    present = [(c, lab) for c, lab in axes if c in cog and cog[c].notna().sum() >= 5]
    cen = within_center(cog, [c for c, _ in present])
    print("\n=== binding vs score/metric (Pearson r) -- POOLED is confounded; WITHIN-Ab is honest ===")
    print(f"  {'axis':16s} {'n':>5s} {'pooled_r':>9s} {'within_r':>9s} {'within_rho':>11s}  expect")
    expect = {"composite": "+", "rank_in_group": "-", "epitope_rmsd": "-", "overall_rmsd": "-",
              "mean_pae": "-", "cylinder_clashes": "-", "af3_clashes": "-"}
    for c, lab in present:
        pooled = corr_pair(cog[c], cog[Y])
        within = corr_pair(cen[c + "_c"], cen["_yc"])
        if pooled and within:
            print(f"  {c:16s} {pooled[3]:5d} {pooled[0]:+9.3f} {within[0]:+9.3f} "
                  f"{within[1]:+11.3f}  {expect.get(c, '?')}")
    print("  (composite +, ranks/costs - would mean the design score tracks real binding.)")

    # Q4: 8VDL accessible (epitope) vs blocked (hotspots)
    v = cog[cog["category"] == "scaffolded8VDL"] if "category" in cog else cog.iloc[0:0]
    if len(v):
        print("\nQ4  8VDL accessible-vs-blocked (median cognate enrichment):")
        for tgt, g in v.groupby("target"):
            print(f"    {tgt:16s} n={len(g):3d}  median={g[Y].median():+.3f}")
    else:
        print("\nQ4  8VDL: no cognate 8VDL rows (its antibody is C7; re-run join with --ab-map "
              "mapping 8vdl/C7 to confirm the accessible-vs-blocked contrast).")

    # Q5: scaffolded (C1) vs linear tiled (C4), matched by target
    if "category" in cog:
        c1 = cog[cog["category"] == "scaffoldedAbEpitope"]
        c4 = cog[cog["category"] == "tiled30mer"]
        shared = sorted(set(c1["target"]) & set(c4["target"])) if len(c1) and len(c4) else []
        if shared:
            print(f"\nQ5  scaffolded(C1) vs linear(C4), {len(shared)} shared targets "
                  f"(median enrichment):")
            print(f"    {'target':12s} {'C1':>8s} {'C4':>8s}")
            for t in shared[:15]:
                print(f"    {t:12s} {c1[c1.target == t][Y].median():+8.3f} "
                      f"{c4[c4.target == t][Y].median():+8.3f}")
        else:
            print("\nQ5  scaffolded-vs-linear: no shared targets with cognate enrichment yet "
                  "(depends on which antibodies were assayed).")
    return cog


def fig_scatter(cog: pd.DataFrame) -> None:
    if len(cog) == 0 or "cognate_log_noab" not in cog:
        print("[fig] scatter skipped (no cognate rows).")
        return
    abs_ = cog["cognate_ab"].value_counts().index.tolist()
    n = len(abs_)
    ncol = min(4, n)
    nrow = int(np.ceil(n / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(3.3 * ncol, 3.2 * nrow), squeeze=False)
    for ax, ab in zip(axes.flat, abs_):
        g = cog[cog["cognate_ab"] == ab]
        x, y = g["cognate_log_noab"], g["cognate_log_ab"]
        ax.scatter(x, y, s=16, c="crimson", edgecolors="darkred", linewidths=0.3)
        hi = max(float(np.nanmax(x)) if len(x) else 1, float(np.nanmax(y)) if len(y) else 1, 1) * 1.05
        ax.plot([0, hi], [0, hi], "--", c="0.4", lw=0.8)
        ax.set_xlim(0, hi); ax.set_ylim(0, hi)
        ax.set_title(f"{ab} (n={len(g)})")
        ax.set_xlabel(r"$\log_{10}(1+\mathrm{NoAb})$")
        ax.set_ylabel(r"$\log_{10}(1+\mathrm{Ab})$")
    for ax in axes.flat[n:]:
        ax.axis("off")
    fig.suptitle("DP4 cognate designs: binders sit above the NoAb diagonal", fontsize=16)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    FIG_SCATTER.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG_SCATTER, dpi=150)
    print(f"[fig] wrote {FIG_SCATTER}")


def fig_metric(cog: pd.DataFrame) -> None:
    axes_spec = [(c, lab) for c, lab in SCORE_AXES + METRIC_AXES
                 if c in cog and cog[c].notna().sum() >= 5]
    if len(cog) == 0 or not axes_spec:
        print("[fig] metric-binding skipped (no scored cognate rows).")
        return
    cen = within_center(cog, [c for c, _ in axes_spec])
    ncol = len(axes_spec)
    fig, axes = plt.subplots(1, ncol, figsize=(3.2 * ncol, 3.4), squeeze=False)
    for j, (c, lab) in enumerate(axes_spec):
        ax = axes[0, j]
        m = cen[c + "_c"].notna() & cen["_yc"].notna()
        ax.scatter(cen.loc[m, c + "_c"], cen.loc[m, "_yc"], s=14, c="0.5", linewidths=0)
        res = corr_pair(cen[c + "_c"], cen["_yc"])
        if res:
            r = res[0]
            xs = np.linspace(cen.loc[m, c + "_c"].min(), cen.loc[m, c + "_c"].max(), 50)
            b1, b0 = np.polyfit(cen.loc[m, c + "_c"], cen.loc[m, "_yc"], 1)
            ax.plot(xs, b0 + b1 * xs, "k", lw=2)
            ax.set_title(f"{lab}\nwithin-Ab r={r:+.2f}")
        ax.set_xlabel(f"{lab} (centered)")
        if j == 0:
            ax.set_ylabel("within-Ab\nlog-enrichment")
    fig.suptitle("DP4: does the design score/metric predict binding? (within-antibody)", fontsize=15)
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    FIG_METRIC.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG_METRIC, dpi=150)
    print(f"[fig] wrote {FIG_METRIC}")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--metrics-csv", type=Path, default=CSV,
                    help="the join output to analyze (default results/dp4_binding_metrics.csv)")
    ap.add_argument("--no-figures", action="store_true", help="print the report only")
    args = ap.parse_args()

    if not args.metrics_csv.exists():
        sys.exit(f"ERROR: {args.metrics_csv} not found. Build it first:\n"
                 f"    python scripts/build_dp4_binding_join.py --assay data/dp4_binding/<run>.csv")
    d = pd.read_csv(args.metrics_csv, low_memory=False)
    cog = report(d)
    if not args.no_figures:
        fig_scatter(cog)
        fig_metric(cog)


if __name__ == "__main__":
    main()
