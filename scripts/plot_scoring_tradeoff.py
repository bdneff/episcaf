#!/usr/bin/env python3
"""plot_scoring_tradeoff.py -- the accessibility-vs-fold tradeoff figure, comprehensive.

One figure that answers "how does the scorer decide, and what did switching to the soft-gate change?"
Top: what each metric measures, the soft-gate score itself, and the selection policy. Bottom: a 2-D
scatter (epitope RMSD x, real Fab clash y; lower-left = best) showing where the top-k selection lands
under three scorers -- the old percentile, a naive clash-up-weighted sigmoid, and the adopted soft-gate
(the CURRENT `antibody_softgate`, epitope-PAE midpoint 2.5) -- pooled and for example targets.

Data: C1 needs its full pool (local `metrics_whole_epitope_103.csv`); C2's pool is the cluster
`metrics_dual_island.parquet`; 8VDL's candidates are `results/dp4_8vdl_top10_allmetrics.csv`. Any pool
that is absent is skipped (the figure still renders the rest).

  # local (C1 + 8VDL):
  python scripts/plot_scoring_tradeoff.py \
     --c1 $D/known_antigen/analysis/data/metrics_whole_epitope_103.csv \
     --vdl results/dp4_8vdl_top10_allmetrics.csv --out scratch_figs/scoring_tradeoff.png
  # full (on the cluster, add C2):
  #   --c2 $WS/runs/dual_island_rfd3/05_analysis/metrics_dual_island.parquet
"""
from __future__ import annotations
import argparse, sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
from episcaf_analysis.score import score            # noqa: E402
from episcaf_analysis.presets import PRESETS         # noqa: E402

CLASH, ERMSD, ORMSD, EPAE, MPAE = "af3_n_clash_res", "epitope_chunk_rmsd", "overall_rmsd", "epitope_pae", "mean_pae"
GREY, RED, BLUE = "#9a9a9a", "#c0392b", "#1f6fb2"

# ---- the three scorers compared (soft-gate is read live from presets, so it tracks the shipped dials) --
def _sig(mid, k, w):  return dict(better="low", transform="sigmoid", midpoint=mid, k=k, weight=w)
def scorers():
    pct = {c: dict(weight=w, better="low", transform="percentile")
           for c, w in [(CLASH, .35), (ERMSD, .35), (ORMSD, .15), (EPAE, .15)]}
    clash50 = {CLASH: _sig(6, .5, .50), ERMSD: _sig(1, 4, .25), ORMSD: _sig(2, 2, .15), EPAE: _sig(5, 1, .10)}
    soft = {k: (v.copy() if isinstance(v, dict) else v) for k, v in PRESETS["antibody_softgate"].items()}
    return {
        "percentile (old)":            (GREY, dict(scope="pooled", gate=None, metrics=pct)),
        "sigmoid, clash weight .50":   (RED,  dict(scope="pooled", gate=None, metrics=clash50)),
        "soft-gate (adopted)":         (BLUE, {k: v for k, v in soft.items() if k != "select"}),
    }

def selected_median(df, preset, group, topk, idval=None):
    """Median (epitope RMSD, clash) of the top-k selected per group under `preset`."""
    p = {k: (v.copy() if isinstance(v, dict) else v) for k, v in preset.items()}
    p["select"] = None
    s = score(df, p)
    s["_r"] = s.groupby(group, sort=False)["composite"].rank(ascending=False, method="first")
    sel = s[s._r <= topk]
    if idval is not None:
        sel = sel[sel[group].astype(str) == str(idval)]
    return (pd.to_numeric(sel[ERMSD], errors="coerce").median(),
            pd.to_numeric(sel[CLASH], errors="coerce").median(), len(sel))

def load(path):
    if path is None or not Path(path).exists():
        return None
    d = pd.read_parquet(path) if str(path).endswith(".parquet") else pd.read_csv(path, low_memory=False)
    if "status" in d.columns:
        d = d[d["status"].astype(str).str.lower().eq("ok")].copy()
    return d


# ---- header text -----------------------------------------------------------------------------------
METRIC_DEFS = [
    ("epitope RMSD",       "epitope-residue backbone RMSD (design vs AF3):\nis the epitope in its native shape?  lower = truer"),
    ("overall RMSD",       "whole-construct backbone RMSD:\nis the entire fold right?  lower = truer"),
    ("epitope PAE",        "AF3 aligned-error within the epitope block:\nconfidence in the epitope's internal geometry"),
    ("scaffold / mean PAE","PAE of the scaffold, and of the whole matrix\n(the four-filter's < 5 gate); diagnostics"),
    ("AF3 clash",          "scaffold residues clashing with the REAL antibody\n(known-Ab arms) -- ACCESSIBILITY; lower = more open"),
    ("cylinder clash",     "scaffold residues in the native-aware cylinder:\nthe accessibility SURROGATE when no antibody exists"),
    ("ptm",                "AF3 predicted TM-score:\nglobal confidence in the model"),
]
SELECTION = [
    r"$\bf{Rank,\ don't\ gate.}$  Score every design; keep the top-$n$ per group. No hard threshold --",
    "   a hard fold floor would empty 3 of 87 C2 islands; the soft gate crushes misfolds toward 0 but",
    "   never TO 0, so every target still ships its best.",
    r"$\bf{Depth\ =\ top}$-$n$ per group:  C1 top-20 / epitope,  C2 top-20 / island,  C3 top-10 / window.",
    r"$\bf{Global\text{-}pass\ promotion.}$  Any design clearing all four Lawson filters (epitope RMSD $\leq$1,",
    "   overall $\leq$2, mean PAE < 5, clash = 0) is floated above every non-passer (a soft AND).",
    r"$\bf{Midpoints}$ sit at the four-filter thresholds; the epitope-PAE midpoint was retuned to 2.5 from",
    "   the data (good epitopes sit ~2 A, not 5). Weights/steepness are a prior, to be fit on DP4 binding.",
]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--c1"); ap.add_argument("--c2"); ap.add_argument("--vdl")
    ap.add_argument("--topk", type=int, default=10)
    ap.add_argument("--out", default="scratch_figs/scoring_tradeoff.png")
    args = ap.parse_args()

    SC = scorers()
    # rows: (title, dataframe, group-col, example-id-or-None, x/y annotation-note-or-None)
    rows = []
    c1 = load(args.c1); c2 = load(args.c2); vdl = load(args.vdl)
    if c1 is not None:
        rows.append(("C1 * whole epitope (pooled median)", c1, "id", None))
        rows.append(("6cyf_0P * whole epitope", c1, "id", "6cyf_0P"))
    if c2 is not None:
        c2 = c2.copy()
        c2["_g"] = c2["id"].astype(str) + "|" + c2.get("island_index", 0).astype(str)
        rows.append(("C2 * single island (pooled median)", c2, "_g", None))
    if vdl is not None:
        vdl = vdl.copy()
        vdl["_g"] = vdl["run"].astype(str)
        rows.append(("8VDL hotspots (John's example)", vdl[vdl.run == "hotspots"], "_g", None))
        rows.append(("8VDL epitope", vdl[vdl.run == "epitope"], "_g", None))
    if not rows:
        sys.exit("no pools given -- pass at least --c1")

    n = len(rows)
    # rows: [0] defs|score header, [1] selection policy (full width), [2] legend strip, [3..] scatters
    fig = plt.figure(figsize=(13, 6.0 + 2.4 * n))
    gs = fig.add_gridspec(3 + n, 2, height_ratios=[3.2, 2.0, 0.35] + [2.3] * n,
                          hspace=0.6, wspace=0.12, left=0.075, right=0.975, top=0.955, bottom=0.045)

    fig.suptitle("Soft-gate cuts clashes while keeping the fold -- and how the selection is made",
                 fontsize=17, fontweight="bold", y=0.99)

    # --- header row 0: metric defs (left) + score table (right) ---
    axL = fig.add_subplot(gs[0, 0]); axL.axis("off")
    axL.text(0, 1.04, "What each metric measures", fontsize=13, fontweight="bold", va="top")
    y = 0.88
    for name, desc in METRIC_DEFS:
        axL.text(0.0, y, name, fontsize=8.8, fontweight="bold", va="top", color="#111")
        axL.text(0.34, y + 0.005, desc, fontsize=7.3, va="top", color="#333", linespacing=1.25)
        y -= 0.142
    axR = fig.add_subplot(gs[0, 1]); axR.axis("off")
    axR.text(0, 1.02, "The soft-gate score", fontsize=13, fontweight="bold", va="top")
    axR.text(0.02, 0.84, r"$S=\sum_m w_m\,\sigma(x_m)$,      $\sigma(x)=\dfrac{1}{1+e^{\,k(x-m)}}$",
             fontsize=12, va="top")
    tbl = [("AF3 clash", ".45", "6", "0.5", "accessibility (ranked)"),
           ("epitope RMSD", ".25", "1", "4", "fidelity (soft gate)"),
           ("overall RMSD", ".20", "2", "4", "fold (soft gate)"),
           ("epitope PAE", ".10", "2.5", "1.2", "confidence (soft gate)")]
    axR.text(0.02, 0.56, f"{'metric':<14}{'w':>5}{'m':>5}{'k':>6}   role", fontsize=9,
             family="monospace", va="top", color="#666")
    yy = 0.46
    for m, w, mid, k, role in tbl:
        col = "#111" if m == "AF3 clash" else BLUE
        axR.text(0.02, yy, f"{m:<14}{w:>5}{mid:>5}{k:>6}   {role}", fontsize=9,
                 family="monospace", va="top", color=col)
        yy -= 0.11
    axR.text(0.02, 0.0, "Each metric squashed by its own sigmoid, then\nweighted-summed. Broad on clash "
             "(ranks accessibility);\nsteep on fold (soft gates). Large k -> a hard gate,\nbut finite, so "
             "nothing is ever zeroed.", fontsize=7.5, va="top", color="#777", linespacing=1.3)

    # --- header row 1: selection policy (full width) ---
    axS = fig.add_subplot(gs[1, :]); axS.axis("off")
    axS.text(0, 1.0, "How we decide what ships (the selection criteria)", fontsize=13,
             fontweight="bold", va="top")
    axS.text(0.0, 0.80, "\n".join(SELECTION), fontsize=9, va="top", color="#222", linespacing=1.5)

    # --- legend strip (row 2) ---
    axLg = fig.add_subplot(gs[2, :]); axLg.axis("off")
    handles = [Line2D([0], [0], marker="o", color="w", markerfacecolor=c, markersize=13, label=l)
               for l, (c, _) in SC.items()]
    axLg.legend(handles=handles, loc="center", ncol=3, frameon=False, fontsize=12,
                handletextpad=0.4, columnspacing=2.5)

    # --- scatter rows ---
    for i, (title, df, group, idval) in enumerate(rows):
        ax = fig.add_subplot(gs[i + 3, :])
        pts = []
        for name, (color, preset) in SC.items():
            x, yv, k = selected_median(df, preset, group, args.topk, idval)
            if np.isfinite(x) and np.isfinite(yv):
                ax.scatter([x], [yv], s=170, color=color, edgecolor="white", linewidth=1.5, zorder=3)
                pts.append((x, yv, color, name))
        ax.axvspan(1.0, ax.get_xlim()[1] if ax.get_xlim()[1] > 1 else 2.2, color="#efefef", zorder=0)
        ax.axvline(1.0, ls="--", color="#bbb", lw=1)
        ax.text(0.99, 0.93, "RMSD > 1 Å", transform=ax.transAxes, ha="right", fontsize=8, color="#999", style="italic")
        ax.set_title(title, fontsize=11, fontweight="bold", loc="left")
        ax.set_xlabel("epitope RMSD (Å)   ·   lower = truer fold", fontsize=9)
        ax.set_ylabel("AF3 clashing residues\n· lower = more accessible", fontsize=9)
        ax.margins(0.25)
        ax.grid(True, color="#eee", lw=0.7)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)

    fig.text(0.5, 0.012, "Pooled top-%d median per group · scripts/plot_scoring_tradeoff.py · "
             "soft-gate read live from presets.py (epitope-PAE midpoint 2.5)" % args.topk,
             ha="center", fontsize=8, color="#888")
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.out, dpi=140, bbox_inches="tight")
    print(f"[tradeoff] wrote {args.out}  ({n} rows: {[r[0] for r in rows]})")


if __name__ == "__main__":
    main()
