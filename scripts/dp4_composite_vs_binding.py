#!/usr/bin/env python3
"""Composite score vs experimental binding, at the DESIGN level, for EVERY tested scaffolded design.

John's helixVsComposite plot shows one point per epitope (the best design's composite). This shows
the whole distribution: every scaffolded design whose antibody was assayed, split into the ones that
bound (John's hitIDs) and the ones that did not. That is the honest "does composite predict binding"
picture, across the full range of composite rather than the top per epitope.

Label: a design is a HIT if its library_member is in any epitope's hitIDs in
data/dp4_binding/john/scaffoldedEpitopeSummary.csv. The tested universe is scaffolded designs
(composite present) whose epitope appears in that summary (i.e. its antibody was assayed).

Run:  /usr/bin/python3 scripts/dp4_composite_vs_binding.py
"""
from pathlib import Path
import numpy as np, pandas as pd
from scipy import stats
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
LIB = ROOT / "data/libraries/dp4_library.csv"
SUMM = ROOT / "data/dp4_binding/john/scaffoldedEpitopeSummary.csv"
OUT = ROOT / "data/dp4_binding/figs/composite_vs_binding.png"

BOUND, MISS = "#2a6f97", "#b0b8c0"   # bound / did-not-bind

def main():
    lib = pd.read_csv(LIB, low_memory=False)
    summ = pd.read_csv(SUMM, encoding="latin-1")
    summ.columns = [c.strip() for c in summ.columns]

    # design-level hit set from hitIDs (DP4_<n>_<m> -> library_member DP4_<n>)
    hits = set()
    for ids in summ["hitIDs"].dropna():
        for tok in str(ids).split(","):
            tok = tok.strip()
            if tok.startswith("DP4_"):
                hits.add("DP4_" + tok.split("_")[1])
    assayed = set(summ["epitope"].astype(str).str.lower())
    print(f"assayed epitopes in summary: {len(assayed)}   design-level hits (hitIDs): {len(hits)}")

    lib["epitope"] = lib["target"].astype(str).str.split("_").str[0].str.lower()
    d = lib[lib["composite"].notna() & lib["epitope"].isin(assayed)].copy()
    d["bound"] = d["library_member"].isin(hits)
    print(f"tested scaffolded designs (composite present, antibody assayed): {len(d):,}")
    print(f"  bound: {int(d['bound'].sum()):,}   did not bind: {int((~d['bound']).sum()):,}"
          f"   overall hit rate: {d['bound'].mean():.3f}")

    b = d.loc[d.bound, "composite"]; m = d.loc[~d.bound, "composite"]
    print(f"  median composite  bound={b.median():.3f}  miss={m.median():.3f}")
    U, p = stats.mannwhitneyu(b, m, alternative="greater")
    auc = U / (len(b) * len(m))   # P(composite_bound > composite_miss)
    print(f"  AUC (composite separates bound vs miss) = {auc:.3f}   Mann-Whitney p = {p:.2g}")

    # ---- figure: distribution (left) + hit-rate vs composite (right) ----
    fig, (axA, axB) = plt.subplots(1, 2, figsize=(12, 5), gridspec_kw={"width_ratios": [1, 1.15]})

    # A: full distribution, every design, split by outcome
    for i, (lab, s, c) in enumerate([("did not bind", m, MISS), ("bound", b, BOUND)]):
        x = np.random.default_rng(0).normal(i, 0.07, len(s))
        axA.scatter(x, s, s=6, color=c, alpha=0.25, edgecolors="none", zorder=2)
        axA.boxplot(s, positions=[i], widths=0.5, showfliers=False,
                    medianprops=dict(color="black", lw=2), boxprops=dict(color="black"),
                    whiskerprops=dict(color="black"), capprops=dict(color="black"), zorder=3)
    axA.set_xticks([0, 1]); axA.set_xticklabels([f"did not bind\n(n={len(m):,})", f"bound\n(n={len(b):,})"])
    axA.set_ylabel("composite score"); axA.set_title("Every tested scaffolded design")
    axA.spines[["top", "right"]].set_visible(False)

    # B: probability of binding vs composite (equal-count bins)
    dd = d.sort_values("composite")
    nbin = 12
    dd["bin"] = pd.qcut(dd["composite"], nbin, labels=False, duplicates="drop")
    g = dd.groupby("bin").agg(x=("composite", "median"), rate=("bound", "mean"), n=("bound", "size"))
    axB.plot(g["x"], g["rate"], "-o", color=BOUND, lw=2, ms=6, zorder=3)
    axB.scatter(dd["composite"], dd["bound"].astype(int) * 1.0, s=5, color="#c9ced6", alpha=0.3, zorder=1)
    for _, r in g.iterrows():
        axB.annotate(f"{int(r['n'])}", (r["x"], r["rate"]), textcoords="offset points",
                     xytext=(0, 7), ha="center", fontsize=7, color="#5a6570")
    axB.set_xlabel("composite score"); axB.set_ylabel("fraction that bound (hit rate)")
    axB.set_title("Binding rate rises with composite"); axB.set_ylim(-0.03, 1.03)
    axB.spines[["top", "right"]].set_visible(False)

    fig.suptitle(f"DP4: composite vs experimental binding  —  {len(d):,} designs, "
                 f"{int(d['bound'].sum()):,} bound  (AUC {auc:.2f}, p={p:.1g})",
                 fontsize=13, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, dpi=150)
    print(f"\nwrote {OUT.relative_to(ROOT)}")

if __name__ == "__main__":
    main()
