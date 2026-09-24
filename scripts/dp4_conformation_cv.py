#!/usr/bin/env python3
"""Conformational selection with a real collective variable from static structures.

We cannot measure a single peptide's conformational distribution without an ensemble (MD). But every
design gives one honest coordinate along a real collective variable: the epitope's RMSD to the native
antibody-binding conformation -- how far the scaffolded epitope sits from the shape the antibody
recognizes. Across the design ensemble, binding rate falls as that distance grows: designs near the
binding conformation bind most. (This is a between-design/epitope trend, not a within-molecule
distribution; the single-molecule version awaits apo/holo MD.)

Run:  /usr/bin/python3 scripts/dp4_conformation_cv.py
"""
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[1]
LIB=ROOT/"data/libraries/dp4_library.csv"; SUMM=ROOT/"data/dp4_binding/john/scaffoldedEpitopeSummary.csv"
OUT=ROOT/"data/dp4_binding/figs/conformation_cv.png"
BLUE="#2456E6"; GREEN="#1f8a4c"; MUTE="#5b6570"

def main():
    lib=pd.read_csv(LIB, low_memory=False)
    summ=pd.read_csv(SUMM, encoding="latin-1"); summ.columns=[c.strip() for c in summ.columns]
    hits=set()
    for ids in summ["hitIDs"].dropna():
        for t in str(ids).split(","):
            t=t.strip()
            if t.startswith("DP4_"): hits.add("DP4_"+t.split("_")[1])
    d=lib[lib.category.isin(["scaffoldedAbEpitope","scaffoldedSingleIsland","metricSpaceTitration"])].copy()
    d["bound"]=d.library_member.isin(hits).astype(int)
    d["r"]=pd.to_numeric(d.epitope_rmsd, errors="coerce"); d=d.dropna(subset=["r"])

    edges=[0,0.5,1.0,1.5,2.0,2.5,3.0,4.0,7.0]
    d["bin"]=pd.cut(d.r, edges)
    g=d.groupby("bin").agg(rate=("bound","mean"), n=("bound","size"), c=("r","mean")).dropna()
    err=np.sqrt(g.rate*(1-g.rate)/g.n)*100

    plt.rcParams.update({"font.family":"DejaVu Sans"})
    fig,ax=plt.subplots(figsize=(9.4,6.2))
    ax.set_xlim(0,4.2); ax.set_ylim(0,17.5)
    ax.axvspan(0,0.5, color=GREEN, alpha=0.07, zorder=0)
    # trend + error bars + points sized by number of designs
    ax.plot(g.c, g.rate*100, color=BLUE, lw=1.6, alpha=0.45, zorder=2)
    ax.errorbar(g.c, g.rate*100, yerr=err, fmt="none", ecolor=BLUE, elinewidth=1.4, alpha=0.5, capsize=3, zorder=2)
    ax.scatter(g.c, g.rate*100, s=np.clip(g.n/5.0,40,300), color=BLUE, edgecolors="white", linewidths=1.3, zorder=3)
    # labels placed in empty regions (top-right is empty; the shaded band's floor is empty)
    ax.text(0.25, 1.1, "near the\nbinding shape", color=GREEN, fontsize=11.5, ha="center", va="bottom", fontweight="bold")
    # size legend: two reference circles for designs-per-bin
    for nv in [300,1000]:
        ax.scatter([],[], s=np.clip(nv/5.0,40,300), color=BLUE, edgecolors="white", linewidths=1.3, label=f"{nv:,} designs")
    ax.legend(loc="upper right", title="point size = designs in the bin", labelspacing=1.7, borderpad=1.0,
              handletextpad=1.4, frameon=True, framealpha=0.95, fontsize=11, title_fontsize=10)
    # which arms this considers
    ax.text(0.015, 0.965, f"C1 + C2 + C5  —  scaffolded designs across the metric range (n = {len(d):,})",
            transform=ax.transAxes, fontsize=11, color=MUTE, va="top")
    ax.set_xlabel("epitope RMSD from the native antibody-binding conformation  (Å)", fontsize=13)
    ax.set_ylabel("binding rate  (%)", fontsize=13)
    ax.set_title("Designs closer to the antibody-binding conformation bind more often",
                 fontsize=14, fontweight="bold", pad=12)
    ax.spines[["top","right"]].set_visible(False); ax.grid(axis="y", alpha=0.14)
    fig.tight_layout(); fig.savefig(OUT, dpi=160)
    print(g.assign(err=err).round(3).to_string()); print(f"wrote {OUT.relative_to(ROOT)}")

if __name__=="__main__":
    main()
