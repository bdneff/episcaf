#!/usr/bin/env python3
"""Non-circular validation: does the composite predict binding on UNSELECTED designs?

C5 (metricSpaceTitration) spans the metric space on purpose; it was NOT selected for high composite, so
it is the honest test of the scorer that C1/C2 cannot give (those are already the top-ranked designs).
We recompute the adopted composite (antibody_softgate preset) on C5's raw metrics and ask whether binding
rate rises with composite -- a prediction, not a selection artifact. C5's own mean is the unselected
base rate.

Run:  /usr/bin/python3 scripts/dp4_c5_titration.py
"""
from pathlib import Path
import sys
import numpy as np, pandas as pd
from scipy import stats
from sklearn.metrics import roc_auc_score
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from episcaf_analysis.score import score
from episcaf_analysis.presets import ANTIBODY_SOFTGATE

LIB = ROOT / "data/libraries/dp4_library.csv"
SUMM = ROOT / "data/dp4_binding/john/scaffoldedEpitopeSummary.csv"
OUT = ROOT / "data/dp4_binding/figs/c5_titration.png"

def main():
    lib = pd.read_csv(LIB, low_memory=False)
    summ = pd.read_csv(SUMM, encoding="latin-1"); summ.columns=[c.strip() for c in summ.columns]
    hits=set()
    for ids in summ["hitIDs"].dropna():
        for t in str(ids).split(","):
            t=t.strip()
            if t.startswith("DP4_"): hits.add("DP4_"+t.split("_")[1])

    c5 = lib[lib.category=="metricSpaceTitration"].copy()
    c5["epitope"] = c5["target"].astype(str).str.split("_").str[0].str.lower()
    c5["bound"] = c5["library_member"].isin(hits).astype(int)
    # recompute the adopted composite on C5's raw metrics (rename to the preset's metric keys)
    c5r = c5.rename(columns={"af3_clashes":"af3_n_clash_res", "epitope_rmsd":"epitope_chunk_rmsd"}).copy()
    c5r["id"] = c5r["library_member"]; c5r["antigen"] = c5r["epitope"]
    c5["composite"] = score(c5r, ANTIBODY_SOFTGATE)["composite"].values

    y = c5["bound"].values; s = c5["composite"].values
    base = y.mean()
    auc = roc_auc_score(y, s)
    # deconfounded (within-epitope) signal: remove per-epitope means, then Spearman
    mc = s - c5.groupby("epitope")["composite"].transform("mean")
    bc = y - c5.groupby("epitope")["bound"].transform("mean")
    rho = stats.spearmanr(mc, bc).correlation
    print(f"C5 (unselected): {len(c5):,} designs, {int(y.sum())} binders, base rate {base:.1%}")
    print(f"composite predicts binding on UNSELECTED designs:")
    print(f"  pooled ROC-AUC              {auc:.3f}")
    print(f"  within-epitope Spearman     {rho:+.3f}  (deconfounded)")

    # titration curve: binding rate across composite quintiles
    c5["bin"] = pd.qcut(s, 5, labels=False, duplicates="drop")
    g = c5.groupby("bin").agg(rate=("bound","mean"), n=("bound","size"), lo=("composite","min"), hi=("composite","max"))
    print("\nbinding rate by composite quintile (unselected C5):")
    for b,row in g.iterrows():
        print(f"  Q{int(b)+1}  composite [{row.lo:.2f},{row.hi:.2f}]  n={int(row.n):4d}  bind {row.rate:.1%}")

    fig,ax=plt.subplots(figsize=(8.5,6))
    x=np.arange(len(g))
    ax.bar(x, g["rate"].values*100, color="#2a6f97", width=0.68, zorder=3)
    ax.axhline(base*100, ls="--", color="#c1502e", lw=1.8, zorder=2,
               label=f"C5 unselected base rate ({base:.1%})")
    for xi,r,n in zip(x, g["rate"].values, g["n"].values):
        ax.text(xi, r*100+0.4, f"{r:.0%}", ha="center", fontsize=10)
    ax.set_xticks(x); ax.set_xticklabels([f"Q{i+1}\n(low→high)" if i in (0,len(g)-1) else f"Q{i+1}" for i in range(len(g))])
    ax.set_xlabel("composite score quintile, on UNSELECTED designs (C5)", fontsize=12)
    ax.set_ylabel("binding rate (%)", fontsize=12)
    ax.set_title("The composite predicts binding on designs it never selected\n"
                 f"(C5 metric-space titration; pooled AUC {auc:.2f}, within-epitope ρ {rho:+.2f})",
                 fontsize=12.5, fontweight="bold")
    ax.legend(fontsize=10); ax.spines[["top","right"]].set_visible(False); ax.grid(axis="y",alpha=0.15)
    fig.tight_layout(); OUT.parent.mkdir(parents=True,exist_ok=True); fig.savefig(OUT,dpi=150)
    print(f"\nwrote {OUT.relative_to(ROOT)}")

if __name__ == "__main__":
    main()
