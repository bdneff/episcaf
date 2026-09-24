#!/usr/bin/env python3
"""The selection problem as real metric space: hundreds of thousands of designs, a tiny passing corner.

Each point is an AF3-scored antibody-arm design (C1/C2/8VDL) from dp4_superset.csv. Axes are the two
questions the filters ask: epitope RMSD (is the fold right?) and antibody clashes (is the epitope
accessible?). The 1,134 four-filter passers (0.45%) hug the low-RMSD / zero-clash corner; everything
else is wrong-fold, blocked, or both. A science-talk alternative to a funnel diagram.

Run:  /usr/bin/python3 scripts/dp4_selection_scatter.py
"""
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
SUP = ROOT/"data/libraries/dp4_superset.csv"; OUT = ROOT/"data/dp4_binding/figs/selection_scatter.png"
GREEN="#1f8a4c"; GREY="#9aa3ad"

def main():
    df = pd.read_csv(SUP, low_memory=False)
    ab = df[df.component.isin(["C1","C2","8VDL"])].copy()
    for c in ["epitope_rmsd","af3_clashes"]: ab[c]=pd.to_numeric(ab[c],errors="coerce")
    ab = ab.dropna(subset=["epitope_rmsd","af3_clashes"])
    gp = ab.is_global_pass.astype(str).str.lower().isin(["true","1","1.0"]).values
    rng = np.random.default_rng(0); jit=lambda v: v+rng.uniform(-0.35,0.35,len(v))
    nonp = ab[~gp].sample(min(22000,int((~gp).sum())), random_state=0); pas = ab[gp]

    fig,ax=plt.subplots(figsize=(9.2,6.4))
    ax.scatter(nonp.epitope_rmsd, np.clip(jit(nonp.af3_clashes),0,26), s=7, c=GREY, alpha=0.13, edgecolors="none", zorder=1)
    ax.scatter(pas.epitope_rmsd, jit(pas.af3_clashes), s=22, c=GREEN, alpha=0.95, edgecolors="white", linewidths=0.4, zorder=3)
    ax.set_xlim(0,4); ax.set_ylim(-0.6,26)
    ax.set_xlabel("epitope RMSD (Å)  —  is the fold right?", fontsize=13)
    ax.set_ylabel("antibody clashes  —  is the epitope accessible?", fontsize=13)
    ax.annotate(f"{len(pas):,} designs pass every filter\n(0.45% of ~253,000)", xy=(0.55,0.4), xytext=(1.7,7),
                fontsize=13, color=GREEN, fontweight="bold", arrowprops=dict(arrowstyle="-|>", color=GREEN, lw=2))
    ax.text(2.6,22,"every other design is\nwrong-fold, blocked, or both", fontsize=12, color="#6b7280", style="italic")
    ax.spines[["top","right"]].set_visible(False)
    fig.tight_layout(); fig.savefig(OUT, dpi=150)
    print(f"passers {len(pas):,} of {len(ab):,} = {gp.mean():.3%}\nwrote {OUT.relative_to(ROOT)}")

if __name__ == "__main__":
    main()
