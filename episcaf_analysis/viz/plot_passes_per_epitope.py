#!/usr/bin/env python3
"""
plot_passes_per_epitope.py  ->  manuscript/figures/passes_per_epitope.png

Four-filter passing designs per epitope target, RFD1+MPNN vs RFD3+MPNN (same DP3 epitopes).
Shows how concentrated the passes are (a few targets dominate) and how the two methods
compare target by target. Filters: epitope RMSD <= 1, overall <= 2, mean PAE < 5, 0 clashes.
"""
import numpy as np, pandas as pd
from pathlib import Path
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

KA   = Path("/Users/bneff/Desktop/projects/episcaf/known_antigen/analysis")
RFD1 = KA/"full_run/metrics_full_rfd1_mpnn_LAWSON.csv"
RFD3 = KA/"data/metrics_cylinder_full.csv"
OUT  = Path(__file__).resolve().parents[2]/"manuscript/figures/passes_per_epitope.png"

def num(s): return pd.to_numeric(s, errors="coerce")
def clash_count(x):
    if pd.isna(x): return np.nan
    return len(str(x).strip().strip("[]").split())

r1 = pd.read_csv(RFD1, usecols=["id","epitope_chunk_rmsd_vs_mpnn","overall_rmsd","mean_pae","af3_clash_resindices"], low_memory=False)
r3 = pd.read_csv(RFD3, usecols=["id","epitope_chunk_rmsd","overall_rmsd","mean_pae","af3_n_clash_res"], low_memory=False)

def passes_per_epitope(d, epi, clash):
    p = (num(d[epi])<=1)&(num(d.overall_rmsd)<=2)&(num(d.mean_pae)<5)&(clash==0)
    return d.assign(_p=p.astype(int)).groupby("id")._p.sum()

p1 = passes_per_epitope(r1, "epitope_chunk_rmsd_vs_mpnn", r1.af3_clash_resindices.map(clash_count))
p3 = passes_per_epitope(r3, "epitope_chunk_rmsd", num(r3.af3_n_clash_res))

tab = pd.concat([p1.rename("RFD1"), p3.rename("RFD3")], axis=1).fillna(0).astype(int)
tab = tab[(tab.RFD1>0)|(tab.RFD3>0)].sort_values(["RFD1","RFD3"], ascending=False)

x = np.arange(len(tab)); w = 0.42
fig, ax = plt.subplots(figsize=(max(9, 0.5*len(tab)), 5))
ax.bar(x-w/2, tab.RFD1, w, color="#1f77b4", label=f"RFD1+MPNN ({int(tab.RFD1.sum()):,} passes)")
ax.bar(x+w/2, tab.RFD3, w, color="#d62728", label=f"RFD3+MPNN ({int(tab.RFD3.sum()):,} passes)")
ax.set_xticks(x); ax.set_xticklabels(tab.index, rotation=60, ha="right", fontsize=8)
ax.set_ylabel("four-filter passing designs"); ax.set_xlabel("epitope target")
ax.set_title(f"Passing designs per epitope target, RFD1 vs RFD3 "
             f"({len(tab)} targets with any passes)", fontsize=12, fontweight="bold")
ax.legend(frameon=False);
for s in ("top","right"): ax.spines[s].set_visible(False)
fig.tight_layout(); fig.savefig(OUT, dpi=140, bbox_inches="tight")
print("wrote", OUT)
print(tab.to_string())
