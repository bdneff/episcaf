#!/usr/bin/env python3
"""
plot_passes_overlay.py  ->  manuscript/figures/passes_overlay.png

All designs (grey fill) vs passing designs (colored outline), per metric, per row.
  DP3 passes = four-filter (epi<=1, overall<=2, mean_pae<5, af3_n_clash_res==0).
  DP4 passes = top-3 per epitope by composite (rank_in_epitope <= 3).
Answers: how do experimentally-relevant *passing* designs look on each metric,
especially the cylinder clash, for the known-antibody DP3 set vs the tiled-12mer DP4 set.

Inputs (local copies off-cluster; on-cluster see configs/paths.py):
  DP3 all     : metrics_cylinder_full.csv
  DP4 all     : metrics_12mer.csv
  DP4 passes  : composite_12mer_top5_allscored.csv  (has rank_in_epitope)
"""
import numpy as np, pandas as pd
from pathlib import Path
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

LOCAL = Path("/Users/bneff/Desktop/projects/episcaf")
DP3   = LOCAL/"known_antigen/analysis/data/metrics_cylinder_full.csv"
DP4   = LOCAL/"12mer_tiling/analysis/data/metrics_12mer.csv"
DP4TOP= LOCAL/"12mer_tiling/analysis/data/composite_12mer_top5_allscored.csv"
OUT   = Path(__file__).resolve().parents[2]/"manuscript/figures/passes_overlay.png"

def num(s): return pd.to_numeric(s, errors="coerce")

dp3 = pd.read_csv(DP3, low_memory=False)
ispass = ((num(dp3.epitope_chunk_rmsd)<=1)&(num(dp3.overall_rmsd)<=2)
          &(num(dp3.mean_pae)<5)&(num(dp3.af3_n_clash_res)==0))
m12 = pd.read_csv(DP4, low_memory=False)
top3 = pd.read_csv(DP4TOP, low_memory=False)
top3 = top3[num(top3.rank_in_epitope)<=3]

COLS=[("epitope_chunk_rmsd","Epitope RMSD (A)"),("epitope_pae","Epitope PAE"),
      ("overall_rmsd","Overall RMSD (A)"),("ptm","pTM"),
      ("af3_n_clash_res","AF3 clashing res\n(real antibody)"),("cylinder_ca_clashes","Cylinder clashes")]
rows=[("DP3 mAb", dp3, dp3[ispass], "#1a1a1a", True)]
for ag,c in [("1d2k","#1f77b4"),("4wat","#ff7f0e"),("6m0j","#2ca02c")]:
    rows.append((f"{ag.upper()} 12mer (DP4)", m12[m12.antigen==ag], top3[top3.antigen==ag], c, False))

xlim={}
for key,_ in COLS:
    v=[num(a[key]).dropna().values for _,a,_,_,hasab in rows if (key in a and not (key=="af3_n_clash_res" and not hasab))]
    v=np.concatenate(v) if v else np.array([0,1])
    lo=max(0.15,np.percentile(v,0.5)) if key=="ptm" else 0
    xlim[key]=(lo,np.percentile(v,99.5))

def edges(key,lo,hi):
    if key in ("cylinder_ca_clashes","af3_n_clash_res"): return np.arange(np.floor(lo),np.ceil(hi)+2,2)
    if key=="ptm": return np.linspace(lo,hi,26)
    return np.linspace(lo,hi,40)

fig,axes=plt.subplots(len(rows),len(COLS),figsize=(20,11.5))
for r,(name,alld,passd,color,hasab) in enumerate(rows):
    for c,(key,label) in enumerate(COLS):
        ax=axes[r,c]
        if key=="af3_n_clash_res" and not hasab:
            ax.text(0.5,0.5,"n/a\n(no antibody)",ha="center",va="center",style="italic",
                    color="#999",transform=ax.transAxes,fontsize=10)
            ax.set_xticks([]); ax.set_yticks([]); [s.set_visible(False) for s in ax.spines.values()]; continue
        lo,hi=xlim[key]; e=edges(key,lo,hi)
        ax.hist(np.clip(num(alld[key]).dropna().values,lo,hi),bins=e,density=True,color="#dddddd",label="all designs")
        xp=np.clip(num(passd[key]).dropna().values,lo,hi) if key in passd else np.array([])
        if len(xp):
            ax.hist(xp,bins=e,density=True,histtype="step",color=color,lw=2,label="passes")
            ax.axvline(np.median(xp),ls="--",color=color,lw=1.1)
        ax.set_xlim(lo,hi); ax.set_yticks([])
        if r==0: ax.set_title(label,fontsize=12,fontweight="bold")
        if c==0:
            ax.set_ylabel(f"{name}\nall n={len(alld):,}\npass n={len(passd):,}",color=color,
                          fontweight="bold",rotation=0,ha="right",va="center",labelpad=46,fontsize=10)
        if r==0 and c==len(COLS)-1: ax.legend(frameon=False,fontsize=9,loc="upper right")
fig.suptitle("All designs (grey) vs passing designs (color): DP3 four-filter passes, DP4 top-3 per epitope\n"
             "density-normalized; dashed line marks the median of passes",fontsize=13,fontweight="bold")
fig.tight_layout(rect=[0.03,0,1,0.95])
fig.savefig(OUT,dpi=130,bbox_inches="tight")
print("wrote",OUT)
print("DP3 all",len(dp3),"pass",int(ispass.sum()),
      "| DP4 top3:",{ag:int((top3.antigen==ag).sum()) for ag in ["1d2k","4wat","6m0j"]})
