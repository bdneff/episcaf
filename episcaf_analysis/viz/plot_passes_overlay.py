#!/usr/bin/env python3
"""
plot_passes_overlay.py  ->  manuscript/figures/passes_overlay.png

Per metric, per row: the full design population in actual COUNTS (grey), with a dashed
line marking where the passing designs center, plus a zoomed inset showing the accepted
(passing) distribution's own shape and counts.
  DP3 passes = four-filter (epi<=1, overall<=2, mean_pae<5, af3_n_clash_res==0).
  DP4 passes = top-3 per epitope by composite (rank_in_epitope <= 3).

Counts (not density) so you can see how small and where the selected slice is; the inset
recovers the accepted shape, which would otherwise be invisible next to the full population.

Inputs (local copies off-cluster; on-cluster see configs/paths.py):
  DP3 all : metrics_cylinder_full.csv | DP4 all : metrics_12mer.csv
  DP4 passes : composite_12mer_top5_allscored.csv (rank_in_epitope)
"""
import numpy as np, pandas as pd
from pathlib import Path
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
plt.rcParams.update({"font.size": 14, "axes.titlesize": 16, "axes.labelsize": 15, "xtick.labelsize": 12, "ytick.labelsize": 12, "legend.fontsize": 12, "figure.titlesize": 18})  # paper-legible fonts

LOCAL = Path("/Users/bneff/Desktop/projects/episcaf")
DP3   = LOCAL/"known_antigen/analysis/data/metrics_native_cyl_full.csv"  # full table + native-aware carve (dp3_native_cylinder.sbatch)
DP4   = LOCAL/"12mer_tiling/analysis/data/metrics_12mer.csv"
DP4TOP= LOCAL/"12mer_tiling/analysis/data/composite_12mer_top5_allscored.csv"
OUT   = Path(__file__).resolve().parents[2]/"manuscript/figures/passes_overlay.png"

def num(s): return pd.to_numeric(s, errors="coerce")

dp3 = pd.read_csv(DP3, low_memory=False)
ispass = ((num(dp3.epitope_chunk_rmsd)<=1)&(num(dp3.overall_rmsd)<=2)
          &(num(dp3.mean_pae)<5)&(num(dp3.af3_n_clash_res)==0))
m12 = pd.read_csv(DP4, low_memory=False)
top3 = pd.read_csv(DP4TOP, low_memory=False); top3 = top3[num(top3.rank_in_epitope)<=3]

COLS=[("epitope_chunk_rmsd","Epitope\nRMSD (A)"),("mean_pae","Mean PAE\n(global)"),
      ("epitope_pae","Epitope\nPAE"),
      ("overall_rmsd","Overall\nRMSD (A)"),("ptm","pTM"),
      ("af3_n_clash_res","AF3 clash\n(real Ab)"),
      ("cylinder_ca_clashes","Cylinder\n(plain)"),
      ("cylinder_native_aware","Cylinder\n(native-aware)")]
rows=[("DP3 mAb", dp3, dp3[ispass], "#c0392b", True)]
for ag,c in [("1d2k","#1f77b4"),("4wat","#ff7f0e"),("6m0j","#2ca02c")]:
    rows.append((f"{ag.upper()} 12mer (DP4)", m12[m12.antigen==ag], top3[top3.antigen==ag], c, False))

# passing-set medians for the two cylinder columns (cited in the fig:overlay caption / sec:cylinterp)
for name,_,passd,_,_ in rows:
    p=num(passd.get("cylinder_ca_clashes")).median(); na=num(passd.get("cylinder_native_aware")).median()
    print(f"  {name:18s} passes: cylinder plain median={p:.1f}  native-aware median={na:.1f}")

xlim={}
for key,_ in COLS:
    v=[num(a[key]).dropna().values for _,a,_,_,hasab in rows if (key in a and not (key=="af3_n_clash_res" and not hasab))]
    v=np.concatenate(v) if v else np.array([0,1])
    lo=max(0.15,np.percentile(v,0.5)) if key=="ptm" else 0
    xlim[key]=(lo,np.percentile(v,99.5))

def edges(key,lo,hi):
    if key in ("cylinder_ca_clashes","cylinder_native_aware","af3_n_clash_res"): return np.arange(np.floor(lo),np.ceil(hi)+2,2)
    if key=="ptm": return np.linspace(lo,hi,26)
    return np.linspace(lo,hi,40)

fig,axes=plt.subplots(len(rows),len(COLS),figsize=(2.2*len(COLS),9.5))
for r,(name,alld,passd,color,hasab) in enumerate(rows):
    for c,(key,label) in enumerate(COLS):
        ax=axes[r,c]
        if key=="af3_n_clash_res" and not hasab:
            ax.text(0.5,0.5,"n/a\n(no antibody)",ha="center",va="center",style="italic",
                    color="#999",transform=ax.transAxes,fontsize=14)
            ax.set_xticks([]); ax.set_yticks([]); [s.set_visible(False) for s in ax.spines.values()]; continue
        lo,hi=xlim[key]; e=edges(key,lo,hi)
        xa=np.clip(num(alld[key]).dropna().values,lo,hi)
        xp=np.clip(num(passd[key]).dropna().values,lo,hi) if key in passd else np.array([])
        ax.hist(xa,bins=e,color="#cfcfcf")                       # all designs, LEFT axis (counts)
        ax.set_xlim(lo,hi); ax.tick_params(axis='x',labelsize=7); ax.tick_params(axis='y',labelsize=6)
        ax.locator_params(axis='y',nbins=3)
        if len(xp):
            ax2=ax.twinx()                                       # passes, RIGHT axis (counts), colored
            ax2.hist(xp,bins=e,color=color,alpha=0.55)
            ax2.hist(xp,bins=e,histtype="step",color=color,lw=1.4)
            ax2.axvline(np.median(xp),ls="--",color=color,lw=1.1)
            ax2.set_xlim(lo,hi); ax2.locator_params(axis='y',nbins=3)
            ax2.tick_params(axis='y',colors=color,labelsize=6)
            ax2.spines['right'].set_color(color)
        if r==0: ax.set_title(label,fontsize=17,fontweight="bold")
        if c==0:
            ax.set_ylabel(f"{name}\nall n={len(alld):,}\npass n={len(passd):,}",color=color,
                          fontweight="bold",rotation=0,ha="right",va="center",labelpad=46,fontsize=14)
fig.suptitle("All designs (grey, LEFT axis) vs passing designs (color, RIGHT colored axis); "
             "both true counts, dashed line = median of passes\n"
             "DP3 passes = four-filter; DP4 passes = top-3 per epitope",
             fontsize=18,fontweight="bold")
fig.tight_layout(rect=[0.03,0,1,0.95])
fig.savefig(OUT,dpi=130,bbox_inches="tight")
print("wrote",OUT)
print("DP3 all",len(dp3),"pass",int(ispass.sum()),
      "| DP4 top3:",{ag:int((top3.antigen==ag).sum()) for ag in ["1d2k","4wat","6m0j"]})
