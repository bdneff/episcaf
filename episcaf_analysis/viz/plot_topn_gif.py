#!/usr/bin/env python3
"""
plot_topn_gif.py  ->  scratch_figs/topn_sweep.gif

Animated version of the top-n overlay: one moving outline per DP4 row that sweeps the
per-epitope selection depth n = 1 .. 192 (the full budget) and closes onto the grey population,
so you watch selection relax instead of reading a stack of static traces. DP3 successes (red)
are the fixed reference in the top row. Same 8 metric columns as plot_topn_overlay.py.

Ranking is the UNGATED composite (rank, don't gate -- selection policy), rank per (antigen, id)
over all 192 designs, so at n=192 the outline lands exactly on grey.
  ranked input: scratch_figs/_scored_full_12mer_ranked.csv   (rank_full 1..192)
  built by: python -m episcaf_analysis.score --preset twelvemer --topk 999 ... then rank per (antigen,id)

NOTE: anim.save renders the figure canvas verbatim (no per-frame tight bbox), so left margin and
title width are reserved explicitly here -- do not rely on bbox_inches='tight'.
"""
import numpy as np, pandas as pd
from pathlib import Path
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
plt.rcParams.update({"font.size": 13, "axes.titlesize": 15, "xtick.labelsize": 8})

LOCAL = Path("/Users/bneff/Desktop/projects/episcaf")
DP3   = LOCAL/"known_antigen/analysis/data/metrics_native_cyl_full.csv"
RANK  = Path("scratch_figs/_scored_full_12mer_ranked.csv")
OUT   = Path("scratch_figs/topn_sweep.gif")

BLUE="#08519c"; DP3_C="#c0392b"; POP_C="#d9d9d9"
COLS=[("epitope_chunk_rmsd","Epitope\nRMSD (A)"),("mean_pae","Mean PAE\n(global)"),
      ("epitope_pae","Epitope\nPAE"),("overall_rmsd","Overall\nRMSD (A)"),("ptm","pTM"),
      ("af3_n_clash_res","AF3 clash\n(real Ab)"),
      ("cylinder_ca_clashes","Cylinder\n(plain)"),
      ("cylinder_native_aware","Cylinder\n(native-aware)")]

def num(s): return pd.to_numeric(s, errors="coerce")

dp3=pd.read_csv(DP3, low_memory=False)
ispass=((num(dp3.epitope_chunk_rmsd)<=1)&(num(dp3.overall_rmsd)<=2)
        &(num(dp3.mean_pae)<5)&(num(dp3.af3_n_clash_res)==0))
dp3p=dp3[ispass]
rk=pd.read_csv(RANK, low_memory=False)

# (label, animated-population df, static-pass df, is_dp3_reference)
rows=[("DP3 mAb\n(succeeded)", None, dp3p, True)]
for ag in ["1d2k","4wat","6m0j"]:
    rows.append((f"{ag.upper()}\n12mer", rk[rk.antigen==ag], None, False))

def edges(key,lo,hi):
    if key in ("cylinder_ca_clashes","cylinder_native_aware","af3_n_clash_res"):
        return np.arange(np.floor(lo),np.ceil(hi)+2,2)
    if key=="ptm": return np.linspace(lo,hi,26)
    return np.linspace(lo,hi,34)

xlim={}; EDG={}
for key,_ in COLS:
    parts=[num(rk[key]).dropna().values] if key in rk else []
    if key in dp3: parts.append(num(dp3[key]).dropna().values)
    v=np.concatenate(parts) if parts else np.array([0,1])
    lo=max(0.15,np.percentile(v,0.5)) if key=="ptm" else 0
    hi=np.percentile(v,99.5); xlim[key]=(lo,hi); EDG[key]=edges(key,lo,hi)

NS=sorted(set(list(range(1,21))+list(range(20,193,4))+[192]))
FRAMES=[1]*6+NS+[192]*10

fig,axes=plt.subplots(len(rows),len(COLS),figsize=(2.5*len(COLS),2.05*len(rows)+1.3))
# precompute grey population (static) and a FIXED per-panel y-max so the axis never rescales:
# grey stays put and the blue outline visibly SHRINKS onto it (peak is at the tightest n=1).
POP={}; YMAX={}
uniq_n=sorted(set(FRAMES))
for r,(name,alld,passd,ref) in enumerate(rows):
    for c,(key,_) in enumerate(COLS):
        if key=="af3_n_clash_res" and not ref: continue
        lo,hi=xlim[key]; e=EDG[key]; hmax=0.0
        if ref:
            hp,_=np.histogram(np.clip(num(passd[key]).dropna().values,lo,hi),bins=e,density=True)
            hmax=hp.max() if len(hp) else 1.0
        else:
            if key not in alld: continue
            POP[(r,c)]=np.clip(num(alld[key]).dropna().values,lo,hi)
            hg,_=np.histogram(POP[(r,c)],bins=e,density=True); hmax=max(hmax,hg.max())
            for n in uniq_n:                                        # tallest blue = tightest n
                xs=np.clip(num(alld[num(alld.rank_full)<=n][key]).dropna().values,lo,hi)
                if len(xs):
                    h,_=np.histogram(xs,bins=e,density=True); hmax=max(hmax,h.max())
        YMAX[(r,c)]=1.08*hmax if hmax>0 else 1.0

def draw(n):
    for r,(name,alld,passd,ref) in enumerate(rows):
        for c,(key,label) in enumerate(COLS):
            ax=axes[r,c]; ax.clear()
            if key=="af3_n_clash_res" and not ref:
                ax.text(0.5,0.5,"n/a\n(no antibody)",ha="center",va="center",style="italic",
                        color="#999",transform=ax.transAxes,fontsize=12)
                ax.set_xticks([]); ax.set_yticks([])
                [s.set_visible(False) for s in ax.spines.values()]
            else:
                lo,hi=xlim[key]; e=EDG[key]
                if ref:
                    xp=np.clip(num(passd[key]).dropna().values,lo,hi)
                    ax.hist(xp,bins=e,density=True,color=DP3_C,alpha=0.20,zorder=1)
                    ax.hist(xp,bins=e,density=True,histtype="step",color=DP3_C,lw=2.2,zorder=3)
                else:
                    ax.hist(POP[(r,c)],bins=e,density=True,color=POP_C,zorder=1)   # grey target
                    sel=alld[num(alld.rank_full)<=n]
                    xs=np.clip(num(sel[key]).dropna().values,lo,hi)
                    if len(xs):
                        ax.hist(xs,bins=e,density=True,histtype="step",color=BLUE,lw=2.4,zorder=4)
                ax.set_xlim(lo,hi); ax.set_ylim(0,YMAX[(r,c)])   # FIXED scale: grey fixed, blue shrinks
                ax.tick_params(axis='x',labelsize=8)
                for s in ("top","right","left"): ax.spines[s].set_visible(False)
            ax.set_yticks([])
            if r==0 and not (key=="af3_n_clash_res" and not ref):
                ax.set_title(label,fontweight="bold",fontsize=13)
            elif r==0:
                ax.set_title(label,fontweight="bold",fontsize=13)
            if c==0:
                ax.set_ylabel(name,color=(DP3_C if ref else "#111"),fontweight="bold",
                              rotation=0,ha="right",va="center",labelpad=12,fontsize=12)
    fig.suptitle(f"Selecting the top {n:>3d} designs / epitope\n"
                 "blue outline relaxes toward the full population (grey) as n grows to 192; "
                 "red = DP3 designs that succeeded",
                 fontsize=15,fontweight="bold")

handles=[Patch(fc=POP_C,ec="none",label="full population (all 192/epitope)"),
         Line2D([0],[0],color=BLUE,lw=2.4,label="DP4 selected: top n / epitope"),
         Line2D([0],[0],color=DP3_C,lw=2.2,label="DP3: designs that succeeded")]
fig.legend(handles=handles,loc="lower center",ncol=3,frameon=False,bbox_to_anchor=(0.5,0.005))

anim=FuncAnimation(fig,lambda f: draw(f),frames=FRAMES,interval=120)
fig.subplots_adjust(left=0.075,right=0.995,top=0.86,bottom=0.10,wspace=0.18,hspace=0.28)
anim.save(OUT,writer=PillowWriter(fps=9))
print("wrote",OUT,"frames",len(FRAMES),"cols",len(COLS))
