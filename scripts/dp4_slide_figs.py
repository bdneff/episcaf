#!/usr/bin/env python3
"""Slide-ready figures for the DP4 talk (clean, large-font, white background).

Produces:
  arm_binding_rates.png   per-arm hit rate, C5 as the unselected baseline
  epitope_context.png     per-epitope binding rate vs helix content (binding is epitope-specific)
  selection_funnel.png    250k designs -> 0.45% pass filter -> tested -> binders
  conformational_selection.png  concept: affinity vs the epitope's conformational distribution

Run:  /usr/bin/python3 scripts/dp4_slide_figs.py
"""
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Polygon

ROOT = Path(__file__).resolve().parents[1]
LIB = ROOT / "data/libraries/dp4_library.csv"
SUMM = ROOT / "data/dp4_binding/john/scaffoldedEpitopeSummary.csv"
STRUCTF = ROOT / "data/dp4_binding/epitope_struct_features.csv"
FIG = ROOT / "data/dp4_binding/figs"
BLUE, ORANGE, TEAL, GREY = "#1f6feb", "#c1502e", "#0e8f83", "#b8c0c9"
plt.rcParams.update({"font.size":13, "axes.spines.top":False, "axes.spines.right":False})

def load():
    lib = pd.read_csv(LIB, low_memory=False)
    summ = pd.read_csv(SUMM, encoding="latin-1"); summ.columns=[c.strip() for c in summ.columns]
    summ["epitope"]=summ.epitope.astype(str).str.lower()
    hits=set()
    for ids in summ["hitIDs"].dropna():
        for t in str(ids).split(","):
            t=t.strip()
            if t.startswith("DP4_"): hits.add("DP4_"+t.split("_")[1])
    lib["epitope"]=lib["target"].astype(str).str.split("_").str[0].str.lower()
    lib["bound"]=lib["library_member"].isin(hits).astype(int)
    return lib, summ

def fig_arms(lib):
    order=[("scaffoldedAbEpitope","C1\nwhole-epitope\n(top-ranked)"),
           ("metricSpaceTitration","C5\nmetric titration\n(unselected)"),
           ("scaffoldedEpitopeControl","C6\ncontrols"),
           ("scaffoldedSingleIsland","C2\nsingle-island")]
    rates=[lib[lib.category==k].bound.mean()*100 for k,_ in order]
    labs=[l for _,l in order]
    base=lib[lib.category=="metricSpaceTitration"].bound.mean()*100
    fig,ax=plt.subplots(figsize=(9,6))
    cols=[BLUE, TEAL, GREY, GREY]
    ax.bar(range(4),rates,color=cols,width=0.66,zorder=3)
    ax.axhline(base,ls="--",color=TEAL,lw=1.8,zorder=2)
    ax.text(3.4,base+0.3,f"unselected baseline {base:.1f}%",color=TEAL,ha="right",fontsize=11)
    for i,r in enumerate(rates): ax.text(i,r+0.3,f"{r:.1f}%",ha="center",fontsize=13,fontweight="bold")
    ax.set_xticks(range(4)); ax.set_xticklabels(labs,fontsize=11)
    ax.set_ylabel("binding rate (%)")
    ax.set_title("Selecting good designs beats the unselected background",fontsize=14,fontweight="bold")
    ax.grid(axis="y",alpha=0.15); fig.tight_layout(); fig.savefig(FIG/"arm_binding_rates.png",dpi=150); plt.close(fig)

def fig_epitope_context(lib, summ):
    d=lib[lib.composite.notna()].copy()
    st=pd.read_csv(STRUCTF)
    d=d[d.epitope.isin(set(st.epitope))]
    g=d.groupby("epitope").agg(rate=("bound","mean"),n=("bound","size")).reset_index()
    g=g.merge(summ[["epitope","helix"]].drop_duplicates("epitope"),on="epitope",how="left")
    fig,ax=plt.subplots(figsize=(9,6))
    ax.scatter(g.helix,g.rate*100,s=np.clip(g.n,20,300),c=BLUE,alpha=0.6,edgecolors="white",zorder=3)
    ax.set_xlabel("epitope helix content (native structure)")
    ax.set_ylabel("binding rate of its designs (%)")
    ax.set_title("Binding is epitope-specific: helical epitopes scaffold into binders",fontsize=13.5,fontweight="bold")
    ax.grid(alpha=0.15); fig.tight_layout(); fig.savefig(FIG/"epitope_context.png",dpi=150); plt.close(fig)

def fig_funnel():
    stages=[("candidate scaffolds\n(AF3-scored)",253318),("pass in-silico filter\n(0.45%)",1134),
            ("tested in the assay\n(antibody arms)",2720),("actually bind",217)]
    fig,ax=plt.subplots(figsize=(9,6))
    maxw=1.0
    for i,(lab,n) in enumerate(stages):
        w=maxw*(np.log10(n)/np.log10(253318))
        y=len(stages)-1-i
        ax.add_patch(FancyBboxPatch((0.5-w/2,y-0.34),w,0.68,boxstyle="round,pad=0.01",
                     fc=[BLUE,TEAL,ORANGE,"#2a9d3a"][i],ec="white",zorder=3))
        ax.text(0.5,y,f"{n:,}",ha="center",va="center",color="white",fontweight="bold",fontsize=15,zorder=4)
        ax.text(1.02,y,lab,ha="left",va="center",fontsize=12)
    ax.set_xlim(0,1.7); ax.set_ylim(-0.6,len(stages)-0.4); ax.axis("off")
    ax.set_title("The selection problem: which designs do we make?",fontsize=14,fontweight="bold",loc="left")
    fig.tight_layout(); fig.savefig(FIG/"selection_funnel.png",dpi=150); plt.close(fig)

def fig_conformational():
    fig,axes=plt.subplots(1,3,figsize=(13,4.6))
    titles=["unstructured peptide","scaffolded, wrong shape","scaffolded, native shape"]
    aff=["low affinity","no affinity","high affinity"]
    x=np.linspace(0,10,400); target=7.0
    curves=[0.5+0.05*np.sin(x), 0.15+1.6*np.exp(-((x-3)**2)/0.6), 0.15+1.9*np.exp(-((x-target)**2)/0.5)]
    for ax,t,a,c in zip(axes,titles,aff,curves):
        ax.plot(x,c,color=ORANGE,lw=3)
        ax.axvline(target,ls=":",color=BLUE,lw=2)
        ax.fill_between(x,c,where=(np.abs(x-target)<0.25),color=BLUE,alpha=0.5)
        ax.set_title(t,fontsize=13,fontweight="bold"); ax.text(0.5,-0.22,a,transform=ax.transAxes,ha="center",fontsize=12,color="0.3")
        ax.set_xticks([]); ax.set_yticks([]); ax.set_ylim(0,2.2)
        ax.set_xlabel("epitope conformation",fontsize=10)
    axes[0].set_ylabel("relative occupancy",fontsize=11)
    fig.suptitle("Affinity depends on how often the epitope is in its antibody-binding shape (dotted line)",
                 fontsize=13.5,fontweight="bold")
    fig.tight_layout(rect=[0,0.02,1,0.94]); fig.savefig(FIG/"conformational_selection.png",dpi=150); plt.close(fig)

def main():
    lib,summ=load()
    fig_arms(lib); fig_epitope_context(lib,summ); fig_funnel(); fig_conformational()
    for f in ["arm_binding_rates","epitope_context","selection_funnel","conformational_selection"]:
        print("wrote", (FIG/(f+".png")).relative_to(ROOT))

if __name__ == "__main__":
    main()
