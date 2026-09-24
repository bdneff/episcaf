#!/usr/bin/env python3
"""Selected known-antibody arm (C1+C2): bind rate vs pass rate for the three selectors.

Here Lawson's four-filter is a real comparator (it passes ~13% of the tested pool). We plot, for the
composite and the learned model, the binding rate among the designs kept (y) as a function of how many
we keep (x = pass rate), and place the four-filter as its single pass/fail point. A selector is better
if it reaches a higher bind rate at the same pass rate.

Model = logistic, out-of-fold under grouped CV by epitope, on design metrics + epitope native context.

Run:  /usr/bin/python3 scripts/dp4_arm_selected.py
"""
from pathlib import Path
import numpy as np, pandas as pd
from sklearn.model_selection import GroupKFold, cross_val_predict
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[1]
LIB=ROOT/"data/libraries/dp4_library.csv"; SUMM=ROOT/"data/dp4_binding/john/scaffoldedEpitopeSummary.csv"
STRUCTF=ROOT/"data/dp4_binding/epitope_struct_features.csv"; OUT=ROOT/"data/dp4_binding/figs/arm_c1c2_headtohead.png"
DESIGN=["cylinder_clashes","epitope_rmsd","overall_rmsd","epitope_pae","scaffold_pae","mean_pae","ptm","af3_clashes","island_index"]
EPI=["helix","strand","n_islands","epi_bfac_z","epi_gnm_dfit","epi_dsasa"]
BLUE="#2456E6"; ORANGE="#c1502e"; GREEN="#2C8C4A"; MUTE="#5b6570"; PURPLE="#7b3fa0"

def main():
    lib=pd.read_csv(LIB, low_memory=False)
    summ=pd.read_csv(SUMM, encoding="latin-1"); summ.columns=[c.strip() for c in summ.columns]
    summ=summ.rename(columns={"#islands":"n_islands"}); summ["epitope"]=summ.epitope.astype(str).str.lower()
    st=pd.read_csv(STRUCTF)
    epi=summ[["epitope","helix","strand","n_islands"]].drop_duplicates("epitope").merge(
        st[["epitope","epi_bfac_z","epi_gnm_dfit","epi_dsasa"]], on="epitope", how="inner")
    hits=set()
    for ids in summ["hitIDs"].dropna():
        for t in str(ids).split(","):
            t=t.strip()
            if t.startswith("DP4_"): hits.add("DP4_"+t.split("_")[1])
    lib["epitope"]=lib["target"].astype(str).str.split("_").str[0].str.lower()
    d=lib[lib.category.isin(["scaffoldedAbEpitope","scaffoldedSingleIsland"]) &
          lib.composite.notna() & lib.epitope.isin(set(epi.epitope))].merge(epi,on="epitope",how="left").copy()
    d["bound"]=d.library_member.isin(hits).astype(int); y=d.bound.values; base=y.mean()

    oof=cross_val_predict(make_pipeline(SimpleImputer(strategy="median"),StandardScaler(),
        LogisticRegression(max_iter=3000,class_weight="balanced")),
        d[DESIGN+EPI].values,y,groups=d.epitope.values,cv=GroupKFold(5),method="predict_proba")[:,1]
    comp=d.composite.values
    gp=(d.is_global_pass.astype(str).str.lower().isin(["true","1","1.0"])).values
    fp=gp.mean(); prec_ff=y[gp].mean()

    def prec(s, minn=30):
        o=np.argsort(-s); yy=y[o]; k=np.arange(1,len(yy)+1); pr=np.cumsum(yy)/k; m=k>=minn
        return k[m]/len(yy)*100, pr[m]*100
    print(f"C1/C2 selected: {len(d):,} designs, {int(y.sum())} binders, base {base:.1%}")
    print(f"four-filter: pass {fp:.0%}, bind {prec_ff:.0%}")
    xm,ym=prec(oof); xc,yc=prec(comp)
    fig,ax=plt.subplots(figsize=(9.4,6))
    ax.axhline(base*100,ls="--",color=MUTE,lw=1.8,label=f"random pick ({base:.0%})")
    ax.plot(xc,yc,color=ORANGE,lw=2.4,label="hand-set composite")
    ax.plot(xm,ym,color=BLUE,lw=2.9,label="logistic model (design + context)")
    ax.plot(fp*100,prec_ff*100,marker="*",ms=22,color=PURPLE,markeredgecolor="white",zorder=6,label="Lawson four-filter")
    ax.annotate(f"four-filter:\npass {fp:.0%}, bind {prec_ff:.0%}",(fp*100,prec_ff*100),
                xytext=(fp*100+8,prec_ff*100-9),fontsize=10.5,color=PURPLE,
                arrowprops=dict(arrowstyle="->",color=PURPLE))
    ax.set_xlabel("% of designs kept  (pass rate)",fontsize=13)
    ax.set_ylabel("binding rate among those kept  (%)",fontsize=13)
    ax.set_title("Known-antibody arm (C1+C2): the learned model beats\nthe composite and Lawson's four-filter at the same pass rate",fontsize=13,fontweight="bold")
    ax.legend(fontsize=11.5,loc="upper right"); ax.spines[["top","right"]].set_visible(False); ax.grid(alpha=0.15)
    ax.set_xlim(0,100); ax.set_ylim(0,None)
    fig.tight_layout(); fig.savefig(OUT,dpi=150); print(f"wrote {OUT.relative_to(ROOT)}")

if __name__=="__main__":
    main()
