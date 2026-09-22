#!/usr/bin/env python3
"""Improved binder-model evaluation, folding in the day's lessons.

Changes from dp4_binder_ml.py:
 - add the native-context structural features from the RCSB complexes (epitope size, native exposure,
   radius of gyration, composition) on top of design metrics + secondary structure;
 - feature ABLATION (design only -> + secondary structure -> + native context) to see what each adds;
 - REPEATED stratified grouped CV (multiple shuffled epitope partitions) for stable mean +/- sd,
   since a single held-out split was noisy;
 - lean architectures only (composite baseline, logistic, gradient boosting, right-sized MLP); the
   oversized MLP is dropped after we confirmed it overfits.

Honest evaluation stays grouped by epitope (test antibodies never seen in training).

Run:  /usr/bin/python3 scripts/dp4_binder_ml_v2.py
"""
from pathlib import Path
import numpy as np, pandas as pd
from sklearn.model_selection import StratifiedGroupKFold
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import HistGradientBoostingClassifier, RandomForestClassifier
from sklearn.neural_network import MLPClassifier
from sklearn.metrics import roc_auc_score, average_precision_score
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
LIB = ROOT / "data/libraries/dp4_library.csv"
SUMM = ROOT / "data/dp4_binding/john/scaffoldedEpitopeSummary.csv"
STRUCTF = ROOT / "data/dp4_binding/epitope_struct_features.csv"
OUT = ROOT / "data/dp4_binding/figs/binder_ml_v2.png"

DESIGN = ["cylinder_clashes","epitope_rmsd","overall_rmsd","epitope_pae","scaffold_pae",
          "mean_pae","ptm","af3_clashes","island_index"]
SS = ["helix","strand","n_islands"]                       # drop loop (helix+strand+loop=1, redundant)
NATIVE = ["epi_size","epi_rel_sasa","epi_rg","epi_frac_hydrophobic","epi_frac_charged","epi_frac_aromatic"]
DYN = ["epi_bfac_z","epi_gnm_dfit","epi_dsasa"]   # holo rigidity, Δflex, ΔSASA (the new features that separated)
SETS = {"design metrics": DESIGN, "+ secondary structure": DESIGN+SS, "+ dynamics": DESIGN+SS+DYN}

def load():
    lib = pd.read_csv(LIB, low_memory=False)
    summ = pd.read_csv(SUMM, encoding="latin-1"); summ.columns=[c.strip() for c in summ.columns]
    summ = summ.rename(columns={"#islands":"n_islands"}); summ["epitope"]=summ.epitope.astype(str).str.lower()
    st = pd.read_csv(STRUCTF)
    epi = summ[["epitope","helix","strand","n_islands"]].merge(st[["epitope"]+DYN], on="epitope", how="inner")
    hits=set()
    for ids in summ["hitIDs"].dropna():
        for t in str(ids).split(","):
            t=t.strip()
            if t.startswith("DP4_"): hits.add("DP4_"+t.split("_")[1])
    lib["epitope"]=lib["target"].astype(str).str.split("_").str[0].str.lower()
    d=lib[lib.composite.notna() & lib.epitope.isin(set(epi.epitope))].copy()
    d["bound"]=d.library_member.isin(hits).astype(int)
    return d.merge(epi, on="epitope", how="left")

def repeated_cv(make, X, y, groups, repeats=3):
    au, ap = [], []
    for s in range(repeats):
        oof=np.zeros(len(y))
        for tr,te in StratifiedGroupKFold(5, shuffle=True, random_state=s).split(X,y,groups):
            m=make(); m.fit(X[tr],y[tr]); oof[te]=m.predict_proba(X[te])[:,1]
        au.append(roc_auc_score(y,oof)); ap.append(average_precision_score(y,oof))
    return np.array(au), np.array(ap)

def main():
    d=load(); y=d.bound.values; g=d.epitope.values; base=y.mean()
    print(f"designs={len(d):,} bound={int(y.sum())} epitopes={d.epitope.nunique()} base rate={base:.3f}")
    pipe = lambda est: make_pipeline(SimpleImputer(strategy="median"), StandardScaler(), est)
    factories = {
        "logistic":   lambda: pipe(LogisticRegression(max_iter=3000, class_weight="balanced")),
        "grad-boost": lambda: make_pipeline(SimpleImputer(strategy="median"),
                        HistGradientBoostingClassifier(class_weight="balanced", max_depth=3,
                        learning_rate=0.05, max_iter=300, random_state=0)),
        "random-forest": lambda: make_pipeline(SimpleImputer(strategy="median"),
                        RandomForestClassifier(400, min_samples_leaf=20, class_weight="balanced", random_state=0)),
        "MLP (32,16)":lambda: pipe(MLPClassifier(hidden_layer_sizes=(32,16), alpha=1e-3,
                        early_stopping=True, max_iter=800, random_state=0)),
    }
    ac = roc_auc_score(y, d.composite); pc = average_precision_score(y, d.composite)
    print(f"\nBASELINE composite: ROC-AUC {ac:.3f}  PR-AUC {pc:.3f}\n")
    res = {}   # (model,set) -> (auc_mean,auc_sd,ap_mean,ap_sd)
    print(f"{'model':14s} {'feature set':16s} {'ROC-AUC':>13s} {'PR-AUC':>13s}")
    for name, make in factories.items():
        for sname, feats in SETS.items():
            au, ap = repeated_cv(make, d[feats].values, y, g)
            res[(name,sname)] = (au.mean(),au.std(),ap.mean(),ap.std())
            print(f"{name:14s} {sname:16s}  {au.mean():.3f}±{au.std():.3f}   {ap.mean():.3f}±{ap.std():.3f}")

    # ---- figure: AUC by feature set, grouped by model ----
    fig,(a1,a2)=plt.subplots(1,2,figsize=(13,5))
    setnames=list(SETS); models=list(factories); x=np.arange(len(setnames)); w=0.8/len(models)
    cols={"logistic":"#4a8bb5","grad-boost":"#1b3a4b","random-forest":"#0e8f83","MLP (32,16)":"#c1502e"}
    off=lambda i:(i-(len(models)-1)/2)*w
    for i,m in enumerate(models):
        a1.bar(x+off(i),[res[(m,s)][0] for s in setnames],w,yerr=[res[(m,s)][1] for s in setnames],
               capsize=3,color=cols[m],label=m)
    a1.axhline(ac,ls="--",c="0.5",lw=1.2); a1.text(len(setnames)-1,ac+0.004,f"composite {ac:.2f}",fontsize=8,ha="right",color="0.4")
    a1.set_xticks(x); a1.set_xticklabels(setnames); a1.set_ylim(0.50,0.86)
    a1.set_ylabel("ROC-AUC (repeated grouped CV)"); a1.set_title("Does each feature block help?")
    a1.legend(fontsize=8); a1.spines[["top","right"]].set_visible(False)
    for i,m in enumerate(models):
        a2.bar(x+off(i),[res[(m,s)][2] for s in setnames],w,yerr=[res[(m,s)][3] for s in setnames],
               capsize=3,color=cols[m])
    a2.axhline(pc,ls="--",c="0.5",lw=1.2); a2.axhline(base,ls=":",c="0.7",lw=1); a2.text(0,base+0.005,f"chance {base:.2f}",fontsize=8,color="0.5")
    a2.set_xticks(x); a2.set_xticklabels(setnames); a2.set_ylabel("PR-AUC"); a2.set_title("Precision-recall (8% base rate)")
    a2.spines[["top","right"]].set_visible(False)
    fig.suptitle(f"DP4 binder model v2 — feature ablation × architecture (grouped CV, {d.epitope.nunique()} epitopes)",
                 fontsize=12, fontweight="bold")
    fig.tight_layout(rect=[0,0,1,0.95]); OUT.parent.mkdir(parents=True,exist_ok=True); fig.savefig(OUT,dpi=150)
    print(f"\nwrote {OUT.relative_to(ROOT)}")

if __name__ == "__main__":
    main()
