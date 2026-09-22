#!/usr/bin/env python3
"""Order-independent contribution of each feature BLOCK, for gradient boosting.

The ablation adds blocks in one fixed order (design -> secondary structure -> dynamics), so each block's
apparent value depends on what is already in. The order-free answer is the Shapley value: average each
block's marginal ROC-AUC gain over ALL orderings of the three blocks. Three blocks means only 2^3 = 8
subsets to score, so this is cheap and exhaustive. Feature COLUMN order inside the model is not searched,
because gradient boosting is invariant to it (splits are chosen by gain over all features).

Scoring: repeated stratified grouped CV by epitope, same as the ablation table. Empty set = 0.5.

Run:  /usr/bin/python3 scripts/dp4_feature_shapley.py
"""
from itertools import combinations
from math import factorial
from pathlib import Path
import numpy as np, pandas as pd
from sklearn.model_selection import StratifiedGroupKFold
from sklearn.pipeline import make_pipeline
from sklearn.impute import SimpleImputer
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.metrics import roc_auc_score

ROOT = Path(__file__).resolve().parents[1]
LIB = ROOT / "data/libraries/dp4_library.csv"
SUMM = ROOT / "data/dp4_binding/john/scaffoldedEpitopeSummary.csv"
STRUCTF = ROOT / "data/dp4_binding/epitope_struct_features.csv"

DESIGN = ["cylinder_clashes","epitope_rmsd","overall_rmsd","epitope_pae","scaffold_pae",
          "mean_pae","ptm","af3_clashes","island_index"]
SS  = ["helix","strand","n_islands"]
DYN = ["epi_bfac_z","epi_gnm_dfit","epi_dsasa"]
BLOCKS = {"design": DESIGN, "2° structure": SS, "dynamics": DYN}

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

def auc_of(feats, d, y, g, repeats=3):
    if not feats: return 0.5
    X = d[feats].values
    aus=[]
    for s in range(repeats):
        oof=np.zeros(len(y))
        for tr,te in StratifiedGroupKFold(5, shuffle=True, random_state=s).split(X,y,g):
            m=make_pipeline(SimpleImputer(strategy="median"),
                            HistGradientBoostingClassifier(class_weight="balanced", max_depth=3,
                            learning_rate=0.05, max_iter=300, random_state=0))
            m.fit(X[tr],y[tr]); oof[te]=m.predict_proba(X[te])[:,1]
        aus.append(roc_auc_score(y,oof))
    return float(np.mean(aus))

def main():
    d=load(); y=d.bound.values; g=d.epitope.values
    names=list(BLOCKS)
    # score all 8 subsets once
    val={}
    for r in range(len(names)+1):
        for combo in combinations(names, r):
            feats=[c for b in combo for c in BLOCKS[b]]
            val[frozenset(combo)]=auc_of(feats, d, y, g)
    print("AUC of every feature-block subset:")
    for k in sorted(val, key=lambda s:(len(s), sorted(s))):
        print(f"  {('{'+', '.join(sorted(k))+'}') if k else '{} (none)':32s} {val[k]:.3f}")
    # Shapley value of each block
    n=len(names); shap={}
    for b in names:
        s=0.0
        for r in range(n):
            for combo in combinations([x for x in names if x!=b], r):
                w=factorial(len(combo))*factorial(n-len(combo)-1)/factorial(n)
                s+=w*(val[frozenset(combo)|{b}]-val[frozenset(combo)])
        shap[b]=s
    print("\nShapley value (average marginal ROC-AUC gain over ALL orderings):")
    for b in names:
        print(f"  {b:16s} {shap[b]:+.3f}")
    print(f"  {'sum + 0.5':16s} {0.5+sum(shap.values()):.3f}   (= full-model AUC {val[frozenset(names)]:.3f})")

if __name__ == "__main__":
    main()
