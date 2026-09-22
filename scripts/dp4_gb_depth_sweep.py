#!/usr/bin/env python3
"""Do feature INTERACTIONS help the gradient booster? A direct test of whether split ordering matters.

Greedy tree induction is myopic, so in principle a non-greedy split ordering could do better -- but only
if the features actually interact. We test that head-on by sweeping tree depth. Depth 1 = decision stumps:
one feature per tree, NO interactions and no within-tree ordering at all (a pure additive model). Deeper
trees allow interactions and an internal split order. If depth 2/4/8 do not beat depth 1, then interactions
(and therefore any cleverer ordering) buy nothing for this data.

We also print plain logistic regression (zero interactions) as the reference: if the trees cannot beat it,
there is no interaction signal to exploit.

Scoring: repeated stratified grouped CV by epitope, same features as the ablation (design + 2 structure +
dynamics). Read the gap against the CV noise (the +/- column), not the third decimal.

Run:  /usr/bin/python3 scripts/dp4_gb_depth_sweep.py
"""
from pathlib import Path
import numpy as np, pandas as pd
from sklearn.model_selection import StratifiedGroupKFold
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.metrics import roc_auc_score

ROOT = Path(__file__).resolve().parents[1]
LIB = ROOT / "data/libraries/dp4_library.csv"
SUMM = ROOT / "data/dp4_binding/john/scaffoldedEpitopeSummary.csv"
STRUCTF = ROOT / "data/dp4_binding/epitope_struct_features.csv"

DESIGN = ["cylinder_clashes","epitope_rmsd","overall_rmsd","epitope_pae","scaffold_pae",
          "mean_pae","ptm","af3_clashes","island_index"]
EPI = ["helix","strand","n_islands","epi_bfac_z","epi_gnm_dfit","epi_dsasa"]
F = DESIGN + EPI

def load():
    lib = pd.read_csv(LIB, low_memory=False)
    summ = pd.read_csv(SUMM, encoding="latin-1"); summ.columns=[c.strip() for c in summ.columns]
    summ = summ.rename(columns={"#islands":"n_islands"}); summ["epitope"]=summ.epitope.astype(str).str.lower()
    st = pd.read_csv(STRUCTF)
    epi = summ[["epitope","helix","strand","n_islands"]].merge(
        st[["epitope","epi_bfac_z","epi_gnm_dfit","epi_dsasa"]], on="epitope", how="inner")
    hits=set()
    for ids in summ["hitIDs"].dropna():
        for t in str(ids).split(","):
            t=t.strip()
            if t.startswith("DP4_"): hits.add("DP4_"+t.split("_")[1])
    lib["epitope"]=lib["target"].astype(str).str.split("_").str[0].str.lower()
    d=lib[lib.composite.notna() & lib.epitope.isin(set(epi.epitope))].copy()
    d["bound"]=d.library_member.isin(hits).astype(int)
    return d.merge(epi, on="epitope", how="left")

def cv_auc(make, X, y, g, repeats=5):
    au=[]
    for s in range(repeats):
        oof=np.zeros(len(y))
        for tr,te in StratifiedGroupKFold(5, shuffle=True, random_state=s).split(X,y,g):
            m=make(); m.fit(X[tr],y[tr]); oof[te]=m.predict_proba(X[te])[:,1]
        au.append(roc_auc_score(y,oof))
    return np.mean(au), np.std(au)

def main():
    d=load(); y=d.bound.values; g=d.epitope.values; X=d[F].values
    print(f"designs {len(d):,}  epitopes {d.epitope.nunique()}  bound {int(y.sum())}\n")
    print(f"{'model':38s} {'ROC-AUC':>14s}")
    logit = lambda: make_pipeline(SimpleImputer(strategy="median"), StandardScaler(),
                                  LogisticRegression(max_iter=3000, class_weight="balanced"))
    a,s = cv_auc(logit, X, y, g)
    print(f"{'logistic regression (no interactions)':38s}  {a:.3f}±{s:.3f}")
    for depth in [1,2,3,4,6,8]:
        make = (lambda dep: lambda: make_pipeline(SimpleImputer(strategy="median"),
                HistGradientBoostingClassifier(class_weight="balanced", max_depth=dep,
                learning_rate=0.05, max_iter=300, random_state=0)))(depth)
        a,s = cv_auc(make, X, y, g)
        tag = "  (stumps: no interactions)" if depth==1 else ""
        print(f"{'grad-boost, max_depth='+str(depth):38s}  {a:.3f}±{s:.3f}{tag}")

if __name__ == "__main__":
    main()
