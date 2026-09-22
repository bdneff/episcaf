#!/usr/bin/env python3
"""Do feature interactions help at all? A clean, dependency-free test with interaction constraints.

The depth sweep conflates two things: a deeper tree is both more nonlinear per feature AND allows more
cross-feature interaction. Here we separate them. sklearn's HistGradientBoostingClassifier takes
interaction_cst='no_interactions', which forces every tree to split on a single feature, so the model is
purely ADDITIVE (a gradient-boosted GAM) yet can still be arbitrarily nonlinear per feature via depth.
We sweep depth with interactions FORBIDDEN vs ALLOWED, at matched flexibility.

Reading it: if 'additive' stays high while 'interactions allowed' falls as depth grows, then interactions
are not signal here, they are overfitting -- and no cleverer split ordering could help, because ordering
only matters once features interact. logistic regression (linear, additive) is the floor.

Also writes per-feature partial-dependence shapes for the additive model, so we can see the nonlinear
response the model learns for each feature (why an additive-nonlinear model beats the linear one).

Scoring: repeated stratified grouped CV by epitope. Same features as the ablation.

Run:  /usr/bin/python3 scripts/dp4_interaction_test.py
"""
from pathlib import Path
import numpy as np, pandas as pd
from sklearn.model_selection import StratifiedGroupKFold
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.inspection import PartialDependenceDisplay
from sklearn.metrics import roc_auc_score
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
LIB = ROOT / "data/libraries/dp4_library.csv"
SUMM = ROOT / "data/dp4_binding/john/scaffoldedEpitopeSummary.csv"
STRUCTF = ROOT / "data/dp4_binding/epitope_struct_features.csv"
OUT_PD = ROOT / "data/dp4_binding/figs/additive_shapes.png"

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

def hgb(depth, additive):
    kw = dict(class_weight="balanced", max_depth=depth, learning_rate=0.05, max_iter=300, random_state=0)
    if additive: kw["interaction_cst"] = "no_interactions"
    return make_pipeline(SimpleImputer(strategy="median"), HistGradientBoostingClassifier(**kw))

def main():
    d=load(); y=d.bound.values; g=d.epitope.values; X=d[F].values
    print(f"designs {len(d):,}  epitopes {d.epitope.nunique()}  bound {int(y.sum())}\n")
    a,s = cv_auc(lambda: make_pipeline(SimpleImputer(strategy="median"), StandardScaler(),
                 LogisticRegression(max_iter=3000, class_weight="balanced")), X, y, g)
    print(f"logistic regression (linear, additive):  {a:.3f}±{s:.3f}\n")
    print(f"{'max_depth':>10s} {'additive (no interactions)':>28s} {'interactions allowed':>22s}")
    for depth in [1,2,3,4,6,8]:
        aa,sa = cv_auc(lambda: hgb(depth, True),  X, y, g)
        af,sf = cv_auc(lambda: hgb(depth, False), X, y, g)
        print(f"{depth:>10d} {aa:>19.3f}±{sa:.3f} {af:>15.3f}±{sf:.3f}")

    # per-feature shapes from the additive model (fit on all data, for visualization only)
    best = hgb(4, True); best.fit(X, y)
    fig,axes = plt.subplots(3,5, figsize=(16,8))
    PartialDependenceDisplay.from_estimator(best, X, features=list(range(len(F))),
        feature_names=F, ax=axes.ravel()[:len(F)], line_kw={"color":"#2a6f97","lw":2})
    for ax in axes.ravel()[len(F):]: ax.set_visible(False)
    fig.suptitle("Per-feature response of the additive model (interactions forbidden)\n"
                 "each panel: how that one feature shifts predicted binding, holding it alone",
                 fontsize=12, fontweight="bold")
    fig.tight_layout(rect=[0,0,1,0.94]); OUT_PD.parent.mkdir(parents=True,exist_ok=True); fig.savefig(OUT_PD,dpi=140)
    print(f"\nwrote {OUT_PD.relative_to(ROOT)}")

if __name__ == "__main__":
    main()
