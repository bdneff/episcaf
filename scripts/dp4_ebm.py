#!/usr/bin/env python3
"""Explainable Boosting Machine: the glassbox test of whether ANY feature interactions help.

An EBM is a gradient-boosted Generalized Additive Model. With interactions=0 it is purely additive (one
learned shape per feature). With interactions=k it also fits the k pairwise terms that most reduce error
and reports each one's importance. So it answers our question two ways at once: does adding interaction
terms raise held-out AUC, and how large are the interaction terms' importances next to the main effects.

Needs the `interpret` package (throwaway venv), so run with that venv's python:
  <venv>/bin/python scripts/dp4_ebm.py

Scoring: repeated stratified grouped CV by epitope. Same features as the ablation.
"""
from pathlib import Path
import numpy as np, pandas as pd
from sklearn.model_selection import StratifiedGroupKFold
from sklearn.metrics import roc_auc_score
from interpret.glassbox import ExplainableBoostingClassifier

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

def cv_auc(n_inter, Xdf, y, g, repeats=3):
    au=[]
    for s in range(repeats):
        oof=np.zeros(len(y))
        for tr,te in StratifiedGroupKFold(5, shuffle=True, random_state=s).split(Xdf, y, g):
            ebm = ExplainableBoostingClassifier(interactions=n_inter, feature_names=F, random_state=0)
            ebm.fit(Xdf.iloc[tr], y[tr]); oof[te]=ebm.predict_proba(Xdf.iloc[te])[:,1]
        au.append(roc_auc_score(y,oof))
    return np.mean(au), np.std(au)

def main():
    d=load(); y=d.bound.values; g=d.epitope.values
    Xdf = d[F].apply(pd.to_numeric, errors="coerce")
    print(f"designs {len(d):,}  epitopes {d.epitope.nunique()}  bound {int(y.sum())}\n")
    for n_inter, tag in [(0,"additive only (pure GAM)"), (10,"+ up to 10 pairwise interactions")]:
        a,s = cv_auc(n_inter, Xdf, y, g)
        print(f"EBM, {tag:34s}  ROC-AUC {a:.3f}±{s:.3f}")

    # full-data fit to read the importances: how big are interaction terms vs main effects?
    ebm = ExplainableBoostingClassifier(interactions=10, feature_names=F, random_state=0).fit(Xdf, y)
    imp = dict(zip(ebm.term_names_, ebm.term_importances()))
    mains = {k:v for k,v in imp.items() if " & " not in k}
    inters = {k:v for k,v in imp.items() if " & " in k}
    print(f"\ntotal importance -- main effects {sum(mains.values()):.3f}   "
          f"interaction terms {sum(inters.values()):.3f}   "
          f"(interactions are {100*sum(inters.values())/sum(imp.values()):.1f}% of the model)")
    print("\ntop main effects:")
    for k,v in sorted(mains.items(), key=lambda kv:-kv[1])[:6]: print(f"  {k:20s} {v:.3f}")
    print("\npairwise interaction terms the EBM chose (importance):")
    for k,v in sorted(inters.items(), key=lambda kv:-kv[1]): print(f"  {k:32s} {v:.3f}")

if __name__ == "__main__":
    main()
