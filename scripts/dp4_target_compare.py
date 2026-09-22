#!/usr/bin/env python3
"""Binary label vs continuous pseudo-affinity: which target trains a better binder predictor?

Same designs, same features, grouped CV by epitope. One model is a CLASSIFIER on hit/no-hit; the
other a REGRESSOR on the dilution-series pseudo-affinity (our pseudo-dG). We compare them on the same
footing by ranking designs with each model's output and scoring against the binary hit label (AUC),
and we also report how well the regressor recovers the graded affinity (Spearman).

Universe: cognate designs of the placed antibodies (those with a pseudo-affinity) that also carry
design metrics. Small (~25 epitopes), so read the gap, not the third decimal.

Features: design metrics + secondary structure + the dynamics/binding-change features
(holo B-factor rigidity, Δflex, ΔSASA).

Run:  /usr/bin/python3 scripts/dp4_target_compare.py
"""
from pathlib import Path
import numpy as np, pandas as pd
from scipy import stats
from sklearn.model_selection import GroupKFold, cross_val_predict
from sklearn.ensemble import HistGradientBoostingClassifier, HistGradientBoostingRegressor
from sklearn.metrics import roc_auc_score, average_precision_score, roc_curve
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
LIB = ROOT / "data/libraries/dp4_library.csv"
SUMM = ROOT / "data/dp4_binding/john/scaffoldedEpitopeSummary.csv"
STRUCTF = ROOT / "data/dp4_binding/epitope_struct_features.csv"
AFF = ROOT / "data/dp4_binding/dp4_pseudoaffinity.csv"
OUT = ROOT / "data/dp4_binding/figs/target_compare.png"

DESIGN = ["cylinder_clashes","epitope_rmsd","overall_rmsd","epitope_pae","scaffold_pae",
          "mean_pae","ptm","af3_clashes","island_index"]
EPI = ["helix","strand","n_islands","epi_bfac_z","epi_gnm_dfit","epi_dsasa"]

def main():
    aff = pd.read_csv(AFF)
    lib = pd.read_csv(LIB, low_memory=False)
    summ = pd.read_csv(SUMM, encoding="latin-1"); summ.columns=[c.strip() for c in summ.columns]
    summ = summ.rename(columns={"#islands":"n_islands"}); summ["epitope"]=summ.epitope.str.lower()
    st = pd.read_csv(STRUCTF)
    epi = summ[["epitope","helix","strand","n_islands"]].merge(
        st[["epitope","epi_bfac_z","epi_gnm_dfit","epi_dsasa"]], on="epitope", how="inner")

    d = (aff.merge(lib[["library_member"]+DESIGN], on="library_member", how="inner")
            .merge(epi, on="epitope", how="inner"))
    d = d[d["cylinder_clashes"].notna()].copy()
    y_bin = d["is_hit"].astype(int).values
    y_reg = d["pseudo_affinity"].values
    g = d["epitope"].values
    X = d[DESIGN+EPI].values
    print(f"designs {len(d):,}  epitopes {d.epitope.nunique()}  hits {y_bin.sum()}  base rate {y_bin.mean():.3f}")

    cv = GroupKFold(5)
    clf = HistGradientBoostingClassifier(class_weight="balanced", max_depth=3, learning_rate=0.05, max_iter=300, random_state=0)
    reg = HistGradientBoostingRegressor(max_depth=3, learning_rate=0.05, max_iter=300, random_state=0)
    p_clf = cross_val_predict(clf, X, y_bin, groups=g, cv=cv, method="predict_proba")[:,1]
    p_reg = cross_val_predict(reg, X, y_reg, groups=g, cv=cv)

    auc_clf = roc_auc_score(y_bin, p_clf); ap_clf = average_precision_score(y_bin, p_clf)
    auc_reg = roc_auc_score(y_bin, p_reg); ap_reg = average_precision_score(y_bin, p_reg)
    rho_reg = stats.spearmanr(p_reg, y_reg).correlation
    print("\ntrained on BINARY hit label (classifier):")
    print(f"   ranks binders: ROC-AUC {auc_clf:.3f}   PR-AUC {ap_clf:.3f}")
    print("trained on CONTINUOUS pseudo-affinity (regressor):")
    print(f"   ranks binders: ROC-AUC {auc_reg:.3f}   PR-AUC {ap_reg:.3f}")
    print(f"   recovers graded affinity: Spearman(pred, affinity) {rho_reg:+.2f}")
    print(f"\n=> {'regressor (dG)' if auc_reg>auc_clf else 'classifier (binary)'} ranks binders better "
          f"(ΔAUC {auc_reg-auc_clf:+.3f})")

    fig,(a1,a2)=plt.subplots(1,2,figsize=(12,5))
    for p,lab,c in [(p_clf,f"binary classifier (AUC {auc_clf:.2f})","#2a6f97"),
                    (p_reg,f"pseudo-affinity regressor (AUC {auc_reg:.2f})","#c1502e")]:
        fpr,tpr,_=roc_curve(y_bin,p); a1.plot(fpr,tpr,lw=2,color=c,label=lab)
    a1.plot([0,1],[0,1],ls="--",c="0.7",lw=1); a1.set_xlabel("false positive rate"); a1.set_ylabel("true positive rate")
    a1.set_title("Ranking binders: binary vs dG target"); a1.legend(fontsize=8,loc="lower right"); a1.spines[["top","right"]].set_visible(False)
    a2.scatter(p_reg, y_reg, s=8, c=np.where(y_bin,"#2a6f97","#c9ced6"), alpha=0.5, edgecolors="none")
    a2.set_xlabel("predicted pseudo-affinity (out-of-fold)"); a2.set_ylabel("actual pseudo-affinity")
    a2.set_title(f"Regressor recovers graded binding (ρ={rho_reg:+.2f})"); a2.spines[["top","right"]].set_visible(False)
    fig.suptitle("DP4: binary label vs pseudo-affinity (dG) as the training target", fontsize=12, fontweight="bold")
    fig.tight_layout(rect=[0,0,1,0.95]); OUT.parent.mkdir(parents=True,exist_ok=True); fig.savefig(OUT,dpi=150)
    print(f"\nwrote {OUT.relative_to(ROOT)}")

if __name__ == "__main__":
    main()
