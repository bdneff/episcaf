#!/usr/bin/env python3
"""Model vs the hand-designed composite, as a design-prioritization tool.

Rank every design by a score, screen from the best down, and track what fraction of the real binders
you have found. A better score finds the binders sooner. We compare the composite (hand-set) to two
fitted models, logistic regression and gradient-boosted trees, each shown twice: with design metrics
only and with the full feature set (design + secondary structure + dynamics). This shows both the
model choice and what the added features buy. All fitted scores are out-of-fold (grouped CV by
epitope). Random picking is the diagonal, and Lawson's four-filter is a single pass/fail operating
point.

Run:  /usr/bin/python3 scripts/dp4_model_vs_composite.py
"""
from pathlib import Path
import numpy as np, pandas as pd
from sklearn.model_selection import StratifiedGroupKFold
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.metrics import roc_curve, roc_auc_score
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
LIB = ROOT / "data/libraries/dp4_library.csv"
SUMM = ROOT / "data/dp4_binding/john/scaffoldedEpitopeSummary.csv"
STRUCTF = ROOT / "data/dp4_binding/epitope_struct_features.csv"
OUT = ROOT / "data/dp4_binding/figs/model_vs_composite.png"
OUT_ROC = ROOT / "data/dp4_binding/figs/roc_comparison.png"
OUT_PREC = ROOT / "data/dp4_binding/figs/precision_vs_budget.png"

DESIGN = ["cylinder_clashes","epitope_rmsd","overall_rmsd","epitope_pae","scaffold_pae",
          "mean_pae","ptm","af3_clashes","island_index"]
EPI = ["helix","strand","n_islands","epi_bfac_z","epi_gnm_dfit","epi_dsasa"]
FULL = DESIGN + EPI

# (model, feature set) combos, and their plot styling: color by model, dashed=design only, solid=full.
COMBOS = [("logistic","design only",DESIGN), ("logistic","all features",FULL),
          ("grad-boost","design only",DESIGN), ("grad-boost","all features",FULL)]
COLOR = {"logistic":"#2a6f97", "grad-boost":"#6a4c93"}
DASH  = {"design only":(0,(5,2)), "all features":"solid"}
LW    = {"design only":1.8, "all features":2.7}

def make(model):
    if model == "logistic":
        return make_pipeline(SimpleImputer(strategy="median"), StandardScaler(),
                             LogisticRegression(max_iter=3000, class_weight="balanced"))
    return make_pipeline(SimpleImputer(strategy="median"),
                         HistGradientBoostingClassifier(class_weight="balanced", max_depth=3,
                         learning_rate=0.05, max_iter=300, random_state=0))

def gain(score, y):
    o = np.argsort(-score); yy = y[o]
    return np.concatenate([[0], np.arange(1,len(yy)+1)/len(yy)]), np.concatenate([[0], np.cumsum(yy)/yy.sum()])

def precision_curve(score, y, min_n=25):
    """Precision = bind rate among the top-scoring designs, as a function of how many you keep.
    x is the fraction of designs selected; y is the fraction of those that bind. Starts at min_n
    designs so the very top is not one lucky pick."""
    o = np.argsort(-score); yy = y[o]
    k = np.arange(1, len(yy)+1)
    prec = np.cumsum(yy) / k
    m = k >= min_n
    return k[m]/len(yy), prec[m]

def oof_repeated(model, X, y, groups, repeats=5):
    """Out-of-fold probabilities averaged over several shuffled grouped-CV partitions, so the AUC and
    curves match the repeated-CV ablation table rather than depending on one noisy split."""
    P = np.zeros(len(y))
    for r in range(repeats):
        for tr, te in StratifiedGroupKFold(5, shuffle=True, random_state=r).split(X, y, groups):
            m = make(model); m.fit(X[tr], y[tr]); P[te] += m.predict_proba(X[te])[:,1]
    return P / repeats

def main():
    lib = pd.read_csv(LIB, low_memory=False)
    summ = pd.read_csv(SUMM, encoding="latin-1"); summ.columns=[c.strip() for c in summ.columns]
    summ = summ.rename(columns={"#islands":"n_islands"}); summ["epitope"]=summ.epitope.str.lower()
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
    d=d.merge(epi,on="epitope",how="left")
    y=d.bound.values

    # out-of-fold scores for each (model, feature set), averaged over repeated shuffled grouped CV
    oof = {}
    for model, fs, feats in COMBOS:
        oof[(model,fs)] = oof_repeated(model, d[feats].values, y, d.epitope.values)
    xc,yc = gain(d.composite.values, y)
    def cap(x,yv,frac): return yv[np.searchsorted(x,frac)]
    # Lawson's four-filter: a hard pass/fail set (is_global_pass). Operating point = "make all passers".
    d["gp"]=(d.is_global_pass.astype(str).str.lower().isin(["true","1","1.0"])).astype(int)
    fp = d.gp.mean(); bp = y[d.gp.values==1].sum()/y.sum()

    base = y.mean()
    print(f"designs {len(d):,}, binders {y.sum()}, base bind rate {base:.1%}")
    print(f"{'method':30s} {'ROC-AUC':>8s} {'top10%':>8s} {'top20%':>8s}")
    print(f"{'composite (hand-set)':30s} {roc_auc_score(y,d.composite):8.3f} "
          f"{cap(xc,yc,0.10):8.0%} {cap(xc,yc,0.20):8.0%}")
    for model, fs, feats in COMBOS:
        s = oof[(model,fs)]; xg,yg = gain(s,y)
        print(f"{model+', '+fs:30s} {roc_auc_score(y,s):8.3f} {cap(xg,yg,0.10):8.0%} {cap(xg,yg,0.20):8.0%}")
    print(f"\nfour-filter (Lawson): makes {fp:.0%} of designs, captures {bp:.0%} of binders")
    print(f"\nprecision = bind rate AMONG the selected designs, at the four-filter's budget ({fp:.0%} selected):")
    print(f"  {'random (= base rate)':30s} {base:6.1%}   lift 1.0x")
    print(f"  {'four-filter (Lawson)':30s} {bp*y.sum()/(fp*len(y)):6.1%}   lift {bp/fp:.1f}x")
    for model, fs, feats in COMBOS:
        if fs != "all features": continue
        s = oof[(model,fs)]; xg,yg = gain(s,y); capf = cap(xg,yg,fp)
        print(f"  {model+', '+fs:30s} {capf*y.sum()/(fp*len(y)):6.1%}   lift {capf/fp:.1f}x")

    # ---- gain / cumulative-capture curve ----
    fig,ax=plt.subplots(figsize=(9,6.8))
    ax.plot([0,1],[0,1],ls="--",color="#b0b8c0",lw=1.5,label="random picking")
    ax.plot(xc,yc,color="#c1502e",lw=2.7,label="composite (hand-set)")
    for model, fs, feats in COMBOS:
        xg,yg = gain(oof[(model,fs)], y)
        ax.plot(xg,yg,color=COLOR[model],ls=DASH[fs],lw=LW[fs],label=f"{model}, {fs}")
    ax.plot([0,fp,1],[0,bp,1],color="#0e8f83",lw=2.2,label="four-filter (Lawson)")
    ax.plot([fp],[bp],marker="o",ms=7,color="#0e8f83",markeredgecolor="white",zorder=6)
    ax.axvline(0.20,color="0.6",lw=1,ls=":"); ax.text(0.205,0.02,"top 20%",fontsize=9,color="0.4")
    ax.set_xlabel("fraction of designs screened  (ranked best first)", fontsize=12)
    ax.set_ylabel("fraction of real binders found", fontsize=12)
    ax.set_title("Finding binders faster: logistic vs gradient boosting, before and after features\n"
                 "(higher = better at prioritizing which designs to make)", fontsize=12.5, fontweight="bold")
    ax.set_xlim(0,1); ax.set_ylim(0,1); ax.legend(fontsize=10, loc="lower right")
    ax.spines[["top","right"]].set_visible(False); ax.grid(alpha=0.15)
    fig.tight_layout(); OUT.parent.mkdir(parents=True,exist_ok=True); fig.savefig(OUT,dpi=150)
    print(f"wrote {OUT.relative_to(ROOT)}")

    # ---- ROC: same binary truth (binder/not); scores are curves, the threshold rule is a point ----
    fc,tc,_ = roc_curve(y, d.composite); acomp = roc_auc_score(y, d.composite)
    gp = d.gp.values
    tp=((gp==1)&(y==1)).sum(); fn=((gp==0)&(y==1)).sum()
    fpp=((gp==1)&(y==0)).sum(); tn=((gp==0)&(y==0)).sum()
    tpr_f, fpr_f = tp/(tp+fn), fpp/(fpp+tn)
    fig2,ax2=plt.subplots(figsize=(7.6,7))
    ax2.plot([0,1],[0,1],ls="--",color="#b0b8c0",lw=1.5,label="random (AUC 0.50)")
    ax2.plot(fc,tc,color="#c1502e",lw=2.7,label=f"composite (AUC {acomp:.2f})")
    for model, fs, feats in COMBOS:
        s = oof[(model,fs)]; f,t,_ = roc_curve(y,s)
        ax2.plot(f,t,color=COLOR[model],ls=DASH[fs],lw=LW[fs],label=f"{model}, {fs} (AUC {roc_auc_score(y,s):.2f})")
    ax2.plot([fpr_f],[tpr_f],marker="o",ms=13,color="#0e8f83",markeredgecolor="white",zorder=6,
             label="four-filter (one operating point)")
    ax2.set_xlabel("false positive rate  (non-binders wrongly selected)", fontsize=11)
    ax2.set_ylabel("true positive rate  (binders caught)", fontsize=11)
    ax2.set_title("Predicting binders: logistic vs gradient boosting, before and after features\n"
                  "(a curve above-left of another is better; the threshold rule is one point)",
                  fontsize=11.5, fontweight="bold")
    ax2.set_xlim(0,1); ax2.set_ylim(0,1); ax2.legend(fontsize=9, loc="lower right")
    ax2.spines[["top","right"]].set_visible(False); ax2.grid(alpha=0.15)
    fig2.tight_layout(); fig2.savefig(OUT_ROC,dpi=150)
    print(f"wrote {OUT_ROC.relative_to(ROOT)}")

    # ---- precision: bind rate AMONG the designs a rule selects, vs how many you keep ----
    scores_named = [("composite", d.composite.values, "#c1502e"),
                    ("logistic, all features", oof[("logistic","all features")], COLOR["logistic"]),
                    ("grad-boost, all features", oof[("grad-boost","all features")], COLOR["grad-boost"])]
    print("\nprecision (bind rate among selected) at fixed budgets:")
    print(f"  {'budget':>7s} {'random':>7s} {'composite':>10s} {'logistic':>9s} {'grad-boost':>11s}")
    for f in [0.05,0.10,0.20]:
        vals=[]
        for _,sc,_ in scores_named:
            xg,yg=gain(sc,y); vals.append(cap(xg,yg,f)*y.sum()/(f*len(y)))
        print(f"  {f:7.0%} {base:7.1%} {vals[0]:10.1%} {vals[1]:9.1%} {vals[2]:11.1%}")
    prec_ff = bp*y.sum()/(fp*len(y))

    fig3,ax3=plt.subplots(figsize=(9,6.5))
    ax3.axhline(base,ls="--",color="#b0b8c0",lw=1.6,label=f"random / base rate ({base:.0%})")
    for name,sc,col in scores_named:
        xp,yp=precision_curve(sc,y); ax3.plot(xp,yp,color=col,lw=2.7,label=name)
    ax3.plot([fp],[prec_ff],marker="o",ms=13,color="#0e8f83",markeredgecolor="white",zorder=6,
             label=f"four-filter (Lawson)")
    ax3.annotate(f"four-filter:\nkeep {fp:.0%}, {prec_ff:.0%} of them bind",(fp,prec_ff),
                 xytext=(fp+0.07,prec_ff+0.05),fontsize=9.5,color="#0a6b60",
                 arrowprops=dict(arrowstyle="->",color="#0e8f83"))
    ax3.set_xlabel("fraction of designs you keep  (ranked best first)", fontsize=12)
    ax3.set_ylabel("bind rate among the designs you keep", fontsize=12)
    ax3.set_title("How likely a selected design is to bind\n(precision vs how many you choose to make)",
                  fontsize=13, fontweight="bold")
    ax3.set_xlim(0,1); ax3.set_ylim(0,None); ax3.legend(fontsize=10, loc="upper right")
    ax3.spines[["top","right"]].set_visible(False); ax3.grid(alpha=0.15)
    fig3.tight_layout(); fig3.savefig(OUT_PREC,dpi=150)
    print(f"wrote {OUT_PREC.relative_to(ROOT)}")

if __name__ == "__main__":
    main()
