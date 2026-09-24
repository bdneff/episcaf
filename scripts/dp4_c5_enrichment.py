#!/usr/bin/env python3
"""Non-circular enrichment: on UNSELECTED designs (C5), does the model concentrate binders?

C5 spans the metric space on purpose (it was NOT selected for a high score), so it is the honest arm to
demonstrate enrichment. We train the binder model on the OTHER arms' designs from OTHER epitopes and
predict each held-out epitope's C5 designs (leave-one-epitope-out), so a C5 design is scored by a model
that never saw its epitope's antibody. Then we screen C5 from the top and track the binding rate among
what we keep, against C5's own unselected base rate (6.2%). This is Brandon's "enrichment score in actual
binding" done without circularity, on the arm built for it (manuscript sec:dp4c5).

We compare the fitted model to the recomputed composite (the hand-set scorer) and to random.

Run:  /usr/bin/python3 scripts/dp4_c5_enrichment.py
"""
from pathlib import Path
import sys
import numpy as np, pandas as pd
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]; sys.path.insert(0, str(ROOT))
from episcaf_analysis.score import score
from episcaf_analysis.presets import ANTIBODY_SOFTGATE

LIB = ROOT/"data/libraries/dp4_library.csv"; SUMM = ROOT/"data/dp4_binding/john/scaffoldedEpitopeSummary.csv"
STRUCTF = ROOT/"data/dp4_binding/epitope_struct_features.csv"; OUT = ROOT/"data/dp4_binding/figs/c5_enrichment.png"
DESIGN=["cylinder_clashes","epitope_rmsd","overall_rmsd","epitope_pae","scaffold_pae","mean_pae","ptm","af3_clashes","island_index"]
EPI=["helix","strand","n_islands","epi_bfac_z","epi_gnm_dfit","epi_dsasa"]
BLUE="#2456E6"; ORANGE="#c1502e"; GREEN="#2C8C4A"; MUTE="#5b6570"

def load():
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
    lib["bound"]=lib["library_member"].isin(hits).astype(int)
    lib=lib.merge(epi, on="epitope", how="left")
    return lib, epi

def main():
    lib, epi = load()
    train = lib[lib.category.isin(["scaffoldedAbEpitope","scaffoldedSingleIsland"])].copy()  # C1+C2
    c5 = lib[lib.category=="metricSpaceTitration"].copy()                                     # unselected
    c5 = c5[c5.epitope.isin(set(epi.epitope))].copy()
    y5 = c5["bound"].values; base = y5.mean()

    # leave-one-epitope-out: train on C1/C2 of other epitopes, predict this epitope's C5
    oof = np.full(len(c5), np.nan)
    for e in c5.epitope.unique():
        tr = train[train.epitope != e]
        if tr["bound"].nunique() < 2: continue
        m = make_pipeline(SimpleImputer(strategy="median"), StandardScaler(),
                          LogisticRegression(max_iter=3000, class_weight="balanced"))
        m.fit(tr[DESIGN+EPI].values, tr["bound"].values)
        idx = (c5.epitope==e).values
        oof[idx] = m.predict_proba(c5.loc[idx, DESIGN+EPI].values)[:,1]
    ok = ~np.isnan(oof); c5=c5[ok].copy(); y5=y5[ok]; oof=oof[ok]

    # refit tests: grouped CV WITHIN C5 (train on other C5 epitopes) -- C5 is the arm built for weight-fitting
    from sklearn.model_selection import StratifiedGroupKFold
    def cv_within(feats):
        X=c5[feats].values; g=c5["epitope"].values; p=np.zeros(len(c5))
        for tr,te in StratifiedGroupKFold(5, shuffle=True, random_state=0).split(X,y5,g):
            mm=make_pipeline(SimpleImputer(strategy="median"),StandardScaler(),
                             LogisticRegression(max_iter=3000,class_weight="balanced"))
            mm.fit(X[tr],y5[tr]); p[te]=mm.predict_proba(X[te])[:,1]
        return p
    oof_c5full = cv_within(DESIGN+EPI)      # refit on C5, all features
    oof_c5design = cv_within(DESIGN)        # refit on C5, DESIGN METRICS ONLY (= refit the composite)

    # recompute the hand-set composite on C5 for comparison
    c5r=c5.rename(columns={"af3_clashes":"af3_n_clash_res","epitope_rmsd":"epitope_chunk_rmsd"}).copy()
    c5r["id"]=c5r["library_member"]; c5r["antigen"]=c5r["epitope"]
    comp = score(c5r, ANTIBODY_SOFTGATE)["composite"].values

    def prec_curve(s, minn=30):
        o=np.argsort(-s); yy=y5[o]; k=np.arange(1,len(yy)+1); pr=np.cumsum(yy)/k
        m=k>=minn; return k[m]/len(yy), pr[m]
    def cap(s, frac):
        o=np.argsort(-s); yy=y5[o]; k=int(round(frac*len(yy))); return yy[:k].mean()

    aucs = {"transfer  C1/C2 -> C5":roc_auc_score(y5,oof),
            "refit on C5 (design + context)":roc_auc_score(y5,oof_c5full),
            "refit on C5 (design metrics ONLY)":roc_auc_score(y5,oof_c5design),
            "hand-set composite":roc_auc_score(y5,comp)}
    print(f"C5 unselected: {len(c5):,} designs, {int(y5.sum())} binders, base rate {base:.1%}\n")
    print("AUC on unselected C5 (0.5 = no signal):")
    for k,v in aucs.items(): print(f"  {k:36s} {v:.3f}")
    best = oof_c5full
    print(f"\nenrichment vs the {base:.1%} baseline:")
    print(f"{'keep top':>9s} {'context model':>14s} {'design-only':>12s} {'composite':>10s}   context fold")
    for f in [0.05,0.10,0.20]:
        print(f"{f:>9.0%} {cap(best,f):>14.1%} {cap(oof_c5design,f):>12.1%} {cap(comp,f):>10.1%}   {cap(best,f)/base:>5.1f}x")

    # Lawson's four-filter, evaluated on C5, as a single pass/fail operating point
    num=lambda c: pd.to_numeric(c5[c],errors="coerce")
    ff=((num("epitope_rmsd")<=1)&(num("overall_rmsd")<=2)&(num("mean_pae")<5)&(num("af3_clashes")==0)).values
    fp=ff.mean(); prec_ff=y5[ff].mean() if ff.sum() else np.nan
    print(f"\nLawson four-filter on C5: {int(ff.sum())} pass ({fp:.1%}), binding rate among passers {prec_ff:.1%}  ({prec_ff/base:.1f}x)")
    print(f"hand-set composite AUC on C5: {roc_auc_score(y5,comp):.3f}")

    xm,ym=prec_curve(best); xc,yc=prec_curve(comp)
    fig,ax=plt.subplots(figsize=(9.4,6))
    ax.axhline(base*100, ls="--", color=GREEN, lw=2, label=f"unselected baseline ({base:.1%})")
    ax.plot(xc*100, yc*100, color=ORANGE, lw=2.4, label="hand-set composite score")
    ax.plot(xm*100, ym*100, color=BLUE, lw=2.9, label="logistic (design + epitope context)")
    ax.text(0.985, 0.06, f"Lawson four-filter passes only {int(ff.sum())} of {len(c5):,} here\n"
            "(C5 spans filter-rejects by design, so the hard gate selects nothing)",
            transform=ax.transAxes, ha="right", va="bottom", fontsize=10.5, color="#7b3fa0",
            bbox=dict(boxstyle="round", fc="white", ec="#7b3fa0", alpha=0.9))
    ax.set_xlabel("% of unselected designs kept (ranked best first)", fontsize=13)
    ax.set_ylabel("binding rate among those kept (%)", fontsize=13)
    ax.set_title("Enriching binders among UNSELECTED designs (C5)\nlearned model vs the hand-set composite and "
                 "Lawson's four-filter", fontsize=13.5, fontweight="bold")
    ax.legend(fontsize=11.5, loc="upper right"); ax.spines[["top","right"]].set_visible(False); ax.grid(alpha=0.15)
    ax.set_xlim(0,100); ax.set_ylim(0,None)
    fig.tight_layout(); fig.savefig(OUT, dpi=150); print(f"\nwrote {OUT.relative_to(ROOT)}")

if __name__ == "__main__":
    main()
