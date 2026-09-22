#!/usr/bin/env python3
"""Learn a binder classifier from the raw design metrics, and ask whether a fitted model beats the
hand-set composite score. This is the "fit the dials, don't guess them" analysis DP4's C5 titration
arm was built to enable: C5 and the controls span the parameter space (including designs the filters
reject), so the metrics finally have the variance needed to fit.

Honest evaluation: grouped cross-validation by EPITOPE, so every fold is tested on antibodies it never
trained on (designs cluster by epitope; a design-level split would leak the between-antibody signal).
With ~56 epitopes that is the true sample size for any epitope-level feature, so models are simple and
compared to the composite-alone baseline.

Label: a design is a hit if its library_member is in any epitope's hitIDs
(data/dp4_binding/john/scaffoldedEpitopeSummary.csv). Universe: scaffolded designs with metrics whose
antibody was assayed.

Run:  /usr/bin/python3 scripts/dp4_binder_ml.py
"""
from pathlib import Path
import numpy as np, pandas as pd
from sklearn.model_selection import GroupKFold, cross_val_predict
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.metrics import roc_auc_score, average_precision_score, roc_curve
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
LIB = ROOT / "data/libraries/dp4_library.csv"
SUMM = ROOT / "data/dp4_binding/john/scaffoldedEpitopeSummary.csv"
OUT = ROOT / "data/dp4_binding/figs/binder_ml.png"

DESIGN = ["cylinder_clashes", "epitope_rmsd", "overall_rmsd", "epitope_pae",
          "scaffold_pae", "mean_pae", "ptm", "af3_clashes", "island_index"]
EPI = ["helix", "strand", "loop", "n_islands"]
POS, NEG, BASE, GB = "#2a6f97", "#c1502e", "#9aa4ad", "#1b3a4b"

def load():
    lib = pd.read_csv(LIB, low_memory=False)
    summ = pd.read_csv(SUMM, encoding="latin-1"); summ.columns = [c.strip() for c in summ.columns]
    hits = set()
    for ids in summ["hitIDs"].dropna():
        for t in str(ids).split(","):
            t = t.strip()
            if t.startswith("DP4_"): hits.add("DP4_" + t.split("_")[1])
    epi = summ.rename(columns={"#islands": "n_islands"})[["epitope","helix","strand","loop","n_islands"]].copy()
    epi["epitope"] = epi["epitope"].astype(str).str.lower()
    lib["epitope"] = lib["target"].astype(str).str.split("_").str[0].str.lower()
    d = lib[lib["composite"].notna() & lib["epitope"].isin(set(epi["epitope"]))].copy()
    d["bound"] = d["library_member"].isin(hits).astype(int)
    return d.merge(epi, on="epitope", how="left")

def oof(model, X, y, groups):
    p = cross_val_predict(model, X, y, groups=groups, cv=GroupKFold(5), method="predict_proba")[:, 1]
    return p, roc_auc_score(y, p), average_precision_score(y, p)

def main():
    d = load()
    y = d["bound"].values; g = d["epitope"].values; base = y.mean()
    print(f"designs={len(d):,}  bound={int(y.sum())}  epitopes={d['epitope'].nunique()}  base rate={base:.3f}")

    logit = lambda: make_pipeline(SimpleImputer(strategy="median"), StandardScaler(),
                                  LogisticRegression(max_iter=2000, class_weight="balanced"))
    hgb = lambda: HistGradientBoostingClassifier(class_weight="balanced", max_depth=3,
                                                 learning_rate=0.05, max_iter=300, random_state=0)
    results = {}   # name -> (proba, auc, ap, color)
    ac, ap0 = roc_auc_score(y, d["composite"]), average_precision_score(y, d["composite"])
    results["composite (hand-set)"] = (d["composite"].values, ac, ap0, BASE)
    for name, mdl, feats, col in [
        ("logistic · design",          logit(), DESIGN,       POS),
        ("logistic · design+epitope",  logit(), DESIGN+EPI,   "#4a8bb5"),
        ("grad-boost · design+epitope",hgb(),   DESIGN+EPI,   GB)]:
        p, a, ap = oof(mdl, d[feats].values, y, g)
        results[name] = (p, a, ap, col)
    print(f"\n{'model':30s} {'ROC-AUC':>8s} {'PR-AUC':>8s}")
    for k, (_, a, ap, _) in results.items():
        print(f"{k:30s} {a:8.3f} {ap:8.3f}")

    lr = logit().fit(d[DESIGN+EPI].values, y)
    coef = dict(zip(DESIGN+EPI, lr.named_steps["logisticregression"].coef_[0]))

    # ---------------- figure ----------------
    fig, (a1, a2, a3) = plt.subplots(1, 3, figsize=(15, 4.6), gridspec_kw={"width_ratios":[1,1.1,1.1]})

    # A: AUC lift bars
    names = list(results)
    aucs = [results[n][1] for n in names]; aps = [results[n][2] for n in names]
    x = np.arange(len(names))
    a1.bar(x-0.2, aucs, 0.4, color=[results[n][3] for n in names], label="ROC-AUC")
    a1.bar(x+0.2, aps, 0.4, color=[results[n][3] for n in names], alpha=0.5, label="PR-AUC")
    a1.axhline(base, ls=":", c="0.5", lw=1); a1.text(len(names)-0.5, base+0.01, f"chance PR={base:.2f}", fontsize=8, ha="right", color="0.4")
    a1.set_xticks(x); a1.set_xticklabels([n.replace(" · ", "\n") for n in names], fontsize=8, rotation=0)
    a1.set_ylim(0, 1); a1.set_title("Fitted models beat the hand-set score"); a1.legend(fontsize=8, loc="upper left")
    a1.spines[["top","right"]].set_visible(False)

    # B: ROC curves (out-of-fold)
    for n in names:
        p, auc, _, col = results[n]
        fpr, tpr, _ = roc_curve(y, p)
        a2.plot(fpr, tpr, color=col, lw=2, label=f"{n}  (AUC {auc:.2f})")
    a2.plot([0,1],[0,1], ls="--", c="0.7", lw=1)
    a2.set_xlabel("false positive rate"); a2.set_ylabel("true positive rate")
    a2.set_title("ROC, held-out epitopes"); a2.legend(fontsize=7.5, loc="lower right")
    a2.spines[["top","right"]].set_visible(False)

    # C: fitted weights
    items = sorted(coef.items(), key=lambda t: t[1])
    labels = [k for k,_ in items]; vals = [v for _,v in items]
    cols = [POS if v>0 else NEG for v in vals]
    a3.barh(range(len(labels)), vals, color=cols)
    a3.set_yticks(range(len(labels))); a3.set_yticklabels(labels, fontsize=8)
    a3.axvline(0, c="0.6", lw=1)
    a3.set_xlabel("standardized weight  (+ favors binding)")
    a3.set_title("What the fit leans on"); a3.spines[["top","right"]].set_visible(False)

    fig.suptitle("DP4 binder model — learn 'good binder' from design + epitope features "
                 f"(grouped CV by epitope, {len(d):,} designs / {d['epitope'].nunique()} epitopes)",
                 fontsize=12, fontweight="bold")
    fig.tight_layout(rect=[0,0,1,0.94])
    OUT.parent.mkdir(parents=True, exist_ok=True); fig.savefig(OUT, dpi=150)
    print(f"\nwrote {OUT.relative_to(ROOT)}")

if __name__ == "__main__":
    main()
