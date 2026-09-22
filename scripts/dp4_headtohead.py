#!/usr/bin/env python3
"""Head-to-head of binder models on a held-out set of DIVERSE epitopes.

Test set: farthest-point sampling picks a spread-out set of binder and non-binder epitopes (both
classes, spanning the feature space), and ALL designs of those epitopes are held out. Training uses
the remaining epitopes. Scoring designs by epitope this way is the honest generalization test:
nothing about a test antibody is seen in training.

Models: composite (hand-set baseline), logistic, gradient-boosted trees, random forest, a right-sized
MLP, and a deliberately oversized MLP (to confirm it overfits: high train AUC, low test AUC).

Features: design metrics + epitope structure (John's helix/strand/loop + our native-context features
from the crystal complexes: size, native exposure, radius of gyration, composition).

Caveat: one held-out split of ~12 epitopes is a noisy AUC estimate (few groups). The robust,
qualitative result is the train-vs-test gap, especially for the oversized MLP.

Run:  /usr/bin/python3 scripts/dp4_headtohead.py
"""
from pathlib import Path
import numpy as np, pandas as pd
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import HistGradientBoostingClassifier, RandomForestClassifier
from sklearn.neural_network import MLPClassifier
from sklearn.metrics import roc_auc_score, average_precision_score, roc_curve
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
LIB = ROOT / "data/libraries/dp4_library.csv"
SUMM = ROOT / "data/dp4_binding/john/scaffoldedEpitopeSummary.csv"
STRUCTF = ROOT / "data/dp4_binding/epitope_struct_features.csv"
OUT = ROOT / "data/dp4_binding/figs/headtohead.png"

DESIGN = ["cylinder_clashes","epitope_rmsd","overall_rmsd","epitope_pae","scaffold_pae",
          "mean_pae","ptm","af3_clashes","island_index"]
EPI = ["helix","strand","loop","n_islands","epi_size","epi_rel_sasa","epi_rg",
       "epi_frac_hydrophobic","epi_frac_charged","epi_frac_aromatic"]
SEL = ["helix","epi_size","epi_rel_sasa","n_islands"]   # axes for diverse test selection

def fps(X, k):
    idx = [int(np.argmin(np.linalg.norm(X - X.mean(0), axis=1)))]
    d = np.full(len(X), np.inf)
    for _ in range(k - 1):
        d = np.minimum(d, np.linalg.norm(X - X[idx[-1]], axis=1))
        d[idx] = -1
        idx.append(int(np.argmax(d)))
    return idx

def load():
    lib = pd.read_csv(LIB, low_memory=False)
    summ = pd.read_csv(SUMM, encoding="latin-1"); summ.columns = [c.strip() for c in summ.columns]
    summ = summ.rename(columns={"#islands": "n_islands"})
    summ["epitope"] = summ["epitope"].astype(str).str.lower()
    st = pd.read_csv(STRUCTF)
    epi = summ[["epitope","helix","strand","loop","n_islands","scafBinding?"]].merge(
        st[["epitope","epi_size","epi_rel_sasa","epi_rg","epi_frac_hydrophobic",
            "epi_frac_charged","epi_frac_aromatic"]], on="epitope", how="inner")
    epi["yes"] = epi["scafBinding?"].astype(str).str.strip().str.upper().eq("YES")
    hits = set()
    for ids in summ["hitIDs"].dropna():
        for t in str(ids).split(","):
            t = t.strip()
            if t.startswith("DP4_"): hits.add("DP4_" + t.split("_")[1])
    lib["epitope"] = lib["target"].astype(str).str.split("_").str[0].str.lower()
    d = lib[lib["composite"].notna() & lib["epitope"].isin(set(epi["epitope"]))].copy()
    d["bound"] = d["library_member"].isin(hits).astype(int)
    d = d.merge(epi.drop(columns=["scafBinding?"]), on="epitope", how="left")
    return d, epi

def choose_test(epi, per_class=6):
    from sklearn.preprocessing import StandardScaler as SS
    test = []
    for cls in [True, False]:
        sub = epi[epi.yes == cls].reset_index(drop=True)
        X = SS().fit_transform(sub[SEL].fillna(sub[SEL].median()).values)
        pick = fps(X, min(per_class, len(sub)))
        test += sub.loc[pick, "epitope"].tolist()
    return set(test)

def scores(model, Xtr, ytr, Xte, yte):
    model.fit(Xtr, ytr)
    ptr, pte = model.predict_proba(Xtr)[:,1], model.predict_proba(Xte)[:,1]
    return (roc_auc_score(ytr,ptr), roc_auc_score(yte,pte),
            average_precision_score(ytr,ptr), average_precision_score(yte,pte), pte)

def main():
    d, epi = load()
    test_epi = choose_test(epi)
    te = d[d.epitope.isin(test_epi)]; tr = d[~d.epitope.isin(test_epi)]
    print(f"held-out test epitopes ({len(test_epi)}): {sorted(test_epi)}")
    print(f"  test binders/non among epitopes: "
          f"{epi[epi.epitope.isin(test_epi)].yes.sum()} YES / {(~epi[epi.epitope.isin(test_epi)].yes).sum()} NO")
    print(f"train designs {len(tr):,} ({tr.bound.sum()} bound) | test designs {len(te):,} ({te.bound.sum()} bound)")

    F = DESIGN + EPI
    Xtr, ytr = tr[F].values, tr.bound.values
    Xte, yte = te[F].values, te.bound.values

    pipe = lambda est: make_pipeline(SimpleImputer(strategy="median"), StandardScaler(), est)
    models = {
        "logistic":        pipe(LogisticRegression(max_iter=3000, class_weight="balanced")),
        "grad-boost":      HistGradientBoostingClassifier(class_weight="balanced", max_depth=3,
                             learning_rate=0.05, max_iter=300, random_state=0),
        "random-forest":   RandomForestClassifier(400, min_samples_leaf=20, class_weight="balanced", random_state=0),
        "MLP (right-size)":pipe(MLPClassifier(hidden_layer_sizes=(32,16), alpha=1e-3,
                             early_stopping=True, max_iter=800, random_state=0)),
        "MLP (oversized)": pipe(MLPClassifier(hidden_layer_sizes=(512,512,512,512), alpha=1e-6,
                             early_stopping=False, max_iter=1500, random_state=0)),
    }
    # baseline: composite, no fit
    res = {"composite (baseline)": (roc_auc_score(ytr, tr.composite), roc_auc_score(yte, te.composite),
                                    average_precision_score(ytr, tr.composite),
                                    average_precision_score(yte, te.composite), te.composite.values)}
    for imp in [SimpleImputer(strategy="median")]:  # fill NaNs for tree models
        Xtr_f = imp.fit_transform(Xtr); Xte_f = imp.transform(Xte)
    for name, mdl in models.items():
        X1, X2 = (Xtr, Xte) if "pipeline" in str(type(mdl)).lower() else (Xtr_f, Xte_f)
        res[name] = scores(mdl, X1, ytr, X2, yte)

    print(f"\n{'model':22s} {'train AUC':>9s} {'test AUC':>9s} {'test PR':>8s}   gap")
    for k,(atr,ate,ptr,pte,_) in res.items():
        print(f"{k:22s} {atr:9.3f} {ate:9.3f} {pte:8.3f}   {atr-ate:+.3f}")

    # ---- figure ----
    names = list(res)
    fig,(a1,a2) = plt.subplots(1,2, figsize=(13,5), gridspec_kw={"width_ratios":[1.15,1]})
    x = np.arange(len(names))
    a1.bar(x-0.2,[res[n][0] for n in names],0.4,color="#c9ced6",label="train AUC")
    a1.bar(x+0.2,[res[n][1] for n in names],0.4,color="#2a6f97",label="test AUC (held-out epitopes)")
    a1.set_xticks(x); a1.set_xticklabels([n.replace(" (","\n(") for n in names],fontsize=8)
    a1.set_ylim(0.4,1.02); a1.axhline(0.5,ls=":",c="0.6",lw=1)
    a1.set_ylabel("ROC-AUC"); a1.set_title("Train vs held-out test (watch the oversized MLP gap)")
    a1.legend(fontsize=8); a1.spines[["top","right"]].set_visible(False)
    for n in names:
        _,ate,_,_,p = res[n]; fpr,tpr,_ = roc_curve(yte,p)
        a2.plot(fpr,tpr,lw=2,label=f"{n} ({ate:.2f})")
    a2.plot([0,1],[0,1],ls="--",c="0.7",lw=1)
    a2.set_xlabel("false positive rate"); a2.set_ylabel("true positive rate")
    a2.set_title("ROC on held-out epitopes"); a2.legend(fontsize=7.5,loc="lower right")
    a2.spines[["top","right"]].set_visible(False)
    fig.suptitle(f"DP4 binder models — head to head on {len(test_epi)} held-out diverse epitopes "
                 f"({len(te):,} test designs)", fontsize=12, fontweight="bold")
    fig.tight_layout(rect=[0,0,1,0.95]); OUT.parent.mkdir(parents=True,exist_ok=True); fig.savefig(OUT,dpi=150)
    print(f"\nwrote {OUT.relative_to(ROOT)}")

if __name__ == "__main__":
    main()
