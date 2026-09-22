#!/usr/bin/env python3
"""Head-to-head of binder models on held-out DIVERSE epitopes, repeated over many random splits.

A single held-out split of ~12 epitopes is a noisy AUC estimate (few groups, few test binders), so a
one-off ranking of the good models is not trustworthy. Here we repeat the split many times: each repeat
starts farthest-point sampling from a random seed epitope, so it picks a different spread-out set of 6
binder + 6 non-binder epitopes, holds ALL their designs out, trains on the rest, and scores train and
held-out AUC. We report mean +/- sd over repeats. This separates two questions the single split
confounds: (1) does a model overfit (train-vs-held-out gap), and (2) which model generalizes best
(held-out mean across repeats, with its spread).

Models: composite (hand-set baseline), logistic, gradient-boosted trees, random forest, a right-sized
MLP, and a deliberately oversized MLP (to confirm it overfits: high train AUC, low held-out AUC).

Features: design metrics + epitope structure (John's helix/strand/loop + our native-context features
from the crystal complexes: size, native exposure, radius of gyration, composition).

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
from sklearn.metrics import roc_auc_score, average_precision_score
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
LIB = ROOT / "data/libraries/dp4_library.csv"
SUMM = ROOT / "data/dp4_binding/john/scaffoldedEpitopeSummary.csv"
STRUCTF = ROOT / "data/dp4_binding/epitope_struct_features.csv"
OUT = ROOT / "data/dp4_binding/figs/headtohead.png"

DESIGN = ["cylinder_clashes","epitope_rmsd","overall_rmsd","epitope_pae","scaffold_pae",
          "mean_pae","ptm","af3_clashes","island_index"]
# Same feature set as the modeling section (binder_ml_v2 "+ all new features"): design metrics +
# secondary structure + dynamics, so Table 3 and Table 2 compare the same inputs.
EPI = ["helix","strand","n_islands","epi_bfac_z","epi_gnm_dfit","epi_dsasa"]
SEL = ["helix","epi_size","epi_rel_sasa","n_islands"]   # axes for diverse test selection (not model inputs)
REPEATS = 20
PER_CLASS = 6

def fps(X, k, start):
    """Farthest-point sampling from a given start index."""
    idx = [int(start)]
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
        st[["epitope","epi_size","epi_rel_sasa","epi_rg","epi_bfac_z","epi_gnm_dfit","epi_dsasa"]],
        on="epitope", how="inner")
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

def choose_test(epi, rng, per_class=PER_CLASS):
    from sklearn.preprocessing import StandardScaler as SS
    test = []
    for cls in [True, False]:
        sub = epi[epi.yes == cls].reset_index(drop=True)
        X = SS().fit_transform(sub[SEL].fillna(sub[SEL].median()).values)
        start = rng.integers(len(sub))
        pick = fps(X, min(per_class, len(sub)), start)
        test += sub.loc[pick, "epitope"].tolist()
    return set(test)

def build_models():
    pipe = lambda est: make_pipeline(SimpleImputer(strategy="median"), StandardScaler(), est)
    return {
        "logistic":         lambda: pipe(LogisticRegression(max_iter=3000, class_weight="balanced")),
        "grad-boost":       lambda: make_pipeline(SimpleImputer(strategy="median"),
                              HistGradientBoostingClassifier(class_weight="balanced", max_depth=3,
                              learning_rate=0.05, max_iter=300, random_state=0)),
        "random-forest":    lambda: make_pipeline(SimpleImputer(strategy="median"),
                              RandomForestClassifier(400, min_samples_leaf=20, class_weight="balanced", random_state=0)),
        "MLP (right-size)": lambda: pipe(MLPClassifier(hidden_layer_sizes=(32,16), alpha=1e-3,
                              early_stopping=True, max_iter=800, random_state=0)),
        "MLP (oversized)":  lambda: pipe(MLPClassifier(hidden_layer_sizes=(512,512,512,512), alpha=1e-6,
                              early_stopping=False, max_iter=800, random_state=0)),
    }

def main():
    d, epi = load()
    F = DESIGN + EPI
    factories = build_models()
    names = ["composite (baseline)"] + list(factories)
    tr_auc = {n: [] for n in names}; te_auc = {n: [] for n in names}
    n_test_binders = []
    for r in range(REPEATS):
        rng = np.random.default_rng(r)
        test_epi = choose_test(epi, rng)
        te = d[d.epitope.isin(test_epi)]; tr = d[~d.epitope.isin(test_epi)]
        ytr, yte = tr.bound.values, te.bound.values
        if yte.sum() == 0 or yte.sum() == len(yte): continue      # degenerate split, skip
        n_test_binders.append(int(yte.sum()))
        # composite baseline, no fit
        tr_auc["composite (baseline)"].append(roc_auc_score(ytr, tr.composite))
        te_auc["composite (baseline)"].append(roc_auc_score(yte, te.composite))
        Xtr, Xte = tr[F].values, te[F].values
        for name, make in factories.items():
            m = make(); m.fit(Xtr, ytr)
            tr_auc[name].append(roc_auc_score(ytr, m.predict_proba(Xtr)[:,1]))
            te_auc[name].append(roc_auc_score(yte, m.predict_proba(Xte)[:,1]))

    nrep = len(te_auc["composite (baseline)"])
    print(f"{REPEATS} random diverse held-out splits ({PER_CLASS}+{PER_CLASS} epitopes each), "
          f"{nrep} usable; median {int(np.median(n_test_binders))} test binders per split\n")
    print(f"{'model':22s} {'train AUC':>14s} {'held-out AUC':>16s} {'gap':>8s}")
    stats = {}
    for n in names:
        a, b = np.array(tr_auc[n]), np.array(te_auc[n])
        stats[n] = (a.mean(), a.std(), b.mean(), b.std())
        print(f"{n:22s}  {a.mean():.3f}+/-{a.std():.3f}   {b.mean():.3f}+/-{b.std():.3f}   {a.mean()-b.mean():+.3f}")

    # ---- figure: (A) train vs held-out mean+/-sd; (B) held-out spread across repeats ----
    fig,(a1,a2) = plt.subplots(1,2, figsize=(13,5.4), gridspec_kw={"width_ratios":[1.1,1]})
    x = np.arange(len(names)); lbl = [n.replace(" (","\n(") for n in names]
    a1.bar(x-0.2,[stats[n][0] for n in names],0.4,yerr=[stats[n][1] for n in names],capsize=3,
           color="#c9ced6",label="train AUC")
    a1.bar(x+0.2,[stats[n][2] for n in names],0.4,yerr=[stats[n][3] for n in names],capsize=3,
           color="#2a6f97",label="held-out AUC")
    a1.set_xticks(x); a1.set_xticklabels(lbl,fontsize=8)
    a1.set_ylim(0.4,1.02); a1.axhline(0.5,ls=":",c="0.6",lw=1)
    a1.set_ylabel("ROC-AUC"); a1.set_title(f"Train vs held-out over {nrep} diverse splits\n(big gap = overfitting)")
    a1.legend(fontsize=8); a1.spines[["top","right"]].set_visible(False)
    data = [np.array(te_auc[n]) for n in names]
    bp = a2.boxplot(data, vert=True, widths=0.6, patch_artist=True, showfliers=False)
    for patch in bp["boxes"]: patch.set_facecolor("#dbe6ef"); patch.set_edgecolor("#2a6f97")
    for med in bp["medians"]: med.set_color("#c1502e"); med.set_linewidth(2)
    for i,n in enumerate(names):
        a2.scatter(np.full(len(data[i]), i+1)+np.random.default_rng(1).uniform(-0.12,0.12,len(data[i])),
                   data[i], s=10, color="#2a6f97", alpha=0.5, edgecolors="none", zorder=3)
    a2.set_xticks(x+1); a2.set_xticklabels(lbl,fontsize=8)
    a2.axhline(0.5,ls=":",c="0.6",lw=1)
    a2.set_ylabel("held-out ROC-AUC"); a2.set_title("Held-out spread across splits\n(is a model reliably ahead?)")
    a2.spines[["top","right"]].set_visible(False)
    fig.suptitle(f"DP4 binder models — repeated held-out on diverse epitopes ({nrep} random splits)",
                 fontsize=12, fontweight="bold")
    fig.tight_layout(rect=[0,0,1,0.94]); OUT.parent.mkdir(parents=True,exist_ok=True); fig.savefig(OUT,dpi=150)
    print(f"\nwrote {OUT.relative_to(ROOT)}")

if __name__ == "__main__":
    main()
