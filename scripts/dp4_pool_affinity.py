#!/usr/bin/env python3
"""Reconstruct antibody->pool, then derive a continuous pseudo-affinity per design from the dilution
series. This replaces the binary hit/no-hit label with a graded binding strength.

Step 1 — placement. For each epitope John called a binder, its HIT designs (hitIDs) should enrich in
whichever pool contains that antibody. Using the hits (strong signal) instead of the median over all
cognate designs sharpens the assignment. Validated: 5fhx should land in Pool8 (John's thread says so).
8pww came individually, so it is its own condition.

Step 2 — pseudo-affinity. For each cognate design of a placed antibody, take its condition's five
dilutions (10X..0.1X), compute log-enrichment vs the NoAbNoSer baseline at each, and integrate over
log-concentration (trapezoid). A design that stays enriched down to low antibody concentration scores
high: an EC50-like relative affinity, not a calibrated Kd.

Run:  /usr/bin/python3 scripts/dp4_pool_affinity.py
"""
from pathlib import Path
import numpy as np, pandas as pd
from scipy import stats
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
COUNTS = ROOT / "data/dp4_binding/john/counts_annotated.csv"
SUMM = ROOT / "data/dp4_binding/john/scaffoldedEpitopeSummary.csv"
OUTCSV = ROOT / "data/dp4_binding/dp4_pseudoaffinity.csv"
OUTPNG = ROOT / "data/dp4_binding/figs/pseudoaffinity.png"

DIL = ["10X", "3X", "1X", "0.3X", "0.1X"]
LOGC = {"10X": 1.0, "3X": np.log10(3), "1X": 0.0, "0.3X": np.log10(0.3), "0.1X": -1.0}
CONDS = [f"DP4_mAbs_Pool{k}" for k in range(1, 9)] + ["DP4_8pww"]

def main():
    c = pd.read_csv(COUNTS, low_memory=False)
    summ = pd.read_csv(SUMM, encoding="latin-1"); summ.columns = [x.strip() for x in summ.columns]
    num = lambda s: pd.to_numeric(s.astype(str).str.replace('"', '', regex=False).str.strip(), errors="coerce")
    l1p = lambda s: np.log10(1.0 + num(s))
    noab = c[[x for x in c.columns if x.startswith("DP4_NoAbNoSer")]].apply(num)
    base = np.log10(1.0 + noab).mean(axis=1)
    c["epi"] = c["target"].astype(str).str.split("_").str[0].str.lower()
    c["lm_short"] = "DP4_" + c["library_member"].astype(str).str.split("_").str[1]

    # per-condition enrichment at every dilution (only where the column exists)
    enr = {}   # (cond, dil) -> series
    for cond in CONDS:
        for dil in DIL:
            col = f"{cond}_{dil}"
            if col in c.columns:
                e = l1p(c[col]) - base
                enr[(cond, dil)] = e - e.median()   # remove per-channel background offset (within-channel centering)

    # ---- step 1: place each binder-antibody in a condition using its hit designs ----
    hit_by_epi = {}
    for _, r in summ.iterrows():
        e = str(r["epitope"]).lower(); ids = r["hitIDs"]
        if pd.isna(ids): continue
        hit_by_epi[e] = {"DP4_" + t.strip().split("_")[1] for t in str(ids).split(",") if t.strip().startswith("DP4_")}
    place = {}
    for e, hits in hit_by_epi.items():
        m = c["lm_short"].isin(hits)
        if m.sum() == 0: continue
        score = {cond: enr[(cond, "1X")][m].mean() for cond in CONDS if (cond, "1X") in enr}
        s = sorted(score.items(), key=lambda t: -t[1])
        place[e] = (s[0][0], s[0][1], s[0][1] - s[1][1])   # cond, top, margin
    print("placement (epitope -> condition, margin), binders only:")
    pool_map = {}
    for e, (cond, top, marg) in sorted(place.items(), key=lambda t: t[1][0]):
        tag = "  <-- 5fhx check" if e == "5fhx" else ""
        pool_map.setdefault(cond, []).append(e)
        print(f"  {e:8s} {cond.replace('DP4_mAbs_',''):12s} margin {marg:5.2f}{tag}")

    # ---- step 2: pseudo-affinity for every cognate design of a placed antibody ----
    rows = []
    for e, (cond, top, marg) in place.items():
        dils = [d for d in DIL if (cond, d) in enr]
        xs = np.array([LOGC[d] for d in dils])
        order = np.argsort(xs); xs = xs[order]; dils = [dils[i] for i in order]
        sub = c[c["epi"] == e]
        Y = np.vstack([enr[(cond, d)][sub.index].values for d in dils]).T   # designs x dilutions
        aff = np.trapz(Y, xs, axis=1)                                        # integrate over log-conc
        for lm, comp, a, y_low, hit in zip(sub["lm_short"], sub["composite"], aff,
                                           enr[(cond, dils[0])][sub.index].values,
                                           sub["lm_short"].isin(hit_by_epi.get(e, set()))):
            rows.append((lm, e, cond.replace("DP4_mAbs_", ""), a, y_low, comp, int(hit)))
    df = pd.DataFrame(rows, columns=["library_member", "epitope", "condition", "pseudo_affinity",
                                     "enr_lowdil", "composite", "is_hit"])
    df.to_csv(OUTCSV, index=False)
    print(f"\nplaced antibodies: {len(place)}   cognate designs with a pseudo-affinity: {len(df):,}"
          f"   (hits: {df.is_hit.sum()})")
    d2 = df[df.composite.notna()]
    rho, p = stats.spearmanr(d2.composite, d2.pseudo_affinity)
    print(f"composite vs pseudo-affinity: Spearman {rho:+.2f} (p={p:.1g}, n={len(d2):,})")
    print(f"wrote {OUTCSV.relative_to(ROOT)}")

    # ---- figure ----
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(12, 5))
    a1.scatter(d2.composite, d2.pseudo_affinity, s=8, c=np.where(d2.is_hit, "#2a6f97", "#c9ced6"),
               alpha=0.5, edgecolors="none")
    a1.set_xlabel("composite score"); a1.set_ylabel("pseudo-affinity (titration AUC)")
    a1.set_title(f"Composite vs graded binding  (ρ={rho:+.2f})"); a1.spines[["top","right"]].set_visible(False)
    # example titration curves: strongest, median, weak among hits
    hits = df[df.is_hit == 1].sort_values("pseudo_affinity")
    picks = [hits.iloc[-1], hits.iloc[len(hits)//2], df[df.is_hit==0].sort_values("pseudo_affinity").iloc[len(df[df.is_hit==0])//2]]
    labels = ["strong hit", "median hit", "typical non-hit"]
    for row, lab, col in zip(picks, labels, ["#1b3a4b", "#2a6f97", "#c1502e"]):
        cond = "DP4_mAbs_" + row.condition if not row.condition.startswith("DP4") else row.condition
        cond = cond if cond in [c2 for c2 in CONDS] else ("DP4_8pww" if row.condition=="8pww" else cond)
        idx = c.index[c.lm_short == row.library_member][0]
        dils = [d for d in DIL if (cond, d) in enr]
        xs = [LOGC[d] for d in dils]; ys = [enr[(cond, d)].loc[idx] for d in dils]
        o = np.argsort(xs)
        a2.plot(np.array(xs)[o], np.array(ys)[o], "-o", color=col, label=f"{lab} ({row.epitope})")
    a2.axhline(0, ls=":", c="0.6"); a2.set_xlabel("log10 antibody concentration (10X … 0.1X)")
    a2.set_ylabel("log-enrichment vs baseline"); a2.set_title("Titration curves = the affinity signal")
    a2.legend(fontsize=8); a2.spines[["top","right"]].set_visible(False)
    fig.suptitle("DP4 pseudo-affinity from the dilution series", fontsize=13, fontweight="bold")
    fig.tight_layout(rect=[0,0,1,0.95]); OUTPNG.parent.mkdir(parents=True, exist_ok=True); fig.savefig(OUTPNG, dpi=150)
    print(f"wrote {OUTPNG.relative_to(ROOT)}")

if __name__ == "__main__":
    main()
