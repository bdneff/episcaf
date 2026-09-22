#!/usr/bin/env python3
"""Context-dependent scaffolding rules: within a given epitope context, which design metric predicts
binding — and does the answer change across contexts?

This tests Brandon's hypothesis (the right way to scaffold depends on the native context) directly,
instead of forcing one global model. For each design metric we compute its WITHIN-antibody correlation
with binding (fixed-effect centered, so it is the honest design-level signal), separately inside each
context stratum (helical vs not, small vs large interface, one vs two islands). If a metric predicts
in one context and not another, that is a context-specific design rule.

Sign: negative means "lower metric -> more binding" (the expected direction for clash / RMSD / PAE).

Run:  /usr/bin/python3 scripts/dp4_context_rules.py
"""
from pathlib import Path
import numpy as np, pandas as pd
from scipy import stats
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
LIB = ROOT / "data/libraries/dp4_library.csv"
SUMM = ROOT / "data/dp4_binding/john/scaffoldedEpitopeSummary.csv"
STRUCTF = ROOT / "data/dp4_binding/epitope_struct_features.csv"
OUT = ROOT / "data/dp4_binding/figs/context_rules.png"

METRICS = ["cylinder_clashes","epitope_rmsd","overall_rmsd","epitope_pae","scaffold_pae",
           "mean_pae","ptm","af3_clashes","composite"]

def within_corr(d, metric):
    d = d.dropna(subset=[metric])
    if d.epitope.nunique() < 3: return np.nan
    mc = d[metric] - d.groupby("epitope")[metric].transform("mean")
    bc = d["bound"] - d.groupby("epitope")["bound"].transform("mean")
    if mc.std() == 0 or bc.std() == 0: return np.nan
    return stats.spearmanr(mc, bc).correlation

def main():
    lib = pd.read_csv(LIB, low_memory=False)
    summ = pd.read_csv(SUMM, encoding="latin-1"); summ.columns=[c.strip() for c in summ.columns]
    summ = summ.rename(columns={"#islands":"n_islands"}); summ["epitope"]=summ.epitope.str.lower()
    st = pd.read_csv(STRUCTF)
    epi = summ[["epitope","helix","n_islands"]].merge(st[["epitope","epi_dsasa","epi_bfac_z"]], on="epitope", how="inner")
    hits=set()
    for ids in summ["hitIDs"].dropna():
        for t in str(ids).split(","):
            t=t.strip()
            if t.startswith("DP4_"): hits.add("DP4_"+t.split("_")[1])
    lib["epitope"]=lib["target"].astype(str).str.split("_").str[0].str.lower()
    d=lib[lib.composite.notna() & lib.epitope.isin(set(epi.epitope))].copy()
    d["bound"]=d.library_member.isin(hits).astype(int)
    d=d.merge(epi,on="epitope",how="left")

    hmed=d.helix.median(); smed=d.epi_dsasa.median(); rmed=d.epi_bfac_z.median()
    strata = {
        "all": d,
        "loopy\n(helix<med)":  d[d.helix<hmed],
        "helical\n(helix≥med)": d[d.helix>=hmed],
        "small iface\n(ΔSASA<med)": d[d.epi_dsasa<smed],
        "big iface\n(ΔSASA≥med)":   d[d.epi_dsasa>=smed],
        "rigid\n(Bfac<med)":   d[d.epi_bfac_z<rmed],
        "flexible\n(Bfac≥med)":d[d.epi_bfac_z>=rmed],
        "1 island":            d[d.n_islands==1],
        "2 islands":           d[d.n_islands==2],
    }
    M = np.full((len(METRICS), len(strata)), np.nan)
    for j,(sname,sub) in enumerate(strata.items()):
        for i,met in enumerate(METRICS):
            M[i,j]=within_corr(sub, met)
    tab = pd.DataFrame(M, index=METRICS, columns=list(strata))
    pd.set_option("display.width",200,"display.max_columns",20)
    print("within-antibody corr of each design metric with binding, by context (neg = lower→binds):\n")
    print(tab.round(2).to_string())

    fig,ax=plt.subplots(figsize=(12,5.5))
    im=ax.imshow(M, cmap="RdBu_r", vmin=-0.4, vmax=0.4, aspect="auto")
    ax.set_xticks(range(len(strata))); ax.set_xticklabels(list(strata), fontsize=8)
    ax.set_yticks(range(len(METRICS))); ax.set_yticklabels(METRICS, fontsize=9)
    for i in range(len(METRICS)):
        for j in range(len(strata)):
            if not np.isnan(M[i,j]): ax.text(j,i,f"{M[i,j]:+.2f}",ha="center",va="center",fontsize=7,
                                             color="white" if abs(M[i,j])>0.22 else "0.3")
    ax.axvline(0.5,color="k",lw=1)
    cb=fig.colorbar(im,ax=ax,shrink=0.8); cb.set_label("within-antibody Spearman with binding")
    ax.set_title("Which design metric predicts binding, by epitope context\n(does the rule change across columns?)",
                 fontsize=12, fontweight="bold")
    fig.tight_layout(); OUT.parent.mkdir(parents=True,exist_ok=True); fig.savefig(OUT,dpi=150)
    print(f"\nwrote {OUT.relative_to(ROOT)}")

if __name__ == "__main__":
    main()
