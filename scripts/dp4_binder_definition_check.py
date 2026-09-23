#!/usr/bin/env python3
"""Provenance check: exactly what defines the 217 binders and the 8.0% base rate.

John asked what criteria produced "217 bind, base rate 8.0%". This script traces it end to end so the
answer is exact, not remembered. The binder label is NOT a threshold we invented: a design is a binder
iff its library_member appears in that epitope's hitIDs column in John's scaffoldedEpitopeSummary.csv
(his own PepSeq hit calls). The 8% denominator is the analysis universe = scaffolded designs that carry a
composite score AND belong to an epitope for which we extracted native-context structural features.

Run:  /usr/bin/python3 scripts/dp4_binder_definition_check.py
"""
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
LIB = ROOT / "data/libraries/dp4_library.csv"
SUMM = ROOT / "data/dp4_binding/john/scaffoldedEpitopeSummary.csv"
STRUCTF = ROOT / "data/dp4_binding/epitope_struct_features.csv"

def parse_hits(summ):
    hits = set()
    for ids in summ["hitIDs"].dropna():
        for t in str(ids).split(","):
            t = t.strip()
            if t.startswith("DP4_"):
                hits.add("DP4_" + t.split("_")[1])
    return hits

def main():
    lib = pd.read_csv(LIB, low_memory=False)
    summ = pd.read_csv(SUMM, encoding="latin-1"); summ.columns = [c.strip() for c in summ.columns]
    summ["epitope"] = summ["epitope"].astype(str).str.lower()
    st = pd.read_csv(STRUCTF)
    lib["epitope"] = lib["target"].astype(str).str.split("_").str[0].str.lower()

    hits = parse_hits(summ)
    print(f"John's summary: {len(summ)} epitope rows")
    if "#hits" in summ.columns:
        print(f"  sum of the '#hits' column           = {pd.to_numeric(summ['#hits'],errors='coerce').sum():.0f}")
    print(f"  distinct hit library_members parsed  = {len(hits)}   (a design is a 'binder' iff it is in this set)")

    scaff = lib[lib["composite"].notna()].copy()               # scaffolded arms carry a composite
    print(f"\nScaffolded designs with a composite score (metrics universe): {len(scaff):,}")
    print(f"  of those, in John's hitIDs (binders)                       : {scaff['library_member'].isin(hits).sum()}"
          f"  ({scaff['library_member'].isin(hits).mean():.1%})")

    epi56 = set(st["epitope"])
    sub = scaff[scaff["epitope"].isin(epi56)].copy()
    sub["bound"] = sub["library_member"].isin(hits).astype(int)
    print(f"\nAnalysis universe used in the ML (also restricted to the {len(epi56)} epitopes with "
          f"structural features):")
    print(f"  designs = {len(sub):,}   binders = {int(sub['bound'].sum())}   base rate = {sub['bound'].mean():.1%}")

    # where do John's hits sit relative to our universe?
    in_scaff = scaff["library_member"].isin(hits).sum()
    print(f"\nAccounting for all {len(hits)} of John's distinct hits:")
    print(f"  in the scaffolded+composite universe : {in_scaff}")
    print(f"  of those, in the {len(epi56)}-epitope struct subset : {int(sub['bound'].sum())}")
    print(f"  hits NOT in the scaffolded+composite universe (other arms / no metrics): {len(hits)-in_scaff}")

    # per-epitope binder counts within the analysis universe (sanity vs John's per-epitope #hits)
    by = sub.groupby("epitope")["bound"].sum().sort_values(ascending=False)
    print(f"\ntop epitopes by binders in the analysis universe:")
    print(by.head(8).to_string())
    print(f"\nepitopes with >=1 binder: {(by>0).sum()}/{len(by)}")

if __name__ == "__main__":
    main()
