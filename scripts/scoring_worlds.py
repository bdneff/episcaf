#!/usr/bin/env python3
"""Compare scoring 'worlds' on one metrics pool: top-k median tables, per-target and pooled.

Built for John's question (2026-07-15): does saturating RMSD/PAE and penalizing clash surface
lower-clash mAb designs? Four worlds:

  percentile   current scorer (rank-based, scale-blind)
  sigmoid      naive: sigmoid transforms, OLD weights (.35/.35/.15/.15) -- known to backfire
  mixed        sigmoid transforms + clash up-weighted to .50 (preset `antibody_sigmoid`)
  gated        Lawson-style fold-quality FLOOR (gate overall_rmsd + epitope_pae) then rank
               clash + epitope RMSD -- rejects the misfolded low-clash designs `mixed` can grab

Transforms are absolute (sigmoid) or rank-based (percentile); see episcaf_analysis/score.py.
Runs on any pool. The shipped C1 is the native-103 redo: local
`known_antigen/analysis/data/metrics_whole_epitope_103.csv` (140,716) -- NOT metrics_native_cyl_full.csv
(that is the older 104-mer pool). C2's pool is the cluster `metrics_dual_island.parquet`.

  python scripts/scoring_worlds.py --metrics-csv <pool.csv|.parquet> --topk 10 --ids 6cyf_0P
"""
from __future__ import annotations

import argparse
import sys

import pandas as pd

sys.path.insert(0, __file__.rsplit("/scripts/", 1)[0])
from episcaf_analysis.score import score

CLASH, ERMSD, ORMSD, EPAE = "af3_n_clash_res", "epitope_chunk_rmsd", "overall_rmsd", "epitope_pae"
COLS = [ERMSD, ORMSD, EPAE, CLASH]


def sig(mid, k):
    return dict(better="low", transform="sigmoid", midpoint=mid, k=k)


def world_presets(group, topk, gate_overall, gate_pae):
    pct = {c: dict(weight=w, better="low", transform="percentile")
           for c, w in [(CLASH, .35), (ERMSD, .35), (ORMSD, .15), (EPAE, .15)]}
    sigm = {CLASH: {**sig(6, .5), "weight": .35}, ERMSD: {**sig(1, 4), "weight": .35},
            ORMSD: {**sig(2, 2), "weight": .15}, EPAE: {**sig(5, 1), "weight": .15}}
    mixed = {CLASH: {**sig(6, .5), "weight": .50}, ERMSD: {**sig(1, 4), "weight": .25},
             ORMSD: {**sig(2, 2), "weight": .15}, EPAE: {**sig(5, 1), "weight": .10}}
    # gated: fold-quality floor handled by pre-filter (below); rank only clash + epitope RMSD
    gated = {CLASH: {**sig(6, .5), "weight": .60}, ERMSD: {**sig(1, 4), "weight": .40}}
    base = dict(scope="pooled", antigen_col="antigen", select=dict(group=group, topk=topk))
    return {
        "percentile (current)":            dict(base, gate=None, metrics=pct),
        "sigmoid naive (.35/.35/.15/.15)":  dict(base, gate=None, metrics=sigm),
        "mixed (sig + clash .50)":          dict(base, gate=None, metrics=mixed),
        f"gated (fold floor o<={gate_overall} pae<{gate_pae}, rank clash+rmsd)":
                                            dict(base, gate=None, metrics=gated, _prefilter=True),
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--metrics-csv", required=True)
    ap.add_argument("--group", default="id")
    ap.add_argument("--topk", type=int, default=10)
    ap.add_argument("--ids", nargs="*", default=None)
    ap.add_argument("--id-col", default="id")
    ap.add_argument("--gate-overall", type=float, default=3.0, help="gated world: overall_rmsd <= this")
    ap.add_argument("--gate-pae", type=float, default=8.0, help="gated world: epitope_pae < this")
    args = ap.parse_args()

    df = pd.read_parquet(args.metrics_csv) if args.metrics_csv.endswith(".parquet") \
        else pd.read_csv(args.metrics_csv, low_memory=False)
    if "status" in df.columns:
        df = df[df["status"].astype(str).str.lower().eq("ok")].copy()
    # multi-column grouping (e.g. C2 selects per (id, island_index)): build a combined key column,
    # since score.py groups by a single column. Single-column groups pass through unchanged.
    gcols = [g for g in args.group.split(",") if g in df.columns]
    if not gcols:
        sys.exit(f"[worlds] none of the group columns {args.group!r} are in the CSV: {list(df.columns)}")
    if len(gcols) > 1:
        group = "_grp"
        df[group] = df[gcols].astype(str).agg("|".join, axis=1)
    else:
        group = gcols[0]
    print(f"loaded {len(df)} designs from {args.metrics_csv}  (group by {gcols} -> {group})")

    worlds = world_presets(group, args.topk, args.gate_overall, args.gate_pae)
    out = {}
    for name, p in worlds.items():
        d = df
        if p.pop("_prefilter", False):
            d = df[(pd.to_numeric(df[ORMSD], errors="coerce") <= args.gate_overall)
                   & (pd.to_numeric(df[EPAE], errors="coerce") < args.gate_pae)].copy()
            print(f"  [{name}] fold-quality floor kept {len(d)}/{len(df)} designs")
        out[name] = score(d.copy(), p)

    scopes = [(f"{i} only", i) for i in (args.ids or [])] + [("ALL targets pooled", None)]
    for label, idf in scopes:
        print(f"\n===== top-{args.topk} medians — {label} =====")
        print(f"{'world':48} {'epi_rmsd':>9} {'ovr_rmsd':>9} {'epi_pae':>8} {'clash':>7} {'n':>5}")
        for name, o in out.items():
            oo = o if idf is None else o[o[args.id_col].astype(str) == str(idf)]
            m = {c: pd.to_numeric(oo[c], errors="coerce").median() for c in COLS}
            print(f"{name:48} {m[ERMSD]:>9.2f} {m[ORMSD]:>9.2f} {m[EPAE]:>8.2f} {m[CLASH]:>7.1f} {len(oo):>5}")


if __name__ == "__main__":
    main()
