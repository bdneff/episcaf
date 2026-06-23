#!/usr/bin/env python3
"""
stage05_summarize.py -- report on the stage05 metrics table (named, committed, so the
numbers are reproducible instead of pasted one-liners). Reads the parquet, prints a
status/metric/feasibility report, and optionally writes a per-island gate summary to
results/ for the manuscript to cite.

Usage:
  python3 scripts/stage05_summarize.py \
      --metrics runs/dual_island_rfd3/05_analysis/metrics_dual_island.parquet \
      --out_island_csv results/dual_island_gate_summary.csv

The "gate" is the composite scorer's gate (epitope_chunk_rmsd <= GATE, default 2.5); the
"four-filter" is the DP3 ground-truth pass (overall<=2, epitope<=1, mean_pae<5, clash==0).
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

METRICS = ["overall_rmsd", "epitope_chunk_rmsd", "mean_pae", "epitope_pae", "scaffold_pae",
           "ptm", "af3_n_clash_res", "cylinder_ca_clashes", "cylinder_native_aware",
           "native_in_cylinder"]


def num(s):
    return pd.to_numeric(s, errors="coerce")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--metrics", required=True, help="stage05 .parquet (or .csv)")
    ap.add_argument("--gate", type=float, default=2.5, help="scorer gate on epitope_chunk_rmsd")
    ap.add_argument("--topk", type=int, default=5, help="designs needed per island for selection")
    ap.add_argument("--out_island_csv", default=None, help="write per-island gate summary here")
    args = ap.parse_args()

    p = Path(args.metrics)
    df = pd.read_parquet(p) if p.suffix == ".parquet" else pd.read_csv(p)

    n_isl = df.groupby(["id", "island_index"]).ngroups
    print(f"rows: {len(df):,}   epitopes: {df.id.nunique()}   islands: {n_isl}")

    print("\n=== status ===")
    print(df.status.value_counts().to_string())
    ok = df[df.status == "ok"].copy()
    print(f"\nstatus==ok: {len(ok):,}/{len(df):,}")

    # Why clash/cylinder may be null -- the accessibility path has its own status.
    if "af3_clash_status" in df.columns:
        print("\n=== af3_clash_status (explains null clash/cylinder) ===")
        print(df.af3_clash_status.value_counts().to_string())
        bad = df[df.af3_clash_status != "ok"]
        if len(bad):
            print("\nepitopes with any non-ok clash status (id: n_designs):")
            print(bad.groupby("id").size().sort_values(ascending=False).head(12).to_string())

    print("\n=== metric summary (status==ok) ===")
    print(f"{'metric':22s} {'median':>8s} {'min':>8s} {'max':>8s} {'nonnull':>10s}")
    for c in METRICS:
        if c in ok.columns:
            v = num(ok[c])
            print(f"{c:22s} {v.median():8.2f} {v.min():8.2f} {v.max():8.2f} {v.notna().sum():10,}")

    print("\n=== per-filter individual pass counts (status==ok) ===")
    e = num(ok.epitope_chunk_rmsd); o = num(ok.overall_rmsd)
    mp = num(ok.mean_pae); cl = num(ok.af3_n_clash_res)
    print(f"  scorer gate epitope_rmsd<={args.gate:<4g} : {(e <= args.gate).sum():>8,}  "
          f"({100*(e <= args.gate).mean():.1f}%)")
    print(f"  epitope_rmsd<=1                : {(e <= 1).sum():>8,}")
    print(f"  overall_rmsd<=2               : {(o <= 2).sum():>8,}")
    print(f"  mean_pae<5                    : {(mp < 5).sum():>8,}")
    print(f"  clash==0                      : {(cl == 0).sum():>8,}  (of {cl.notna().sum():,} with clash)")
    four = (o <= 2) & (e <= 1) & (mp < 5) & (cl == 0)
    print(f"  DP3 four-filter (real clash)  : {four.sum():>8,}  ({100*four.mean():.2f}%)")

    # Top-k feasibility: how many islands have >= topk designs past the gate.
    print(f"\n=== top-{args.topk}/island feasibility (gate epitope_rmsd<={args.gate}) ===")
    ok["pass_gate"] = e <= args.gate
    g = ok.groupby(["id", "island_index"])
    isl = g.agg(n_designs=("pass_gate", "size"),
                n_gated=("pass_gate", "sum"),
                island_size=("island_size", "first")).reset_index()
    enough = (isl.n_gated >= args.topk).sum()
    print(f"  islands with >= {args.topk} gated designs : {enough} / {len(isl)}")
    print(f"  islands with 0 gated designs       : {(isl.n_gated == 0).sum()}")
    short = isl[isl.n_gated < args.topk].sort_values("n_gated")
    if len(short):
        print(f"\n  islands SHORT of {args.topk} (will yield fewer; note honestly):")
        for _, r in short.iterrows():
            print(f"    {r.id} island {int(r.island_index)} (size {int(r.island_size)}): "
                  f"{int(r.n_gated)} gated of {int(r.n_designs)}")

    if args.out_island_csv:
        out = Path(args.out_island_csv)
        out.parent.mkdir(parents=True, exist_ok=True)
        isl.sort_values(["id", "island_index"]).to_csv(out, index=False)
        print(f"\nwrote per-island gate summary -> {out}")


if __name__ == "__main__":
    main()
