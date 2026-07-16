#!/usr/bin/env python3
"""Compare two scoring presets on the same metrics file, per target.

Built to answer John's question (2026-07-15): does switching RMSD/PAE to saturating sigmoids and
penalizing clash absolutely (preset `antibody_softgate`) surface lower-clash designs than the current
percentile scorer (preset `antibody`)? Shows the top-k per target under each preset side by side, so
the reordering — and any designs the new preset pulls in from outside the current top-k — is visible.

Runs wherever the full candidate metrics file lives. NOTE the shipped C1 is the native-103 redo, whose
pool is on the cluster (`metrics_whole_epitope.csv`), NOT the local 104-mer `metrics_native_cyl_full.csv`
(that is the older RFD-comparison pool, 150,948 rows, and would demo on stale designs). C2's pool
(`metrics_dual_island.parquet`) is on the cluster too.

Examples (on the cluster, from the repo root):
  python scripts/compare_scoring.py \
    --metrics-csv runs/whole_epitope_rfd3/05_analysis/metrics_whole_epitope.csv \
    --preset-a antibody --preset-b antibody_softgate --group id --topk 10 --ids 6cyf_0P
  python scripts/compare_scoring.py \
    --metrics-csv runs/dual_island_rfd3/05_analysis/metrics_dual_island.parquet \
    --group id,island_index --topk 5 --ids 6cyf_0P
"""
from __future__ import annotations

import argparse
import sys

import pandas as pd

sys.path.insert(0, __file__.rsplit("/scripts/", 1)[0])
from episcaf_analysis.presets import PRESETS
from episcaf_analysis.score import score

SHOW = ["epitope_chunk_rmsd", "overall_rmsd", "epitope_pae", "af3_n_clash_res"]


def load(path: str) -> pd.DataFrame:
    df = pd.read_parquet(path) if path.endswith(".parquet") else pd.read_csv(path, low_memory=False)
    if "status" in df.columns:
        df = df[df["status"].astype(str).str.lower().eq("ok")].copy()
    return df


def preset_with(name: str, group: str, topk: int) -> dict:
    p = {k: (v.copy() if isinstance(v, dict) else v) for k, v in PRESETS[name].items()}
    p["select"] = dict(group=group, topk=topk)
    return p


def show(df: pd.DataFrame, title: str) -> None:
    print(f"\n{title}")
    hdr = f"  {'rk':>3} " + " ".join(f"{c.replace('epitope_chunk_','epi_').replace('af3_n_clash_res','clash'):>10}" for c in SHOW) + f" {'composite':>9}"
    print(hdr)
    for i, r in enumerate(df.itertuples(), 1):
        vals = " ".join(f"{getattr(r, c):>10.2f}" for c in SHOW)
        print(f"  {i:>3} {vals} {r.composite:>9.3f}")
    cl = df["af3_n_clash_res"]
    print(f"  --> clash  min {cl.min():.0f}  median {cl.median():.0f}  max {cl.max():.0f}")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--metrics-csv", required=True, help="full candidate pool (.csv or .parquet)")
    ap.add_argument("--preset-a", default="antibody", choices=sorted(PRESETS))
    ap.add_argument("--preset-b", default="antibody_softgate", choices=sorted(PRESETS))
    ap.add_argument("--group", default="id", help="selection group column(s); comma-joined uses the first present")
    ap.add_argument("--topk", type=int, default=10)
    ap.add_argument("--ids", nargs="*", default=None, help="restrict the display to these id values")
    ap.add_argument("--id-col", default="id")
    args = ap.parse_args()

    df = load(args.metrics_csv)
    print(f"loaded {len(df)} designs from {args.metrics_csv}")
    group = next((g for g in args.group.split(",") if g in df.columns), args.group.split(",")[0])

    a = score(df.copy(), preset_with(args.preset_a, group, args.topk))
    b = score(df.copy(), preset_with(args.preset_b, group, args.topk))

    ids = args.ids or sorted(df[args.id_col].astype(str).unique())
    for tid in ids:
        ta = a[a[args.id_col].astype(str) == str(tid)]
        tb = b[b[args.id_col].astype(str) == str(tid)]
        if not len(ta) and not len(tb):
            print(f"\n[{tid}] no rows"); continue
        print(f"\n{'='*72}\nTARGET {tid}\n{'='*72}")
        show(ta, f"[{tid}]  A = {args.preset_a}  (current)  top-{args.topk}")
        show(tb, f"[{tid}]  B = {args.preset_b}  (candidate)  top-{args.topk}")
        # which of B's picks were NOT in A's top-k (true re-selection)?
        key = "predID" if "predID" in df.columns else None
        if key:
            gained = set(tb[key]) - set(ta[key])
            print(f"\n  B pulled in {len(gained)}/{len(tb)} designs that were outside A's top-{args.topk}")


if __name__ == "__main__":
    main()
