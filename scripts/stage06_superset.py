#!/usr/bin/env python3
"""Build the ALL-DESIGNS superset of the DP4 scaffolded arms (John, 2026-07-16).

The shipped dp4_library.csv holds only the SELECTED designs; this emits every candidate
design in the scaffolded pools with its metrics, so the overall distributions can be poked at.
Metrics only -- the 103-mer sequences were case-encoded only for the selected subset, so unselected
designs have no sequence (that would need folding-out ~335k constructs). One row per AF3 design:

  component            C1 / C2 / C3 (scaffoldedAbEpitope / SingleIsland / Polyclonal)
  target, id, island_index, predID                 identity
  epitope_chunk_rmsd, overall_rmsd, epitope_pae, mean_pae, ptm,
  af3_n_clash_res, cylinder_native_aware           the metrics (blank where a component lacks one)
  composite, rank_in_group                         shipped-preset composite + within-group rank
  is_global_pass                                   clears ALL four Lawson filters (C1/C2 only)

C1 and C3 pools are local; C2 (metrics_dual_island.parquet) and 8VDL live on the cluster -- run this
there for those and concat. Usage (one component at a time; --append to add to an existing file):

  python scripts/stage06_superset.py --component C1 \
    --metrics-csv $D/known_antigen/analysis/data/metrics_whole_epitope_103.csv \
    --out data/libraries/dp4_superset_metrics.csv
  python scripts/stage06_superset.py --component C3 \
    --metrics-csv $D/12mer_tiling/analysis/data/metrics_12mer.csv --append \
    --out data/libraries/dp4_superset_metrics.csv
"""
from __future__ import annotations

import argparse
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from episcaf_analysis.presets import PRESETS
from episcaf_analysis.score import score

# component -> (category label, scoring preset, group cols, has real antibody clash)
COMPONENTS = {
    "C1": ("scaffoldedAbEpitope",   "antibody",  ["id"],                True),
    "C2": ("scaffoldedSingleIsland","antibody",  ["id", "island_index"], True),
    "C3": ("scaffoldedPolyclonal",  "twelvemer", ["antigen", "id"],      False),
}
OUT_COLS = ["component", "category", "target", "id", "island_index", "predID",
            "epitope_chunk_rmsd", "overall_rmsd", "epitope_pae", "mean_pae", "ptm",
            "af3_n_clash_res", "cylinder_native_aware",
            "composite", "rank_in_group", "is_global_pass"]
# the four-filter (global mean PAE, exactly as sec:fourfilter)
FOUR_FILTER = {"epitope_chunk_rmsd": 1.0, "overall_rmsd": 2.0, "mean_pae": 5.0, "af3_n_clash_res": 0.0}


def build(component: str, path: str) -> pd.DataFrame:
    cat, preset_name, gcols, has_ab = COMPONENTS[component]
    df = pd.read_parquet(path) if path.endswith(".parquet") else pd.read_csv(path, low_memory=False)
    if "status" in df.columns:
        df = df[df["status"].astype(str).str.lower().eq("ok")].copy()
    print(f"[{component}] {len(df)} designs from {path}")

    # composite under the shipped preset (rank ALL rows; no top-k cut here)
    preset = {k: (v.copy() if isinstance(v, dict) else v) for k, v in PRESETS[preset_name].items()}
    preset["select"] = None
    scored = score(df.copy(), preset)

    # rank within the shipped selection group
    gcols = [g for g in gcols if g in scored.columns] or [gcols[0]]
    if len(gcols) > 1:
        scored["_grp"] = scored[gcols].astype(str).agg("|".join, axis=1)
        gkey = "_grp"
    else:
        gkey = gcols[0]
    scored["rank_in_group"] = scored.groupby(gkey)["composite"].rank(ascending=False, method="first")

    # global four-filter pass (known-antibody arms only)
    if has_ab and all(c in scored.columns for c in FOUR_FILTER):
        ok = pd.Series(True, index=scored.index)
        for c, t in FOUR_FILTER.items():
            x = pd.to_numeric(scored[c], errors="coerce")
            ok &= (x == 0) if t == 0.0 else (x <= t)
        scored["is_global_pass"] = ok
    else:
        scored["is_global_pass"] = pd.NA

    scored["component"] = component
    scored["category"] = cat
    scored["target"] = scored["id"] if "id" in scored.columns else pd.NA
    for c in OUT_COLS:
        if c not in scored.columns:
            scored[c] = pd.NA
    return scored[OUT_COLS]


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--component", required=True, choices=sorted(COMPONENTS))
    ap.add_argument("--metrics-csv", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--append", action="store_true", help="append to --out instead of overwriting")
    args = ap.parse_args()

    out = build(args.component, args.metrics_csv)
    header = not (args.append and os.path.exists(args.out))
    out.to_csv(args.out, mode="a" if args.append else "w", header=header, index=False)
    gp = out["is_global_pass"]
    print(f"[{args.component}] wrote {len(out)} rows -> {args.out}"
          f"  (global-pass: {int(gp.sum()) if gp.notna().any() else 'n/a'})")


if __name__ == "__main__":
    main()
