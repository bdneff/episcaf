#!/usr/bin/env python3
"""
stage06_select.py -- composite selection for a DP4 component.

Scores a metrics table with a presets.py preset (episcaf_analysis.score.score) and emits the
FULL ranking: every design keeps its composite plus rank_in_group (1 = best in its group). We do
NOT hard-cut to top-k here on purpose -- selection depth is the elastic buffer decided last, once
the DP4 budget (minibinder count) is fixed (memory: dp4-budget). Apply the top-n cut downstream at
library assembly (filter rank_in_group <= n). Pass --topk to also write a pre-cut convenience copy.

Grouping supports MULTIPLE columns (score.py groups on one), which C2 needs: the single-island
deliverable is top-n per (id, island_index), not per id.

Usage:
  # C1 whole-epitope known-Ab (local): top-n per epitope
  python scripts/stage06_select.py --preset antibody \
      --metrics-csv known_antigen/analysis/data/metrics_native_cyl_full.csv \
      --group id --out results/dp4_C1_whole_epitope_ranked.csv

  # C2 single-island known-Ab (RUN ON GEMINI against the cluster metrics): top-n per island
  python scripts/stage06_select.py --preset antibody \
      --metrics-csv <run>/05_analysis/metrics_dual_island.parquet \
      --group id,island_index --out results/dp4_C2_single_island_ranked.csv

  # C3 scaffolded 12mer tiles (local): top-n per epitope within antigen
  python scripts/stage06_select.py --preset twelvemer \
      --metrics-csv 12mer_tiling/analysis/data/metrics_12mer.csv \
      --group antigen,id --out results/dp4_C3_12mer_ranked.csv
"""
from __future__ import annotations
import argparse, sys
from pathlib import Path
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
from episcaf_analysis.score import score            # noqa: E402
from episcaf_analysis.presets import PRESETS         # noqa: E402


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--preset", required=True, choices=sorted(PRESETS))
    ap.add_argument("--metrics-csv", required=True, help=".csv or .parquet metrics table")
    ap.add_argument("--group", required=True, help="comma list of grouping columns (e.g. id or id,island_index)")
    ap.add_argument("--out", required=True)
    ap.add_argument("--topk", type=int, default=0, help="if >0, also write <out>.top<k>.csv pre-cut")
    args = ap.parse_args()

    p = Path(args.metrics_csv)
    df = pd.read_parquet(p) if p.suffix == ".parquet" else pd.read_csv(p, low_memory=False)
    group = [g.strip() for g in args.group.split(",")]
    missing = [g for g in group if g not in df.columns]
    if missing:
        sys.exit(f"[select] grouping columns absent from metrics: {missing}\n  have: {list(df.columns)}")

    # score with the preset's transforms/weights but DISABLE its single-col top-k select,
    # so we get the composite on every row and do our own multi-col grouping + ranking.
    preset = {k: (v.copy() if isinstance(v, dict) else v) for k, v in PRESETS[args.preset].items()}
    preset["select"] = None
    scored = score(df, preset)

    scored["rank_in_group"] = (scored.groupby(group, sort=False)["composite"]
                                     .rank(ascending=False, method="first").astype(int))
    scored = scored.sort_values(group + ["rank_in_group"]).reset_index(drop=True)

    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    scored.to_csv(args.out, index=False)
    ngroups = scored.groupby(group, sort=False).ngroups
    print(f"[select] preset={args.preset} group={group}  {len(scored):,} rows, {ngroups} groups "
          f"-> {args.out}")
    print(f"[select] group size (designs/group): min={scored.groupby(group).size().min()} "
          f"median={int(scored.groupby(group).size().median())} max={scored.groupby(group).size().max()}")
    if args.topk > 0:
        cut = scored[scored.rank_in_group <= args.topk]
        outk = f"{args.out.rsplit('.',1)[0]}.top{args.topk}.csv"
        cut.to_csv(outk, index=False)
        print(f"[select] pre-cut top-{args.topk}: {len(cut):,} rows -> {outk}")


if __name__ == "__main__":
    main()
