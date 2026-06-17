#!/usr/bin/env python3
"""
02_compile_contigs.py

Reads the DP2 "superset" parquet snapshot and produces a compact contig table for RFD3.

Design goals:
- Deterministic, readable output
- Uses contig_string already present in the parquet (no regeneration)
- Adds seeds / repeats so you can fan out RFD3 generations cleanly
"""

import argparse
import os
import pandas as pd


KEEP_COLS = [
    "id",
    "assay_scaffolded_epitope_id",
    "assay_scaffolded_epitope_seq",
    "contig_id",
    "contig_string",
    "contig_length",
    "epitope_resindices",
    "epitope_boolmask",
    "scaffolded_epitope_resindices",
    "assay_scaffolded_epitope_resindices",
    # useful for comparing to Lawson later:
    "overall_rmsd",
    "epitope_chunk_rmsd_vs_mpnn",
    "mean_pae",
    "mpnn_clash_resindices",
    "af3_clash_resindices",
]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in_parquet", required=True, help="DP2 parquet snapshot inside run/00_inputs/")
    ap.add_argument("--out_parquet", required=True, help="Output contigs parquet inside run/01_contigs/")
    ap.add_argument("--max_rows", type=int, default=0, help="If >0, take only first N rows (for quick tests)")
    ap.add_argument("--seeds", type=str, default="0,1,2,3", help="Comma-separated seeds to expand per contig")
    ap.add_argument("--designs_per_seed", type=int, default=1, help="Replicates per (contig, seed)")
    args = ap.parse_args()

    df = pd.read_parquet(args.in_parquet)

    # basic sanity checks
    required = ["contig_string", "contig_length", "id"]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns: {missing}")

    # keep only a compact set of columns that exist
    keep = [c for c in KEEP_COLS if c in df.columns]
    df = df[keep].copy()

    # optionally limit size for a quick test
    if args.max_rows and args.max_rows > 0:
        df = df.head(args.max_rows).copy()

    # expand by seeds and replicate count
    seeds = [int(x.strip()) for x in args.seeds.split(",") if x.strip() != ""]
    rows = []
    for _, r in df.iterrows():
        base_id = r.get("assay_scaffolded_epitope_id", None) or r["id"]
        contig_id = r.get("contig_id", None)
        for seed in seeds:
            for rep in range(args.designs_per_seed):
                design_id = f"{base_id}__contig{contig_id}__seed{seed}__rep{rep}"
                out = r.to_dict()
                out["design_id"] = design_id
                out["seed"] = seed
                out["rep"] = rep
                rows.append(out)

    out_df = pd.DataFrame(rows)

    # ensure output directory exists
    out_dir = os.path.dirname(args.out_parquet)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    out_df.to_parquet(args.out_parquet, index=False)

    # also write a lightweight CSV for eyeballing
    csv_path = args.out_parquet.replace(".parquet", ".csv")
    out_df[["design_id", "contig_string", "contig_length", "seed", "rep"]].to_csv(csv_path, index=False)

    print(f"Wrote: {args.out_parquet}")
    print(f"Wrote: {csv_path}")
    print(f"Rows: {len(out_df)} (from {len(df)} base contigs, seeds={seeds}, designs_per_seed={args.designs_per_seed})")


if __name__ == "__main__":
    main()
