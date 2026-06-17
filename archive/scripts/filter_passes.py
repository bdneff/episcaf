#!/usr/bin/env python3
"""
filter_passes.py

Read metrics_by_pred.csv and write out all rows that pass user-specified filters.
Also reports pass rate among "finished" rows (rows that have all required metrics present).

Example:
  python scripts/filter_passes.py \
    --csv runs/run_test_rfd3_nompmn/04_filter/metrics_by_pred.csv \
    --out runs/run_test_rfd3_nompmn/04_filter/passes_strict.csv \
    --pae_max 5 --rmsd_all_max 2 --rmsd_epi_max 1
"""

import argparse
from pathlib import Path

import pandas as pd


REQUIRED_METRICS = ["pae_mean_all", "rmsd_ca_all", "rmsd_ca_epitope"]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv", required=True, help="Path to metrics_by_pred.csv")
    ap.add_argument("--out", required=True, help="Output CSV path for passing rows")
    ap.add_argument("--pae_max", type=float, default=5.0, help="Require pae_mean_all < pae_max")
    ap.add_argument("--rmsd_all_max", type=float, default=2.0, help="Require rmsd_ca_all < rmsd_all_max")
    ap.add_argument("--rmsd_epi_max", type=float, default=1.0, help="Require rmsd_ca_epitope < rmsd_epi_max")
    ap.add_argument("--sort_by", default="pae_mean_all", help="Column to sort passes by (ascending)")
    ap.add_argument("--top", type=int, default=25, help="Print top N passing rows to stdout")
    args = ap.parse_args()

    df = pd.read_csv(args.csv)

    missing_cols = [c for c in REQUIRED_METRICS if c not in df.columns]
    if missing_cols:
        raise SystemExit(f"ERROR: missing required columns in CSV: {missing_cols}\n"
                         f"Columns present: {list(df.columns)}")

    # "Finished" = has all required metrics (i.e., AF3/RFD3 present and RMSDs computed)
    finished_mask = df[REQUIRED_METRICS].notna().all(axis=1)
    finished = df[finished_mask].copy()

    # Apply strict filters ONLY within finished rows
    pass_mask_finished = (
        (finished["pae_mean_all"] < args.pae_max) &
        (finished["rmsd_ca_all"] < args.rmsd_all_max) &
        (finished["rmsd_ca_epitope"] < args.rmsd_epi_max)
    )
    passed = finished[pass_mask_finished].copy()

    # Sort (if column exists)
    if args.sort_by in passed.columns:
        passed = passed.sort_values(args.sort_by, ascending=True)

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    passed.to_csv(out_path, index=False)

    n_total = len(df)
    n_finished = len(finished)
    n_pass = len(passed)

    overall_rate = 100.0 * n_pass / n_total if n_total else 0.0
    finished_rate = 100.0 * n_pass / n_finished if n_finished else 0.0

    print(f"Total rows (design×pred): {n_total:,}")
    print(f"Finished rows (have PAE + RMSDs): {n_finished:,}")
    print(f"Passed filters: {n_pass:,} / {n_finished:,} ({finished_rate:.4f}% among finished)")
    print(f"Overall: {n_pass:,} / {n_total:,} ({overall_rate:.6f}% of all rows)")
    print(f"Output written to: {out_path}")

    if n_pass and args.top > 0:
        cols = [
            "assay_scaffolded_epitope_id",
            "pred",
            "pae_mean_all",
            "rmsd_ca_all",
            "rmsd_ca_epitope",
            "n_ca_epitope_common",
            "ptm",
            "ranking_score",
            "rfd3_path",
            "af3_path",
        ]
        cols = [c for c in cols if c in passed.columns]
        #print(f"\nTop {min(args.top, n_pass)} passing rows:")
        #print(passed[cols].head(args.top).to_string(index=False))


if __name__ == "__main__":
    main()
