#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path
import pandas as pd
import numpy as np

from lawson_rmsd import overall_rmsd_mpnn_vs_af3_window, epitope_chunk_rmsd

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dp2", default="datasets/dp2.parquet", help="Path to dp2.parquet")
    ap.add_argument("--root", required=True, help="ROOT data dir (sourced_antibody_v1/no_antibody)")
    ap.add_argument("--n", type=int, default=0, help="If >0, limit rows (debug)")
    ap.add_argument("--seed", type=int, default=1, help="Random seed for sampling (only used if --sample >0)")
    ap.add_argument("--sample", type=int, default=0, help="If >0, sample N rows from dp2 instead of full")
    ap.add_argument("--out_csv", default="runs/dp2_lawson_validation.csv", help="Output CSV")
    ap.add_argument("--progress_every", type=int, default=200, help="Print progress every N rows")
    args = ap.parse_args()

    dp2_path = Path(args.dp2).resolve()
    root = Path(args.root).resolve()
    out_csv = Path(args.out_csv).resolve()
    out_csv.parent.mkdir(parents=True, exist_ok=True)

    dp2 = pd.read_parquet(dp2_path)

    required = [
        "id","contig_id","rfd_id","mpnn_id","assay_scaffolded_epitope_id",
        "scaffolded_epitope_chunk_resindices","assay_scaffolded_epitope_chunk_resindices",
    ]
    missing = [c for c in required if c not in dp2.columns]
    if missing:
        raise SystemExit(f"dp2 missing required cols: {missing}")

    # Identify stored metric column names (dp2 uses these names in your printout)
    if "overall_rmsd" not in dp2.columns:
        raise SystemExit("dp2 missing 'overall_rmsd' column (stored overall rmsd)")
    if "epitope_chunk_rmsd_vs_mpnn" not in dp2.columns:
        raise SystemExit("dp2 missing 'epitope_chunk_rmsd_vs_mpnn' column (stored epitope rmsd)")

    # Optional: limit/sample
    if args.sample and args.sample > 0:
        dp2 = dp2.sample(n=min(args.sample, len(dp2)), random_state=args.seed).reset_index(drop=True)
    if args.n and args.n > 0:
        dp2 = dp2.head(args.n).reset_index(drop=True)

    # Normalize token
    dp2["assay_scaffolded_epitope_id"] = dp2["assay_scaffolded_epitope_id"].astype(str).str.lower()

    rows = []
    file_missing = 0
    seqmatch_fail = 0
    ok = 0

    for i, r in dp2.iterrows():
        pid = str(r["id"])
        contig = int(r["contig_id"])
        rfd_id = int(r["rfd_id"])
        mpnn = int(r["mpnn_id"])
        tok = str(r["assay_scaffolded_epitope_id"]).lower()

        mpnn_pdb = root / "proteinmpnn" / pid / str(contig) / f"{pid}_{rfd_id}_fixed_dldesign_{mpnn}.pdb"
        af3_cif = root / "af3_predictions" / tok / "seed-1_sample-0" / "model.cif.gz"
        if not af3_cif.exists():
            af3_cif = root / "af3_predictions" / tok / f"{tok}_model.cif.gz"

        out = {
            "id": pid,
            "contig_id": contig,
            "rfd_id": rfd_id,
            "mpnn_id": mpnn,
            "tok": tok,
            "mpnn_pdb": str(mpnn_pdb),
            "af3_cif": str(af3_cif),
            "pdb_ok": bool(mpnn_pdb.exists()),
            "cif_ok": bool(af3_cif.exists()),
            "dp2_overall_rmsd": float(r["overall_rmsd"]) if pd.notna(r["overall_rmsd"]) else np.nan,
            "dp2_epi_rmsd": float(r["epitope_chunk_rmsd_vs_mpnn"]) if pd.notna(r["epitope_chunk_rmsd_vs_mpnn"]) else np.nan,
            "calc_overall_rmsd": np.nan,
            "calc_epi_rmsd": np.nan,
            "d_overall": np.nan,
            "d_epi": np.nan,
            "af3_window_start": np.nan,
            "af3_window_end": np.nan,
            "status": "ok",
        }

        if not (mpnn_pdb.exists() and af3_cif.exists()):
            out["status"] = "file_missing"
            file_missing += 1
            rows.append(out)
            continue

        try:
            overall, start, end = overall_rmsd_mpnn_vs_af3_window(mpnn_pdb, af3_cif)
            out["calc_overall_rmsd"] = float(overall)
            out["af3_window_start"] = int(start)
            out["af3_window_end"] = int(end)
        except Exception:
            out["status"] = "seqmatch_fail"
            seqmatch_fail += 1
            rows.append(out)
            continue

        try:
            mpnn_ris = list(r["scaffolded_epitope_chunk_resindices"])
            af3_ris  = list(r["assay_scaffolded_epitope_chunk_resindices"])
            epi = epitope_chunk_rmsd(mpnn_pdb, af3_cif, mpnn_ris, af3_ris)
            out["calc_epi_rmsd"] = float(epi)
        except Exception:
            out["status"] = "epi_fail"
            rows.append(out)
            continue

        # diffs only meaningful if dp2 stored values are finite
        if np.isfinite(out["dp2_overall_rmsd"]) and np.isfinite(out["calc_overall_rmsd"]):
            out["d_overall"] = float(abs(out["dp2_overall_rmsd"] - out["calc_overall_rmsd"]))
        if np.isfinite(out["dp2_epi_rmsd"]) and np.isfinite(out["calc_epi_rmsd"]):
            out["d_epi"] = float(abs(out["dp2_epi_rmsd"] - out["calc_epi_rmsd"]))

        ok += 1
        rows.append(out)

        if args.progress_every > 0 and (i + 1) % args.progress_every == 0:
            print(f"[progress] {i+1}/{len(dp2)} ok={ok} file_missing={file_missing} seqmatch_fail={seqmatch_fail}")

    df = pd.DataFrame(rows)
    df.to_csv(out_csv, index=False)

    print("\n=== DONE ===")
    print("tested rows:", len(df))
    print("file_missing:", int((df['status'] == 'file_missing').sum()))
    print("seqmatch_fail:", int((df['status'] == 'seqmatch_fail').sum()))
    print("epi_fail:", int((df['status'] == 'epi_fail').sum()))
    # diffs summary (only where diffs are present)
    dd = df.loc[np.isfinite(df["d_overall"]) & np.isfinite(df["d_epi"]), ["d_overall","d_epi"]]
    print("\nabs-diff summary:")
    if len(dd) == 0:
        print("No rows with finite diffs (check dp2 stored cols).")
    else:
        print(dd.describe().to_string())
    print("\nwrote:", out_csv)

if __name__ == "__main__":
    main()
