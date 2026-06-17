#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import re
from pathlib import Path
from typing import Any, Optional, Tuple

import numpy as np
import pandas as pd


def find_square_pae_matrix(obj: Any, path: str = "") -> Optional[Tuple[str, np.ndarray]]:
    """
    Recursively look for a square numeric matrix under any key/path containing 'pae'.
    Returns (path, matrix) for the first hit, else None.
    """
    if isinstance(obj, dict):
        for k, v in obj.items():
            hit = find_square_pae_matrix(v, f"{path}.{k}" if path else k)
            if hit is not None:
                return hit
        return None

    if isinstance(obj, list):
        # matrix-like list-of-lists?
        if obj and isinstance(obj[0], list) and obj[0] and all(isinstance(x, (int, float)) for x in obj[0]):
            if "pae" in path.lower():
                try:
                    A = np.asarray(obj, dtype=float)
                    if A.ndim == 2 and A.shape[0] == A.shape[1] and A.shape[0] > 0:
                        return (path, A)
                except Exception:
                    return None
        # otherwise recurse
        for i, v in enumerate(obj):
            hit = find_square_pae_matrix(v, f"{path}[{i}]")
            if hit is not None:
                return hit
    return None


def mean_pae_from_af3_dir(af3_dir: Path) -> Tuple[float, str]:
    """
    Prefer *_confidences.json, then *_data.json, then *_summary_confidences.json.
    Returns (mean_pae, source_string).
    Raises if no PAE matrix found.
    """
    candidates = []
    candidates += sorted(af3_dir.glob("*_confidences.json"))
    candidates += sorted(af3_dir.glob("*_data.json"))
    candidates += sorted(af3_dir.glob("*_summary_confidences.json"))

    if not candidates:
        raise FileNotFoundError(f"No JSON candidates in {af3_dir}")

    last_err = None
    for jf in candidates:
        try:
            obj = json.loads(jf.read_text())
            hit = find_square_pae_matrix(obj)
            if hit is None:
                continue
            key_path, A = hit
            return float(A.mean()), f"{jf.name}:{key_path}"
        except Exception as e:
            last_err = e

    raise ValueError(f"No PAE matrix found in {af3_dir} (last_err={last_err})")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in_csv", required=True, help="Input CSV (e.g., metrics_my_run_COMBINED.csv)")
    ap.add_argument("--run_dir", required=True, help="Run dir (e.g., runs/run_test_rfd3_nompmn)")
    ap.add_argument("--out_csv", required=True, help="Output CSV with mean_pae column added")
    ap.add_argument("--progress_every", type=int, default=500)
    args = ap.parse_args()

    in_csv = Path(args.in_csv)
    run_dir = Path(args.run_dir)
    out_csv = Path(args.out_csv)
    out_csv.parent.mkdir(parents=True, exist_ok=True)

    af3_root = run_dir / "03_af3" / "outputs"
    if not af3_root.exists():
        raise SystemExit(f"AF3 outputs root not found: {af3_root}")

    df = pd.read_csv(in_csv)
    if "tok" not in df.columns and "assay_scaffolded_epitope_id" in df.columns:
        df = df.rename(columns={"assay_scaffolded_epitope_id": "tok"})
    if "tok" not in df.columns:
        raise SystemExit("Input CSV must have column 'tok' (or 'assay_scaffolded_epitope_id').")
    if "pred" not in df.columns:
        raise SystemExit("Input CSV must have column 'pred'.")

    # Only need compute once per unique (tok,pred)
    pairs = df[["tok", "pred"]].drop_duplicates().copy()
    pairs["tok"] = pairs["tok"].astype(str).str.lower()
    pairs["pred"] = pairs["pred"].astype(int)

    # Build AF3 dir by globbing: <tok>...0_model_<pred> (directory)
    # This matches your AF3 layout.
    def find_af3_dir(tok: str, pred: int) -> Optional[Path]:
        pat = f"{tok}*0_model_{pred}"
        hit = next(af3_root.glob(pat), None)
        return hit if (hit is not None and hit.is_dir()) else None

    mean_vals = []
    n_missing_dir = 0
    n_no_pae = 0

    for i, r in pairs.iterrows():
        tok = r["tok"]
        pred = int(r["pred"])
        af3_dir = find_af3_dir(tok, pred)

        if af3_dir is None:
            mean_vals.append((tok, pred, np.nan, "missing_af3_dir"))
            n_missing_dir += 1
        else:
            try:
                mp, src = mean_pae_from_af3_dir(af3_dir)
                mean_vals.append((tok, pred, mp, src))
            except Exception as e:
                mean_vals.append((tok, pred, np.nan, f"no_pae:{type(e).__name__}"))
                n_no_pae += 1

        if args.progress_every and (len(mean_vals) % args.progress_every == 0):
            print(f"[progress] computed {len(mean_vals)}/{len(pairs)} | missing_dir={n_missing_dir} no_pae={n_no_pae}")

    pae_df = pd.DataFrame(mean_vals, columns=["tok","pred","mean_pae","mean_pae_src"])

    # Merge back to full df
    df2 = df.copy()
    df2["tok"] = df2["tok"].astype(str).str.lower()
    df2["pred"] = df2["pred"].astype(int)

    out = df2.merge(pae_df, on=["tok","pred"], how="left")
    out.to_csv(out_csv, index=False)

    print("DONE")
    print("input rows:", len(df))
    print("unique (tok,pred):", len(pairs))
    print("missing_af3_dir:", n_missing_dir)
    print("no_pae:", n_no_pae)
    print("wrote:", out_csv)


if __name__ == "__main__":
    main()
