#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
import MDAnalysis as mda
from MDAnalysis.analysis.align import alignto


def _as_list(x) -> List[int]:
    """dp2 parquet list columns sometimes come in as numpy arrays, lists, or objects."""
    if x is None or (isinstance(x, float) and np.isnan(x)):
        return []
    if isinstance(x, (list, tuple)):
        return [int(i) for i in x]
    if hasattr(x, "tolist"):
        return [int(i) for i in x.tolist()]
    # last resort: try iter
    try:
        return [int(i) for i in list(x)]
    except Exception:
        return []


def load_dp2_maps(dp2_parquet: Path):
    dp2 = pd.read_parquet(dp2_parquet)

    required = [
        "assay_scaffolded_epitope_id",
        "id",
        "assay_scaffolded_epitope_chunk_resindices",
        "epitope_chunk_resindices",
    ]
    miss = [c for c in required if c not in dp2.columns]
    if miss:
        raise SystemExit(f"dp2 parquet missing required cols: {miss}")

    dp2["assay_scaffolded_epitope_id"] = dp2["assay_scaffolded_epitope_id"].astype(str).str.lower()
    dp2["id"] = dp2["id"].astype(str).str.strip()

    tok_to_af3 = {}
    g = dp2.groupby("assay_scaffolded_epitope_id", sort=False)["assay_scaffolded_epitope_chunk_resindices"]
    for tok, series in g:
        v = None
        for x in series:
            lst = _as_list(x)
            if lst:
                v = lst
                break
        if v is None:
            v = _as_list(series.iloc[0])
        tok_to_af3[tok] = v

    id_to_true = {}
    g2 = dp2.groupby("id", sort=False)["epitope_chunk_resindices"]
    for rid, series in g2:
        v = None
        for x in series:
            lst = _as_list(x)
            if lst:
                v = lst
                break
        if v is None:
            v = _as_list(series.iloc[0])
        id_to_true[rid] = v

    idtok_to_truth = {}
    if "af3_clash_resindices" in dp2.columns:
        for rid, tok, arr in dp2[["id", "assay_scaffolded_epitope_id", "af3_clash_resindices"]].itertuples(index=False):
            idtok_to_truth[(rid, tok)] = _as_list(arr)

    return tok_to_af3, id_to_true, idtok_to_truth
def find_af3_model_file(af3_root: Path, tok: str, pred: Optional[int], layout: str) -> Optional[Path]:
    """
    layout:
      - "mine":   af3_root/<tok>*0_model_<pred>/<dir>_model.cif
      - "lawson": af3_root/<tok>/<tok>_model.cif.gz  (or any *_model.cif.gz)
                 and fallback into seed-* subdirs
    """
    tok = tok.lower()

    if layout == "mine":
        if pred is None:
            raise ValueError("layout=mine requires pred")
        d = next(af3_root.glob(f"{tok}*0_model_{int(pred)}"), None)
        if d is None:
            return None
        if d.is_dir():
            p = d / f"{d.name}_model.cif"
            if p.exists():
                return p
            # sometimes gz
            pgz = d / f"{d.name}_model.cif.gz"
            if pgz.exists():
                return pgz
        return None

    if layout == "lawson":
        d = af3_root / tok
        if d.is_dir():
            # common case: <tok>/<tok>_model.cif.gz
            p1 = d / f"{tok}_model.cif.gz"
            if p1.exists():
                return p1
            p2 = d / f"{tok}_model.cif"
            if p2.exists():
                return p2
            # otherwise find any *_model.cif(.gz) at top-level
            hit = next(d.glob("*_model.cif.gz"), None)
            if hit is not None:
                return hit
            hit = next(d.glob("*_model.cif"), None)
            if hit is not None:
                return hit

            # fallback: seed-* dirs
            hit = next(d.glob("seed-*/**/*_model.cif.gz"), None)
            if hit is not None:
                return hit
            hit = next(d.glob("seed-*/**/*_model.cif"), None)
            if hit is not None:
                return hit
        return None

    raise ValueError(f"Unknown layout: {layout}")


def ensure_segids(u: mda.Universe) -> None:
    """
    Lawson assumes segid A for the scaffold model and segid B/C for antibody.
    AF3 CIFs sometimes come in with blank segids. If AF3 has only one segment, set it to 'A'.
    """
    try:
        segids = [s.segid for s in u.segments]
        if len(u.segments) == 1 and (segids[0] is None or str(segids[0]).strip() == ""):
            u.segments[0].segid = "A"
    except Exception:
        pass


def compute_af3_clash_resindices(
    af3_model_file: Path,
    true_complex_pdb: Path,
    af3_epitope_ris: List[int],
    true_epitope_ris: List[int],
    cutoff_A: float = 4.0,
) -> List[int]:
    """
    Replicates Lawson AF3 clash logic:

    1) load AF3 predicted scaffold-only structure as segid A
    2) load true complex, take antibody atoms segid B or segid C
    3) align AF3 epitope backbone onto true epitope backbone
    4) merge (AF3 scaffold) + (true antibody)
    5) find antibody residues (B/C, heavy atoms) within cutoff of any *non-epitope* scaffold residues
    6) return those antibody residue resindices (in merged universe, like Lawson)
    """
    af3_u = mda.Universe(str(af3_model_file))
    ensure_segids(af3_u)

    true_u = mda.Universe(str(true_complex_pdb))
    ab_atoms = true_u.select_atoms("segid B or segid C")
    if len(ab_atoms) == 0:
        raise RuntimeError(f"No antibody atoms found with segid B/C in {true_complex_pdb}")

    # align AF3 epitope chunk to true epitope chunk
    alignto(
        af3_u.residues[af3_epitope_ris].atoms.select_atoms("backbone"),
        true_u.residues[true_epitope_ris].atoms.select_atoms("backbone"),
    )

    combined = mda.Merge(af3_u.atoms, ab_atoms)

    # unintended scaffold residues = segid A residues not in epitope chunk
    tmp = np.full(len(af3_u.residues), False)
    tmp[af3_epitope_ris] = True

    unintended_res = combined.select_atoms("segid A").residues[~tmp]

    # antibody residues (B/C) within cutoff of unintended scaffold residues, heavy atoms only
    sel = combined.atoms.select_atoms(
        f"((segid B or segid C) and (not name H*)) and around {cutoff_A:g} group unintended",
        unintended=unintended_res.atoms,
    )
    return list(sel.residues.resindices)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in_csv", required=True, help="Your combined metrics CSV (must include columns: id, tok, pred if layout=mine)")
    ap.add_argument("--dp2_parquet", required=True, help="datasets/dp2.parquet (contains epitope indices and Lawson af3_clash_resindices)")
    ap.add_argument("--true_dir", required=True, help="Directory of true complexes (id.pdb)")
    ap.add_argument("--af3_root", required=True, help="AF3 predictions root")
    ap.add_argument("--layout", choices=["mine", "lawson"], required=True)
    ap.add_argument("--out_csv", required=True)
    ap.add_argument("--cutoff_A", type=float, default=4.0)
    ap.add_argument("--limit", type=int, default=0, help="If >0, only process first N rows (debug)")
    ap.add_argument("--progress_every", type=int, default=200)
    args = ap.parse_args()

    in_csv = Path(args.in_csv)
    dp2_parquet = Path(args.dp2_parquet)
    true_dir = Path(args.true_dir)
    af3_root = Path(args.af3_root)

    df = pd.read_csv(in_csv)
    # normalize columns
    if "tok" not in df.columns and "assay_scaffolded_epitope_id" in df.columns:
        df = df.rename(columns={"assay_scaffolded_epitope_id": "tok"})
    if "tok" not in df.columns or "id" not in df.columns:
        raise SystemExit("Input CSV must have columns: id and tok (or assay_scaffolded_epitope_id).")
    if args.layout == "mine" and "pred" not in df.columns:
        raise SystemExit("layout=mine requires column pred in input CSV.")

    df["tok"] = df["tok"].astype(str).str.lower()
    df["id"] = df["id"].astype(str).str.strip()
    if "pred" in df.columns:
        df["pred"] = df["pred"].astype(int)

    tok_to_af3, id_to_true, idtok_to_truth = load_dp2_maps(dp2_parquet)

    rows = []
    n_missing_af3 = 0
    n_missing_true = 0
    n_missing_maps = 0
    n_fail = 0

    it = df.itertuples(index=False)
    if args.limit and args.limit > 0:
        it = list(it)[: args.limit]

    for i, r in enumerate(it, start=1):
        rid = getattr(r, "id")
        tok = getattr(r, "tok")
        pred = getattr(r, "pred", None) if args.layout == "mine" else None

        af3_ep = tok_to_af3.get(tok, None)
        true_ep = id_to_true.get(rid, None)
        if af3_ep is None or true_ep is None:
            n_missing_maps += 1
            rows.append({"id": rid, "tok": tok, "pred": pred, "af3_clash_resindices": np.nan, "af3_n_clash_res": np.nan,
                         "truth_af3_n_clash_res": np.nan, "truth_match": np.nan, "status": "missing_dp2_maps"})
            continue

        true_pdb = true_dir / f"{rid}.pdb"
        if not true_pdb.exists():
            n_missing_true += 1
            rows.append({"id": rid, "tok": tok, "pred": pred, "af3_clash_resindices": np.nan, "af3_n_clash_res": np.nan,
                         "truth_af3_n_clash_res": np.nan, "truth_match": np.nan, "status": "missing_true_pdb"})
            continue

        af3_model = find_af3_model_file(af3_root, tok, pred, args.layout)
        if af3_model is None or not af3_model.exists():
            n_missing_af3 += 1
            rows.append({"id": rid, "tok": tok, "pred": pred, "af3_clash_resindices": np.nan, "af3_n_clash_res": np.nan,
                         "truth_af3_n_clash_res": np.nan, "truth_match": np.nan, "status": "missing_af3_model"})
            continue

        try:
            clashes = compute_af3_clash_resindices(
                af3_model, true_pdb, af3_ep, true_ep, cutoff_A=args.cutoff_A
            )
            truth = idtok_to_truth.get((rid, tok), [])
            # Lawson stored list could be empty/None; compare lengths as the main sanity check
            truth_len = len(truth) if truth is not None else np.nan
            ok = (len(clashes) == truth_len) if isinstance(truth_len, int) else np.nan

            rows.append({
                "id": rid,
                "tok": tok,
                "pred": pred,
                "af3_model_file": str(af3_model),
                "af3_n_clash_res": len(clashes),
                "af3_clash_resindices": clashes,
                "truth_af3_n_clash_res": truth_len,
                "truth_match": ok,
                "status": "ok",
            })
        except Exception as e:
            import traceback
            print("\n=== FIRST FAIL DEBUG ===")
            print("id:", rid, "tok:", tok, "pred:", pred)
            print("true_pdb:", true_pdb)
            print("af3_model:", af3_model)
            traceback.print_exc()
            raise

        if args.progress_every and (i % args.progress_every == 0):
            print(f"[progress] {i} rows | missing_af3={n_missing_af3} missing_true={n_missing_true} missing_maps={n_missing_maps} fail={n_fail}")

    out = pd.DataFrame(rows)
    out.to_csv(args.out_csv, index=False)

    print("DONE")
    print("wrote:", args.out_csv)
    print("missing_af3_model:", n_missing_af3)
    print("missing_true_pdb:", n_missing_true)
    print("missing_dp2_maps:", n_missing_maps)
    print("fail:", n_fail)
    if "truth_match" in out.columns:
        try:
            m = out["truth_match"].dropna()
            if len(m) > 0:
                print("truth_match_rate:", float(m.mean()))
        except Exception:
            pass


if __name__ == "__main__":
    main()
