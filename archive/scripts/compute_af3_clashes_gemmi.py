#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path
from typing import List, Optional, Tuple, Any, Dict
import gzip

import numpy as np
import pandas as pd
import gemmi

try:
    from scipy.spatial import cKDTree as KDTree
except Exception:
    KDTree = None


def _as_list(x) -> List[int]:
    if x is None or (isinstance(x, float) and np.isnan(x)):
        return []
    if isinstance(x, (list, tuple)):
        return [int(i) for i in x]
    if hasattr(x, "tolist"):
        return [int(i) for i in x.tolist()]
    try:
        return [int(i) for i in list(x)]
    except Exception:
        return []


def read_structure_gemmi(p: Path) -> gemmi.Structure:
    s = str(p)
    if s.endswith(".gz"):
        with gzip.open(p, "rt") as f:
            doc = gemmi.cif.read_string(f.read())
        st = gemmi.make_structure_from_block(doc.sole_block())
        st.setup_entities()
        return st
    if p.suffix.lower() == ".cif":
        doc = gemmi.cif.read(str(p))
        st = gemmi.make_structure_from_block(doc.sole_block())
        st.setup_entities()
        return st
    st = gemmi.read_structure(str(p))
    st.setup_entities()
    return st


def get_chain(model: gemmi.Model, name: str) -> gemmi.Chain:
    for ch in model:
        if ch.name == name:
            return ch
    raise KeyError(f"Chain {name} not found")


def is_h(atom: gemmi.Atom) -> bool:
    try:
        if atom.element.name == "H":
            return True
    except Exception:
        pass
    nm = (atom.name or "").strip().upper()
    return nm.startswith("H")


def atom_pos(res: gemmi.Residue, atom_name: str) -> Optional[np.ndarray]:
    a = res.find_atom(atom_name, altloc="*")
    if not a:
        return None
    p = a.pos
    return np.array([p.x, p.y, p.z], dtype=float)


def kabsch(P: np.ndarray, Q: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """
    Return (R, t) such that (R @ x + t) maps P onto Q.
    P,Q: (N,3)
    """
    if P.shape != Q.shape or P.ndim != 2 or P.shape[1] != 3:
        raise ValueError("kabsch requires P and Q as (N,3) arrays of same shape")

    cP = P.mean(axis=0)
    cQ = Q.mean(axis=0)
    X = P - cP
    Y = Q - cQ

    C = X.T @ Y
    V, S, Wt = np.linalg.svd(C)
    d = np.sign(np.linalg.det(V @ Wt))
    D = np.diag([1.0, 1.0, d])
    R = V @ D @ Wt
    t = cQ - (R @ cP)
    return R, t


def apply_rt(xyz: np.ndarray, R: np.ndarray, t: np.ndarray) -> np.ndarray:
    return (xyz @ R.T) + t


def compute_clashes_gemmi(
    af3_cif: Path,
    true_pdb: Path,
    af3_epitope_ris: List[int],
    true_epitope_ris: List[int],
    cutoff_A: float = 4.0,
) -> List[int]:
    """
    Lawson-style clash metric using gemmi only.

    True complex convention:
      chain A = antigen
      chain B = antibody heavy
      chain C = antibody light

    Returns merged-style residue indices:
      merged_resindex = len(scaffold_residues) + antibody_res_counter (B then C)
    """
    st_af3 = read_structure_gemmi(af3_cif)
    st_true = read_structure_gemmi(true_pdb)

    m_af3 = st_af3[0]
    m_true = st_true[0]

    ch_scaf = get_chain(m_af3, "A")
    ch_ag   = get_chain(m_true, "A")
    ch_H    = get_chain(m_true, "B")
    ch_L    = get_chain(m_true, "C")

    scaf_res = list(ch_scaf)
    ag_res   = list(ch_ag)
    H_res    = list(ch_H)
    L_res    = list(ch_L)

    if max(af3_epitope_ris, default=-1) >= len(scaf_res):
        raise IndexError(f"af3_epitope_ris out of range: max={max(af3_epitope_ris)} len(scaf)={len(scaf_res)}")
    if max(true_epitope_ris, default=-1) >= len(ag_res):
        raise IndexError(f"true_epitope_ris out of range: max={max(true_epitope_ris)} len(ag)={len(ag_res)}")

    # alignment sets from backbone atoms
    names = ("N", "CA", "C", "O")
    P = []
    Q = []
    for i_scaf, i_true in zip(af3_epitope_ris, true_epitope_ris):
        rS = scaf_res[i_scaf]
        rT = ag_res[i_true]
        for nm in names:
            pS = atom_pos(rS, nm)
            pT = atom_pos(rT, nm)
            if pS is None or pT is None:
                continue
            P.append(pS)
            Q.append(pT)

    if len(P) < 3:
        raise RuntimeError(f"Not enough backbone atoms for alignment: n={len(P)}")

    P = np.asarray(P, dtype=float)
    Q = np.asarray(Q, dtype=float)
    R, t = kabsch(P, Q)

    # unintended scaffold atom coords transformed into true frame
    epi_mask = np.zeros(len(scaf_res), dtype=bool)
    epi_mask[af3_epitope_ris] = True

    unintended = []
    for i, r in enumerate(scaf_res):
        if epi_mask[i]:
            continue
        for a in r:
            if is_h(a):
                continue
            p = a.pos
            unintended.append([p.x, p.y, p.z])

    if not unintended:
        return []

    unintended_xyz = apply_rt(np.asarray(unintended, dtype=float), R, t)
    cutoff = float(cutoff_A)

    # KDTree / fallback
    if KDTree is not None:
        tree = KDTree(unintended_xyz)

        def residue_hits(res_atoms_xyz: np.ndarray) -> bool:
            d, _ = tree.query(res_atoms_xyz, k=1, workers=-1)
            return bool(np.any(d <= cutoff))

    else:
        U = unintended_xyz
        cutoff2 = cutoff * cutoff

        def residue_hits(res_atoms_xyz: np.ndarray) -> bool:
            # brute force min distance^2
            for i0 in range(0, len(res_atoms_xyz), 256):
                Pchunk = res_atoms_xyz[i0:i0+256]
                diff = Pchunk[:, None, :] - U[None, :, :]
                d2 = np.sum(diff * diff, axis=2)
                if np.any(np.min(d2, axis=1) <= cutoff2):
                    return True
            return False

    clashing = set()
    merged_offset = len(scaf_res)

    ab_counter = 0
    for chain_res in (H_res, L_res):
        for r in chain_res:
            heavy = []
            for a in r:
                if is_h(a):
                    continue
                p = a.pos
                heavy.append([p.x, p.y, p.z])
            if heavy:
                heavy_xyz = np.asarray(heavy, dtype=float)
                if residue_hits(heavy_xyz):
                    clashing.add(merged_offset + ab_counter)
            ab_counter += 1

    return sorted(clashing)


def find_af3_cif_mine(af3_root: Path, tok: str, pred: int) -> Optional[Path]:
    # your layout: <tok>*0_model_<pred>/<dir>_model.cif (or .cif.gz)
    pat = f"{tok}*0_model_{int(pred)}"
    d = next(af3_root.glob(pat), None)
    if d is None:
        return None
    if d.is_dir():
        p = d / f"{d.name}_model.cif"
        if p.exists():
            return p
        pgz = d / f"{d.name}_model.cif.gz"
        if pgz.exists():
            return pgz
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in_csv", required=True)
    ap.add_argument("--dp2_parquet", required=True)
    ap.add_argument("--true_dir", required=True)
    ap.add_argument("--af3_root", required=True)
    ap.add_argument("--cutoff_A", type=float, default=4.0)
    ap.add_argument("--out_csv", required=True)
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--progress_every", type=int, default=200)
    args = ap.parse_args()

    df = pd.read_csv(args.in_csv)
    if "tok" not in df.columns and "assay_scaffolded_epitope_id" in df.columns:
        df = df.rename(columns={"assay_scaffolded_epitope_id": "tok"})
    if "tok" not in df.columns or "id" not in df.columns or "pred" not in df.columns:
        raise SystemExit("in_csv must contain columns: id, tok, pred")

    # keep exact id case (files are case sensitive)
    df["id"] = df["id"].astype(str).str.strip()
    df["tok"] = df["tok"].astype(str).str.strip().str.lower()
    df["pred"] = df["pred"].astype(int)

    dp2 = pd.read_parquet(args.dp2_parquet)
    dp2["id"] = dp2["id"].astype(str).str.strip()
    dp2["assay_scaffolded_epitope_id"] = dp2["assay_scaffolded_epitope_id"].astype(str).str.strip().str.lower()

    # tok -> af3 epitope chunk indices
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

    # id -> true epitope chunk indices
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

    true_dir = Path(args.true_dir)
    af3_root = Path(args.af3_root)

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
        pred = int(getattr(r, "pred"))

        af3_ep = tok_to_af3.get(tok)
        true_ep = id_to_true.get(rid)
        if af3_ep is None or true_ep is None:
            n_missing_maps += 1
            rows.append({"id": rid, "tok": tok, "pred": pred, "af3_n_clash_res": np.nan, "af3_clash_resindices": np.nan, "status": "missing_dp2_maps"})
            continue

        true_pdb = true_dir / f"{rid}.pdb"
        if not true_pdb.exists():
            n_missing_true += 1
            rows.append({"id": rid, "tok": tok, "pred": pred, "af3_n_clash_res": np.nan, "af3_clash_resindices": np.nan, "status": "missing_true_pdb"})
            continue

        af3_cif = find_af3_cif_mine(af3_root, tok, pred)
        if af3_cif is None or not af3_cif.exists():
            n_missing_af3 += 1
            rows.append({"id": rid, "tok": tok, "pred": pred, "af3_n_clash_res": np.nan, "af3_clash_resindices": np.nan, "status": "missing_af3_model"})
            continue

        try:
            clashes = compute_clashes_gemmi(af3_cif, true_pdb, af3_ep, true_ep, cutoff_A=args.cutoff_A)
            rows.append({
                "id": rid,
                "tok": tok,
                "pred": pred,
                "af3_model_file": str(af3_cif),
                "af3_n_clash_res": len(clashes),
                "af3_clash_resindices": clashes,
                "status": "ok",
            })
        except Exception as e:
            n_fail += 1
            rows.append({"id": rid, "tok": tok, "pred": pred, "af3_n_clash_res": np.nan, "af3_clash_resindices": np.nan, "status": f"fail:{type(e).__name__}"})

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


if __name__ == "__main__":
    main()
