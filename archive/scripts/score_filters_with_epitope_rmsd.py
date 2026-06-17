#!/usr/bin/env python3
from __future__ import annotations

import argparse, gzip, json, math, re
from pathlib import Path
from typing import Optional, Tuple, Dict, Any

import numpy as np
import pandas as pd
import gemmi


def kabsch_rmsd(P: np.ndarray, Q: np.ndarray) -> float:
    if P.shape != Q.shape or P.shape[0] < 3:
        return float("nan")
    Pc = P - P.mean(axis=0)
    Qc = Q - Q.mean(axis=0)
    C = Pc.T @ Qc
    V, _, Wt = np.linalg.svd(C)
    d = np.sign(np.linalg.det(V @ Wt))
    D = np.diag([1.0, 1.0, d])
    U = V @ D @ Wt
    diff = (Pc @ U) - Qc
    return float(np.sqrt((diff * diff).sum() / P.shape[0]))


def parse_index_list(x):
    if x is None or (isinstance(x, float) and math.isnan(x)):
        return []
    if isinstance(x, (list, tuple, np.ndarray)):
        return [int(i) for i in x]
    s = str(x).strip().replace("[", "").replace("]", "").replace(",", " ")
    out = []
    for tok in s.split():
        try:
            out.append(int(tok))
        except ValueError:
            pass
    return out


def read_structure(path: Path) -> gemmi.Structure:
    if str(path).endswith(".gz"):
        with gzip.open(path, "rt") as f:
            doc = gemmi.cif.read_string(f.read())
        return gemmi.make_structure_from_block(doc.sole_block())
    if path.suffix.lower() == ".cif":
        doc = gemmi.cif.read(str(path))
        return gemmi.make_structure_from_block(doc.sole_block())
    return gemmi.read_structure(str(path))


def ca_coords(struct: gemmi.Structure, chain_id: str = "A") -> np.ndarray:
    model = struct[0]
    chain_names = [c.name for c in model]
    chain = model[chain_id] if chain_id in chain_names else model[0]
    coords = []
    for res in chain:
        ca = res.find_atom("CA", altloc="*")
        if ca:
            p = ca.pos
            coords.append([p.x, p.y, p.z])
    return np.array(coords, dtype=float)


def ca_coords_resids(struct: gemmi.Structure, resids_1based: set[int], chain_id: str = "A") -> np.ndarray:
    model = struct[0]
    chain_names = [c.name for c in model]
    chain = model[chain_id] if chain_id in chain_names else model[0]
    coords = []
    for res in chain:
        if res.seqid.num in resids_1based:
            ca = res.find_atom("CA", altloc="*")
            if ca:
                p = ca.pos
                coords.append([p.x, p.y, p.z])
    return np.array(coords, dtype=float)

def ca_coord_map(struct: gemmi.Structure, chain_id: str = "A") -> dict[int, np.ndarray]:
    model = struct[0]
    chain_names = [c.name for c in model]
    chain = model[chain_id] if chain_id in chain_names else model[0]
    out: dict[int, np.ndarray] = {}
    for res in chain:
        ca = res.find_atom("CA", altloc="*")
        if ca:
            p = ca.pos
            out[res.seqid.num] = np.array([p.x, p.y, p.z], dtype=float)
    return out


def coords_for_resids(struct: gemmi.Structure, resids_1based: set[int], chain_id: str = "A") -> np.ndarray:
    m = ca_coord_map(struct, chain_id)
    common = sorted([rid for rid in resids_1based if rid in m])
    return np.array([m[rid] for rid in common], dtype=float)


def _json_load(p: Optional[Path]) -> Optional[dict]:
    if p is None:
        return None
    try:
        return json.loads(p.read_text())
    except Exception:
        return None


def pae_mean_from_conf(conf_json: Optional[Path]) -> Optional[float]:
    d = _json_load(conf_json)
    if not isinstance(d, dict):
        return None
    for k in ("pae", "predicted_aligned_error", "predicted_aligned_error_matrix"):
        if k in d:
            try:
                a = np.array(d[k], dtype=float)
                return float(np.nanmean(a))
            except Exception:
                return None
    return None


def summary_scalars(summary_json: Optional[Path]) -> dict[str, Optional[float]]:
    d = _json_load(summary_json)
    out = {"ptm": None, "ranking_score": None, "chain_pair_pae_min": None, "fraction_disordered": None, "has_clash": None}
    if not isinstance(d, dict):
        return out

    def to_float(x):
        try:
            if x is None:
                return None
            return float(x)
        except Exception:
            return None

    out["ptm"] = to_float(d.get("ptm"))
    out["ranking_score"] = to_float(d.get("ranking_score"))
    out["fraction_disordered"] = to_float(d.get("fraction_disordered"))
    out["has_clash"] = to_float(d.get("has_clash"))

    # chain_pair_pae_min is often [[0.76]] etc.
    cppm = d.get("chain_pair_pae_min")
    try:
        if isinstance(cppm, list) and cppm and isinstance(cppm[0], list) and cppm[0]:
            out["chain_pair_pae_min"] = float(cppm[0][0])
        elif isinstance(cppm, (int, float)):
            out["chain_pair_pae_min"] = float(cppm)
    except Exception:
        pass

    return out


def find_af3_files_in_dir(d: Path) -> Tuple[Optional[Path], Optional[Path], Optional[Path]]:
    """
    Returns: (af3_cif, conf_json_with_matrix, summary_conf_json)
    """
    af3_cif = next(iter(d.glob("*_model.cif")), None)

    # Prefer top-level *_confidences.json (often contains 'pae' matrix)
    conf = next(iter(d.glob("*_confidences.json")), None)
    # Fallback to nested seed-* / confidences.json
    if conf is None:
        conf = next(iter(d.rglob("confidences.json")), None)

    summary = next(iter(d.glob("*_summary_confidences.json")), None)
    if summary is None:
        summary = next(iter(d.rglob("summary_confidences.json")), None)

    # last resort for cif
    if af3_cif is None:
        af3_cif = next(iter(d.rglob("model.cif")), None)

    return af3_cif, conf, summary


def index_af3_outputs(af3_root: Path) -> Dict[tuple[str, int], Path]:
    idx: Dict[tuple[str, int], Path] = {}
    if not af3_root.exists():
        return idx
    tok_re = re.compile(r"^([0-9a-fA-F]{32})")
    pred_re = re.compile(r"_0_model_([0-7])\b", re.IGNORECASE)
    for d in af3_root.iterdir():
        if not d.is_dir():
            continue
        m_tok = tok_re.match(d.name)
        m_pred = pred_re.search(d.name)
        if not (m_tok and m_pred):
            continue
        tok = m_tok.group(1).lower()
        pred = int(m_pred.group(1))
        idx.setdefault((tok, pred), d)
    return idx


def index_rfd3_outputs(rfd_root: Path) -> Dict[tuple[str, int], Path]:
    idx: Dict[tuple[str, int], Path] = {}
    if not rfd_root.exists():
        return idx
    tok_re = re.compile(r"^([0-9a-fA-F]{32})")
    pred_re = re.compile(r"_0_model_([0-7])\.cif\.gz$", re.IGNORECASE)
    for p in rfd_root.iterdir():
        if not p.is_file():
            continue
        m_tok = tok_re.match(p.name)
        m_pred = pred_re.search(p.name)
        if not (m_tok and m_pred):
            continue
        tok = m_tok.group(1).lower()
        pred = int(m_pred.group(1))
        idx.setdefault((tok, pred), p)
    return idx


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--run_dir", required=True)
    ap.add_argument("--contigs_parquet", default="datasets/dp2.parquet")
    ap.add_argument("--out_csv", default=None)
    ap.add_argument("--pred_max", type=int, default=7)
    ap.add_argument("--limit", type=int, default=0)
    args = ap.parse_args()

    run_dir = Path(args.run_dir).resolve()
    contigs_pq = Path(args.contigs_parquet).resolve()
    out_csv = Path(args.out_csv).resolve() if args.out_csv else (run_dir / "04_filter" / "metrics_by_pred.csv")
    out_csv.parent.mkdir(parents=True, exist_ok=True)

    df = pd.read_parquet(contigs_pq)
    if "assay_scaffolded_epitope_id" not in df.columns:
        raise KeyError(f"'assay_scaffolded_epitope_id' not found. Columns: {list(df.columns)}")

    df["assay_scaffolded_epitope_id"] = df["assay_scaffolded_epitope_id"].astype(str).str.lower()
    if args.limit and args.limit > 0:
        df = df.head(args.limit)

    af3_root = run_dir / "03_af3" / "outputs"
    rfd_root = run_dir / "02_rfd3" / "outputs"

    af3_index = index_af3_outputs(af3_root)
    rfd3_index = index_rfd3_outputs(rfd_root)

    rows = []
    missing = 0

    for _, r in df.iterrows():
        tok = r["assay_scaffolded_epitope_id"]
        ep0 = parse_index_list(r.get("assay_scaffolded_epitope_resindices", r.get("epitope_resindices")))
        ep_resids = set(i + 1 for i in ep0)

        for pred in range(args.pred_max + 1):
            af3_dir = af3_index.get((tok, pred))
            af3_cif, conf_json, summary_json = (None, None, None)
            if af3_dir:
                af3_cif, conf_json, summary_json = find_af3_files_in_dir(af3_dir)

            rfd3_file = rfd3_index.get((tok, pred))

            summ = summary_scalars(summary_json)
            row = {
                "assay_scaffolded_epitope_id": tok,
                "pred": pred,
                "id": r.get("id"),
                "contig_id": r.get("contig_id"),
                "rfd3_path": str(rfd3_file) if rfd3_file else None,
                "af3_dir": str(af3_dir) if af3_dir else None,
                "af3_path": str(af3_cif) if af3_cif else None,
                "af3_conf_path": str(conf_json) if conf_json else None,
                "af3_summary_path": str(summary_json) if summary_json else None,
                "pae_mean_all": pae_mean_from_conf(conf_json) if conf_json else None,
                "ptm": summ["ptm"],
                "ranking_score": summ["ranking_score"],
                "chain_pair_pae_min": summ["chain_pair_pae_min"],
                "fraction_disordered": summ["fraction_disordered"],
                "has_clash": summ["has_clash"],
                "rmsd_ca_all": None,
                "rmsd_ca_epitope": None,
                "n_ca_epitope_common": None,
            }

            if (rfd3_file is None) or (af3_cif is None):
                missing += 1
                rows.append(row)
                continue

            try:
                st_rfd = read_structure(Path(rfd3_file))
                st_af3 = read_structure(Path(af3_cif))

                P = ca_coords(st_rfd, "A")
                Q = ca_coords(st_af3, "A")
                row["rmsd_ca_all"] = kabsch_rmsd(P, Q)

                Pe = coords_for_resids(st_rfd, ep_resids, "A")
                Qe = coords_for_resids(st_af3, ep_resids, "A")
                if Pe.shape == Qe.shape and Pe.shape[0] >= 3:
                    row["rmsd_ca_epitope"] = kabsch_rmsd(Pe, Qe)
                else:
                    row["rmsd_ca_epitope"] = float("nan")
                row["n_ca_epitope_common"] = int(Pe.shape[0])
            except Exception:
                pass

            rows.append(row)

    out = pd.DataFrame(rows)
    out.to_csv(out_csv, index=False)
    print(f"Wrote {len(out)} rows -> {out_csv}")
    print(f"Missing pairs (no rfd3 or no af3 cif): {missing}")
    print(f"AF3 index size (tok,pred)->dir: {len(af3_index)}")
    print(f"RFD3 index size (tok,pred)->cif.gz: {len(rfd3_index)}")


if __name__ == "__main__":
    main()
