#!/usr/bin/env python3
from __future__ import annotations

import argparse
import gzip
import json
import math
import re
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
import gemmi
import MDAnalysis as mda
from MDAnalysis.lib.distances import distance_array


# ------------------------ helpers ------------------------

def parse_index_list(x: Any) -> List[int]:
    if x is None:
        return []
    if isinstance(x, float) and math.isnan(x):
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
    P_rot = Pc @ U
    diff = P_rot - Qc
    return float(np.sqrt((diff * diff).sum() / P.shape[0]))


def kabsch_fit(P: np.ndarray, Q: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """Return (R,t) such that P_fit = (P @ R) + t best aligns to Q."""
    Pc = P - P.mean(axis=0)
    Qc = Q - Q.mean(axis=0)
    C = Pc.T @ Qc
    V, _, Wt = np.linalg.svd(C)
    d = np.sign(np.linalg.det(V @ Wt))
    D = np.diag([1.0, 1.0, d])
    R = V @ D @ Wt
    t = Q.mean(axis=0) - (P.mean(axis=0) @ R)
    return R, t


def read_structure(path: Path) -> gemmi.Structure:
    s = str(path)
    if s.endswith(".gz"):
        with gzip.open(path, "rt") as f:
            doc = gemmi.cif.read_string(f.read())
        return gemmi.make_structure_from_block(doc.sole_block())
    if path.suffix.lower() == ".cif":
        doc = gemmi.cif.read(str(path))
        return gemmi.make_structure_from_block(doc.sole_block())
    return gemmi.read_structure(str(path))


def get_chain(struct: gemmi.Structure, chain_id: str = "A") -> gemmi.Chain:
    model = struct[0]
    for ch in model:
        if ch.name == chain_id:
            return ch
    return model[0]


def chain_residues(chain: gemmi.Chain) -> List[gemmi.Residue]:
    return [r for r in chain]


def res_atom_pos(res: gemmi.Residue, name: str) -> Optional[np.ndarray]:
    a = res.find_atom(name, altloc="*")
    if not a:
        return None
    p = a.pos
    return np.array([p.x, p.y, p.z], dtype=float)


def ca_coords(struct: gemmi.Structure, chain_id: str = "A") -> np.ndarray:
    ch = get_chain(struct, chain_id)
    coords = []
    for res in chain_residues(ch):
        p = res_atom_pos(res, "CA")
        if p is not None:
            coords.append(p)
    return np.array(coords, dtype=float)


def ca_coords_residx(struct: gemmi.Structure, residx_0based: List[int], chain_id: str = "A") -> np.ndarray:
    ch = get_chain(struct, chain_id)
    res = chain_residues(ch)
    coords = []
    for i in residx_0based:
        if 0 <= i < len(res):
            p = res_atom_pos(res[i], "CA")
            if p is not None:
                coords.append(p)
    return np.array(coords, dtype=float)


BB_ATOMS = ("N", "CA", "C", "O")


def bb_coords(struct: gemmi.Structure, chain_id: str = "A") -> np.ndarray:
    ch = get_chain(struct, chain_id)
    coords = []
    for res in chain_residues(ch):
        for an in BB_ATOMS:
            p = res_atom_pos(res, an)
            if p is not None:
                coords.append(p)
    return np.array(coords, dtype=float)


def bb_coords_residx(struct: gemmi.Structure, residx_0based: List[int], chain_id: str = "A") -> np.ndarray:
    ch = get_chain(struct, chain_id)
    res = chain_residues(ch)
    coords = []
    for i in residx_0based:
        if 0 <= i < len(res):
            for an in BB_ATOMS:
                p = res_atom_pos(res[i], an)
                if p is not None:
                    coords.append(p)
    return np.array(coords, dtype=float)


PAE_KEYS = ("pae", "predicted_aligned_error", "predicted_aligned_error_matrix")


def load_json(path: Path) -> Optional[dict]:
    try:
        return json.loads(path.read_text())
    except Exception:
        return None


def mean_pae_from_conf(conf_json: Path) -> Tuple[Optional[float], str]:
    d = load_json(conf_json)
    if d is None:
        return None, "conf_json_unreadable"
    for k in PAE_KEYS:
        if k in d:
            arr = np.array(d[k], dtype=float)
            return float(np.nanmean(arr)), "has_pae_key"
    return None, "no_pae_key"


def extract_af3_scalar_metrics(summary_or_conf: dict) -> Dict[str, Any]:
    out = {}
    for k in ("ptm", "ranking_score", "fraction_disordered", "has_clash"):
        if k in summary_or_conf:
            out[k] = summary_or_conf.get(k)
    if "chain_pair_pae_min" in summary_or_conf:
        try:
            out["chain_pair_pae_min"] = float(np.array(summary_or_conf["chain_pair_pae_min"], dtype=float).min())
        except Exception:
            out["chain_pair_pae_min"] = None
    return out


def find_best_conf_file(af3_dir: Path, stem: str) -> Tuple[Optional[Path], Optional[Path]]:
    """
    Prefer:
      1) <stem>_confidences.json (top-level)
      2) seed-*/confidences.json
      3) <stem>_summary_confidences.json (fallback, scalars only)
    Returns (conf_path, summary_path)
    """
    conf_top = af3_dir / f"{stem}_confidences.json"
    summ_top = af3_dir / f"{stem}_summary_confidences.json"

    conf_seed = None
    for p in sorted(af3_dir.glob("seed-*/confidences.json")):
        conf_seed = p
        break

    conf_path = conf_top if conf_top.exists() else (conf_seed if conf_seed and conf_seed.exists() else None)
    summary_path = summ_top if summ_top.exists() else None
    return conf_path, summary_path


TOKEN_RE = re.compile(r"^([0-9a-f]{32})__contig", re.IGNORECASE)
PRED_RE = re.compile(r"_0_model_(\d+)$")


def parse_token_and_pred(dirname: str) -> Tuple[Optional[str], Optional[int]]:
    m = TOKEN_RE.search(dirname)
    token = m.group(1).lower() if m else None
    m2 = PRED_RE.search(dirname)
    pred = int(m2.group(1)) if m2 else None
    return token, pred


def sel_chain_or_segid(u: mda.Universe, chain_letters: List[str]) -> mda.core.groups.AtomGroup:
    seg_sel = " or ".join([f"segid {c}" for c in chain_letters])
    ag = u.select_atoms(seg_sel)
    if len(ag) > 0:
        return ag
    ch_sel = " or ".join([f"chainid {c}" for c in chain_letters])
    return u.select_atoms(ch_sel)


def pick_antibody_atoms(true_u: mda.Universe) -> mda.core.groups.AtomGroup:
    ab = sel_chain_or_segid(true_u, ["B", "C"]).select_atoms("protein and not name H*")
    if len(ab) > 0:
        return ab
    return true_u.select_atoms("protein and not (segid A or chainid A) and not name H*")


def gemmi_res_heavy_coords(res: gemmi.Residue) -> np.ndarray:
    coords = []
    for a in res:
        if a.element.name == "H":
            continue
        p = a.pos
        coords.append([p.x, p.y, p.z])
    return np.array(coords, dtype=float) if coords else np.zeros((0, 3), dtype=float)


def compute_af3_clash_resindices(
    af3_cif: Path,
    true_pdb: Path,
    af3_chunk_residx: List[int],    # AF3 chain-A residue index (0-based)
    true_chunk_residx: List[int],   # TRUE chain-A residue index (0-based)
    cutoff: float = 4.0,
) -> Tuple[Optional[List[int]], str]:
    """
    Returns (clashing_residx_list, status)
    status: ok / missing_chainA / bad_indices / too_few_pairs / no_ab
    """
    true_u = mda.Universe(str(true_pdb))
    true_ag = sel_chain_or_segid(true_u, ["A"])
    if len(true_ag) == 0:
        return None, "missing_true_chainA"
    true_res = true_ag.residues

    ab_atoms = pick_antibody_atoms(true_u)
    if len(ab_atoms) == 0:
        return None, "no_ab_atoms"
    ab_pos = ab_atoms.positions

    st = read_structure(af3_cif)
    chA = get_chain(st, "A")
    af3_res = chain_residues(chA)

    if not af3_chunk_residx or not true_chunk_residx:
        return None, "bad_indices"
    if max(af3_chunk_residx) >= len(af3_res) or max(true_chunk_residx) >= len(true_res):
        return None, "bad_indices"

    # Align using CA pairs (robust). (Backbone alignment is possible too but CA is stable.)
    P_list, Q_list = [], []
    for i_af3, i_true in zip(af3_chunk_residx, true_chunk_residx):
        caP = res_atom_pos(af3_res[i_af3], "CA")
        caQ = true_res[i_true].atoms.select_atoms("name CA and not name H*")
        if caP is None or len(caQ) != 1:
            continue
        P_list.append(caP)
        Q_list.append(caQ.positions[0])

    if len(P_list) < 3:
        return None, "too_few_pairs"

    P = np.vstack(P_list)
    Q = np.vstack(Q_list)
    R, t = kabsch_fit(P, Q)

    # unintended residues = all AF3 chain-A residues excluding chunk
    mask = np.zeros(len(af3_res), dtype=bool)
    mask[np.array(af3_chunk_residx, dtype=int)] = True
    unintended_idx = np.where(~mask)[0]

    clashing = []
    for i in unintended_idx:
        heavy = gemmi_res_heavy_coords(af3_res[i])
        if heavy.shape[0] == 0:
            continue
        heavy_fit = (heavy @ R) + t
        d = distance_array(heavy_fit, ab_pos)
        if np.any(d < cutoff):
            clashing.append(int(i))

    return clashing, "ok"


# ------------------------ main ------------------------

def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--run_dir", required=True)
    ap.add_argument("--dp2_parquet", default="datasets/dp2.parquet")
    ap.add_argument("--true_dir", required=True)
    ap.add_argument("--out_csv", default=None)
    ap.add_argument("--limit", type=int, default=0, help="Limit number of AF3 dirs (for testing). 0 = no limit.")
    args = ap.parse_args()

    run_dir = Path(args.run_dir).resolve()
    af3_root = run_dir / "03_af3" / "outputs"
    rfd3_root = run_dir / "02_rfd3" / "outputs"
    true_dir = Path(args.true_dir).resolve()

    out_csv = Path(args.out_csv).resolve() if args.out_csv else (run_dir / "04_filter" / "metrics_all.csv")
    out_csv.parent.mkdir(parents=True, exist_ok=True)

    # dp2 parquet
    dp2 = pd.read_parquet(Path(args.dp2_parquet).resolve())
    dp2["assay_scaffolded_epitope_id"] = dp2["assay_scaffolded_epitope_id"].astype(str).str.lower()

    # Choose dp2 columns to merge (keep broad)
    dp2_keep = [
        "assay_scaffolded_epitope_id",
        "id",
        "contig_id",
        "rfd_id",
        "mpnn_id",
        "epitope_chunk_resindices",
        "assay_scaffolded_epitope_chunk_resindices",
        "epitope_resindices",
        "assay_scaffolded_epitope_resindices",
        "epitope_chunk_rmsd_vs_mpnn",
        "overall_rmsd",
        "mean_pae",
        "scaffolded_epitope_seq",
        "assay_scaffolded_epitope_seq",
    ]
    dp2_keep = [c for c in dp2_keep if c in dp2.columns]
    dp2s = dp2[dp2_keep].copy()

    # Index AF3 dirs once
    af3_dirs = sorted([p for p in af3_root.iterdir() if p.is_dir()])
    if args.limit and args.limit > 0:
        af3_dirs = af3_dirs[: args.limit]

    rows: List[Dict[str, Any]] = []
    n = 0

    for d in af3_dirs:
        n += 1
        token, pred = parse_token_and_pred(d.name)

        # AF3 files
        stem = d.name
        af3_cif = d / f"{stem}_model.cif"
        if not af3_cif.exists():
            # some layouts also include seed-*/model.cif, but your top-level exists; we keep strict
            af3_cif = None

        conf_path, summary_path = find_best_conf_file(d, stem)
        af3_done = (d / "_DONE").exists()

        # Basic row
        row: Dict[str, Any] = {
            "assay_scaffolded_epitope_id": token,
            "pred": pred,
            "af3_dir": str(d),
            "af3_done": bool(af3_done),
            "af3_path": str(af3_cif) if af3_cif else None,
            "af3_conf_path": str(conf_path) if conf_path else None,
            "af3_summary_path": str(summary_path) if summary_path else None,
        }

        # Merge dp2 info by token (left join later is expensive row-by-row; quick lookup now)
        # We'll do a full merge at end; but we want per-row chunk indices now, so do small merge by token+pred later.
        rows.append(row)

    base = pd.DataFrame(rows)

    # Merge dp2 on token
    merged = base.merge(dp2s, on="assay_scaffolded_epitope_id", how="left", suffixes=("", "_dp2"))

    # Now compute metrics per row
    out_rows: List[Dict[str, Any]] = []
    missing_pair = 0

    for _, r in merged.iterrows():
        token = r.get("assay_scaffolded_epitope_id")
        pred = r.get("pred")

        af3_path = Path(r["af3_path"]) if isinstance(r.get("af3_path"), str) else None
        conf_path = Path(r["af3_conf_path"]) if isinstance(r.get("af3_conf_path"), str) else None
        summary_path = Path(r["af3_summary_path"]) if isinstance(r.get("af3_summary_path"), str) else None

        # RFD3 path uses AF3 dir name convention: token__contig..__0_model_pred.cif.gz
        # Your 02_rfd3 outputs are named exactly like AF3 dir + .cif.gz
        rfd3_cif_gz = rfd3_root / f"{Path(r['af3_dir']).name}.cif.gz"
        if not rfd3_cif_gz.exists():
            # fallback: try scan by token+pred (rare)
            missing_pair += 1
            rfd3_cif_gz = None

        row = dict(r)

        row["rfd3_path"] = str(rfd3_cif_gz) if rfd3_cif_gz else None

        # AF3 scalar + PAE
        pae_mean = None
        pae_reason = None
        scalars: Dict[str, Any] = {}

        if conf_path and conf_path.exists():
            pae_mean, pae_reason = mean_pae_from_conf(conf_path)
            dconf = load_json(conf_path)
            if dconf:
                scalars.update(extract_af3_scalar_metrics(dconf))
        elif summary_path and summary_path.exists():
            pae_mean, pae_reason = None, "no_conf_file"
            dsum = load_json(summary_path)
            if dsum:
                scalars.update(extract_af3_scalar_metrics(dsum))
        else:
            pae_mean, pae_reason = None, "no_conf_or_summary"

        row["pae_mean_all"] = pae_mean
        row["pae_reason"] = pae_reason
        for k, v in scalars.items():
            row[k] = v

        # RMSDs between RFD3 and AF3 (CA + backbone)
        if af3_path and af3_path.exists() and rfd3_cif_gz:
            try:
                st_rfd = read_structure(Path(rfd3_cif_gz))
                st_af3 = read_structure(af3_path)

                # all CA
                P = ca_coords(st_rfd, "A")
                Q = ca_coords(st_af3, "A")
                row["rmsd_ca_all"] = kabsch_rmsd(P, Q)

                # all backbone
                Pb = bb_coords(st_rfd, "A")
                Qb = bb_coords(st_af3, "A")
                row["rmsd_bb_all"] = kabsch_rmsd(Pb, Qb)

                # epitope indices: prefer dp2 assay_scaffolded_epitope_resindices if present else epitope_resindices
                ep0 = parse_index_list(row.get("assay_scaffolded_epitope_resindices")) or parse_index_list(row.get("epitope_resindices"))
                row["n_epitope_idx"] = len(ep0)

                Pe = ca_coords_residx(st_rfd, ep0, "A")
                Qe = ca_coords_residx(st_af3, ep0, "A")
                row["rmsd_ca_epitope"] = kabsch_rmsd(Pe, Qe) if (Pe.shape == Qe.shape and Pe.shape[0] >= 3) else float("nan")

                Pbe = bb_coords_residx(st_rfd, ep0, "A")
                Qbe = bb_coords_residx(st_af3, ep0, "A")
                row["rmsd_bb_epitope"] = kabsch_rmsd(Pbe, Qbe) if (Pbe.shape == Qbe.shape and Pbe.shape[0] >= 3) else float("nan")

            except Exception as e:
                row["rmsd_ca_all"] = None
                row["rmsd_bb_all"] = None
                row["rmsd_ca_epitope"] = None
                row["rmsd_bb_epitope"] = None
                row["rmsd_error"] = repr(e)
        else:
            row["rmsd_ca_all"] = None
            row["rmsd_bb_all"] = None
            row["rmsd_ca_epitope"] = None
            row["rmsd_bb_epitope"] = None

        # Clash metric (AF3 vs true antibody)
        af3_chunk = parse_index_list(row.get("assay_scaffolded_epitope_chunk_resindices"))
        true_chunk = parse_index_list(row.get("epitope_chunk_resindices"))
        pid = row.get("id")
        true_pdb = true_dir / f"{pid}.pdb" if isinstance(pid, str) else None

        if af3_path and af3_path.exists() and true_pdb and true_pdb.exists() and len(af3_chunk) >= 3 and len(true_chunk) >= 3:
            l, st = compute_af3_clash_resindices(af3_path, true_pdb, af3_chunk, true_chunk)
            row["af3_clash_status"] = st
            row["af3_clash_resindices"] = l
            row["af3_n_clash_res"] = (len(l) if isinstance(l, list) else None)
            row["af3_has_clash"] = (len(l) > 0) if isinstance(l, list) else None
        else:
            row["af3_clash_status"] = "not_eligible"
            row["af3_clash_resindices"] = None
            row["af3_n_clash_res"] = None
            row["af3_has_clash"] = None

        out_rows.append(row)

    out = pd.DataFrame(out_rows)
    out.to_csv(out_csv, index=False)

    print(f"Wrote {len(out):,} rows -> {out_csv}")
    if missing_pair:
        print(f"Missing RFD3 pairs: {missing_pair:,}")
    print("PAE reason counts:")
    print(out["pae_reason"].value_counts(dropna=False).to_string())
    print("AF3 clash status counts:")
    print(out["af3_clash_status"].value_counts(dropna=False).to_string())


if __name__ == "__main__":
    main()
