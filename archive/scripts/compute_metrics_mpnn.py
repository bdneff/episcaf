#!/usr/bin/env python3
"""
compute_metrics_mpnn.py  —  RFD3+MPNN → AF3 metric extraction

Computes four filters per (token, pred, mpnn_id) triple, matching Lawson's definitions:
  1. overall_rmsd          : backbone RMSD, MPNN scaffold vs AF3 prediction    (threshold ≤ 2 Å)
  2. epitope_chunk_rmsd    : backbone RMSD, epitope chunk only, MPNN vs AF3    (threshold ≤ 1 Å)
  3. mean_pae              : mean predicted aligned error from confidences.json (threshold < 5)
  4. af3_n_clash_res       : antibody residues within 4 Å of non-epitope atoms (threshold == 0)

Key difference from compute_metrics.py (RFD3-only):
  - Designed structure is the MPNN all-atom PDB, not the RFD3 CIF
  - pred_id format: {token}_pred{N}_fixed_dldesign_{M}
  - The MPNN sequence is used as the reference for RMSD (not the RFD3 sequence)

Usage:
    python scripts/compute_metrics_mpnn.py \
        --run_dir     runs/run_rfd3_mpnn \
        --dp2_parquet datasets/dp2.parquet \
        --true_dir    /tgen_labs/altin/alphafold3/workspace/episcaf-experiments/data/abdb/complex_pdbfiles/cleaned \
        --out_csv     runs/run_rfd3_mpnn/04_filter/metrics.csv

    # test on a small subset first
    python scripts/compute_metrics_mpnn.py \
        --run_dir     runs/run_rfd3_mpnn \
        --dp2_parquet datasets/dp2.parquet \
        --true_dir    /tgen_labs/altin/.../cleaned \
        --out_csv     runs/run_rfd3_mpnn/04_filter/metrics_test.csv \
        --limit       500
"""

from __future__ import annotations

import argparse
import gzip
import json
import math
import re
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
import gemmi
import MDAnalysis as mda
from MDAnalysis.analysis import rms as mda_rms
from MDAnalysis.lib.distances import distance_array


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

AA3TO1 = {
    "ALA": "A", "CYS": "C", "ASP": "D", "GLU": "E", "PHE": "F",
    "GLY": "G", "HIS": "H", "ILE": "I", "LYS": "K", "LEU": "L",
    "MET": "M", "ASN": "N", "PRO": "P", "GLN": "Q", "ARG": "R",
    "SER": "S", "THR": "T", "VAL": "V", "TRP": "W", "TYR": "Y",
}

BB_ATOMS = ("N", "CA", "C", "O")

# pred_id format: {32-hex-token}_pred{N}_fixed_dldesign_{M}
_PRED_ID_RE = re.compile(
    r"^([0-9a-fA-F]{32})_pred(\d+)_fixed_dldesign_(\d+)$",
    re.IGNORECASE,
)


# ---------------------------------------------------------------------------
# Structure I/O helpers  (identical to compute_metrics.py)
# ---------------------------------------------------------------------------

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


def get_chain(st: gemmi.Structure, chain_id: str = "A") -> gemmi.Chain:
    m = st[0]
    for ch in m:
        if ch.name == chain_id:
            return ch
    return m[0]


def chain_seq(chain: gemmi.Chain) -> str:
    return "".join(AA3TO1.get(r.name.upper(), "X") for r in chain)


def bb_coords_range(chain: gemmi.Chain, start: int, end: int) -> np.ndarray:
    residues = list(chain)[start:end]
    coords = []
    for r in residues:
        for atom_name in BB_ATOMS:
            a = r.find_atom(atom_name, altloc="*")
            if a:
                coords.append([a.pos.x, a.pos.y, a.pos.z])
    return np.array(coords, dtype=float)


def bb_coords_resindices(chain: gemmi.Chain, ris: List[int]) -> np.ndarray:
    residues = list(chain)
    coords = []
    for i in ris:
        r = residues[i]
        for atom_name in BB_ATOMS:
            a = r.find_atom(atom_name, altloc="*")
            if a:
                coords.append([a.pos.x, a.pos.y, a.pos.z])
    return np.array(coords, dtype=float)


def heavy_coords_residue(res: gemmi.Residue) -> np.ndarray:
    coords = []
    for a in res:
        if a.element.name != "H":
            coords.append([a.pos.x, a.pos.y, a.pos.z])
    return np.array(coords, dtype=float) if coords else np.zeros((0, 3), dtype=float)


def ca_coord(res: gemmi.Residue) -> Optional[np.ndarray]:
    a = res.find_atom("CA", altloc="*")
    if not a:
        return None
    return np.array([a.pos.x, a.pos.y, a.pos.z], dtype=float)


# ---------------------------------------------------------------------------
# RMSD / alignment helpers
# ---------------------------------------------------------------------------

def rmsd_superpose(P: np.ndarray, Q: np.ndarray) -> float:
    n = min(len(P), len(Q))
    if n < 3:
        return float("nan")
    return float(mda_rms.rmsd(P[:n], Q[:n], superposition=True))


def find_subseq(haystack: str, needle: str) -> int:
    i = haystack.find(needle)
    if i >= 0:
        return i
    Lh, Ln = len(haystack), len(needle)
    for i in range(Lh - Ln + 1):
        if all(h == "X" or h == n for h, n in zip(haystack[i:i + Ln], needle)):
            return i
    return -1


def kabsch_fit(P: np.ndarray, Q: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    Pc = P - P.mean(axis=0)
    Qc = Q - Q.mean(axis=0)
    V, _, Wt = np.linalg.svd(Pc.T @ Qc)
    d = np.sign(np.linalg.det(V @ Wt))
    R = V @ np.diag([1.0, 1.0, d]) @ Wt
    t = Q.mean(axis=0) - (P.mean(axis=0) @ R)
    return R, t


# ---------------------------------------------------------------------------
# PAE / summary helpers
# ---------------------------------------------------------------------------

def mean_pae_from_conf(conf_json: Path) -> Optional[float]:
    try:
        d = json.loads(conf_json.read_text())
    except Exception:
        return None
    for key in ("pae", "predicted_aligned_error", "predicted_aligned_error_matrix"):
        if key in d:
            try:
                return float(np.nanmean(np.array(d[key], dtype=float)))
            except Exception:
                pass
    return None


def summary_scalars(summary_json: Optional[Path]) -> Dict[str, Any]:
    out: Dict[str, Any] = {
        "ptm": None, "ranking_score": None,
        "fraction_disordered": None, "has_clash_af3": None,
    }
    if summary_json is None or not summary_json.exists():
        return out
    try:
        d = json.loads(summary_json.read_text())
    except Exception:
        return out
    for k in ("ptm", "ranking_score", "fraction_disordered"):
        try:
            out[k] = float(d[k]) if d.get(k) is not None else None
        except Exception:
            pass
    try:
        out["has_clash_af3"] = float(d["has_clash"])
    except Exception:
        pass
    return out


def find_af3_files(af3_dir: Path) -> Tuple[Optional[Path], Optional[Path], Optional[Path]]:
    cif = next(af3_dir.glob("*_model.cif"), None) or next(af3_dir.rglob("model.cif"), None)
    conf = None
    for p in af3_dir.glob("*_confidences.json"):
        if not p.name.endswith("_summary_confidences.json"):
            conf = p
            break
    if conf is None:
        for p in af3_dir.rglob("confidences.json"):
            conf = p
            break
    summary = next(af3_dir.glob("*_summary_confidences.json"), None) or \
              next(af3_dir.rglob("summary_confidences.json"), None)
    return cif, conf, summary


# ---------------------------------------------------------------------------
# Clash detection  (identical logic to compute_metrics.py)
# ---------------------------------------------------------------------------

def compute_af3_clash_resindices(
    af3_cif: Path,
    true_pdb: Path,
    af3_epitope_ris: List[int],
    true_epitope_ris: List[int],
    cutoff: float = 4.0,
) -> Tuple[Optional[List[int]], str]:
    true_u = mda.Universe(str(true_pdb))

    true_ag = true_u.select_atoms("segid A")
    if len(true_ag) == 0:
        true_ag = true_u.select_atoms("chainid A")
    if len(true_ag) == 0:
        return None, "no_true_chainA"

    ab_atoms = true_u.select_atoms("(segid B or segid C) and not name H*")
    if len(ab_atoms) == 0:
        ab_atoms = true_u.select_atoms("(chainid B or chainid C) and not name H*")
    if len(ab_atoms) == 0:
        return None, "no_antibody_atoms"

    true_res = true_ag.residues

    st = read_structure(af3_cif)
    chA = get_chain(st, "A")
    af3_res = list(chA)

    if not af3_epitope_ris or not true_epitope_ris:
        return None, "empty_epitope_indices"
    if max(af3_epitope_ris) >= len(af3_res):
        return None, f"af3_epitope_idx_oob (max={max(af3_epitope_ris)}, nres={len(af3_res)})"
    if max(true_epitope_ris) >= len(true_res):
        return None, f"true_epitope_idx_oob (max={max(true_epitope_ris)}, nres={len(true_res)})"

    P_list, Q_list = [], []
    for i_af3, i_true in zip(af3_epitope_ris, true_epitope_ris):
        p = ca_coord(af3_res[i_af3])
        q_sel = true_res[i_true].atoms.select_atoms("name CA and not name H*")
        if p is not None and len(q_sel) == 1:
            P_list.append(p)
            Q_list.append(q_sel.positions[0])

    if len(P_list) < 3:
        return None, f"too_few_ca_pairs ({len(P_list)})"

    R, t = kabsch_fit(np.vstack(P_list), np.vstack(Q_list))

    mask = np.zeros(len(af3_res), dtype=bool)
    mask[np.array(af3_epitope_ris, dtype=int)] = True
    unintended_idx = np.where(~mask)[0]

    ab_pos = ab_atoms.positions
    clashing = []
    for i in unintended_idx:
        heavy = heavy_coords_residue(af3_res[i])
        if heavy.shape[0] == 0:
            continue
        heavy_fit = (heavy @ R) + t
        if np.any(distance_array(heavy_fit, ab_pos) < cutoff):
            clashing.append(int(i))

    return clashing, "ok"


# ---------------------------------------------------------------------------
# Index AF3 outputs
# ---------------------------------------------------------------------------

def index_af3_outputs(af3_root: Path) -> Dict[Tuple[str, int, int], Path]:
    """Map (token, pred_idx, mpnn_id) -> AF3 output directory."""
    idx: Dict[Tuple[str, int, int], Path] = {}
    if not af3_root.exists():
        return idx
    for d in af3_root.iterdir():
        if not d.is_dir() or d.name.startswith("_"):
            continue
        m = _PRED_ID_RE.match(d.name)
        if m:
            tok     = m.group(1).lower()
            pred    = int(m.group(2))
            mpnn_id = int(m.group(3))
            idx[(tok, pred, mpnn_id)] = d
    return idx


def index_mpnn_pdbs(mpnn_root: Path) -> Dict[Tuple[str, int, int], Path]:
    """Map (token, pred_idx, mpnn_id) -> MPNN all-atom PDB path."""
    idx: Dict[Tuple[str, int, int], Path] = {}
    if not mpnn_root.exists():
        return idx
    for p in mpnn_root.rglob("*_fixed_dldesign_*.pdb"):
        m = _PRED_ID_RE.match(p.stem)
        if m:
            tok     = m.group(1).lower()
            pred    = int(m.group(2))
            mpnn_id = int(m.group(3))
            idx[(tok, pred, mpnn_id)] = p
    return idx


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


# ---------------------------------------------------------------------------
# Core per-triple metric computation
# ---------------------------------------------------------------------------

def compute_triple_metrics(
    mpnn_pdb:         Path,
    af3_cif:          Path,
    conf_json:        Optional[Path],
    summary_json:     Optional[Path],
    mpnn_epitope_ris: List[int],   # 0-based indices into MPNN chain A (= RFD3 epitope indices)
    true_pdb:         Optional[Path],
    true_epitope_ris: List[int],
    clash_cutoff:     float = 4.0,
) -> Dict[str, Any]:
    result: Dict[str, Any] = {
        "mean_pae":             None,
        "overall_rmsd":         None,
        "epitope_chunk_rmsd":   None,
        "af3_window_start":     None,
        "af3_window_end":       None,
        "af3_clash_resindices": None,
        "af3_n_clash_res":      None,
        "af3_has_clash":        None,
        "af3_clash_status":     "not_run",
        "status":               "ok",
    }

    # PAE
    if conf_json and conf_json.exists():
        result["mean_pae"] = mean_pae_from_conf(conf_json)

    result.update(summary_scalars(summary_json))

    # Load structures
    try:
        st_mpnn = read_structure(mpnn_pdb)
        st_af3  = read_structure(af3_cif)
        ch_mpnn = get_chain(st_mpnn, "A")
        ch_af3  = get_chain(st_af3,  "A")
    except Exception as e:
        result["status"] = f"structure_load_fail: {e}"
        return result

    mpnn_seq = chain_seq(ch_mpnn)
    af3_seq  = chain_seq(ch_af3)

    # Locate MPNN sequence inside AF3 output
    ws = find_subseq(af3_seq, mpnn_seq)
    if ws < 0:
        result["status"] = "seqmatch_fail"
        return result

    we = ws + len(mpnn_seq)
    result["af3_window_start"] = ws
    result["af3_window_end"]   = we

    # Overall RMSD: full scaffold backbone, MPNN vs AF3
    try:
        P = bb_coords_range(ch_mpnn, 0, len(mpnn_seq))
        Q = bb_coords_range(ch_af3,  ws, we)
        result["overall_rmsd"] = rmsd_superpose(P, Q)
    except Exception as e:
        result["status"] = f"overall_rmsd_fail: {e}"
        return result

    # Epitope chunk RMSD
    if mpnn_epitope_ris:
        try:
            af3_epitope_ris = [ws + i for i in mpnn_epitope_ris]
            P_epi = bb_coords_resindices(ch_mpnn, mpnn_epitope_ris)
            Q_epi = bb_coords_resindices(ch_af3,  af3_epitope_ris)
            result["epitope_chunk_rmsd"] = rmsd_superpose(P_epi, Q_epi)
        except Exception as e:
            result["status"] = f"epitope_rmsd_fail: {e}"

    # Clash detection
    if true_pdb and true_pdb.exists() and mpnn_epitope_ris and true_epitope_ris:
        af3_epitope_ris = [ws + i for i in mpnn_epitope_ris]
        clash_list, clash_status = compute_af3_clash_resindices(
            af3_cif, true_pdb,
            af3_epitope_ris, true_epitope_ris,
            cutoff=clash_cutoff,
        )
        result["af3_clash_status"]     = clash_status
        result["af3_clash_resindices"] = clash_list
        if clash_list is not None:
            result["af3_n_clash_res"] = len(clash_list)
            result["af3_has_clash"]   = len(clash_list) > 0

    return result


# ---------------------------------------------------------------------------
# Main run
# ---------------------------------------------------------------------------

def run_metrics(args: argparse.Namespace) -> None:
    run_dir   = Path(args.run_dir).resolve()
    mpnn_root = run_dir / "02_mpnn_pdbs"
    af3_root  = run_dir / "03_af3" / "outputs"
    true_dir  = Path(args.true_dir).resolve()
    out_csv   = Path(args.out_csv).resolve()

    dp2 = pd.read_parquet(Path(args.dp2_parquet).resolve())
    dp2["assay_scaffolded_epitope_id"] = dp2["assay_scaffolded_epitope_id"].astype(str).str.lower()

    required = ["assay_scaffolded_epitope_id", "id",
                "epitope_chunk_resindices", "scaffolded_epitope_chunk_resindices"]
    missing_cols = [c for c in required if c not in dp2.columns]
    if missing_cols:
        sys.exit(f"dp2 missing required columns: {missing_cols}")

    # Build per-token lookup for epitope indices and protein id
    dp2_tok = dp2.drop_duplicates("assay_scaffolded_epitope_id").set_index(
        "assay_scaffolded_epitope_id"
    )

    print(f"[metrics] Indexing MPNN PDBs from {mpnn_root} ...")
    mpnn_idx = index_mpnn_pdbs(mpnn_root)
    print(f"[metrics] MPNN index size: {len(mpnn_idx):,}")

    print(f"[metrics] Indexing AF3 outputs from {af3_root} ...")
    af3_idx = index_af3_outputs(af3_root)
    print(f"[metrics] AF3 index size:  {len(af3_idx):,}")

    # Union of all (tok, pred, mpnn_id) keys present in either index
    all_keys = sorted(set(mpnn_idx.keys()) | set(af3_idx.keys()))
    if args.limit and args.limit > 0:
        all_keys = all_keys[:args.limit]

    print(f"[metrics] Total (tok, pred, mpnn_id) triples: {len(all_keys):,}")

    rows = []
    n_ok = n_missing = n_fail = 0

    for i, (tok, pred, mpnn_id) in enumerate(all_keys):
        mpnn_pdb = mpnn_idx.get((tok, pred, mpnn_id))
        af3_dir  = af3_idx.get((tok, pred, mpnn_id))

        af3_cif = conf_json = summary_json = None
        if af3_dir:
            af3_cif, conf_json, summary_json = find_af3_files(af3_dir)

        # Look up dp2 metadata for this token
        dp2_row = dp2_tok.loc[tok] if tok in dp2_tok.index else None
        pid            = dp2_row["id"] if dp2_row is not None else None
        mpnn_epi_ris   = parse_index_list(dp2_row["scaffolded_epitope_chunk_resindices"]) \
                         if dp2_row is not None else []
        true_epi_ris   = parse_index_list(dp2_row["epitope_chunk_resindices"]) \
                         if dp2_row is not None else []
        true_pdb       = (true_dir / f"{pid}.pdb") if pid else None

        base_row: Dict[str, Any] = {
            "token":         tok,
            "pred":          pred,
            "mpnn_id":       mpnn_id,
            "id":            pid,
            "mpnn_pdb":      str(mpnn_pdb) if mpnn_pdb else None,
            "af3_dir":       str(af3_dir)  if af3_dir  else None,
            "mean_pae":              None,
            "overall_rmsd":          None,
            "epitope_chunk_rmsd":    None,
            "af3_window_start":      None,
            "af3_window_end":        None,
            "af3_clash_resindices":  None,
            "af3_n_clash_res":       None,
            "af3_has_clash":         None,
            "af3_clash_status":      "not_run",
            "ptm":                   None,
            "ranking_score":         None,
            "fraction_disordered":   None,
            "has_clash_af3":         None,
            "status":                "missing_pair",
        }

        if mpnn_pdb is None or af3_cif is None:
            n_missing += 1
            rows.append(base_row)
            continue

        metrics = compute_triple_metrics(
            mpnn_pdb         = mpnn_pdb,
            af3_cif          = af3_cif,
            conf_json        = conf_json,
            summary_json     = summary_json,
            mpnn_epitope_ris = mpnn_epi_ris,
            true_pdb         = true_pdb,
            true_epitope_ris = true_epi_ris,
            clash_cutoff     = args.clash_cutoff,
        )
        base_row.update(metrics)

        if metrics["status"] == "ok":
            n_ok += 1
        else:
            n_fail += 1

        rows.append(base_row)

        if (i + 1) % 5000 == 0:
            print(f"[metrics] {i+1:,}/{len(all_keys):,}  ok={n_ok:,}  missing={n_missing:,}  fail={n_fail:,}")

    out = pd.DataFrame(rows)
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(out_csv, index=False)

    # Apply filters
    passing = out.loc[
        out["mean_pae"].notna() &
        (out["mean_pae"] < 5.0) &
        out["overall_rmsd"].notna() &
        (out["overall_rmsd"] <= 2.0) &
        out["epitope_chunk_rmsd"].notna() &
        (out["epitope_chunk_rmsd"] <= 1.0) &
        (out["af3_n_clash_res"].fillna(1) == 0)
    ]

    n_total = len(out)
    n_with_af3 = out["overall_rmsd"].notna().sum()

    print(f"\n[metrics] Wrote {n_total:,} rows -> {out_csv}")
    print(f"[metrics] ok={n_ok:,}  missing={n_missing:,}  fail={n_fail:,}")
    print(f"[metrics] Rows with AF3 output: {n_with_af3:,}")
    print(f"\n=== FILTER SUMMARY (Lawson thresholds) ===")
    print(f"  mean_pae < 5       : {(out['mean_pae'] < 5.0).sum():,}")
    print(f"  overall_rmsd ≤ 2Å  : {(out['overall_rmsd'] <= 2.0).sum():,}")
    print(f"  epitope_rmsd ≤ 1Å  : {(out['epitope_chunk_rmsd'] <= 1.0).sum():,}")
    print(f"  no clashes         : {(out['af3_n_clash_res'].fillna(1) == 0).sum():,}")
    print(f"  ALL FOUR PASS      : {len(passing):,}  ({100*len(passing)/max(n_with_af3,1):.2f}% of designs with AF3 output)")

    print(f"\n=== PASS RATE COMPARISON ===")
    print(f"  Lawson RFD1+MPNN   : ~0.55%  (840 / ~150,000)")
    print(f"  RFD3-only          : ~0.19%  (36 / 18,880)")
    print(f"  RFD3+MPNN (this)   : {100*len(passing)/max(n_with_af3,1):.2f}%  ({len(passing):,} / {n_with_af3:,})")

    if len(passing) > 0:
        print(f"\nTop 20 passing designs:")
        print(
            passing.sort_values("epitope_chunk_rmsd")
            [["id", "token", "pred", "mpnn_id", "overall_rmsd", "epitope_chunk_rmsd", "mean_pae", "af3_n_clash_res"]]
            .head(20).to_string(index=False)
        )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--run_dir",      required=True,
                        help="Run directory (contains 02_mpnn_pdbs/ and 03_af3/outputs/)")
    parser.add_argument("--dp2_parquet",  required=True,
                        help="Path to dp2.parquet for epitope index lookup")
    parser.add_argument("--true_dir",     required=True,
                        help="Directory with true complex PDBs named <id>.pdb")
    parser.add_argument("--out_csv",      required=True,
                        help="Output CSV path")
    parser.add_argument("--clash_cutoff", type=float, default=4.0,
                        help="Clash distance cutoff in Å (default: 4.0)")
    parser.add_argument("--limit",        type=int, default=0,
                        help="Limit to first N triples for testing (0=all)")
    args = parser.parse_args()
    run_metrics(args)


if __name__ == "__main__":
    main()
