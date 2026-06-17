#!/usr/bin/env python3
"""
compute_metrics.py  —  RFD3 → AF3 metric extraction (no-MPNN pipeline)

Computes four filters per (tok, pred) pair, matching Lawson's definitions:
  1. overall_rmsd          : backbone RMSD, full RFD3 scaffold vs AF3 prediction  (threshold ≤ 2 Å)
  2. epitope_chunk_rmsd    : backbone RMSD, epitope chunk only, RFD3 vs AF3        (threshold ≤ 1 Å)
  3. mean_pae              : mean predicted aligned error from AF3 confidences.json (threshold < 5)
  4. af3_clash_resindices  : antibody residue indices within 4 Å of non-epitope    (threshold: list len == 0)

Validation mode (--validate):
  Reads dp2.parquet (Lawson's ground truth), computes overall_rmsd and
  epitope_chunk_rmsd using his MPNN PDBs + AF3 CIFs, and diffs against
  stored dp2 values. Use this FIRST before running on your RFD3 data.

Usage — validation against Lawson dp2:
  python compute_metrics.py validate \\
    --dp2_parquet  datasets/dp2.parquet \\
    --lawson_root  /path/to/sourced_antibody_v1/no_antibody \\
    --out_csv      runs/validation_result.csv \\
    --sample       500

Usage — compute metrics on your RFD3 run:
  python compute_metrics.py run \\
    --run_dir      runs/run_test_rfd3_nompmn \\
    --dp2_parquet  datasets/dp2.parquet \\
    --true_dir     /tgen_labs/altin/.../complex_pdbfiles/cleaned \\
    --out_csv      runs/run_test_rfd3_nompmn/04_filter/metrics.csv \\
    --clash_cutoff 4.0
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

# Regex to extract (token, pred_index) from RFD3/AF3 filenames/dirnames.
# Token: 32 hex chars at the start.
# Pred:  _0_model_<N> anywhere in the name.
_TOK_RE  = re.compile(r"^([0-9a-fA-F]{32})", re.IGNORECASE)
_PRED_RE = re.compile(r"_0_model_(\d+)", re.IGNORECASE)


# ---------------------------------------------------------------------------
# Structure I/O helpers
# ---------------------------------------------------------------------------

def read_structure(path: Path) -> gemmi.Structure:
    """Read a CIF, CIF.GZ, or PDB file into a gemmi Structure."""
    s = str(path)
    if s.endswith(".gz"):
        with gzip.open(path, "rt") as f:
            doc = gemmi.cif.read_string(f.read())
        return gemmi.make_structure_from_block(doc.sole_block())
    if path.suffix.lower() == ".cif":
        doc = gemmi.cif.read(str(path))
        return gemmi.make_structure_from_block(doc.sole_block())
    return gemmi.read_structure(str(path))  # PDB


def get_chain(st: gemmi.Structure, chain_id: str = "A") -> gemmi.Chain:
    m = st[0]
    for ch in m:
        if ch.name == chain_id:
            return ch
    return m[0]  # fallback: first chain


def chain_seq(chain: gemmi.Chain) -> str:
    return "".join(AA3TO1.get(r.name.upper(), "X") for r in chain)


def bb_coords_range(chain: gemmi.Chain, start: int, end: int) -> np.ndarray:
    """Backbone (N,CA,C,O) coordinates for residues [start, end) by 0-based index."""
    residues = list(chain)[start:end]
    coords = []
    for r in residues:
        for atom_name in BB_ATOMS:
            a = r.find_atom(atom_name, altloc="*")
            if a:
                coords.append([a.pos.x, a.pos.y, a.pos.z])
    return np.array(coords, dtype=float)


def bb_coords_resindices(chain: gemmi.Chain, ris: List[int]) -> np.ndarray:
    """Backbone coordinates for a list of 0-based residue indices."""
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
    """All non-hydrogen atom coordinates for a residue."""
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
# RMSD helpers
# ---------------------------------------------------------------------------

def rmsd_superpose(P: np.ndarray, Q: np.ndarray) -> float:
    """
    RMSD after optimal superposition (Kabsch via MDAnalysis).
    Matches Lawson's rms.rmsd(..., superposition=True).
    """
    n = min(len(P), len(Q))
    if n < 3:
        return float("nan")
    return float(mda_rms.rmsd(P[:n], Q[:n], superposition=True))


def find_subseq(haystack: str, needle: str) -> int:
    """
    Find needle inside haystack, allowing 'X' in haystack to match any residue.
    Returns start index (0-based) or -1.
    """
    # fast exact path
    i = haystack.find(needle)
    if i >= 0:
        return i
    Lh, Ln = len(haystack), len(needle)
    for i in range(Lh - Ln + 1):
        if all(h == "X" or h == n for h, n in zip(haystack[i:i + Ln], needle)):
            return i
    return -1


# ---------------------------------------------------------------------------
# Kabsch fit (for clash alignment)
# ---------------------------------------------------------------------------

def kabsch_fit(P: np.ndarray, Q: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """
    Compute rotation R and translation t such that (P @ R) + t ≈ Q.
    Used to superimpose AF3 epitope onto true epitope before clash detection.
    """
    Pc = P - P.mean(axis=0)
    Qc = Q - Q.mean(axis=0)
    V, _, Wt = np.linalg.svd(Pc.T @ Qc)
    d = np.sign(np.linalg.det(V @ Wt))
    R = V @ np.diag([1.0, 1.0, d]) @ Wt
    t = Q.mean(axis=0) - (P.mean(axis=0) @ R)
    return R, t


# ---------------------------------------------------------------------------
# PAE helper
# ---------------------------------------------------------------------------

def mean_pae_from_conf(conf_json: Path) -> Optional[float]:
    """
    Read mean PAE from AF3 confidences.json.
    Tries keys: 'pae', 'predicted_aligned_error', 'predicted_aligned_error_matrix'.
    """
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
        "chain_pair_pae_min": None,
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
    try:
        cppm = d.get("chain_pair_pae_min")
        if isinstance(cppm, list) and cppm and isinstance(cppm[0], list) and cppm[0]:
            out["chain_pair_pae_min"] = float(cppm[0][0])
        elif isinstance(cppm, (int, float)):
            out["chain_pair_pae_min"] = float(cppm)
    except Exception:
        pass
    return out


# ---------------------------------------------------------------------------
# Clash detection  (Lawson-style)
# ---------------------------------------------------------------------------

def compute_af3_clash_resindices(
    af3_cif: Path,
    true_pdb: Path,
    af3_epitope_ris: List[int],    # 0-based indices into AF3 chain A
    true_epitope_ris: List[int],   # 0-based indices into true chain A
    cutoff: float = 4.0,
) -> Tuple[Optional[List[int]], str]:
    """
    Lawson-style clash detection:
      1. Superimpose AF3 epitope chunk CA onto true epitope chunk CA (Kabsch)
      2. Apply that transform to all AF3 non-epitope (scaffold) heavy atoms
      3. Find antibody (chain B/C) heavy atoms within `cutoff` Å
      4. Return list of antibody residue resindices that are clashing

    Returns (clash_resindices_list, status_string).
    """
    # --- Load true complex ---
    true_u = mda.Universe(str(true_pdb))

    # True antigen (chain A) residues
    true_ag = true_u.select_atoms("segid A")
    if len(true_ag) == 0:
        true_ag = true_u.select_atoms("chainid A")
    if len(true_ag) == 0:
        return None, "no_true_chainA"

    # True antibody atoms (chains B + C), no hydrogens
    ab_atoms = true_u.select_atoms("(segid B or segid C) and not name H*")
    if len(ab_atoms) == 0:
        ab_atoms = true_u.select_atoms("(chainid B or chainid C) and not name H*")
    if len(ab_atoms) == 0:
        return None, "no_antibody_atoms"

    true_res = true_ag.residues

    # --- Load AF3 structure (gemmi, avoids MDA CIF issues) ---
    st = read_structure(af3_cif)
    chA = get_chain(st, "A")
    af3_res = list(chA)

    if not af3_epitope_ris or not true_epitope_ris:
        return None, "empty_epitope_indices"
    if max(af3_epitope_ris) >= len(af3_res):
        return None, f"af3_epitope_idx_oob (max={max(af3_epitope_ris)}, nres={len(af3_res)})"
    if max(true_epitope_ris) >= len(true_res):
        return None, f"true_epitope_idx_oob (max={max(true_epitope_ris)}, nres={len(true_res)})"

    # --- Build CA alignment pairs ---
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

    # --- Find unintended (non-epitope) AF3 residues ---
    mask = np.zeros(len(af3_res), dtype=bool)
    mask[np.array(af3_epitope_ris, dtype=int)] = True
    unintended_idx = np.where(~mask)[0]

    ab_pos = ab_atoms.positions

    # --- Check heavy atom distance of each unintended residue to antibody ---
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
# File indexing
# ---------------------------------------------------------------------------

def index_rfd3_outputs(rfd3_root: Path) -> Dict[Tuple[str, int], Path]:
    """Map (token, pred) -> path to RFD3 .cif.gz file."""
    idx: Dict[Tuple[str, int], Path] = {}
    if not rfd3_root.exists():
        return idx
    for p in rfd3_root.iterdir():
        if not p.is_file():
            continue
        m_tok  = _TOK_RE.match(p.name)
        m_pred = _PRED_RE.search(p.name)
        if m_tok and m_pred and (p.name.endswith(".cif.gz") or p.name.endswith(".cif")):
            idx[(m_tok.group(1).lower(), int(m_pred.group(1)))] = p
    return idx


def index_af3_outputs(af3_root: Path) -> Dict[Tuple[str, int], Path]:
    """Map (token, pred) -> path to AF3 output directory."""
    idx: Dict[Tuple[str, int], Path] = {}
    if not af3_root.exists():
        return idx
    for d in af3_root.iterdir():
        if not d.is_dir():
            continue
        m_tok  = _TOK_RE.match(d.name)
        m_pred = _PRED_RE.search(d.name)
        if m_tok and m_pred:
            idx[(m_tok.group(1).lower(), int(m_pred.group(1)))] = d
    return idx


def find_af3_files(af3_dir: Path) -> Tuple[Optional[Path], Optional[Path], Optional[Path]]:
    """Return (model_cif, confidences_json, summary_confidences_json) from an AF3 output dir."""
    cif     = next(af3_dir.glob("*_model.cif"), None) or next(af3_dir.rglob("model.cif"), None)
    conf    = next(af3_dir.glob("*_confidences.json"), None) or next(af3_dir.rglob("confidences.json"), None)
    summary = next(af3_dir.glob("*_summary_confidences.json"), None) or next(af3_dir.rglob("summary_confidences.json"), None)
    return cif, conf, summary


# ---------------------------------------------------------------------------
# Index list parsing  (handles list, ndarray, or stringified "[1 2 3]")
# ---------------------------------------------------------------------------

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
# Core per-pair metric computation
# ---------------------------------------------------------------------------

def compute_pair_metrics(
    rfd3_cif:        Path,
    af3_cif:         Path,
    conf_json:       Optional[Path],
    summary_json:    Optional[Path],
    rfd3_epitope_ris: List[int],   # 0-based indices into RFD3 chain A
    true_pdb:        Optional[Path],
    true_epitope_ris: List[int],   # 0-based indices into true chain A
    clash_cutoff:    float = 4.0,
) -> Dict[str, Any]:
    """
    Compute all four metrics for one (RFD3, AF3) pair.

    Epitope RMSD:
      - The AF3 chain A contains the same sequence as RFD3 chain A (since AF3
        received exactly the RFD3 all-atom sequence as input).
      - We locate the RFD3 sequence inside AF3 chain A via substring match
        (allowing X for unknown residues in AF3).
      - The AF3 epitope resindices = window_start + rfd3_epitope_ris.
      - This is the key difference from Lawson's pipeline (he used pre-computed
        assay_scaffolded_epitope_chunk_resindices; we derive them from the match).
    """
    result: Dict[str, Any] = {
        "mean_pae":              None,
        "overall_rmsd":          None,
        "epitope_chunk_rmsd":    None,
        "af3_window_start":      None,
        "af3_window_end":        None,
        "af3_clash_resindices":  None,
        "af3_n_clash_res":       None,
        "af3_has_clash":         None,
        "af3_clash_status":      "not_run",
        "status":                "ok",
    }

    # --- PAE ---
    if conf_json and conf_json.exists():
        result["mean_pae"] = mean_pae_from_conf(conf_json)

    # --- AF3 scalar metrics ---
    result.update(summary_scalars(summary_json))

    # --- Load structures once ---
    try:
        st_rfd = read_structure(rfd3_cif)
        st_af3 = read_structure(af3_cif)
        ch_rfd = get_chain(st_rfd, "A")
        ch_af3 = get_chain(st_af3, "A")
    except Exception as e:
        result["status"] = f"structure_load_fail: {e}"
        return result

    rfd3_seq = chain_seq(ch_rfd)
    af3_seq  = chain_seq(ch_af3)

    # --- Sequence window match: locate RFD3 inside AF3 ---
    ws = find_subseq(af3_seq, rfd3_seq)
    if ws < 0:
        result["status"] = "seqmatch_fail"
        return result

    we = ws + len(rfd3_seq)
    result["af3_window_start"] = ws
    result["af3_window_end"]   = we

    # --- Overall RMSD (full scaffold backbone, Lawson threshold ≤ 2 Å) ---
    try:
        P = bb_coords_range(ch_rfd, 0, len(rfd3_seq))
        Q = bb_coords_range(ch_af3, ws, we)
        result["overall_rmsd"] = rmsd_superpose(P, Q)
    except Exception as e:
        result["status"] = f"overall_rmsd_fail: {e}"
        return result

    # --- Epitope chunk RMSD (Lawson threshold ≤ 1 Å) ---
    # Key fix: AF3 epitope indices = window_start + RFD3 epitope indices
    if rfd3_epitope_ris:
        try:
            af3_epitope_ris = [ws + i for i in rfd3_epitope_ris]
            P_epi = bb_coords_resindices(ch_rfd, rfd3_epitope_ris)
            Q_epi = bb_coords_resindices(ch_af3, af3_epitope_ris)
            result["epitope_chunk_rmsd"] = rmsd_superpose(P_epi, Q_epi)
        except Exception as e:
            result["status"] = f"epitope_rmsd_fail: {e}"

    # --- Clash detection (Lawson threshold: 0 clashing antibody residues) ---
    if true_pdb and true_pdb.exists() and rfd3_epitope_ris and true_epitope_ris:
        af3_epitope_ris = [ws + i for i in rfd3_epitope_ris]
        clash_list, clash_status = compute_af3_clash_resindices(
            af3_cif, true_pdb,
            af3_epitope_ris, true_epitope_ris,
            cutoff=clash_cutoff,
        )
        result["af3_clash_status"]    = clash_status
        result["af3_clash_resindices"] = clash_list
        if clash_list is not None:
            result["af3_n_clash_res"] = len(clash_list)
            result["af3_has_clash"]   = len(clash_list) > 0

    return result


# ---------------------------------------------------------------------------
# Validation mode: reproduce Lawson's MPNN→AF3 metrics and diff against dp2
# ---------------------------------------------------------------------------

def run_validate(args: argparse.Namespace) -> None:
    """
    Validate metric calculations against Lawson's stored dp2 values.
    Uses his MPNN PDB files + AF3 CIF files as inputs, computes overall_rmsd
    and epitope_chunk_rmsd_vs_mpnn, and reports abs differences.
    """
    dp2 = pd.read_parquet(args.dp2_parquet)
    dp2["assay_scaffolded_epitope_id"] = dp2["assay_scaffolded_epitope_id"].astype(str).str.lower()

    required = ["id", "contig_id", "rfd_id", "mpnn_id", "assay_scaffolded_epitope_id",
                "scaffolded_epitope_chunk_resindices", "assay_scaffolded_epitope_chunk_resindices",
                "overall_rmsd", "epitope_chunk_rmsd_vs_mpnn"]
    missing = [c for c in required if c not in dp2.columns]
    if missing:
        sys.exit(f"dp2 missing columns: {missing}")

    if args.sample and args.sample > 0:
        dp2 = dp2.sample(n=min(args.sample, len(dp2)), random_state=42).reset_index(drop=True)
    if args.n and args.n > 0:
        dp2 = dp2.head(args.n).reset_index(drop=True)

    root = Path(args.lawson_root)
    mpnn_root = root / "proteinmpnn"
    af3_root  = root / "af3_predictions"

    rows = []
    n_ok = n_missing = n_seqfail = n_epifail = 0

    for i, r in dp2.iterrows():
        pid    = str(r["id"])
        contig = int(r["contig_id"])
        rfd_id = int(r["rfd_id"])
        mpnn   = int(r["mpnn_id"])
        tok    = str(r["assay_scaffolded_epitope_id"])

        mpnn_pdb = mpnn_root / pid / str(contig) / f"{pid}_{rfd_id}_fixed_dldesign_{mpnn}.pdb"

        # AF3 CIF: try a few known layouts Lawson used
        af3_cif = af3_root / tok / "seed-1_sample-0" / "model.cif.gz"
        if not af3_cif.exists():
            af3_cif = af3_root / tok / f"{tok}_model.cif.gz"
        if not af3_cif.exists():
            af3_cif = af3_root / tok / f"{tok}_model.cif"

        out = {
            "id": pid, "contig_id": contig, "rfd_id": rfd_id, "mpnn_id": mpnn, "tok": tok,
            "pdb_exists": mpnn_pdb.exists(), "cif_exists": af3_cif.exists(),
            "dp2_overall_rmsd": r.get("overall_rmsd"),
            "dp2_epitope_rmsd": r.get("epitope_chunk_rmsd_vs_mpnn"),
            "calc_overall_rmsd": np.nan, "calc_epitope_rmsd": np.nan,
            "d_overall": np.nan, "d_epitope": np.nan,
            "af3_window_start": np.nan, "status": "ok",
        }

        if not (mpnn_pdb.exists() and af3_cif.exists()):
            out["status"] = "file_missing"
            n_missing += 1
            rows.append(out)
            continue

        # Load structures
        try:
            st_af3 = read_structure(af3_cif)
            ch_af3 = get_chain(st_af3, "A")
            af3_seq = chain_seq(ch_af3)

            # MPNN side: use MDAnalysis (PDB file, segid A)
            u = mda.Universe(str(mpnn_pdb))
            selA = u.select_atoms("segid A")
            if len(selA) == 0:
                selA = u.select_atoms("chainid A")
            mpnn_seq = selA.residues.sequence(format="string")
        except Exception as e:
            out["status"] = f"load_fail: {e}"
            n_seqfail += 1
            rows.append(out)
            continue

        # Sequence match
        ws = find_subseq(af3_seq, mpnn_seq)
        if ws < 0:
            out["status"] = "seqmatch_fail"
            n_seqfail += 1
            rows.append(out)
            continue

        we = ws + len(mpnn_seq)
        out["af3_window_start"] = ws

        # Overall RMSD (Lawson: MPNN chain A backbone vs AF3 window backbone)
        try:
            P = selA.select_atoms("backbone").positions
            Q = bb_coords_range(ch_af3, ws, we)
            out["calc_overall_rmsd"] = rmsd_superpose(P, Q)
        except Exception as e:
            out["status"] = f"overall_rmsd_fail: {e}"
            n_seqfail += 1
            rows.append(out)
            continue

        # Epitope RMSD (Lawson: uses pre-computed assay_scaffolded_epitope_chunk_resindices)
        try:
            mpnn_ris = parse_index_list(r["scaffolded_epitope_chunk_resindices"])
            af3_ris  = parse_index_list(r["assay_scaffolded_epitope_chunk_resindices"])
            P_epi = selA.residues[mpnn_ris].atoms.select_atoms("backbone").positions
            Q_epi = bb_coords_resindices(ch_af3, af3_ris)
            out["calc_epitope_rmsd"] = rmsd_superpose(P_epi, Q_epi)
        except Exception as e:
            out["status"] = f"epitope_rmsd_fail: {e}"
            n_epifail += 1
            rows.append(out)
            continue

        if pd.notna(out["dp2_overall_rmsd"]) and np.isfinite(out["calc_overall_rmsd"]):
            out["d_overall"] = abs(float(out["dp2_overall_rmsd"]) - out["calc_overall_rmsd"])
        if pd.notna(out["dp2_epitope_rmsd"]) and np.isfinite(out["calc_epitope_rmsd"]):
            out["d_epitope"] = abs(float(out["dp2_epitope_rmsd"]) - out["calc_epitope_rmsd"])

        n_ok += 1
        rows.append(out)

        if (i + 1) % (args.progress_every or 200) == 0:
            print(f"[validate] {i+1}/{len(dp2)}  ok={n_ok}  missing={n_missing}  fail={n_seqfail+n_epifail}")

    df = pd.DataFrame(rows)
    Path(args.out_csv).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.out_csv, index=False)

    print("\n=== VALIDATION RESULTS ===")
    print(f"Total rows tested : {len(df)}")
    print(f"ok                : {n_ok}")
    print(f"file_missing      : {n_missing}")
    print(f"seqmatch/load fail: {n_seqfail}")
    print(f"epitope fail      : {n_epifail}")
    finite = df.loc[df["d_overall"].notna() & df["d_epitope"].notna(), ["d_overall", "d_epitope"]]
    if len(finite) > 0:
        print("\nAbs-difference summary (should be ~0 for both columns if correct):")
        print(finite.describe().to_string())
        bad_overall = (finite["d_overall"] > 0.01).sum()
        bad_epitope = (finite["d_epitope"] > 0.01).sum()
        print(f"\nRows with |d_overall| > 0.01 Å : {bad_overall}")
        print(f"Rows with |d_epitope| > 0.01 Å : {bad_epitope}")
        if bad_overall == 0 and bad_epitope == 0:
            print("\n✓ PASS: metric calculations match Lawson's stored values.")
        else:
            print("\n✗ FAIL: discrepancies found. Check the rows above.")
    else:
        print("\nNo rows with finite diffs to compare.")
    print(f"\nWrote: {args.out_csv}")


# ---------------------------------------------------------------------------
# Run mode: compute metrics on your RFD3 data
# ---------------------------------------------------------------------------

def run_metrics(args: argparse.Namespace) -> None:
    run_dir  = Path(args.run_dir).resolve()
    rfd3_root = run_dir / "02_rfd3" / "outputs"
    af3_root  = run_dir / "03_af3" / "outputs"
    true_dir  = Path(args.true_dir).resolve()

    dp2 = pd.read_parquet(Path(args.dp2_parquet).resolve())
    dp2["assay_scaffolded_epitope_id"] = dp2["assay_scaffolded_epitope_id"].astype(str).str.lower()

    required = ["assay_scaffolded_epitope_id", "id",
                "epitope_chunk_resindices", "scaffolded_epitope_chunk_resindices"]
    missing = [c for c in required if c not in dp2.columns]
    if missing:
        sys.exit(f"dp2 missing required columns: {missing}")

    if args.limit and args.limit > 0:
        dp2 = dp2.head(args.limit).copy()

    rfd3_idx = index_rfd3_outputs(rfd3_root)
    af3_idx  = index_af3_outputs(af3_root)

    print(f"[metrics] dp2 rows      : {len(dp2):,}")
    print(f"[metrics] RFD3 index    : {len(rfd3_idx):,}")
    print(f"[metrics] AF3 index     : {len(af3_idx):,}")

    rows = []
    n_missing = n_seqfail = n_ok = 0

    for _, r in dp2.iterrows():
        tok = str(r["assay_scaffolded_epitope_id"]).lower()
        pid = r.get("id")

        # 0-based resindices into RFD3 chain A for the epitope chunk
        rfd3_epi_ris = parse_index_list(r.get("scaffolded_epitope_chunk_resindices"))
        # 0-based resindices into true chain A for the epitope chunk
        true_epi_ris = parse_index_list(r.get("epitope_chunk_resindices"))

        true_pdb = (true_dir / f"{pid}.pdb") if isinstance(pid, str) else None

        for pred in range(args.pred_max + 1):
            rfd3_file = rfd3_idx.get((tok, pred))
            af3_dir   = af3_idx.get((tok, pred))

            af3_cif = conf_json = summary_json = None
            if af3_dir:
                af3_cif, conf_json, summary_json = find_af3_files(af3_dir)

            base_row: Dict[str, Any] = {
                "assay_scaffolded_epitope_id": tok,
                "pred":           pred,
                "id":             pid,
                "contig_id":      r.get("contig_id"),
                "rfd_id":         r.get("rfd_id"),
                "mpnn_id":        r.get("mpnn_id"),
                "contig_string":  r.get("contig_string"),
                "rfd3_path":      str(rfd3_file) if rfd3_file else None,
                "af3_dir":        str(af3_dir)   if af3_dir   else None,
                "af3_path":       str(af3_cif)   if af3_cif   else None,
                "af3_conf_path":  str(conf_json) if conf_json else None,
                # Pass-fail filter columns (populated below)
                "mean_pae":              None,
                "overall_rmsd":          None,
                "epitope_chunk_rmsd":    None,
                "af3_window_start":      None,
                "af3_window_end":        None,
                "af3_clash_resindices":  None,
                "af3_n_clash_res":       None,
                "af3_has_clash":         None,
                "af3_clash_status":      "not_run",
                "status":                "missing_pair",
                # AF3 scalar metrics
                "ptm": None, "ranking_score": None,
                "fraction_disordered": None, "has_clash_af3": None,
                "chain_pair_pae_min": None,
            }

            if rfd3_file is None or af3_cif is None:
                n_missing += 1
                rows.append(base_row)
                continue

            metrics = compute_pair_metrics(
                rfd3_cif         = rfd3_file,
                af3_cif          = af3_cif,
                conf_json        = conf_json,
                summary_json     = summary_json,
                rfd3_epitope_ris = rfd3_epi_ris,
                true_pdb         = true_pdb,
                true_epitope_ris = true_epi_ris,
                clash_cutoff     = args.clash_cutoff,
            )
            base_row.update(metrics)

            if metrics["status"] == "ok":
                n_ok += 1
            else:
                n_seqfail += 1

            rows.append(base_row)

    out = pd.DataFrame(rows)
    out_csv = Path(args.out_csv).resolve()
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(out_csv, index=False)

    # Apply Lawson's filters and print summary
    passing = out.loc[
        out["mean_pae"].notna() &
        (out["mean_pae"] < 5.0) &
        out["overall_rmsd"].notna() &
        (out["overall_rmsd"] <= 2.0) &
        out["epitope_chunk_rmsd"].notna() &
        (out["epitope_chunk_rmsd"] <= 1.0) &
        (out["af3_n_clash_res"].fillna(1) == 0)
    ]

    print(f"\n[metrics] Wrote {len(out):,} rows -> {out_csv}")
    print(f"[metrics] ok={n_ok:,}  missing_pair={n_missing:,}  other_fail={n_seqfail:,}")
    print(f"\n=== FILTER SUMMARY (Lawson thresholds) ===")
    print(f"  mean_pae < 5       : {(out['mean_pae'] < 5.0).sum():,}")
    print(f"  overall_rmsd ≤ 2Å  : {(out['overall_rmsd'] <= 2.0).sum():,}")
    print(f"  epitope_rmsd ≤ 1Å  : {(out['epitope_chunk_rmsd'] <= 1.0).sum():,}")
    print(f"  no clashes         : {(out['af3_n_clash_res'].fillna(1) == 0).sum():,}")
    print(f"  ALL FOUR PASS      : {len(passing):,}  ({100*len(passing)/max(len(out),1):.1f}%)")

    if len(passing) > 0:
        print("\nTop 20 passing designs:")
        print(
            passing.sort_values("epitope_chunk_rmsd")
            [["id", "pred", "overall_rmsd", "epitope_chunk_rmsd", "mean_pae", "af3_n_clash_res"]]
            .head(20).to_string(index=False)
        )


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main() -> None:
    top = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    sub = top.add_subparsers(dest="mode", required=True)

    # --- validate subcommand ---
    p_val = sub.add_parser("validate", help="Reproduce Lawson metrics and diff against dp2")
    p_val.add_argument("--dp2_parquet",    required=True)
    p_val.add_argument("--lawson_root",    required=True,
                       help="Path to sourced_antibody_v1/no_antibody (contains proteinmpnn/ and af3_predictions/)")
    p_val.add_argument("--out_csv",        required=True)
    p_val.add_argument("--sample",         type=int, default=0,
                       help="Random sample N rows (0=all)")
    p_val.add_argument("--n",              type=int, default=0,
                       help="Take first N rows (0=all, applied after --sample)")
    p_val.add_argument("--progress_every", type=int, default=200)

    # --- run subcommand ---
    p_run = sub.add_parser("run", help="Compute metrics for your RFD3 run")
    p_run.add_argument("--run_dir",        required=True)
    p_run.add_argument("--dp2_parquet",    required=True)
    p_run.add_argument("--true_dir",       required=True,
                       help="Directory with true complex PDBs named <id>.pdb")
    p_run.add_argument("--out_csv",        required=True)
    p_run.add_argument("--pred_max",       type=int, default=7,
                       help="Max pred index (default 7 → preds 0..7 = 8 total)")
    p_run.add_argument("--clash_cutoff",   type=float, default=4.0)
    p_run.add_argument("--limit",          type=int, default=0,
                       help="Limit dp2 rows for testing (0=all)")

    args = top.parse_args()

    if args.mode == "validate":
        run_validate(args)
    elif args.mode == "run":
        run_metrics(args)


if __name__ == "__main__":
    main()
