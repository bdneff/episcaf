#!/usr/bin/env python3
from __future__ import annotations

import argparse, gzip, json, math, re
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
import gemmi

import MDAnalysis as mda
from MDAnalysis.analysis import rms
from MDAnalysis.lib.distances import distance_array


AA3TO1 = {
 "ALA":"A","CYS":"C","ASP":"D","GLU":"E","PHE":"F","GLY":"G","HIS":"H","ILE":"I","LYS":"K","LEU":"L",
 "MET":"M","ASN":"N","PRO":"P","GLN":"Q","ARG":"R","SER":"S","THR":"T","VAL":"V","TRP":"W","TYR":"Y"
}

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

def read_gemmi_structure(p: Path) -> gemmi.Structure:
    s = str(p)
    if s.endswith(".gz"):
        with gzip.open(p, "rt") as f:
            doc = gemmi.cif.read_string(f.read())
        return gemmi.make_structure_from_block(doc.sole_block())
    if p.suffix.lower() == ".cif":
        doc = gemmi.cif.read(str(p))
        return gemmi.make_structure_from_block(doc.sole_block())
    return gemmi.read_structure(str(p))

def get_chainA(st: gemmi.Structure) -> gemmi.Chain:
    m = st[0]
    for ch in m:
        if ch.name == "A":
            return ch
    return m[0]

def chain_seq_1letter(chain: gemmi.Chain) -> str:
    return "".join(AA3TO1.get(r.name.upper(), "X") for r in chain)

def find_subseq_allowX(hay: str, needle: str) -> int:
    i = hay.find(needle)
    if i >= 0:
        return i
    H, N = hay, needle
    Lh, Ln = len(H), len(N)
    for i in range(Lh - Ln + 1):
        ok = True
        for a, b in zip(H[i:i+Ln], N):
            if a != "X" and a != b:
                ok = False
                break
        if ok:
            return i
    return -1

def gemmi_bb_positions(chain: gemmi.Chain, start: int, end: int) -> np.ndarray:
    res = [r for r in chain][start:end]
    coords = []
    for r in res:
        for name in ("N","CA","C","O"):
            a = r.find_atom(name, altloc="*")
            if a:
                p = a.pos
                coords.append([p.x,p.y,p.z])
    return np.array(coords, float)

def gemmi_bb_positions_for_resindices(chain: gemmi.Chain, ris: List[int]) -> np.ndarray:
    res = [r for r in chain]
    coords = []
    for i in ris:
        r = res[i]
        for name in ("N","CA","C","O"):
            a = r.find_atom(name, altloc="*")
            if a:
                p = a.pos
                coords.append([p.x,p.y,p.z])
    return np.array(coords, float)

def rmsd_superpose(P: np.ndarray, Q: np.ndarray) -> float:
    n = min(len(P), len(Q))
    if n < 3:
        return float("nan")
    return float(rms.rmsd(P[:n], Q[:n], superposition=True))

# --- clash helpers ---

def kabsch_fit(P: np.ndarray, Q: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    Pc = P - P.mean(axis=0)
    Qc = Q - Q.mean(axis=0)
    C = Pc.T @ Qc
    V, _, Wt = np.linalg.svd(C)
    d = np.sign(np.linalg.det(V @ Wt))
    D = np.diag([1.0, 1.0, d])
    R = V @ D @ Wt
    t = Q.mean(axis=0) - (P.mean(axis=0) @ R)
    return R, t

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

def gemmi_res_ca(res: gemmi.Residue) -> Optional[np.ndarray]:
    a = res.find_atom("CA", altloc="*")
    if not a:
        return None
    p = a.pos
    return np.array([p.x, p.y, p.z], dtype=float)

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
    af3_chunk_resindices: List[int],
    true_chunk_resindices: List[int],
    cutoff: float = 4.0,
) -> Optional[List[int]]:
    true_u = mda.Universe(str(true_pdb))
    true_ag = sel_chain_or_segid(true_u, ["A"])
    if len(true_ag) == 0:
        return None
    true_ag_res = true_ag.residues

    ab_atoms = pick_antibody_atoms(true_u)
    if len(ab_atoms) == 0:
        return None
    ab_pos = ab_atoms.positions

    st = read_gemmi_structure(af3_cif)
    chA = get_chainA(st)
    af3_res = [r for r in chA]

    if not af3_chunk_resindices or not true_chunk_resindices:
        return None
    if max(af3_chunk_resindices) >= len(af3_res) or max(true_chunk_resindices) >= len(true_ag_res):
        return None

    P_list = []
    Q_list = []
    for i_af3, i_true in zip(af3_chunk_resindices, true_chunk_resindices):
        caP = gemmi_res_ca(af3_res[i_af3])
        caQ_ag = true_ag_res[i_true].atoms.select_atoms("name CA and not name H*")
        if caP is None or len(caQ_ag) != 1:
            continue
        P_list.append(caP)
        Q_list.append(caQ_ag.positions[0])

    if len(P_list) < 3:
        return None

    P = np.vstack(P_list)
    Q = np.vstack(Q_list)
    R, t = kabsch_fit(P, Q)

    mask = np.zeros(len(af3_res), dtype=bool)
    mask[np.array(af3_chunk_resindices, dtype=int)] = True
    unintended_idx = np.where(~mask)[0]
    if unintended_idx.size == 0:
        return []

    clashing = []
    for i in unintended_idx:
        heavy = gemmi_res_heavy_coords(af3_res[i])
        if heavy.shape[0] == 0:
            continue
        heavy_fit = (heavy @ R) + t
        d = distance_array(heavy_fit, ab_pos)
        if np.any(d < cutoff):
            clashing.append(int(i))
    return clashing

# --- indexing run outputs ---

def find_af3_files_in_dir(d: Path) -> Tuple[Optional[Path], Optional[Path]]:
    af3_cif = next(iter(d.glob("*_model.cif")), None)
    if af3_cif is None:
        af3_cif = next(iter(d.glob("*_model.cif.gz")), None)
    if af3_cif is None:
        af3_cif = next(iter(d.rglob("model.cif")), None)
    if af3_cif is None:
        af3_cif = next(iter(d.rglob("model.cif.gz")), None)

    conf = next(iter(d.glob("*_confidences.json")), None)
    if conf is None:
        conf = next(iter(d.rglob("confidences.json")), None)

    return af3_cif, conf

def index_af3_outputs(af3_root: Path) -> Dict[tuple[str,int], Path]:
    idx: Dict[tuple[str,int], Path] = {}
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
        idx[(tok, pred)] = d
    return idx

def index_rfd3_outputs(rfd_root: Path) -> Dict[tuple[str,int], Path]:
    idx: Dict[tuple[str,int], Path] = {}
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
        idx[(tok, pred)] = p
    return idx

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--run_dir", required=True)
    ap.add_argument("--dp2_parquet", required=True)
    ap.add_argument("--true_dir", required=True)
    ap.add_argument("--out_csv", required=True)
    ap.add_argument("--pred_max", type=int, default=7)
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--clash_cutoff", type=float, default=4.0)
    args = ap.parse_args()

    run_dir = Path(args.run_dir).resolve()
    dp2 = pd.read_parquet(Path(args.dp2_parquet).resolve())

    if "assay_scaffolded_epitope_id" not in dp2.columns:
        raise SystemExit("dp2.parquet must contain assay_scaffolded_epitope_id")
    for col in ("id","assay_scaffolded_epitope_chunk_resindices","epitope_chunk_resindices"):
        if col not in dp2.columns:
            raise SystemExit(f"dp2.parquet missing required column: {col}")

    dp2["assay_scaffolded_epitope_id"] = dp2["assay_scaffolded_epitope_id"].astype(str).str.lower()

    base = dp2.copy()
    if args.limit and args.limit > 0:
        base = base.head(args.limit).copy()

    af3_root = run_dir / "03_af3" / "outputs"
    rfd_root = run_dir / "02_rfd3" / "outputs"
    af3_index = index_af3_outputs(af3_root)
    rfd3_index = index_rfd3_outputs(rfd_root)

    print(f"[metrics] base rows: {len(base):,}")
    print(f"[metrics] AF3 index entries: {len(af3_index):,}")
    print(f"[metrics] RFD3 index entries: {len(rfd3_index):,}")

    true_dir = Path(args.true_dir).resolve()

    rows = []
    missing_pairs = 0
    seqmatch_fail = 0
    epi_fail = 0
    clash_fail = 0

    for _, r in base.iterrows():
        tok = str(r["assay_scaffolded_epitope_id"]).lower()
        pid = r.get("id")

        af3_chunk = parse_index_list(r.get("assay_scaffolded_epitope_chunk_resindices"))
        rfd_chunk = parse_index_list(r.get("epitope_chunk_resindices"))
        true_chunk = parse_index_list(r.get("epitope_chunk_resindices"))

        for pred in range(args.pred_max + 1):
            af3_dir = af3_index.get((tok, pred))
            rfd3_file = rfd3_index.get((tok, pred))

            af3_cif = conf_json = None
            if af3_dir:
                af3_cif, conf_json = find_af3_files_in_dir(af3_dir)

            row = {
                "assay_scaffolded_epitope_id": tok,
                "pred": pred,
                "id": pid,
                "rfd3_path": str(rfd3_file) if rfd3_file else None,
                "af3_dir": str(af3_dir) if af3_dir else None,
                "af3_path": str(af3_cif) if af3_cif else None,
                "af3_conf_path": str(conf_json) if conf_json else None,
                "mean_pae": pae_mean_from_conf(conf_json) if conf_json else None,
                "overall_rmsd": None,
                "epitope_chunk_rmsd_vs_rfd3": None,
                "af3_window_start": None,
                "af3_window_end": None,
                "af3_clash_resindices": None,
                "af3_has_clash": None,
                "af3_n_clash_res": None,
            }

            if (rfd3_file is None) or (af3_cif is None) or (pid is None):
                missing_pairs += 1
                rows.append(row)
                continue

            # --- overall_rmsd (Lawson-style window) ---
            try:
                st_rfd = read_gemmi_structure(Path(rfd3_file))
                st_af3 = read_gemmi_structure(Path(af3_cif))
                ch_rfd = get_chainA(st_rfd)
                ch_af3 = get_chainA(st_af3)

                rfd_seq = chain_seq_1letter(ch_rfd)
                af3_seq = chain_seq_1letter(ch_af3)

                start = find_subseq_allowX(af3_seq, rfd_seq)
                if start < 0:
                    seqmatch_fail += 1
                else:
                    end = start + len(rfd_seq)
                    P = gemmi_bb_positions(ch_rfd, 0, len([rr for rr in ch_rfd]))
                    Q = gemmi_bb_positions(ch_af3, start, end)
                    row["overall_rmsd"] = rmsd_superpose(P, Q)
                    row["af3_window_start"] = int(start)
                    row["af3_window_end"] = int(end)
            except Exception:
                seqmatch_fail += 1

            # --- epitope chunk rmsd (mapped indices) ---
            try:
                if len(rfd_chunk) > 0 and len(af3_chunk) > 0:
                    st_rfd = read_gemmi_structure(Path(rfd3_file))
                    st_af3 = read_gemmi_structure(Path(af3_cif))
                    ch_rfd = get_chainA(st_rfd)
                    ch_af3 = get_chainA(st_af3)
                    P = gemmi_bb_positions_for_resindices(ch_rfd, rfd_chunk)
                    Q = gemmi_bb_positions_for_resindices(ch_af3, af3_chunk)
                    row["epitope_chunk_rmsd_vs_rfd3"] = rmsd_superpose(P, Q)
            except Exception:
                epi_fail += 1

            # --- clashes (AF3 unintended vs true antibody) ---
            try:
                true_pdb = true_dir / f"{pid}.pdb"
                if true_pdb.exists() and len(af3_chunk) >= 3 and len(true_chunk) >= 3:
                    cl = compute_af3_clash_resindices(Path(af3_cif), true_pdb, af3_chunk, true_chunk, cutoff=args.clash_cutoff)
                    row["af3_clash_resindices"] = cl
                    if cl is not None:
                        row["af3_has_clash"] = bool(len(cl) > 0)
                        row["af3_n_clash_res"] = int(len(cl))
            except Exception:
                clash_fail += 1

            rows.append(row)

    out = pd.DataFrame(rows)
    out_csv = Path(args.out_csv).resolve()
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(out_csv, index=False)

    print(f"[metrics] wrote rows: {len(out):,} -> {out_csv}")
    print(f"[metrics] missing pairs: {missing_pairs:,}")
    print(f"[metrics] seqmatch_fail: {seqmatch_fail:,}")
    print(f"[metrics] epi_fail: {epi_fail:,}")
    print(f"[metrics] clash_fail: {clash_fail:,}")

if __name__ == "__main__":
    main()
