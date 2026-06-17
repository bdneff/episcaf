#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Optional

import numpy as np
import pandas as pd

# You likely already use MDAnalysis; it's good for RMSD + selections
import MDAnalysis as mda


def parse_index_list(x) -> list[int]:
    """Parse list-like epitope indices from parquet columns.
    Handles python lists, strings like "[1,2,3]" or "1,2,3" or "1 2 3".
    """
    if x is None or (isinstance(x, float) and math.isnan(x)):
        return []
    if isinstance(x, (list, tuple, np.ndarray)):
        return [int(i) for i in x]
    s = str(x).strip()
    if not s:
        return []
    # strip brackets
    s = s.replace("[", "").replace("]", "").replace("(", "").replace(")", "")
    # split by comma or whitespace
    parts = [p for p in s.replace(",", " ").split() if p]
    out = []
    for p in parts:
        try:
            out.append(int(p))
        except ValueError:
            continue
    return out


def pick_epitope_indices(row: dict) -> list[int]:
    for c in ("assay_scaffolded_epitope_resindices", "scaffolded_epitope_resindices", "epitope_resindices"):
        if c in row and row[c] is not None:
            idx = parse_index_list(row[c])
            if idx:
                return idx
    return []


def find_first(globs: Iterable[Path]) -> Optional[Path]:
    for p in globs:
        if p.exists():
            return p
    return None


def find_af3_structure(out_dir: Path) -> Optional[Path]:
    # common AF3 output structure filenames
    candidates = []
    candidates += sorted(out_dir.glob("**/*model*.cif"))
    candidates += sorted(out_dir.glob("**/*ranked_0*.cif"))
    candidates += sorted(out_dir.glob("**/*.cif"))
    candidates += sorted(out_dir.glob("**/*.pdb"))
    return candidates[0] if candidates else None


def find_af3_pae_json(out_dir: Path) -> Optional[Path]:
    # common AF3 confidence summary locations
    candidates = []
    candidates += sorted(out_dir.glob("**/*summary*.json"))
    candidates += sorted(out_dir.glob("**/*confidence*.json"))
    candidates += sorted(out_dir.glob("**/*ranking*.json"))
    candidates += sorted(out_dir.glob("**/*.json"))
    # prefer files containing "pae" text
    for p in candidates:
        try:
            txt = p.read_text()
            if "pae" in txt.lower():
                return p
        except Exception:
            continue
    return candidates[0] if candidates else None


def extract_mean_pae(pae_json_path: Path) -> Optional[float]:
    """Try to extract mean PAE from AF3 outputs. Returns None if not found."""
    try:
        d = json.loads(pae_json_path.read_text())
    except Exception:
        return None

    # Heuristics: look for arrays named like "pae", "predicted_aligned_error", etc.
    def walk(obj):
        if isinstance(obj, dict):
            for k,v in obj.items():
                lk = str(k).lower()
                if lk in ("pae", "predicted_aligned_error", "predicted_aligned_error_matrix"):
                    return v
                out = walk(v)
                if out is not None:
                    return out
        elif isinstance(obj, list):
            for it in obj:
                out = walk(it)
                if out is not None:
                    return out
        return None

    pae = walk(d)
    if pae is None:
        return None

    arr = np.array(pae, dtype=float)
    if arr.ndim == 0:
        return float(arr)
    return float(np.nanmean(arr))


def kabsch_rmsd(P: np.ndarray, Q: np.ndarray) -> float:
    """RMSD after optimal superposition (Kabsch). P,Q shape (N,3)."""
    if P.shape != Q.shape or P.shape[0] < 3:
        return float("nan")
    Pc = P - P.mean(axis=0)
    Qc = Q - Q.mean(axis=0)
    C = Pc.T @ Qc
    V, S, Wt = np.linalg.svd(C)
    d = np.sign(np.linalg.det(V @ Wt))
    D = np.diag([1,1,d])
    U = V @ D @ Wt
    P_rot = Pc @ U
    diff = P_rot - Qc
    return float(np.sqrt((diff*diff).sum() / P.shape[0]))


def get_ca_coords(u: mda.Universe, chain_id: str = "A") -> np.ndarray:
    sel = u.select_atoms(f"protein and segid {chain_id} and name CA")
    if sel.n_atoms == 0:
        # some files use chainID in 'chainID' not segid; MDAnalysis often maps to segid.
        sel = u.select_atoms("protein and name CA")
    return sel.positions.copy()


def get_ca_coords_for_resids(u: mda.Universe, resids_1based: list[int], chain_id: str = "A") -> np.ndarray:
    # MDAnalysis selection uses resid numbers; we assume AF3/RFD3 structures use standard 1-based resid numbering.
    resid_str = " ".join(str(r) for r in resids_1based)
    sel = u.select_atoms(f"protein and segid {chain_id} and name CA and resid {resid_str}")
    if sel.n_atoms == 0:
        sel = u.select_atoms(f"protein and name CA and resid {resid_str}")
    return sel.positions.copy()


@dataclass
class RowOut:
    design_id: str
    id: str
    contig_id: int
    rfd3_struct: str | None
    af3_struct: str | None
    pae_json: str | None
    rmsd_ca_all: float | None
    rmsd_ca_epitope: float | None
    pae_mean_all: float | None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--run_dir", required=True)
    ap.add_argument("--contigs_parquet", default=None)
    ap.add_argument("--out_parquet", default=None)
    args = ap.parse_args()

    run_dir = Path(args.run_dir).resolve()
    contigs_pq = Path(args.contigs_parquet).resolve() if args.contigs_parquet else (run_dir / "01_design" / "contigs.parquet")
    af3_out_root = run_dir / "03_af3" / "outputs"

    out_pq = Path(args.out_parquet).resolve() if args.out_parquet else (run_dir / "04_filter" / "metrics.parquet")
    out_pq.parent.mkdir(parents=True, exist_ok=True)

    df = pd.read_parquet(contigs_pq)
    df = df.set_index("design_id", drop=False)

    rows: list[RowOut] = []
    missing = 0

    # iterate over AF3 output dirs (only score those that exist)
    out_dirs = sorted([p for p in af3_out_root.iterdir() if p.is_dir()])
    print(f"[score] Found {len(out_dirs)} AF3 output dirs under {af3_out_root}")

    for out_dir in out_dirs:
        design_id = out_dir.name
        if design_id not in df.index:
            continue
        r = df.loc[design_id].to_dict()
        ep_idx0 = pick_epitope_indices(r)   # 0-based
        ep_resids = [i+1 for i in ep_idx0]  # -> 1-based

        af3_struct = find_af3_structure(out_dir)
        pae_json = find_af3_pae_json(out_dir)
        pae_mean = extract_mean_pae(pae_json) if pae_json else None

        # Find RFD3 structure corresponding to design_id.
        # You may need to tweak these patterns to match your outputs.
        rfd3_candidates = []
        rfd3_candidates += list((run_dir / "02_rfd3" / "outputs").glob(f"**/{design_id}*.cif"))
        rfd3_candidates += list((run_dir / "02_rfd3" / "outputs").glob(f"**/{design_id}*.pdb"))
        rfd3_candidates += list((run_dir / "02_rfd3" / "outputs").glob(f"**/{design_id}*/**/*.cif"))
        rfd3_candidates += list((run_dir / "02_rfd3" / "outputs").glob(f"**/{design_id}*/**/*.pdb"))
        rfd3_candidates = sorted(set([p for p in rfd3_candidates if p.is_file()]))
        rfd3_struct = rfd3_candidates[0] if rfd3_candidates else None

        if not af3_struct or not rfd3_struct:
            missing += 1
            rows.append(RowOut(
                design_id=design_id,
                id=str(r.get("id")),
                contig_id=int(r.get("contig_id", -1)),
                rfd3_struct=str(rfd3_struct) if rfd3_struct else None,
                af3_struct=str(af3_struct) if af3_struct else None,
                pae_json=str(pae_json) if pae_json else None,
                rmsd_ca_all=None,
                rmsd_ca_epitope=None,
                pae_mean_all=pae_mean,
            ))
            continue

        # Load universes
        u_rfd = mda.Universe(str(rfd3_struct))
        u_af3 = mda.Universe(str(af3_struct))

        # RMSD on all CA (superpose)
        P = get_ca_coords(u_rfd)
        Q = get_ca_coords(u_af3)
        rmsd_all = kabsch_rmsd(P, Q)

        # RMSD on epitope CA (using same superposition is typical; here we compute epitope RMSD after its own Kabsch)
        # If you want epitope RMSD after global superposition, we can do that too.
        Pe = get_ca_coords_for_resids(u_rfd, ep_resids)
        Qe = get_ca_coords_for_resids(u_af3, ep_resids)
        rmsd_epi = kabsch_rmsd(Pe, Qe) if (Pe.size and Qe.size and Pe.shape == Qe.shape) else float("nan")

        rows.append(RowOut(
            design_id=design_id,
            id=str(r.get("id")),
            contig_id=int(r.get("contig_id", -1)),
            rfd3_struct=str(rfd3_struct),
            af3_struct=str(af3_struct),
            pae_json=str(pae_json) if pae_json else None,
            rmsd_ca_all=rmsd_all,
            rmsd_ca_epitope=rmsd_epi,
            pae_mean_all=pae_mean,
        ))

    out_df = pd.DataFrame([r.__dict__ for r in rows])
    out_df.to_parquet(out_pq, index=False)
    print(f"[score] wrote {len(out_df)} rows -> {out_pq}")
    print(f"[score] missing structures for {missing} designs (see null paths)")

if __name__ == "__main__":
    main()
