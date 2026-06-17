#!/usr/bin/env python3
"""
compute_decomposed_metrics.py

Extends the existing metrics CSV with decomposed epitope/scaffold metrics:

  1. epitope_pae     : mean PAE restricted to epitope residue indices
  2. scaffold_pae    : mean PAE restricted to non-epitope residue indices
  3. scaffold_rmsd   : backbone RMSD over non-epitope residues only,
                       after aligning on non-epitope residues

Together with the existing epitope_chunk_rmsd and mean_pae, these four metrics
allow the four experimental conditions from the filter decomposition proposal:

  Condition 1 (rigid/rigid)       : strict epitope RMSD + strict scaffold RMSD
  Condition 2 (rigid epi/flex sc) : strict epitope RMSD + relaxed scaffold RMSD
  Condition 3 (flex epi/rigid sc) : relaxed epitope RMSD + strict scaffold RMSD
  Condition 4 (flex/flex)         : relaxed on both (negative control)

PAE decomposition uses the full N x N PAE matrix from AF3 confidences.json.
The epitope PAE is the mean of the submatrix rows AND columns corresponding to
epitope residue indices (i.e. inter-epitope PAE). Scaffold PAE is the mean of
the remaining submatrix.

Usage:
    python scripts/compute_decomposed_metrics.py \
        --metrics_csv  runs/run_rfd3_mpnn/04_filter/metrics_partial.csv \
        --dp2_parquet  datasets/dp2.parquet \
        --mpnn_pdb_dir runs/run_rfd3_mpnn/02_mpnn_pdbs \
        --af3_out_dir  runs/run_rfd3_mpnn/03_af3/outputs \
        --out_csv      runs/run_rfd3_mpnn/04_filter/metrics_decomposed.csv \
        --limit        500
"""

import argparse
import gzip
import json
import math
import os
import re
import sys
import tempfile
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import gemmi
import MDAnalysis as mda
from MDAnalysis.analysis import rms as mda_rms
import numpy as np
import pandas as pd


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

AA3TO1 = {
    "ALA": "A", "CYS": "C", "ASP": "D", "GLU": "E", "PHE": "F",
    "GLY": "G", "HIS": "H", "ILE": "I", "LYS": "K", "LEU": "L",
    "MET": "M", "ASN": "N", "PRO": "P", "GLN": "Q", "ARG": "R",
    "SER": "S", "THR": "T", "VAL": "V", "TRP": "W", "TYR": "Y",
}

_PRED_ID_RE = re.compile(
    r"^([0-9a-fA-F]{32})_pred(\d+)_fixed_dldesign_(\d+)$",
    re.IGNORECASE,
)


# ---------------------------------------------------------------------------
# Structure helpers
# ---------------------------------------------------------------------------

def cif_to_tmp_pdb(cif_path: Path) -> str:
    if str(cif_path).endswith(".gz"):
        with gzip.open(cif_path, "rt") as f:
            doc = gemmi.cif.read_string(f.read())
    else:
        doc = gemmi.cif.read(str(cif_path))
    st = gemmi.make_structure_from_block(doc.sole_block())
    tmp = tempfile.NamedTemporaryFile(suffix=".pdb", delete=False)
    st.write_pdb(tmp.name)
    tmp.close()
    return tmp.name


def chain_seq_from_cif(cif_path: Path) -> str:
    if str(cif_path).endswith(".gz"):
        with gzip.open(cif_path, "rt") as f:
            doc = gemmi.cif.read_string(f.read())
    else:
        doc = gemmi.cif.read(str(cif_path))
    st = gemmi.make_structure_from_block(doc.sole_block())
    ch = st[0][0]
    return "".join(AA3TO1.get(r.name.upper(), "X") for r in ch)


def find_subseq(haystack: str, needle: str) -> int:
    i = haystack.find(needle)
    if i >= 0:
        return i
    for i in range(len(haystack) - len(needle) + 1):
        if all(h == "X" or h == n
               for h, n in zip(haystack[i:i + len(needle)], needle)):
            return i
    return -1


def find_af3_files(af3_dir: Path) -> Tuple[Optional[Path], Optional[Path]]:
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
    return cif, conf


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
# Decomposed PAE
# ---------------------------------------------------------------------------

def decomposed_pae(conf_json: Path, epitope_ris: List[int], window_start: int) -> Dict[str, Optional[float]]:
    """
    Compute epitope and scaffold PAE from the full PAE matrix.

    epitope_ris are 0-based indices into the MPNN sequence.
    window_start is the offset of the MPNN sequence within the AF3 output.
    So AF3 epitope indices = window_start + epitope_ris.

    The PAE matrix from AF3 confidences.json is indexed over all AF3 residues.
    epitope_pae  = mean of PAE submatrix [epi_rows, epi_cols]
    scaffold_pae = mean of PAE submatrix [scaf_rows, scaf_cols]
    """
    result = {'epitope_pae': None, 'scaffold_pae': None}

    try:
        d = json.loads(conf_json.read_text())
    except Exception:
        return result

    pae_matrix = None
    for key in ("pae", "predicted_aligned_error", "predicted_aligned_error_matrix"):
        if key in d:
            try:
                pae_matrix = np.array(d[key], dtype=float)
                break
            except Exception:
                pass

    if pae_matrix is None or pae_matrix.ndim != 2:
        return result

    n = pae_matrix.shape[0]

    # AF3 indices for epitope
    af3_epi_ris = [window_start + i for i in epitope_ris]
    af3_epi_ris = [i for i in af3_epi_ris if i < n]

    if not af3_epi_ris:
        return result

    # scaffold = all AF3 residues not in epitope
    all_ris = np.arange(n)
    epi_mask = np.zeros(n, dtype=bool)
    epi_mask[af3_epi_ris] = True
    scaf_ris = list(np.where(~epi_mask)[0])

    if af3_epi_ris:
        epi_sub = pae_matrix[np.ix_(af3_epi_ris, af3_epi_ris)]
        result['epitope_pae'] = float(np.nanmean(epi_sub))

    if scaf_ris:
        scaf_sub = pae_matrix[np.ix_(scaf_ris, scaf_ris)]
        result['scaffold_pae'] = float(np.nanmean(scaf_sub))

    return result


# ---------------------------------------------------------------------------
# Decomposed RMSD
# ---------------------------------------------------------------------------

def scaffold_rmsd(mpnn_pdb: Path, af3_cif: Path,
                  mpnn_epitope_ris: List[int], window_start: int) -> Optional[float]:
    """
    Compute scaffold RMSD:
      1. Align MPNN non-epitope backbone onto AF3 non-epitope backbone
      2. Compute RMSD on the aligned non-epitope residues

    mpnn_epitope_ris: 0-based epitope indices into MPNN chain A
    window_start: offset of MPNN sequence in AF3 output
    """
    tmp_af3 = None
    try:
        tmp_af3 = cif_to_tmp_pdb(af3_cif)
        af3_u  = mda.Universe(tmp_af3)
        mpnn_u = mda.Universe(str(mpnn_pdb))
    except Exception:
        return None
    finally:
        if tmp_af3 and os.path.exists(tmp_af3):
            os.unlink(tmp_af3)

    mpnn_selA = mpnn_u.select_atoms("segid A")
    if len(mpnn_selA) == 0:
        mpnn_selA = mpnn_u.select_atoms("chainid A")

    n_mpnn = len(mpnn_selA.residues)
    af3_window_end = window_start + n_mpnn

    # build scaffold residue index lists (non-epitope)
    all_mpnn_ris = list(range(n_mpnn))
    epi_set = set(mpnn_epitope_ris)
    scaf_mpnn_ris = [i for i in all_mpnn_ris if i not in epi_set]
    scaf_af3_ris  = [window_start + i for i in scaf_mpnn_ris]

    if not scaf_mpnn_ris:
        return None

    try:
        P = mpnn_selA.residues[scaf_mpnn_ris].atoms.select_atoms("backbone").positions
        Q = af3_u.residues[scaf_af3_ris].atoms.select_atoms("backbone").positions
        if len(P) < 3 or len(Q) < 3:
            return None
        return float(mda_rms.rmsd(P, Q, superposition=True))
    except Exception:
        return None


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--metrics_csv",  required=True,
                        help="Existing metrics CSV from compute_metrics_mpnn.py")
    parser.add_argument("--dp2_parquet",  required=True)
    parser.add_argument("--mpnn_pdb_dir", required=True)
    parser.add_argument("--af3_out_dir",  required=True)
    parser.add_argument("--out_csv",      required=True)
    parser.add_argument("--limit",        type=int, default=0,
                        help="Limit to first N rows for testing (0=all)")
    args = parser.parse_args()

    print("Loading metrics CSV ...")
    df = pd.read_csv(args.metrics_csv, low_memory=False)
    if args.limit and args.limit > 0:
        df = df.head(args.limit).copy()
    print(f"  {len(df):,} rows")

    print("Loading dp2 ...")
    dp2 = pd.read_parquet(args.dp2_parquet)
    dp2["assay_scaffolded_epitope_id"] = dp2["assay_scaffolded_epitope_id"].astype(str).str.lower()
    dp2_tok = dp2.drop_duplicates("assay_scaffolded_epitope_id").set_index("assay_scaffolded_epitope_id")

    mpnn_root = Path(args.mpnn_pdb_dir)
    af3_root  = Path(args.af3_out_dir)
    out_csv   = Path(args.out_csv)

    # build AF3 output index
    print("Indexing AF3 outputs ...")
    af3_idx: Dict[Tuple, Path] = {}
    for d in af3_root.iterdir():
        if d.is_dir() and (d / "_DONE").exists():
            m = _PRED_ID_RE.match(d.name)
            if m:
                af3_idx[(m.group(1).lower(), int(m.group(2)), int(m.group(3)))] = d

    # build MPNN PDB index
    print("Indexing MPNN PDBs ...")
    mpnn_idx: Dict[Tuple, Path] = {}
    for p in mpnn_root.rglob("*_fixed_dldesign_*.pdb"):
        m = _PRED_ID_RE.match(p.stem)
        if m:
            mpnn_idx[(m.group(1).lower(), int(m.group(2)), int(m.group(3)))] = p

    # initialize new columns
    df['epitope_pae']   = np.nan
    df['scaffold_pae']  = np.nan
    df['scaffold_rmsd'] = np.nan

    n_ok = n_skip = n_fail = 0

    for i, row in df.iterrows():
        tok     = str(row['token']).lower()
        pred    = int(row['pred'])
        mpnn_id = int(row['mpnn_id'])

        # skip if AF3 didn't complete or metrics were missing
        if pd.isna(row.get('overall_rmsd')):
            n_skip += 1
            continue

        window_start = int(row['af3_window_start']) if pd.notna(row.get('af3_window_start')) else 0

        # get epitope indices
        dp2_row = dp2_tok.loc[tok] if tok in dp2_tok.index else None
        if dp2_row is None:
            n_skip += 1
            continue
        epi_ris = parse_index_list(dp2_row['scaffolded_epitope_chunk_resindices'])
        if not epi_ris:
            n_skip += 1
            continue

        # get file paths
        af3_dir  = af3_idx.get((tok, pred, mpnn_id))
        mpnn_pdb = mpnn_idx.get((tok, pred, mpnn_id))

        if af3_dir is None or mpnn_pdb is None:
            n_skip += 1
            continue

        af3_cif, conf_json = find_af3_files(af3_dir)
        if af3_cif is None:
            n_skip += 1
            continue

        # decomposed PAE
        if conf_json and conf_json.exists():
            pae_result = decomposed_pae(conf_json, epi_ris, window_start)
            df.at[i, 'epitope_pae']  = pae_result['epitope_pae']
            df.at[i, 'scaffold_pae'] = pae_result['scaffold_pae']

        # scaffold RMSD
        scaf_rmsd = scaffold_rmsd(mpnn_pdb, af3_cif, epi_ris, window_start)
        df.at[i, 'scaffold_rmsd'] = scaf_rmsd

        n_ok += 1

        if (i + 1) % 5000 == 0:
            print(f"  {i+1:,} rows processed  ok={n_ok:,}  skip={n_skip:,}  fail={n_fail:,}")

    out_csv.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_csv, index=False)
    print(f"\nWrote {len(df):,} rows -> {out_csv}")
    print(f"ok={n_ok:,}  skip={n_skip:,}  fail={n_fail:,}")

    # quick summary of new metrics
    print(f"\n=== NEW METRIC SUMMARY ===")
    for col in ['epitope_pae', 'scaffold_pae', 'scaffold_rmsd']:
        valid = df[col].dropna()
        print(f"  {col:20s}  n={len(valid):,}  mean={valid.mean():.2f}  "
              f"median={valid.median():.2f}  min={valid.min():.2f}  max={valid.max():.2f}")

    # show pass counts under each condition using your proposed thresholds
    # (you can tune these after seeing the distributions)
    print(f"\n=== CONDITION PASS COUNTS (example thresholds) ===")
    strict_epi_rmsd  = df['epitope_chunk_rmsd'] <= 1.0
    relaxed_epi_rmsd = df['epitope_chunk_rmsd'] <= 3.0
    strict_scaf_rmsd = df['scaffold_rmsd'] <= 2.0
    relaxed_scaf_rmsd= df['scaffold_rmsd'] <= 5.0
    strict_epi_pae   = df['epitope_pae'] < 5.0
    relaxed_epi_pae  = df['epitope_pae'] < 10.0
    strict_scaf_pae  = df['scaffold_pae'] < 5.0
    relaxed_scaf_pae = df['scaffold_pae'] < 10.0
    no_clash         = df['af3_n_clash_res'].fillna(1) == 0

    conditions = {
        'Cond 1 rigid/rigid':       strict_epi_rmsd  & strict_scaf_rmsd  & strict_epi_pae  & strict_scaf_pae  & no_clash,
        'Cond 2 rigid epi/flex sc': strict_epi_rmsd  & relaxed_scaf_rmsd & strict_epi_pae  & relaxed_scaf_pae & no_clash,
        'Cond 3 flex epi/rigid sc': relaxed_epi_rmsd & strict_scaf_rmsd  & relaxed_epi_pae & strict_scaf_pae  & no_clash,
        'Cond 4 flex/flex (ctrl)':  relaxed_epi_rmsd & relaxed_scaf_rmsd & relaxed_epi_pae & relaxed_scaf_pae & no_clash,
    }

    for name, mask in conditions.items():
        n = mask.sum()
        print(f"  {name:35s}  {n:6,}  ({100*n/max(len(df),1):.3f}%)")


if __name__ == "__main__":
    main()
