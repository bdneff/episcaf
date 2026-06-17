#!/usr/bin/env python3
"""
compute_sasa_retention.py

Compute epitope SASA retention for each design in the known-antibody
RFD3+MPNN dataset and add as columns to metrics_decomposed.csv.

Two versions computed:
  mpnn_sasa_retention : SASA of epitope in MPNN scaffold / native crystal SASA
  af3_sasa_retention  : SASA of epitope in AF3 prediction / native crystal SASA

This allows comparison of which better correlates with af3_n_clash_res,
validating SASA retention as a proxy clash metric for the no-antibody case.

Usage:
    # test
    python scripts/compute_sasa_retention.py \
        --metrics_csv   runs/run_rfd3_mpnn/04_filter/metrics_decomposed.csv \
        --dp2_parquet   datasets/dp2.parquet \
        --true_dir      /tgen_labs/altin/.../cleaned \
        --out_csv       runs/run_rfd3_mpnn/04_filter/metrics_sasa.csv \
        --mpnn_pdb_dir  runs/run_rfd3_mpnn/02_mpnn_pdbs \
        --limit         500

    # full run (sbatch)
    python scripts/compute_sasa_retention.py \
        --metrics_csv   runs/run_rfd3_mpnn/04_filter/metrics_decomposed.csv \
        --dp2_parquet   datasets/dp2.parquet \
        --true_dir      /tgen_labs/altin/.../cleaned \
        --out_csv       runs/run_rfd3_mpnn/04_filter/metrics_sasa.csv \
        --mpnn_pdb_dir  runs/run_rfd3_mpnn/02_mpnn_pdbs
"""

import argparse
import gzip
import math
import os
import re
import tempfile
from pathlib import Path
from typing import List, Optional

import gemmi
import MDAnalysis as mda
from MDAnalysis.analysis.align import alignto
import freesasa
import numpy as np
import pandas as pd


_PRED_ID_RE = re.compile(
    r"^([0-9a-fA-F]{32})_pred(\d+)_fixed_dldesign_(\d+)$",
    re.IGNORECASE,
)


def parse_index_list(x) -> List[int]:
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


def index_mpnn_pdbs(mpnn_root: Path) -> dict:
    idx = {}
    for p in mpnn_root.rglob("*_fixed_dldesign_*.pdb"):
        m = _PRED_ID_RE.match(p.stem)
        if m:
            key = (m.group(1).lower(), int(m.group(2)), int(m.group(3)))
            idx[key] = p
    return idx


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


_SASA_TMP = '/tmp/sasa_tmp.pdb'

def compute_sasa(atoms) -> float:
    """Compute SASA using freesasa via a fixed temp PDB path."""
    atoms.write(_SASA_TMP)
    structure = freesasa.Structure(_SASA_TMP)
    result    = freesasa.calc(structure)
    return result.totalArea()


def get_crystal_sasa(true_pdb: Path, true_epi_ris: List[int]) -> Optional[float]:
    """Compute native epitope SASA from crystal structure (cached by caller)."""
    try:
        u = mda.Universe(str(true_pdb))
        true_res = u.select_atoms("segid A and protein").residues
        if len(true_res) == 0:
            true_res = u.select_atoms("chainid A and protein").residues
        if max(true_epi_ris) >= len(true_res):
            return None
        epi_atoms = true_res[true_epi_ris].atoms.select_atoms("not name H*")
        return compute_sasa(epi_atoms)
    except Exception:
        return None


def sasa_retention_from_pdb(
    scaffold_pdb:  Path,
    true_pdb:      Path,
    scaffold_epi_ris: List[int],
    true_epi_ris:     List[int],
    crystal_sasa:     float,
    chain:            str = "A",
    is_tmp:           bool = False,
) -> Optional[float]:
    """
    Align scaffold epitope onto crystal epitope, compute epitope SASA
    in context of scaffold, return retention ratio.
    """
    try:
        scaffold_u = mda.Universe(str(scaffold_pdb))
        true_u     = mda.Universe(str(true_pdb))

        scaf_sel = scaffold_u.select_atoms(f"segid {chain}")
        if len(scaf_sel) == 0:
            scaf_sel = scaffold_u.select_atoms(f"chainid {chain}")
        scaf_res = scaf_sel.residues

        true_res = true_u.select_atoms("segid A and protein").residues
        if len(true_res) == 0:
            true_res = true_u.select_atoms("chainid A and protein").residues

        if max(scaffold_epi_ris) >= len(scaf_res):
            return None
        if max(true_epi_ris) >= len(true_res):
            return None

        # align scaffold epitope CA onto crystal epitope CA
        scaf_epi_ca = scaf_res[scaffold_epi_ris].atoms.select_atoms("name CA")
        true_epi_ca = true_res[true_epi_ris].atoms.select_atoms("name CA")

        if len(scaf_epi_ca) < 3 or len(scaf_epi_ca) != len(true_epi_ca):
            return None

        alignto(scaf_epi_ca, true_epi_ca)

        # SASA of epitope in context of full scaffold
        scaf_epi_atoms = scaf_res[scaffold_epi_ris].atoms.select_atoms("not name H*")
        scaffold_sasa  = compute_sasa(scaf_epi_atoms)

        if crystal_sasa < 1e-6:
            return None

        return scaffold_sasa / crystal_sasa

    except Exception:
        return None
    finally:
        if is_tmp and scaffold_pdb.exists():
            os.unlink(str(scaffold_pdb))


def find_af3_cif(af3_dir: Path) -> Optional[Path]:
    cif = next(af3_dir.glob("*_model.cif"), None)
    return cif


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--metrics_csv",  required=True)
    parser.add_argument("--dp2_parquet",  required=True)
    parser.add_argument("--true_dir",     required=True)
    parser.add_argument("--out_csv",      required=True)
    parser.add_argument("--mpnn_pdb_dir", default=None)
    parser.add_argument("--limit",        type=int, default=0)
    args = parser.parse_args()

    print("Loading metrics CSV ...")
    df = pd.read_csv(args.metrics_csv, low_memory=False)
    if args.limit > 0:
        df = df.head(args.limit).copy()
    print(f"  {len(df):,} rows")

    print("Loading dp2 ...")
    dp2 = pd.read_parquet(args.dp2_parquet)
    dp2["assay_scaffolded_epitope_id"] = dp2["assay_scaffolded_epitope_id"].astype(str).str.lower()
    dp2_tok = dp2.drop_duplicates("assay_scaffolded_epitope_id").set_index("assay_scaffolded_epitope_id")

    true_dir = Path(args.true_dir)

    mpnn_idx = {}
    if args.mpnn_pdb_dir:
        print("Indexing MPNN PDBs ...")
        mpnn_idx = index_mpnn_pdbs(Path(args.mpnn_pdb_dir))
        print(f"  {len(mpnn_idx):,} PDBs indexed")

    # cache crystal SASA per (pid, true_epi_ris) to avoid recomputing
    crystal_sasa_cache = {}

    df['mpnn_sasa_retention'] = np.nan
    df['af3_sasa_retention']  = np.nan
    n_ok = n_fail = n_skip = 0

    out_csv = Path(args.out_csv)
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    write_header = True

    for i, row in df.iterrows():
        tok     = str(row['token']).lower()
        pred    = int(row['pred'])
        mpnn_id = int(row['mpnn_id'])

        if pd.isna(row.get('overall_rmsd')):
            n_skip += 1
            continue

        dp2_row = dp2_tok.loc[tok] if tok in dp2_tok.index else None
        if dp2_row is None:
            n_skip += 1
            continue

        mpnn_epi_ris = parse_index_list(dp2_row['scaffolded_epitope_chunk_resindices'])
        true_epi_ris = parse_index_list(dp2_row['epitope_chunk_resindices'])
        pid          = dp2_row['id']
        true_pdb     = true_dir / f"{pid}.pdb"

        if not true_pdb.exists() or not mpnn_epi_ris or not true_epi_ris:
            n_skip += 1
            continue

        # get/cache crystal SASA
        cache_key = (pid, tuple(true_epi_ris))
        if cache_key not in crystal_sasa_cache:
            crystal_sasa_cache[cache_key] = get_crystal_sasa(true_pdb, true_epi_ris)
        crystal_sasa = crystal_sasa_cache[cache_key]

        if crystal_sasa is None:
            n_skip += 1
            continue

        did_something = False

        # --- MPNN SASA retention ---
        mpnn_pdb = None
        if 'mpnn_pdb' in row and pd.notna(row['mpnn_pdb']):
            mpnn_pdb = Path(str(row['mpnn_pdb']))

        if mpnn_pdb and Path(mpnn_pdb).exists():
            ret = sasa_retention_from_pdb(
                Path(mpnn_pdb), true_pdb,
                mpnn_epi_ris, true_epi_ris,
                crystal_sasa, chain="A"
            )
            if ret is not None:
                df.at[i, 'mpnn_sasa_retention'] = ret
                did_something = True

        # --- AF3 SASA retention ---
        ws = int(row['af3_window_start']) if pd.notna(row.get('af3_window_start')) else 0
        af3_epi_ris = [ws + i for i in mpnn_epi_ris]

        af3_dir_path = row.get('af3_dir')
        if pd.notna(af3_dir_path):
            af3_dir = Path(str(af3_dir_path))
            af3_cif = find_af3_cif(af3_dir)
            if af3_cif and af3_cif.exists():
                tmp_pdb = Path(cif_to_tmp_pdb(af3_cif))
                ret = sasa_retention_from_pdb(
                    tmp_pdb, true_pdb,
                    af3_epi_ris, true_epi_ris,
                    crystal_sasa, chain="A",
                    is_tmp=True
                )
                if ret is not None:
                    df.at[i, 'af3_sasa_retention'] = ret
                    did_something = True
                elif tmp_pdb.exists():
                    os.unlink(str(tmp_pdb))

        if did_something:
            n_ok += 1
        else:
            n_fail += 1

        # write row incrementally to avoid losing data on crash
        row_out = pd.DataFrame([df.loc[i]])
        row_out.to_csv(out_csv, mode="a", header=write_header, index=False)
        write_header = False

        if (n_ok + n_fail + n_skip) % 500 == 0:
            print(f"  {n_ok+n_fail+n_skip:,} processed  ok={n_ok:,}  fail={n_fail:,}  skip={n_skip:,}")

    out_csv = Path(args.out_csv)
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_csv, index=False)

    print(f"\nWrote {len(df):,} rows -> {out_csv}")
    print(f"ok={n_ok:,}  fail={n_fail:,}  skip={n_skip:,}")

    # correlations
    print(f"\n=== CORRELATIONS WITH af3_n_clash_res ===")
    for col in ['mpnn_sasa_retention', 'af3_sasa_retention']:
        valid = df[df[col].notna() & df['af3_n_clash_res'].notna()]
        if len(valid) > 10:
            corr = valid[col].corr(valid['af3_n_clash_res'])
            print(f"  {col:25s}: r = {corr:.3f}  (n={len(valid):,})")

    print(f"\n=== MEAN SASA RETENTION: PASSING vs FAILING ===")
    passing_mask = (
        (df['mean_pae'] < 5.0) &
        (df['overall_rmsd'] <= 2.0) &
        (df['epitope_chunk_rmsd'] <= 1.0) &
        (df['af3_n_clash_res'].fillna(1) == 0)
    )
    for col in ['mpnn_sasa_retention', 'af3_sasa_retention']:
        valid = df[df[col].notna()]
        pass_mean = valid[passing_mask[valid.index]]['_col_'.replace('_col_', col)].mean() if passing_mask.sum() > 0 else float('nan')
        fail_mean = valid[~passing_mask[valid.index]][col].mean()
        pass_vals = valid[passing_mask[valid.index]][col]
        fail_vals = valid[~passing_mask[valid.index]][col]
        print(f"  {col:25s}: passing={pass_vals.mean():.3f}  failing={fail_vals.mean():.3f}")


if __name__ == "__main__":
    main()
