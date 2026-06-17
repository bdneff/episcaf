#!/usr/bin/env python3
"""
compute_cylinder_clash.py

For each design in metrics_decomposed.csv, place a dummy antibody cylinder
at the epitope face of the AF3 predicted structure and count scaffold residues
that clash with it.

Method:
  1. Load AF3 CIF for the design
  2. Extract epitope CA coordinates
  3. Fit a plane to the epitope CAs using PCA
     - PC1, PC2 define the epitope plane
     - PC3 (normal vector) points outward from the antigen
  4. Place a cylinder of radius CYLINDER_RADIUS and height CYLINDER_HEIGHT
     centered on the epitope centroid, extending along the outward normal
  5. Count scaffold residues (non-epitope) with CA or any atom within cylinder
  6. Add as new columns to metrics_decomposed.csv:
     - cylinder_ca_clashes     : number of scaffold CA atoms inside cylinder
     - cylinder_allatom_clashes: number of scaffold heavy atoms inside cylinder

Parameters (soft-coded at top of script):
  CYLINDER_RADIUS : 12.0 Angstroms
  CYLINDER_HEIGHT : 40.0 Angstroms
  CLASH_BUFFER    : 0.0 Angstroms (additional buffer around cylinder)

Usage:
    # test
    python scripts/compute_cylinder_clash.py \
        --metrics_csv   runs/run_rfd3_mpnn/04_filter/metrics_decomposed.csv \
        --dp2_parquet   datasets/dp2.parquet \
        --out_csv       runs/run_rfd3_mpnn/04_filter/metrics_cylinder.csv \
        --limit         500

    # full run
    sbatch --wrap="python scripts/compute_cylinder_clash.py \\
        --metrics_csv   runs/run_rfd3_mpnn/04_filter/metrics_decomposed.csv \\
        --dp2_parquet   datasets/dp2.parquet \\
        --out_csv       runs/run_rfd3_mpnn/04_filter/metrics_cylinder.csv" \\
        -p compute -J cylinder_clash \\
        --cpus-per-task=4 --mem=64G -t 8:00:00 \\
        --chdir=/home/bneff/rfd3/repo_refactored \\
        -o logs/cylinder_clash.out \\
        -e logs/cylinder_clash.err
"""

import argparse
import gzip
import json
import math
import os
import tempfile
from pathlib import Path
from typing import List, Optional, Tuple

import gemmi
import numpy as np
import pandas as pd

# ---------------------------------------------------------------------------
# Parameters — adjust here
# ---------------------------------------------------------------------------

CYLINDER_RADIUS = 16.0   # Angstroms
CYLINDER_HEIGHT = 20.0   # Angstroms
CLASH_BUFFER    = 0.0    # extra buffer around cylinder
CYLINDER_OFFSET = 4.0    # Angstroms past furthest epitope CA to cylinder base


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

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


def load_af3_atoms(af3_dir: Path) -> Optional[Tuple[np.ndarray, np.ndarray, List[str], List[int]]]:
    """
    Load all heavy atoms from AF3 CIF output.
    Returns (coords, ca_coords, atom_names, res_indices) or None.
    coords: (N, 3) all heavy atom positions
    ca_mask: boolean mask for CA atoms
    res_idx: residue index (0-based) for each atom
    """
    cif_path = next(af3_dir.glob("*_model.cif"), None)
    if cif_path is None:
        return None

    try:
        with gzip.open(str(cif_path) + ".gz", "rt") if not cif_path.exists() else open(cif_path) as _:
            pass
    except Exception:
        pass

    # find the CIF (may or may not be gzipped)
    cif_gz = next(af3_dir.glob("*.cif.gz"), None)
    if cif_gz:
        with gzip.open(cif_gz, "rt") as f:
            doc = gemmi.cif.read_string(f.read())
    elif cif_path:
        doc = gemmi.cif.read(str(cif_path))
    else:
        return None

    try:
        st = gemmi.make_structure_from_block(doc.sole_block())
    except Exception:
        return None

    coords    = []
    ca_mask   = []
    res_idx   = []
    res_count = 0

    for chain in st[0]:
        for res in chain:
            for atom in res:
                if atom.element == gemmi.Element('H'):
                    continue
                coords.append([atom.pos.x, atom.pos.y, atom.pos.z])
                ca_mask.append(atom.name == 'CA')
                res_idx.append(res_count)
            res_count += 1

    if not coords:
        return None

    return (np.array(coords, dtype=float),
            np.array(ca_mask, dtype=bool),
            np.array(res_idx, dtype=int),
            res_count)


def fit_epitope_plane(epi_ca_coords: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """
    Fit a plane to epitope CA coordinates using PCA.
    Returns (centroid, normal_vector).
    normal_vector points in direction of PC3 (smallest variance = normal to plane).
    """
    centroid = epi_ca_coords.mean(axis=0)
    centered = epi_ca_coords - centroid
    _, _, Vt = np.linalg.svd(centered)
    normal = Vt[-1]  # last row = direction of least variance = plane normal
    return centroid, normal


def orient_normal_outward(normal: np.ndarray, centroid: np.ndarray,
                          all_ca_coords: np.ndarray) -> np.ndarray:
    """
    Flip normal vector if it points toward the protein center of mass
    (we want it pointing away from the protein, toward where an antibody would bind).
    """
    protein_com = all_ca_coords.mean(axis=0)
    to_com = protein_com - centroid
    if np.dot(normal, to_com) > 0:
        normal = -normal
    return normal


def point_in_cylinder(points: np.ndarray, centroid: np.ndarray,
                      normal: np.ndarray, radius: float, height: float,
                      buffer: float = 0.0) -> np.ndarray:
    """
    Test which points fall inside a cylinder.
    Cylinder axis: centroid -> centroid + height * normal
    Returns boolean array.
    """
    r = radius + buffer
    h = height + buffer

    # vector from centroid to each point
    v = points - centroid

    # projection along cylinder axis
    proj = v @ normal  # scalar projection, shape (N,)

    # distance from axis
    perp = v - np.outer(proj, normal)
    dist_from_axis = np.linalg.norm(perp, axis=1)

    # inside cylinder: projection in [0, h] and radial distance <= r
    inside = (proj >= 0) & (proj <= h) & (dist_from_axis <= r)
    return inside


def count_cylinder_clashes(
    coords:      np.ndarray,
    ca_mask:     np.ndarray,
    res_idx:     np.ndarray,
    epi_ris:     List[int],
    centroid:    np.ndarray,
    normal:      np.ndarray,
) -> Tuple[int, int]:
    """
    Count scaffold residues (non-epitope) with CA or any heavy atom inside cylinder.
    Returns (ca_clashes, allatom_clashes).
    """
    epi_set = set(epi_ris)

    # scaffold mask — exclude epitope residues
    scaf_mask = np.array([ri not in epi_set for ri in res_idx])

    scaf_coords = coords[scaf_mask]
    scaf_ca     = ca_mask[scaf_mask]
    scaf_ri     = res_idx[scaf_mask]

    inside = point_in_cylinder(scaf_coords, centroid, normal,
                               CYLINDER_RADIUS, CYLINDER_HEIGHT, CLASH_BUFFER)

    # CA clashes: unique residues with CA inside cylinder
    ca_inside_ri = set(scaf_ri[inside & scaf_ca])
    ca_clashes   = len(ca_inside_ri)

    # all-atom clashes: unique residues with any heavy atom inside cylinder
    aa_inside_ri = set(scaf_ri[inside])
    aa_clashes   = len(aa_inside_ri)

    return ca_clashes, aa_clashes


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--metrics_csv",  required=True)
    parser.add_argument("--dp2_parquet",  required=True)
    parser.add_argument("--out_csv",      required=True)
    parser.add_argument("--limit",        type=int, default=0)
    args = parser.parse_args()

    print("Loading metrics CSV ...")
    df = pd.read_csv(args.metrics_csv, low_memory=False)
    if args.limit > 0:
        df = df.head(args.limit).copy()
    print(f"  {len(df):,} rows")

    print("Loading dp2 ...")
    dp2 = pd.read_parquet(args.dp2_parquet)
    dp2["assay_scaffolded_epitope_id"] = (dp2["assay_scaffolded_epitope_id"]
                                          .astype(str).str.lower())
    dp2_tok = (dp2.drop_duplicates("assay_scaffolded_epitope_id")
                  .set_index("assay_scaffolded_epitope_id"))

    print(f"Cylinder: radius={CYLINDER_RADIUS}A  height={CYLINDER_HEIGHT}A")
    print()

    df['cylinder_ca_clashes']      = np.nan
    df['cylinder_allatom_clashes'] = np.nan
    n_ok = n_fail = n_skip = 0

    out_csv = Path(args.out_csv)
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    write_header = True

    for i, row in df.iterrows():
        tok = str(row['token']).lower()

        if pd.isna(row.get('overall_rmsd')):
            n_skip += 1
            row_out = pd.DataFrame([df.loc[i]])
            row_out.to_csv(out_csv, mode='a', header=write_header, index=False)
            write_header = False
            continue

        dp2_row = dp2_tok.loc[tok] if tok in dp2_tok.index else None
        if dp2_row is None:
            n_skip += 1
            row_out = pd.DataFrame([df.loc[i]])
            row_out.to_csv(out_csv, mode='a', header=write_header, index=False)
            write_header = False
            continue

        ws           = int(row['af3_window_start']) if pd.notna(row.get('af3_window_start')) else 0
        mpnn_epi_ris = parse_index_list(dp2_row['scaffolded_epitope_chunk_resindices'])
        af3_epi_ris  = [ws + i for i in mpnn_epi_ris]

        af3_dir_path = row.get('af3_dir')
        if pd.isna(af3_dir_path):
            n_skip += 1
            row_out = pd.DataFrame([df.loc[i]])
            row_out.to_csv(out_csv, mode='a', header=write_header, index=False)
            write_header = False
            continue

        af3_dir = Path(str(af3_dir_path))
        result  = load_af3_atoms(af3_dir)

        if result is None:
            n_fail += 1
            row_out = pd.DataFrame([df.loc[i]])
            row_out.to_csv(out_csv, mode='a', header=write_header, index=False)
            write_header = False
            continue

        coords, ca_mask, res_idx, n_res = result

        if not af3_epi_ris or max(af3_epi_ris) >= n_res:
            n_fail += 1
            row_out = pd.DataFrame([df.loc[i]])
            row_out.to_csv(out_csv, mode='a', header=write_header, index=False)
            write_header = False
            continue

        try:
            # get epitope CA coordinates
            epi_ca_mask = np.array([
                (ri in set(af3_epi_ris)) and ca_mask[j]
                for j, ri in enumerate(res_idx)
            ])
            epi_ca_coords = coords[epi_ca_mask]

            if len(epi_ca_coords) < 3:
                raise ValueError(f"Too few epitope CA atoms: {len(epi_ca_coords)}")

            # fit plane and orient normal
            all_ca_coords = coords[ca_mask]
            centroid, normal = fit_epitope_plane(epi_ca_coords)
            normal = orient_normal_outward(normal, centroid, all_ca_coords)

            # offset cylinder base to 4A past furthest epitope atom
            projections = (epi_ca_coords - centroid) @ normal
            offset = float(projections.max()) + CYLINDER_OFFSET
            shifted_centroid = centroid + offset * normal

            # count clashes
            ca_clashes, aa_clashes = count_cylinder_clashes(
                coords, ca_mask, res_idx, af3_epi_ris, shifted_centroid, normal
            )

            df.at[i, 'cylinder_ca_clashes']      = ca_clashes
            df.at[i, 'cylinder_allatom_clashes'] = aa_clashes
            n_ok += 1

        except Exception as e:
            n_fail += 1

        row_out = pd.DataFrame([df.loc[i]])
        row_out.to_csv(out_csv, mode='a', header=write_header, index=False)
        write_header = False

        if (n_ok + n_fail + n_skip) % 1000 == 0:
            print(f"  {n_ok+n_fail+n_skip:,}  ok={n_ok:,}  fail={n_fail:,}  skip={n_skip:,}")

    print(f"\nDone: ok={n_ok:,}  fail={n_fail:,}  skip={n_skip:,}")
    print(f"Wrote: {out_csv}")

    # correlation summary
    out = pd.read_csv(out_csv, low_memory=False)
    for col in ['cylinder_ca_clashes', 'cylinder_allatom_clashes']:
        for ref in ['af3_n_clash_res']:
            valid = out[out[col].notna() & pd.to_numeric(out[ref], errors='coerce').notna()]
            if len(valid) > 10:
                r = pd.to_numeric(valid[col]).corr(pd.to_numeric(valid[ref]))
                print(f"  {col} vs {ref}: r={r:.3f}  n={len(valid):,}")


if __name__ == "__main__":
    main()
