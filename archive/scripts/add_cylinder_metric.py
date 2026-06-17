#!/usr/bin/env python3
"""
add_cylinder_metric.py

Append cylinder clash columns to a per-design metrics CSV, using the locked
official geometry and the centroid offset convention (matches the validated scan).

Two metrics, both computed on chain-A CAs, both EXCLUDING epitope residues:
    cylinder_ca_clashes      : binary count of non-epitope residues inside the cylinder
    cylinder_weighted_clash  : depth-weighted sum -- each inside non-epitope residue
                               contributes (1 - dist_from_axis/radius)**depth_power,
                               so residues near the central axis (deep in the
                               antibody's path) count ~1 and ones grazing the radial
                               wall count ~0.

Official params: radius=16, height=40, offset=-4 (from epitope centroid).

Example
-------
    python scripts/add_cylinder_metric.py \
        --metrics_csv runs/run_rfd3_mpnn/04_filter/metrics_decomposed.csv \
        --dp2_parquet datasets/dp2.parquet \
        --out_csv     runs/run_rfd3_mpnn/04_filter/metrics_cylinder_full.csv \
        --depth_power 1 --limit 1000
"""

import argparse
import gzip
import math
from pathlib import Path
from typing import List, Optional, Tuple

import gemmi
import numpy as np
import pandas as pd

RADIUS = 16.0
HEIGHT = 40.0
OFFSET = -4.0

TOKEN_COLS   = ["token", "assay_scaffolded_epitope_id"]
AF3_DIR_COLS = ["af3_dir"]
DP2_EPI_COLS = ["assay_scaffolded_epitope_chunk_resindices"]
DP2_EPI_MPNN = ["scaffolded_epitope_chunk_resindices"]
WS_COLS      = ["af3_window_start"]
TRUE_CLASH   = ["af3_n_clash_res"]


def parse_index_list(x) -> List[int]:
    if x is None or (isinstance(x, float) and math.isnan(x)):
        return []
    if isinstance(x, (list, tuple, np.ndarray)):
        return [int(i) for i in x]
    s = str(x).strip().replace("[", "").replace("]", "").replace(",", " ")
    out = []
    for t in s.split():
        try:
            out.append(int(t))
        except ValueError:
            pass
    return out


def first_present(cols, cands):
    for c in cands:
        if c in cols:
            return c
    return None


def read_gemmi(p: Path) -> gemmi.Structure:
    s = str(p)
    if s.endswith(".gz"):
        with gzip.open(p, "rt") as f:
            doc = gemmi.cif.read_string(f.read())
        return gemmi.make_structure_from_block(doc.sole_block())
    if p.suffix.lower() == ".cif":
        doc = gemmi.cif.read(str(p))
        return gemmi.make_structure_from_block(doc.sole_block())
    return gemmi.read_structure(str(p))


def get_chainA(st):
    for ch in st[0]:
        if ch.name == "A":
            return ch
    return st[0][0]


def find_af3_cif(af3_dir: Path) -> Optional[Path]:
    for pat in ("*_model.cif", "*_model.cif.gz", "model.cif", "model.cif.gz"):
        h = next(af3_dir.glob(pat), None) or next(af3_dir.rglob(pat), None)
        if h:
            return h
    return None


def load_af3_ca(af3_dir: Path):
    cif = find_af3_cif(Path(af3_dir))
    if cif is None:
        return None
    try:
        st = read_gemmi(cif)
    except Exception:
        return None
    chA = get_chainA(st)
    coords, ridx = [], []
    for i, res in enumerate(chA):           # one CA per residue
        a = res.find_atom("CA", altloc="*")
        if a:
            coords.append([a.pos.x, a.pos.y, a.pos.z])
            ridx.append(i)
    if not coords:
        return None
    return np.asarray(coords, float), np.asarray(ridx, int)


def cylinder_metrics(ca, res_idx, epi_set, centroid, normal,
                     radius, height, offset, depth_power) -> Tuple[int, float]:
    """Returns (binary_count, weighted_sum) over non-epitope residues inside the
    cylinder. depth = 1 at the axis, 0 at the radial wall."""
    base = centroid + offset * normal
    v = ca - base
    proj = v @ normal
    perp = v - np.outer(proj, normal)
    dist = np.linalg.norm(perp, axis=1)
    inside = (proj >= 0.0) & (proj <= height) & (dist <= radius)

    count = 0
    wsum = 0.0
    for r, ins, d in zip(res_idx, inside, dist):
        # *** epitope residues are always excluded from the clash count ***
        if (not ins) or (int(r) in epi_set):
            continue
        count += 1
        depth = max(0.0, 1.0 - d / radius)       # radial depth, [0,1]
        wsum += depth ** depth_power
    return count, float(wsum)


def compute(af3_dir, epi_ris, depth_power):
    loaded = load_af3_ca(Path(af3_dir))
    if loaded is None:
        return None
    ca, res_idx = loaded
    if not epi_ris or max(epi_ris) >= int(res_idx.max()) + 1:
        return None
    epi_set = set(epi_ris)
    epi_mask = np.array([int(r) in epi_set for r in res_idx])
    if epi_mask.sum() < 3:
        return None
    epi_ca = ca[epi_mask]
    centroid = epi_ca.mean(axis=0)
    _, _, Vt = np.linalg.svd(epi_ca - centroid)
    normal = Vt[-1]
    if np.dot(normal, ca.mean(axis=0) - centroid) > 0:
        normal = -normal
    return cylinder_metrics(ca, res_idx, epi_set, centroid, normal,
                            RADIUS, HEIGHT, OFFSET, depth_power)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--metrics_csv", required=True)
    ap.add_argument("--dp2_parquet", required=True)
    ap.add_argument("--out_csv", required=True)
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--depth_power", type=float, default=1.0,
                    help="exponent on radial depth; 1=linear, 2=steeper (deep "
                         "residues dominate). Default 1.")
    args = ap.parse_args()

    print(f"Cylinder: radius={RADIUS} height={HEIGHT} offset={OFFSET} "
          f"(centroid)  depth_power={args.depth_power}")
    df = pd.read_csv(args.metrics_csv, low_memory=False)
    if args.limit > 0:
        df = df.head(args.limit).copy()
    print(f"  {len(df):,} rows")

    tok_col = first_present(df.columns, TOKEN_COLS)
    dir_col = first_present(df.columns, AF3_DIR_COLS)
    ws_col  = first_present(df.columns, WS_COLS)
    if tok_col is None or dir_col is None:
        raise SystemExit("need token + af3_dir columns")

    dp2 = pd.read_parquet(args.dp2_parquet)
    dp2["assay_scaffolded_epitope_id"] = (dp2["assay_scaffolded_epitope_id"]
                                          .astype(str).str.lower())
    epi_col  = first_present(dp2.columns, DP2_EPI_COLS)
    epi_mpnn = first_present(dp2.columns, DP2_EPI_MPNN)
    dp2_lookup = dp2.drop_duplicates("assay_scaffolded_epitope_id").set_index(
        "assay_scaffolded_epitope_id")

    binc = np.full(len(df), np.nan)
    wgt  = np.full(len(df), np.nan)
    n_ok = n_fail = 0
    for pos, (_, row) in enumerate(df.iterrows()):
        tok = str(row[tok_col]).lower()
        if tok not in dp2_lookup.index or pd.isna(row.get(dir_col)):
            n_fail += 1; continue
        drow = dp2_lookup.loc[tok]
        if epi_col is not None:
            epi_ris = parse_index_list(drow[epi_col])
        else:
            ws = int(row[ws_col]) if (ws_col and pd.notna(row.get(ws_col))) else 0
            epi_ris = [ws + i for i in parse_index_list(drow[epi_mpnn])]
        res = compute(str(row[dir_col]), epi_ris, args.depth_power)
        if res is None:
            n_fail += 1
        else:
            binc[pos], wgt[pos] = res
            n_ok += 1
        if (pos + 1) % 1000 == 0:
            print(f"  {pos+1:,}  ok={n_ok:,} fail={n_fail:,}")

    df["cylinder_ca_clashes"] = binc
    df["cylinder_weighted_clash"] = wgt
    out = Path(args.out_csv); out.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out, index=False)
    print(f"\nWrote {out}  (ok={n_ok:,} fail={n_fail:,})")

    tc = first_present(df.columns, TRUE_CLASH)
    if tc is not None:
        truth = pd.to_numeric(df[tc], errors="coerce")
        print(f"\nCorrelation vs true clash ({tc}):")
        print("  proxy                      Pearson   Spearman")
        for c in ("cylinder_ca_clashes", "cylinder_weighted_clash"):
            p = pd.to_numeric(df[c], errors="coerce")
            m = p.notna() & truth.notna()
            if m.sum() < 10:
                continue
            print(f"  {c:24s}  {p[m].corr(truth[m]):+.3f}    "
                  f"{p[m].corr(truth[m], method='spearman'):+.3f}")


if __name__ == "__main__":
    main()
