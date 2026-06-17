#!/usr/bin/env python3
"""
scan_cylinder_params.py

Grid-search the dummy-antibody-cylinder accessibility metric against the
known-antibody clash filter, reproducing the radius x height x offset sweep
from the cylinder_clash_metric writeup (Table 1).

For every design it:
  1. loads the AF3-predicted chain-A CA coordinates
  2. fits a plane to the epitope CAs (PCA) and orients the normal OUTWARD
     (away from the protein center of mass)
  3. for each (radius, height, offset) places a cylinder with base at
        epitope_centroid + offset * normal
     extending +height along the outward normal
        (offset < 0 pushes the base into the epitope; offset > 0 away)
  4. labels scaffold (non-epitope) residues whose CA is inside as cylinder-clash
  5. compares that residue set to the TRUE antibody-clash residue set

Aggregated across designs, for each parameter combo it reports (micro-averaged
over residues):
  - precision : of residues the cylinder flags, fraction that truly clash
  - recall    : of residues that truly clash,  fraction the cylinder catches
  - pearson_r : correlation of per-design cylinder count vs true clash count
  - n_designs : designs that contributed (had a usable AF3 model + true labels)

The cylinder is CA-based (matches the writeup); the true clash set is whatever
your known-antibody filter wrote (heavy-atom based). Both are in AF3 chain-A
0-based residue index space.

Example
-------
    python scripts/scan_cylinder_params.py \
        --metrics_csv runs/run_rfd3_mpnn/04_filter/metrics_decomposed.csv \
        --dp2_parquet datasets/dp2.parquet \
        --out_csv     runs/run_rfd3_mpnn/04_filter/cylinder_param_scan.csv \
        --limit 500
"""

import argparse
import gzip
import math
from pathlib import Path
from typing import List, Optional, Tuple

import gemmi
import numpy as np
import pandas as pd

# ---------------------------------------------------------------------------
# Default grid (matches the writeup). Override on the command line if needed.
# ---------------------------------------------------------------------------
DEFAULT_RADII   = [8.0, 12.0, 16.0, 20.0]
DEFAULT_HEIGHTS = [20.0, 40.0]
DEFAULT_OFFSETS = [-4.0, -2.0, 0.0, 2.0, 4.0]

# Column-name candidates, tried in order, so the script tolerates the slightly
# different naming across your metrics files.
TOKEN_COLS       = ["token", "assay_scaffolded_epitope_id"]
AF3_DIR_COLS     = ["af3_dir"]
TRUE_CLASH_COLS  = ["af3_clash_resindices"]
DP2_EPI_AF3_COLS = ["assay_scaffolded_epitope_chunk_resindices"]
DP2_EPI_MPNN_COLS = ["scaffolded_epitope_chunk_resindices"]
WINDOW_START_COLS = ["af3_window_start"]


# ---------------------------------------------------------------------------
# Parsing / IO helpers
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


def first_present(df_cols, candidates) -> Optional[str]:
    for c in candidates:
        if c in df_cols:
            return c
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


def find_af3_cif(af3_dir: Path) -> Optional[Path]:
    for pat in ("*_model.cif", "*_model.cif.gz"):
        hit = next(af3_dir.glob(pat), None)
        if hit:
            return hit
    for pat in ("model.cif", "model.cif.gz"):
        hit = next(af3_dir.rglob(pat), None)
        if hit:
            return hit
    return None


def load_af3_ca(af3_dir: Path) -> Optional[Tuple[np.ndarray, np.ndarray]]:
    """
    Returns (ca_coords [M,3], res_idx [M]) for chain A, where res_idx is the
    0-based position of the residue within chain A (same frame the known-antibody
    clash filter uses for af3_clash_resindices). Residues lacking a CA are skipped
    but do not shift the indices of later residues.
    """
    cif = find_af3_cif(Path(af3_dir))
    if cif is None:
        return None
    try:
        st = read_gemmi_structure(cif)
    except Exception:
        return None
    chA = get_chainA(st)

    coords, ridx = [], []
    for i, res in enumerate(chA):
        a = res.find_atom("CA", altloc="*")
        if a:
            coords.append([a.pos.x, a.pos.y, a.pos.z])
            ridx.append(i)
    if not coords:
        return None
    return np.asarray(coords, dtype=float), np.asarray(ridx, dtype=int)


# ---------------------------------------------------------------------------
# Geometry
# ---------------------------------------------------------------------------
def fit_plane_normal(epi_ca: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    centroid = epi_ca.mean(axis=0)
    _, _, Vt = np.linalg.svd(epi_ca - centroid)
    return centroid, Vt[-1]


def orient_outward(normal: np.ndarray, centroid: np.ndarray,
                   all_ca: np.ndarray) -> np.ndarray:
    if np.dot(normal, all_ca.mean(axis=0) - centroid) > 0:
        return -normal
    return normal


def cylinder_clash_set(ca_coords: np.ndarray, res_idx: np.ndarray,
                       epi_set: set, centroid: np.ndarray, normal: np.ndarray,
                       radius: float, height: float, offset: float,
                       base_shift: float = 0.0) -> set:
    # base_shift = 0 -> offset measured from the epitope centroid
    # base_shift = max projection of epitope CAs along the normal
    #              -> offset measured from the epitope surface (furthest CA)
    base = centroid + (base_shift + offset) * normal
    v = ca_coords - base
    proj = v @ normal
    perp = v - np.outer(proj, normal)
    dist = np.linalg.norm(perp, axis=1)
    inside = (proj >= 0.0) & (proj <= height) & (dist <= radius)
    return {int(r) for r, ins in zip(res_idx, inside) if ins and int(r) not in epi_set}


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--metrics_csv", required=True,
                    help="per-design metrics that already contain the true clash result")
    ap.add_argument("--dp2_parquet", required=True)
    ap.add_argument("--out_csv", required=True)
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--radii", type=float, nargs="+", default=DEFAULT_RADII)
    ap.add_argument("--heights", type=float, nargs="+", default=DEFAULT_HEIGHTS)
    ap.add_argument("--offsets", type=float, nargs="+", default=DEFAULT_OFFSETS)
    ap.add_argument("--offset_from", choices=["centroid", "surface"],
                    default="centroid",
                    help="reference point for offset: epitope 'centroid' (matches "
                         "the writeup text) or 'surface' (furthest epitope CA along "
                         "the normal, matches compute_cylinder_clash.py)")
    args = ap.parse_args()
    print(f"offset reference: {args.offset_from}")

    print("Loading metrics ...")
    df = pd.read_csv(args.metrics_csv, low_memory=False)
    if "af3_clash_status" in df.columns:
        before = len(df)
        df = df[df["af3_clash_status"] == "ok"].copy()
        print(f"  kept {len(df):,}/{before:,} rows with af3_clash_status == 'ok'")
    if args.limit > 0:
        df = df.head(args.limit).copy()

    tok_col   = first_present(df.columns, TOKEN_COLS)
    dir_col   = first_present(df.columns, AF3_DIR_COLS)
    clash_col = first_present(df.columns, TRUE_CLASH_COLS)
    ws_col    = first_present(df.columns, WINDOW_START_COLS)
    for label, col in (("token", tok_col), ("af3_dir", dir_col),
                       ("true-clash list", clash_col)):
        if col is None:
            raise SystemExit(
                f"Could not find a {label} column in {args.metrics_csv}. "
                f"Columns present: {list(df.columns)[:40]}")
    print(f"  token={tok_col}  af3_dir={dir_col}  true_clash={clash_col}  "
          f"window_start={ws_col}")

    print("Loading dp2 ...")
    dp2 = pd.read_parquet(args.dp2_parquet)
    dp2_tok = "assay_scaffolded_epitope_id"
    if dp2_tok not in dp2.columns:
        raise SystemExit(f"dp2 missing {dp2_tok}")
    dp2[dp2_tok] = dp2[dp2_tok].astype(str).str.lower()

    epi_af3_col  = first_present(dp2.columns, DP2_EPI_AF3_COLS)
    epi_mpnn_col = first_present(dp2.columns, DP2_EPI_MPNN_COLS)
    if epi_af3_col is None and epi_mpnn_col is None:
        raise SystemExit("dp2 has no epitope residue-index column.")
    dp2_lookup = dp2.drop_duplicates(dp2_tok).set_index(dp2_tok)
    print(f"  epitope indices: af3_col={epi_af3_col}  mpnn_col={epi_mpnn_col}")

    grid = [(r, h, o) for r in args.radii for h in args.heights for o in args.offsets]
    # accumulators per param combo
    tp = {g: 0 for g in grid}
    fp = {g: 0 for g in grid}
    fn = {g: 0 for g in grid}
    cyl_counts  = {g: [] for g in grid}
    true_counts = {g: [] for g in grid}

    n_used = n_skip = 0
    for _, row in df.iterrows():
        tok = str(row[tok_col]).lower()
        true_set = set(parse_index_list(row.get(clash_col)))
        # only designs with a real ground-truth clash result contribute
        if pd.isna(row.get(clash_col)) or tok not in dp2_lookup.index:
            n_skip += 1
            continue

        # epitope residue indices in AF3 chain-A frame
        drow = dp2_lookup.loc[tok]
        if epi_af3_col is not None:
            epi_ris = parse_index_list(drow[epi_af3_col])
        else:
            ws = int(row[ws_col]) if (ws_col and pd.notna(row.get(ws_col))) else 0
            epi_ris = [ws + i for i in parse_index_list(drow[epi_mpnn_col])]
        if len(epi_ris) < 3:
            n_skip += 1
            continue

        loaded = load_af3_ca(Path(str(row[dir_col])))
        if loaded is None:
            n_skip += 1
            continue
        ca_coords, res_idx = loaded

        epi_set = set(epi_ris)
        epi_mask = np.array([int(r) in epi_set for r in res_idx])
        if epi_mask.sum() < 3:
            n_skip += 1
            continue

        centroid, normal = fit_plane_normal(ca_coords[epi_mask])
        normal = orient_outward(normal, centroid, ca_coords)

        # for the 'surface' convention, shift the base out to the furthest
        # epitope CA along the (outward) normal before applying offset
        if args.offset_from == "surface":
            base_shift = float(((ca_coords[epi_mask] - centroid) @ normal).max())
        else:
            base_shift = 0.0

        for g in grid:
            r, h, o = g
            cyl = cylinder_clash_set(ca_coords, res_idx, epi_set,
                                     centroid, normal, r, h, o, base_shift)
            tp[g] += len(cyl & true_set)
            fp[g] += len(cyl - true_set)
            fn[g] += len(true_set - cyl)
            cyl_counts[g].append(len(cyl))
            true_counts[g].append(len(true_set))
        n_used += 1
        if n_used % 200 == 0:
            print(f"  {n_used} designs scored ...")

    print(f"\nScored {n_used} designs (skipped {n_skip}).")

    out_rows = []
    for g in grid:
        r, h, o = g
        prec = tp[g] / (tp[g] + fp[g]) if (tp[g] + fp[g]) else float("nan")
        rec  = tp[g] / (tp[g] + fn[g]) if (tp[g] + fn[g]) else float("nan")
        a = np.asarray(cyl_counts[g], float)
        b = np.asarray(true_counts[g], float)
        pear = float(np.corrcoef(a, b)[0, 1]) if len(a) > 2 and a.std() > 0 and b.std() > 0 else float("nan")
        out_rows.append({
            "radius": r, "height": h, "offset": o, "offset_from": args.offset_from,
            "recall": rec, "precision": prec, "pearson_r": pear,
            "tp": tp[g], "fp": fp[g], "fn": fn[g], "n_designs": n_used,
        })

    out = pd.DataFrame(out_rows).sort_values(["radius", "height", "offset"])
    out_csv = Path(args.out_csv)
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(out_csv, index=False)

    print(f"\nWrote {out_csv}\n")
    show = out[["radius", "height", "offset", "recall", "precision", "pearson_r"]]
    print(show.to_string(index=False,
          formatters={"recall": "{:.3f}".format,
                      "precision": "{:.3f}".format,
                      "pearson_r": "{:.3f}".format}))


if __name__ == "__main__":
    main()
