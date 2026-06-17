#!/usr/bin/env python3
"""
scan_weighted_cylinder.py

Grid-search the weighted-cylinder clash proxy over radius x height x offset x
depth_power, scoring each combo by its correlation (Pearson + Spearman) with the
true clash count (af3_n_clash_res) across designs. depth_power=0 reproduces the
plain binary count, so the incumbent is included in the grid.

Per design the structure is loaded once; every param combo is a cheap vector op
on the cached projection/distance, so you can run this on many designs.

Example
-------
    python scripts/scan_weighted_cylinder.py \
        --metrics_csv runs/run_rfd3_mpnn/04_filter/metrics_decomposed.csv \
        --dp2_parquet datasets/dp2.parquet \
        --out_csv     runs/run_rfd3_mpnn/04_filter/weighted_cyl_scan.csv \
        --limit 8000
"""

import argparse
import gzip
import math
from pathlib import Path
from typing import List, Optional

import gemmi
import numpy as np
import pandas as pd

DEF_RADII   = [12.0, 14.0, 16.0, 18.0, 20.0]
DEF_HEIGHTS = [20.0, 40.0]
DEF_OFFSETS = [-6.0, -4.0, -2.0, 0.0]
DEF_POWERS  = [0.0, 1.0, 2.0, 3.0]

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


def read_gemmi(p: Path):
    s = str(p)
    if s.endswith(".gz"):
        with gzip.open(p, "rt") as f:
            doc = gemmi.cif.read_string(f.read())
        return gemmi.make_structure_from_block(doc.sole_block())
    if p.suffix.lower() == ".cif":
        return gemmi.make_structure_from_block(gemmi.cif.read(str(p)).sole_block())
    return gemmi.read_structure(str(p))


def get_chainA(st):
    for ch in st[0]:
        if ch.name == "A":
            return ch
    return st[0][0]


def find_af3_cif(af3_dir: Path):
    for pat in ("*_model.cif", "*_model.cif.gz", "model.cif", "model.cif.gz"):
        h = next(af3_dir.glob(pat), None) or next(af3_dir.rglob(pat), None)
        if h:
            return h
    return None


def load_plane(af3_dir, epi_ris):
    cif = find_af3_cif(Path(af3_dir))
    if cif is None:
        return None
    try:
        st = read_gemmi(cif)
    except Exception:
        return None
    chA = get_chainA(st)
    coords, ridx = [], []
    for i, res in enumerate(chA):
        a = res.find_atom("CA", altloc="*")
        if a:
            coords.append([a.pos.x, a.pos.y, a.pos.z]); ridx.append(i)
    if not coords:
        return None
    ca = np.asarray(coords, float); res_idx = np.asarray(ridx, int)
    if not epi_ris or max(epi_ris) >= int(res_idx.max()) + 1:
        return None
    epi_set = set(epi_ris)
    epi_mask = np.fromiter((int(r) in epi_set for r in res_idx), bool, len(res_idx))
    if epi_mask.sum() < 3:
        return None
    epi_ca = ca[epi_mask]
    centroid = epi_ca.mean(axis=0)
    _, _, Vt = np.linalg.svd(epi_ca - centroid)
    normal = Vt[-1]
    if np.dot(normal, ca.mean(axis=0) - centroid) > 0:
        normal = -normal
    return ca, epi_mask, centroid, normal


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--metrics_csv", required=True)
    ap.add_argument("--dp2_parquet", required=True)
    ap.add_argument("--out_csv", required=True)
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--radii", type=float, nargs="+", default=DEF_RADII)
    ap.add_argument("--heights", type=float, nargs="+", default=DEF_HEIGHTS)
    ap.add_argument("--offsets", type=float, nargs="+", default=DEF_OFFSETS)
    ap.add_argument("--powers", type=float, nargs="+", default=DEF_POWERS)
    ap.add_argument("--top", type=int, default=15)
    args = ap.parse_args()

    grid = [(r, h, o, p) for o in args.offsets for r in args.radii
            for h in args.heights for p in args.powers]
    print(f"{len(grid)} param combos "
          f"(r={args.radii} h={args.heights} off={args.offsets} pow={args.powers})")

    df = pd.read_csv(args.metrics_csv, low_memory=False)
    if args.limit > 0:
        df = df.head(args.limit).copy()
    tok_col = first_present(df.columns, TOKEN_COLS)
    dir_col = first_present(df.columns, AF3_DIR_COLS)
    ws_col  = first_present(df.columns, WS_COLS)
    tc_col  = first_present(df.columns, TRUE_CLASH)
    if tc_col is None:
        raise SystemExit("need af3_n_clash_res to score the sweep")

    dp2 = pd.read_parquet(args.dp2_parquet)
    dp2["assay_scaffolded_epitope_id"] = (dp2["assay_scaffolded_epitope_id"]
                                          .astype(str).str.lower())
    epi_col  = first_present(dp2.columns, DP2_EPI_COLS)
    epi_mpnn = first_present(dp2.columns, DP2_EPI_MPNN)
    dp2_lookup = dp2.drop_duplicates("assay_scaffolded_epitope_id").set_index(
        "assay_scaffolded_epitope_id")

    per_combo = {g: [] for g in grid}        # per-design weighted sum
    truth = []
    n_ok = n_fail = 0
    for pos, (_, row) in enumerate(df.iterrows()):
        tok = str(row[tok_col]).lower()
        if tok not in dp2_lookup.index or pd.isna(row.get(dir_col)) \
           or pd.isna(row.get(tc_col)):
            n_fail += 1; continue
        drow = dp2_lookup.loc[tok]
        if epi_col is not None:
            epi_ris = parse_index_list(drow[epi_col])
        else:
            ws = int(row[ws_col]) if (ws_col and pd.notna(row.get(ws_col))) else 0
            epi_ris = [ws + i for i in parse_index_list(drow[epi_mpnn])]
        plane = load_plane(str(row[dir_col]), epi_ris)
        if plane is None:
            n_fail += 1; continue
        ca, epi_mask, centroid, normal = plane
        scaf = ~epi_mask

        # cache proj/dist per offset, then apply (radius,height,power)
        by_off = {}
        for o in args.offsets:
            base = centroid + o * normal
            v = ca - base
            proj = v @ normal
            perp = v - np.outer(proj, normal)
            dist = np.linalg.norm(perp, axis=1)
            by_off[o] = (proj, dist)

        for (r, h, o, p) in grid:
            proj, dist = by_off[o]
            inside = scaf & (proj >= 0.0) & (proj <= h) & (dist <= r)
            if not inside.any():
                per_combo[(r, h, o, p)].append(0.0)
                continue
            depth = np.clip(1.0 - dist[inside] / r, 0.0, 1.0)
            per_combo[(r, h, o, p)].append(float((depth ** p).sum()))
        truth.append(float(row[tc_col]))
        n_ok += 1
        if (pos + 1) % 1000 == 0:
            print(f"  {pos+1:,}  ok={n_ok:,} fail={n_fail:,}")

    truth = np.asarray(truth)
    print(f"\nScored {n_ok} designs (skipped {n_fail}).")

    rows = []
    for g in grid:
        vals = np.asarray(per_combo[g])
        if len(vals) < 10 or vals.std() == 0:
            continue
        pear = float(np.corrcoef(vals, truth)[0, 1])
        spear = float(pd.Series(vals).corr(pd.Series(truth), method="spearman"))
        r, h, o, p = g
        rows.append(dict(radius=r, height=h, offset=o, depth_power=p,
                         pearson=pear, spearman=spear))
    out = pd.DataFrame(rows).sort_values("spearman", ascending=False)
    Path(args.out_csv).parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.out_csv, index=False)
    print(f"Wrote {args.out_csv}\n")

    inc = out[(out.radius == 16) & (out.height == 40) &
              (out.offset == -4) & (out.depth_power == 1)]
    print(f"Top {args.top} by Spearman:")
    print(out.head(args.top).to_string(index=False,
          formatters={"pearson": "{:+.3f}".format, "spearman": "{:+.3f}".format}))
    if len(inc):
        print("\nIncumbent (r16 h40 off-4 pow1):")
        print(inc.to_string(index=False,
              formatters={"pearson": "{:+.3f}".format, "spearman": "{:+.3f}".format}))


if __name__ == "__main__":
    main()
