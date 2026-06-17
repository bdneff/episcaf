#!/usr/bin/env python3
"""
add_cone_clash_metric.py

A better clash proxy than the single cylinder: instead of one approach direction
(the PCA normal), sweep a fan of directions across the outward hemisphere (within
--max_cone_deg of the normal) and, for each, count scaffold residues inside the
cylinder. Accessibility is about whether SOME approach is clear, so we summarize
the sweep several ways and let the data say which is the best proxy.

Adds columns:
    cone_normal_clash : clash along the PCA normal (== the original cylinder)
    cone_min_clash    : best (least-blocked) approach over all directions
    cone_p10_clash    : 10th-percentile clash (robust "is there a clear-ish way in")
    cone_mean_clash   : mean clash over directions
    cone_clear_frac   : fraction of directions that are fully clear (clash == 0)

Locked cylinder geometry (matches the validated scan): r=16, h=40, offset=-4
(centroid convention). Cone defaults: 60 deg half-angle, 4 rings x 8 azimuths
(+ the normal) = 33 directions.

At the end it prints the correlation of every summary vs the true clash count
(af3_n_clash_res), next to cone_normal_clash, so you can see immediately whether
the sweep beats the plain cylinder (which was ~0.34).

Example
-------
    python scripts/add_cone_clash_metric.py \
        --metrics_csv runs/run_rfd3_mpnn/04_filter/metrics_cylinder_full.csv \
        --dp2_parquet datasets/dp2.parquet \
        --out_csv     runs/run_rfd3_mpnn/04_filter/metrics_cone.csv \
        --limit 1000
"""

import argparse
import gzip
import math
from pathlib import Path
from typing import List, Optional, Tuple

import gemmi
import numpy as np
import pandas as pd

# locked cylinder geometry
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


def get_chainA(st):
    m = st[0]
    for ch in m:
        if ch.name == "A":
            return ch
    return m[0]


def find_af3_cif(af3_dir: Path) -> Optional[Path]:
    for pat in ("*_model.cif", "*_model.cif.gz"):
        h = next(af3_dir.glob(pat), None)
        if h:
            return h
    for pat in ("model.cif", "model.cif.gz"):
        h = next(af3_dir.rglob(pat), None)
        if h:
            return h
    return None


def load_af3_ca(af3_dir: Path):
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
    return np.asarray(coords, float), np.asarray(ridx, int)


def sample_cone_directions(normal: np.ndarray, max_cone_deg: float,
                           rings: int, n_az: int) -> np.ndarray:
    """Unit directions within max_cone_deg of `normal` (outward hemisphere).
    Includes the normal itself, then `rings` concentric rings of `n_az` each."""
    n = normal / np.linalg.norm(normal)
    a = np.array([1.0, 0.0, 0.0]) if abs(n[0]) < 0.9 else np.array([0.0, 1.0, 0.0])
    u = np.cross(n, a); u /= np.linalg.norm(u)
    v = np.cross(n, u)
    dirs = [n.copy()]
    for th in np.linspace(0.0, math.radians(max_cone_deg), rings + 1)[1:]:
        for k in range(n_az):
            phi = 2 * math.pi * k / n_az
            d = (math.cos(th) * n
                 + math.sin(th) * (math.cos(phi) * u + math.sin(phi) * v))
            dirs.append(d / np.linalg.norm(d))
    return np.asarray(dirs)


def clash_count_along(ca, res_idx, epi_set, centroid, direction):
    base = centroid + OFFSET * direction
    w = ca - base
    proj = w @ direction
    perp = w - np.outer(proj, direction)
    dist = np.linalg.norm(perp, axis=1)
    inside = (proj >= 0.0) & (proj <= HEIGHT) & (dist <= RADIUS)
    return sum(1 for r, ins in zip(res_idx, inside)
               if ins and int(r) not in epi_set)


def cone_summaries(af3_dir, epi_ris, max_cone_deg, rings, n_az):
    loaded = load_af3_ca(Path(af3_dir))
    if loaded is None:
        return None
    ca, res_idx = loaded
    if not epi_ris or max(epi_ris) >= (int(res_idx.max()) + 1):
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

    dirs = sample_cone_directions(normal, max_cone_deg, rings, n_az)
    counts = np.array([clash_count_along(ca, res_idx, epi_set, centroid, d)
                       for d in dirs], dtype=float)
    return dict(
        cone_normal_clash=float(counts[0]),         # direction 0 == the normal
        cone_min_clash=float(counts.min()),
        cone_p10_clash=float(np.percentile(counts, 10)),
        cone_mean_clash=float(counts.mean()),
        cone_clear_frac=float((counts == 0).mean()),
    )


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--metrics_csv", required=True)
    ap.add_argument("--dp2_parquet", required=True)
    ap.add_argument("--out_csv", required=True)
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--max_cone_deg", type=float, default=60.0)
    ap.add_argument("--rings", type=int, default=4)
    ap.add_argument("--n_az", type=int, default=8)
    args = ap.parse_args()

    ndir = 1 + args.rings * args.n_az
    print(f"Cone sweep: max_cone={args.max_cone_deg} deg, {ndir} directions "
          f"(geometry r={RADIUS} h={HEIGHT} off={OFFSET})")

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

    cols = ["cone_normal_clash", "cone_min_clash", "cone_p10_clash",
            "cone_mean_clash", "cone_clear_frac"]
    acc = {c: np.full(len(df), np.nan) for c in cols}
    n_ok = n_fail = 0
    for pos, (_, row) in enumerate(df.iterrows()):
        tok = str(row[tok_col]).lower()
        if tok not in dp2_lookup.index or pd.isna(row.get(dir_col)):
            n_fail += 1
            continue
        drow = dp2_lookup.loc[tok]
        if epi_col is not None:
            epi_ris = parse_index_list(drow[epi_col])
        else:
            ws = int(row[ws_col]) if (ws_col and pd.notna(row.get(ws_col))) else 0
            epi_ris = [ws + i for i in parse_index_list(drow[epi_mpnn])]
        res = cone_summaries(str(row[dir_col]), epi_ris,
                             args.max_cone_deg, args.rings, args.n_az)
        if res is None:
            n_fail += 1
        else:
            for c in cols:
                acc[c][pos] = res[c]
            n_ok += 1
        if (pos + 1) % 1000 == 0:
            print(f"  {pos+1:,}  ok={n_ok:,} fail={n_fail:,}")

    for c in cols:
        df[c] = acc[c]
    out = Path(args.out_csv)
    out.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out, index=False)
    print(f"\nWrote {out}  (ok={n_ok:,} fail={n_fail:,})")

    # head-to-head vs the true clash count
    tc = first_present(df.columns, TRUE_CLASH)
    if tc is not None:
        truth = pd.to_numeric(df[tc], errors="coerce")
        print(f"\nCorrelation of each proxy vs true clash ({tc}):")
        print("  proxy                 Pearson   Spearman")
        ref_cols = cols + (["cylinder_ca_clashes"]
                           if "cylinder_ca_clashes" in df.columns else [])
        for c in ref_cols:
            p = pd.to_numeric(df[c], errors="coerce")
            m = p.notna() & truth.notna()
            if m.sum() < 10:
                continue
            pear = p[m].corr(truth[m])
            spear = p[m].corr(truth[m], method="spearman")
            tag = "  <- original" if c == "cylinder_ca_clashes" else ""
            print(f"  {c:20s}  {pear:+.3f}    {spear:+.3f}{tag}")
        print("\n(higher |corr| = better proxy; clear_frac should be NEGATIVELY "
              "correlated since more clear directions = less true clash)")


if __name__ == "__main__":
    main()
