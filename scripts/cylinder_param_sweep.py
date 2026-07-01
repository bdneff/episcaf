#!/usr/bin/env python3
"""Sweep ALL cylinder geometry parameters and see which choices best predict the real clash.

The cylinder has four geometry dials, all currently hard-set with no recorded provenance or
sweep (see docs/CYLINDER_PARAMS.md): base OFFSET (-4 A from the epitope centroid along the
approach normal), RADIUS (16 A), HEIGHT (40 A), and the native-aware carve distance exclude_dist
(1.0 A in production). This recomputes the cylinder count over a grid of (offset, radius, height,
exclude_dist) IN ONE PASS over the DP3 structures -- the epitope frame is computed once per
design; every grid cell is a cheap in/out test on the same scaffold coordinates -- and reports,
per cell, the AUC for predicting clash-free (af3_n_clash_res==0). We then choose the geometry on
the whole DP3 ground truth rather than inheriting it. Runs on Gemini.

Usage (sample first for speed, then full):
  conda activate ~/rfd3/env/rfd3_py312
  python3 scripts/cylinder_param_sweep.py \
      --metrics_csv <run>/04_filter/metrics_native_cyl_full.csv \
      --dp2_parquet <datasets>/dp2.parquet \
      --native_dir  <abdb cleaned complex dir> \
      --offsets -6,-4,-2,0,2,4 --radii 12,14,16,18,20 --heights 30,40,50 \
      --exclude_dists 1.0 --limit 20000 \
      --out_csv results/cylinder_param_sweep.csv
"""
from __future__ import annotations
import argparse
from collections import defaultdict
from pathlib import Path
import sys
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from scipy.stats import rankdata

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE)); sys.path.insert(0, str(HERE.parent))
from dp3_native_cylinder import (  # identical geometry; only the dials vary  # noqa: E402
    load_af3, load_native_antigen, inside_cylinder, match_epitope, parse_index_list,
    first_present, TOKEN_COLS, AF3_DIR_COLS, ID_COLS, DP2_EPI_COLS, TRUE_CLASH)


def counts_grid(af3_dir, native_pdb, epi_ris, offsets, radii, heights, exclude_dists):
    a = load_af3(af3_dir)
    if a is None:
        return None
    ca, res_idx, seq = a
    epi_set = set(epi_ris)
    if not epi_ris or max(epi_ris) >= len(ca):
        return None
    epi_pos = [p for p, ri in enumerate(res_idx) if ri in epi_set]
    if len(epi_pos) < 3:
        return None
    epi_ca = ca[epi_pos]; epi_seq = "".join(seq[p] for p in epi_pos)
    centroid = epi_ca.mean(0)
    _, _, Vt = np.linalg.svd(epi_ca - centroid); normal = Vt[-1]
    if np.dot(normal, ca.mean(0) - centroid) > 0:
        normal = -normal
    scaf_ca = ca[[p for p in range(len(ca)) if res_idx[p] not in epi_set]]

    # native antigen heavy atoms (aligned once) for the native-aware carve
    tree = None
    nat = load_native_antigen(Path(native_pdb)) if native_pdb else None
    if nat is not None:
        nca, nseq, nheavy, nresi = nat
        m = match_epitope(epi_seq, epi_ca, nseq, nca)
        if m is not None:
            R, t, nepi_set, _ = m
            nheavy_al = (R @ nheavy.T).T + t
            nonepi = np.ones(len(nca), bool); nonepi[list(nepi_set)] = False
            hv = nheavy_al[nonepi[nresi]]
            if len(hv):
                tree = cKDTree(hv)

    out = {}
    for rad in radii:
        for hgt in heights:
            for off in offsets:
                base = centroid + off * normal
                ins = inside_cylinder(scaf_ca, base, normal, r=rad, h=hgt)
                n_plain = int(ins.sum())
                out[(off, rad, hgt, "plain")] = n_plain
                if tree is not None and n_plain > 0:
                    d, _ = tree.query(scaf_ca[ins], k=1)
                    for exd in exclude_dists:
                        out[(off, rad, hgt, f"aware{exd:g}")] = int((d > exd).sum())
                else:
                    for exd in exclude_dists:
                        out[(off, rad, hgt, f"aware{exd:g}")] = n_plain
    return out


def auc_low_is_pos(score, y):
    s = np.asarray(score, float)
    m = ~np.isnan(s)
    r = rankdata(-s[m]); yy = y[m]; npos = yy.sum(); nneg = len(yy) - npos
    if npos == 0 or nneg == 0:
        return float("nan")
    return (r[yy == 1].sum() - npos * (npos + 1) / 2) / (npos * nneg)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--metrics_csv", required=True)
    ap.add_argument("--dp2_parquet", required=True)
    ap.add_argument("--native_dir", required=True)
    ap.add_argument("--offsets", default="-6,-4,-2,0,2,4")
    ap.add_argument("--radii", default="12,14,16,18,20")
    ap.add_argument("--heights", default="30,40,50")
    ap.add_argument("--exclude_dists", default="1.0")
    ap.add_argument("--limit", type=int, default=0, help="random sample this many designs (0=all)")
    ap.add_argument("--out_csv", required=True)
    args = ap.parse_args()
    offsets = [float(x) for x in args.offsets.split(",")]
    radii = [float(x) for x in args.radii.split(",")]
    heights = [float(x) for x in args.heights.split(",")]
    exclude_dists = [float(x) for x in args.exclude_dists.split(",")]

    df = pd.read_csv(args.metrics_csv, low_memory=False)
    if args.limit > 0 and args.limit < len(df):
        df = df.sample(args.limit, random_state=0).reset_index(drop=True)
    tok = first_present(df.columns, TOKEN_COLS); dirc = first_present(df.columns, AF3_DIR_COLS)
    idc = first_present(df.columns, ID_COLS); tcc = first_present(df.columns, TRUE_CLASH)

    dp2 = pd.read_parquet(args.dp2_parquet)
    dp2["assay_scaffolded_epitope_id"] = dp2.assay_scaffolded_epitope_id.astype(str).str.lower()
    epic = first_present(dp2.columns, DP2_EPI_COLS)
    look = (dp2.drop_duplicates("assay_scaffolded_epitope_id")
               .set_index("assay_scaffolded_epitope_id")[epic])
    native_index = {p.stem.lower(): p for p in Path(args.native_dir).glob("*.pdb")}

    acc = defaultdict(list); yl = []; n_ok = n_fail = 0
    for pos, row in enumerate(df.itertuples(index=False)):
        t = str(getattr(row, tok)).lower()
        if t not in look.index or pd.isna(getattr(row, dirc)):
            n_fail += 1; continue
        g = counts_grid(str(getattr(row, dirc)), native_index.get(str(getattr(row, idc)).lower()),
                        parse_index_list(look.get(t)), offsets, radii, heights, exclude_dists)
        if g is None:
            n_fail += 1; continue
        for k, v in g.items():
            acc[k].append(v)
        yl.append(1 if pd.to_numeric(getattr(row, tcc), errors="coerce") == 0 else 0)
        n_ok += 1
        if (pos + 1) % 2000 == 0:
            print(f"  {pos+1}/{len(df)} ok={n_ok} fail={n_fail}", flush=True)

    y = np.asarray(yl)
    print(f"\nn={len(y):,}  clash-free {int(y.sum()):,} ({100*y.mean():.1f}%)  grid cells = {len(acc)}")
    rows = []
    for (off, rad, hgt, variant), counts in acc.items():
        c = np.asarray(counts, float)
        rows.append(dict(offset=off, radius=rad, height=hgt, variant=variant,
                         auc=auc_low_is_pos(c, y),
                         med_free=float(np.median(c[y == 1])) if (y == 1).any() else np.nan,
                         med_clash=float(np.median(c[y == 0])) if (y == 0).any() else np.nan))
    res = pd.DataFrame(rows).sort_values("auc", ascending=False)
    Path(args.out_csv).parent.mkdir(parents=True, exist_ok=True)
    res.to_csv(args.out_csv, index=False)
    print(f"wrote {args.out_csv}\n")
    cur = res[(res.offset == -4) & (res.radius == 16) & (res.height == 40)]
    print("current default (offset -4, R 16, H 40):")
    print(cur.to_string(index=False))
    print("\ntop 12 cells by AUC:")
    print(res.head(12).to_string(index=False))


if __name__ == "__main__":
    main()
