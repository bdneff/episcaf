#!/usr/bin/env python3
"""Sweep the cylinder base OFFSET and see which value best predicts the real antibody clash.

The cylinder base sits OFFSET A from the epitope centroid along the approach normal (current
value -4, i.e. 4 A *below* the epitope). The 8pww probe (scripts/cylinder_fp_probe.py) showed
that at -4 the cylinder scoops up scaffold hugging the epitope BELOW where the paratope sits,
inflating the count without real steric clashes. This recomputes the cylinder count at a range
of offsets IN ONE PASS over the DP3 structures (the epitope frame is computed once per design;
only the base shifts) and reports, per offset, the AUC for predicting clash-free
(af3_n_clash_res==0). The best offset is the one that most cleanly separates real clashes --
chosen on the whole DP3 set, not on 8pww. Runs on Gemini.

Usage (sample first for speed, then full):
  conda activate ~/rfd3/env/rfd3_py312
  python3 scripts/cylinder_offset_sweep.py \
      --metrics_csv <run>/04_filter/metrics_native_cyl_full.csv \
      --dp2_parquet <datasets>/dp2.parquet \
      --native_dir  <abdb cleaned complex dir> \
      --offsets -6,-4,-2,0,2,4,6 --limit 20000 \
      --out_csv results/cylinder_offset_sweep.csv
"""
from __future__ import annotations
import argparse
from pathlib import Path
import sys
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from scipy.stats import rankdata

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE)); sys.path.insert(0, str(HERE.parent))
from dp3_native_cylinder import (  # identical geometry; only the base offset varies  # noqa: E402
    load_af3, load_native_antigen, inside_cylinder, match_epitope, parse_index_list,
    first_present, TOKEN_COLS, AF3_DIR_COLS, ID_COLS, DP2_EPI_COLS, TRUE_CLASH, RADIUS, HEIGHT)


def counts_at_offsets(af3_dir, native_pdb, epi_ris, offsets, exclude_dist):
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

    # native antigen heavy atoms for the native-aware carve (frame-independent of offset)
    heavy_nonepi = None
    nat = load_native_antigen(Path(native_pdb)) if native_pdb else None
    if nat is not None:
        nca, nseq, nheavy, nresi = nat
        m = match_epitope(epi_seq, epi_ca, nseq, nca)
        if m is not None:
            R, t, nepi_set, _ = m
            nheavy_al = (R @ nheavy.T).T + t
            nonepi = np.ones(len(nca), bool); nonepi[list(nepi_set)] = False
            heavy_nonepi = nheavy_al[nonepi[nresi]]
    tree = cKDTree(heavy_nonepi) if (heavy_nonepi is not None and len(heavy_nonepi)) else None

    out = {}
    for off in offsets:
        base = centroid + off * normal
        ins = inside_cylinder(scaf_ca, base, normal, r=RADIUS, h=HEIGHT)
        n_plain = int(ins.sum())
        n_aware = n_plain
        if tree is not None and n_plain > 0:
            d, _ = tree.query(scaf_ca[ins], k=1)
            n_aware = int((d > exclude_dist).sum())
        out[f"cyl_plain_{off:g}"] = n_plain
        out[f"cyl_aware_{off:g}"] = n_aware
    return out


def auc_low_is_pos(score, y):
    s = pd.to_numeric(pd.Series(score), errors="coerce").to_numpy()
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
    ap.add_argument("--offsets", default="-6,-4,-2,0,2,4,6")
    ap.add_argument("--exclude_dist", type=float, default=1.0)
    ap.add_argument("--limit", type=int, default=0, help="random sample this many designs (0=all)")
    ap.add_argument("--out_csv", required=True)
    args = ap.parse_args()
    offsets = [float(x) for x in args.offsets.split(",")]

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

    rows = []; n_ok = n_fail = 0
    for pos, row in enumerate(df.itertuples(index=False)):
        t = str(getattr(row, tok)).lower()
        if t not in look.index or pd.isna(getattr(row, dirc)):
            n_fail += 1; continue
        epi = parse_index_list(look.get(t))
        npdb = native_index.get(str(getattr(row, idc)).lower())
        res = counts_at_offsets(str(getattr(row, dirc)), npdb, epi, offsets, args.exclude_dist)
        if res is None:
            n_fail += 1; continue
        res["af3_n_clash_res"] = pd.to_numeric(getattr(row, tcc), errors="coerce")
        res["id"] = getattr(row, idc)
        rows.append(res); n_ok += 1
        if (pos + 1) % 2000 == 0:
            print(f"  {pos+1}/{len(df)} ok={n_ok} fail={n_fail}", flush=True)

    out = pd.DataFrame(rows)
    Path(args.out_csv).parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.out_csv, index=False)
    print(f"\nwrote {args.out_csv}  ok={n_ok} fail={n_fail}")

    y = (pd.to_numeric(out.af3_n_clash_res, errors="coerce") == 0).astype(int).to_numpy()
    print(f"clash-free {int(y.sum()):,}/{len(y):,} ({100*y.mean():.1f}%)")
    print(f"\n{'offset':>7}  {'AUC plain':>10}  {'AUC native-aware':>16}   (current default is -4)")
    for off in offsets:
        ap_ = auc_low_is_pos(out[f"cyl_plain_{off:g}"], y)
        aa_ = auc_low_is_pos(out[f"cyl_aware_{off:g}"], y)
        mark = "  <- current" if off == -4 else ""
        print(f"{off:+7.0f}  {ap_:10.3f}  {aa_:16.3f}{mark}")


if __name__ == "__main__":
    main()
