#!/usr/bin/env python3
"""
scan_native_cylinder.py

Sweep the native-carve distance (exclude_dist) for the native-aware cylinder and
check whether carving the native-antigen footprint out of the cylinder improves it
as a clash proxy.

For each design we compute, once, the distance from every flagged scaffold CA (CA
inside the cylinder) to the nearest native-antigen heavy atom. The native-aware
count at a given exclude_dist is then just (distance > exclude_dist).sum(), so the
whole sweep is free after one pass over the structures. exclude_dist = 0 reproduces
the plain cylinder (nothing carved), so it is the baseline row.

For each exclude_dist we report, vs the true clash (af3_n_clash_res):
    Pearson, Spearman          correlation of the proxy with the true clash count
    AUC                        ranking power for "clash-free" (af3_n_clash_res==0)
    prec, rec                  precision/recall of clash-free designs if you keep the
                               lowest --op_frac of designs by the proxy
Reported twice: on all processed designs, and on the RMSD-gated subset (the
population the composite actually scores), where the cylinder is meaningful.

Depends on add_native_cylinder.py (same scripts/ dir) for the loaders + geometry.

    python scripts/scan_native_cylinder.py \
        --metrics_csv runs/run_rfd3_mpnn/04_filter/metrics_cylinder_full.csv \
        --dp2_parquet datasets/dp2.parquet \
        --native_dir  /tgen_labs/.../abdb/complex_pdbfiles/cleaned \
        --n_all 500 --n_gated 50 --gate 2.5 \
        --out_csv runs/run_rfd3_mpnn/04_filter/native_cyl_sweep.csv
"""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from scipy.stats import rankdata

import add_native_cylinder as N   # loaders + geometry (cylinder_frame, inside_cylinder, kabsch, ...)

EPI_RMSD_COLS = ["epitope_chunk_rmsd"]
STATUS_COLS   = ["status"]
CLASHST_COLS  = ["af3_clash_status"]


def design_distances(af3_dir, nat, epi_ris):
    """Return (plain_count, dist_to_nearest_native[inside], native_in) or None.
    `nat` is the preloaded native antigen tuple (or None). dist array has one entry
    per flagged (inside) scaffold residue."""
    a = N.load_af3(af3_dir)
    if a is None:
        return None
    ca, res_idx, seq = a
    if not epi_ris or max(epi_ris) >= len(ca):
        return None
    epi_set = set(epi_ris)
    epi_pos = [p for p, ri in enumerate(res_idx) if ri in epi_set]
    if len(epi_pos) < 3:
        return None
    epi_ca = ca[epi_pos]
    epi_seq = "".join(seq[p] for p in epi_pos)
    _, normal, base = N.cylinder_frame(epi_ca, ca)

    scaf_ca = ca[[p for p in range(len(ca)) if res_idx[p] not in epi_set]]
    inside_xyz = scaf_ca[N.inside_cylinder(scaf_ca, base, normal)]
    plain = len(inside_xyz)

    if nat is None:
        return None
    nca, nseq, nheavy, nresi = nat
    m = N.match_epitope(epi_seq, epi_ca, nseq, nca)   # gap-tolerant locate + align
    if m is None:
        return None
    R, t, nepi_set, _ = m
    nca_al = (R @ nca.T).T + t
    nheavy_al = (R @ nheavy.T).T + t
    nonepi = np.ones(len(nca), bool); nonepi[list(nepi_set)] = False
    native_in = int((N.inside_cylinder(nca_al, base, normal) & nonepi).sum())

    heavy_nonepi = nheavy_al[nonepi[nresi]]
    if plain == 0 or len(heavy_nonepi) == 0:
        return plain, np.empty(0), native_in
    d, _ = cKDTree(heavy_nonepi).query(inside_xyz, k=1)
    return plain, d, native_in


def auc_clashfree(metric, clashfree):
    """AUC for the proxy ranking clash-free designs low. Mann-Whitney; tie-safe."""
    pos = clashfree; neg = ~clashfree
    if pos.sum() == 0 or neg.sum() == 0:
        return float("nan")
    r = rankdata(metric)                          # low metric -> low rank (good for clash-free)
    # AUC = P(blocked ranks above clash-free) = how well HIGH metric flags blocked
    auc_blocked = (r[neg].sum() - neg.sum() * (neg.sum() + 1) / 2) / (pos.sum() * neg.sum())
    return auc_blocked

def pr_at(metric, clashfree, frac):
    """Keep the lowest `frac` of designs by metric; precision/recall of clash-free."""
    n = len(metric); k = max(1, int(round(frac * n)))
    keep = np.zeros(n, bool); keep[np.argsort(metric, kind="stable")[:k]] = True
    tp = int((keep & clashfree).sum())
    prec = tp / keep.sum() if keep.sum() else float("nan")
    rec = tp / clashfree.sum() if clashfree.sum() else float("nan")
    return prec, rec


def evaluate(rows, eds, op_frac):
    """rows: list of (plain, d_array, true_clash). Returns DataFrame over eds."""
    truth = np.array([r[2] for r in rows], float)
    clashfree = truth == 0
    out = []
    for ed in eds:
        metric = np.array([int((r[1] > ed).sum()) if len(r[1]) else r[0] for r in rows],
                          float)
        if ed == 0:                               # baseline = plain count
            metric = np.array([r[0] for r in rows], float)
        m = np.isfinite(metric) & np.isfinite(truth)
        if m.sum() < 5 or np.std(metric[m]) == 0:
            continue
        pe = float(np.corrcoef(metric[m], truth[m])[0, 1])
        sp = float(pd.Series(metric[m]).corr(pd.Series(truth[m]), method="spearman"))
        auc = auc_clashfree(metric[m], clashfree[m])
        prec, rec = pr_at(metric[m], clashfree[m], op_frac)
        out.append(dict(exclude_dist=ed, n=int(m.sum()), pearson=pe, spearman=sp,
                        auc=auc, precision=prec, recall=rec))
    return pd.DataFrame(out)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--metrics_csv", required=True)
    ap.add_argument("--dp2_parquet", required=True)
    ap.add_argument("--native_dir", required=True)
    ap.add_argument("--out_csv", required=True)
    ap.add_argument("--n_all", type=int, default=500)
    ap.add_argument("--n_gated", type=int, default=50)
    ap.add_argument("--gate", type=float, default=2.5, help="epitope_chunk_rmsd <")
    ap.add_argument("--exclude_dists", type=float, nargs="+",
                    default=[0, 1, 2, 3, 4, 5, 6, 7, 8])
    ap.add_argument("--op_frac", type=float, default=0.25,
                    help="keep this lowest fraction by proxy for precision/recall")
    ap.add_argument("--shuffle", action="store_true",
                    help="shuffle across antigens before sampling (the metrics file is id-sorted)")
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()

    df = pd.read_csv(args.metrics_csv, low_memory=False)
    sc = N.first_present(df.columns, STATUS_COLS)
    cc = N.first_present(df.columns, CLASHST_COLS)
    if sc: df = df[df[sc].astype(str) == "ok"]
    if cc: df = df[df[cc].astype(str) == "ok"]
    if args.shuffle:
        df = df.sample(frac=1.0, random_state=args.seed).reset_index(drop=True)
    tok = N.first_present(df.columns, N.TOKEN_COLS)
    dirc = N.first_present(df.columns, N.AF3_DIR_COLS)
    idc = N.first_present(df.columns, N.ID_COLS)
    wsc = N.first_present(df.columns, N.WS_COLS)
    tcc = N.first_present(df.columns, N.TRUE_CLASH)
    erc = N.first_present(df.columns, EPI_RMSD_COLS)
    if tcc is None or erc is None:
        raise SystemExit("need af3_n_clash_res and epitope_chunk_rmsd columns")

    dp2 = pd.read_parquet(args.dp2_parquet)
    dp2["assay_scaffolded_epitope_id"] = dp2["assay_scaffolded_epitope_id"].astype(str).str.lower()
    epic = N.first_present(dp2.columns, N.DP2_EPI_COLS)
    epim = N.first_present(dp2.columns, N.DP2_EPI_MPNN)
    look = dp2.drop_duplicates("assay_scaffolded_epitope_id").set_index("assay_scaffolded_epitope_id")

    # pick the two evaluation sets (gated subset can overlap the full set)
    all_rows = df.head(args.n_all)
    gated_rows = df[pd.to_numeric(df[erc], errors="coerce") < args.gate].head(args.n_gated)
    wanted = pd.concat([all_rows, gated_rows]).drop_duplicates(subset=[dirc])
    print(f"processing {len(wanted)} unique designs "
          f"(all={len(all_rows)}, gated<{args.gate}={len(gated_rows)})")

    native_index = {p.stem.lower(): p for p in Path(args.native_dir).glob("*.pdb")}
    nat_cache = {}
    cache = {}; native_ins = []; n_ok = n_fail = 0
    for _, row in wanted.iterrows():
        t = str(row[tok]).lower()
        if t not in look.index or pd.isna(row.get(dirc)):
            n_fail += 1; continue
        dr = look.loc[t]
        if epic is not None:
            epi = N.parse_index_list(dr[epic])
        else:
            ws = int(row[wsc]) if (wsc and pd.notna(row.get(wsc))) else 0
            epi = [ws + i for i in N.parse_index_list(dr[epim])]
        npdb = native_index.get(str(row[idc]).lower()) if idc else None
        if npdb is None:
            nat = None
        else:
            key = str(npdb)
            if key not in nat_cache:
                nat_cache[key] = N.load_native_antigen(npdb)
            nat = nat_cache[key]
        res = design_distances(str(row[dirc]), nat, epi)
        if res is None:
            n_fail += 1; continue
        plain, d, native_in = res
        cache[str(row[dirc])] = (plain, d,
                                 float(pd.to_numeric(row[tcc], errors="coerce")))
        native_ins.append(native_in)
        n_ok += 1
        if n_ok % 100 == 0:
            print(f"  processed {n_ok}  (fail={n_fail})")

    if native_ins:
        nv = np.array(native_ins, float)
        print(f"\nnative_in_cylinder over {len(nv)} designs: "
              f"mean {nv.mean():.1f}  median {np.median(nv):.0f}  "
              f"frac>0 {(nv>0).mean():.2f}  max {nv.max():.0f}")

    def rows_for(sub):
        out = []
        for _, row in sub.iterrows():
            key = str(row.get(dirc))
            if key in cache and np.isfinite(cache[key][2]):
                out.append(cache[key])
        return out

    results = []
    for name, sub in [("ALL", all_rows), (f"GATED(<{args.gate})", gated_rows)]:
        rows = rows_for(sub)
        if len(rows) < 5:
            print(f"\n[{name}] too few usable designs ({len(rows)})"); continue
        tab = evaluate(rows, args.exclude_dists, args.op_frac)
        tab.insert(0, "set", name)
        results.append(tab)
        base = tab[tab.exclude_dist == 0]
        print(f"\n=== {name}  (n={len(rows)}) ===")
        print(tab.to_string(index=False, formatters={
            "pearson": "{:+.3f}".format, "spearman": "{:+.3f}".format,
            "auc": "{:.3f}".format, "precision": "{:.3f}".format,
            "recall": "{:.3f}".format}))
        if len(base):
            b = base.iloc[0]
            best = tab.loc[tab.spearman.idxmax()]
            if best.exclude_dist != 0 and best.spearman > b.spearman + 1e-9:
                print(f"  -> best Spearman at exclude_dist={best.exclude_dist:g} "
                      f"({best.spearman:+.3f} vs plain {b.spearman:+.3f})")
            else:
                print(f"  -> no exclude_dist beats the plain cylinder "
                      f"(plain Spearman {b.spearman:+.3f})")

    if results:
        out = pd.concat(results, ignore_index=True)
        Path(args.out_csv).parent.mkdir(parents=True, exist_ok=True)
        out.to_csv(args.out_csv, index=False)
        print(f"\nWrote {args.out_csv}")


if __name__ == "__main__":
    main()
