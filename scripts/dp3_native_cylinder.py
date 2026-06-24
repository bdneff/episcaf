#!/usr/bin/env python3
"""
dp3_native_cylinder.py -- compute cylinder_native_aware (the native-antigen carve) for the
FULL DP3 set, so Figure 7 / the overlays can show plain AND native-aware side by side for
both pools (the DP4 12-mer table already has native_aware; DP3 did not, at full coverage).

This is a faithful, sped-up promotion of archive/scripts/add_native_cylinder.py: the geometry
and matching are IDENTICAL (it reproduces the existing metrics_native_cyl.csv numbers), with
two changes -- (1) the native complex is loaded ONCE per epitope id (lru_cache) instead of once
per design (the archived script reloaded ~150k times), and (2) the viz code is dropped (not
needed for a batch recompute). Run WITHOUT --gate so every AF3-resolved design gets a value;
the original 44,862-design file was produced with --gate 2.5 (only designs below the gate),
which is why DP3 native_aware was ~30% populated.

For each design: build the locked cylinder over the AF3 epitope; load the native complex
(<id>.pdb, antigen = chain A); locate the epitope in the antigen by sequence and Kabsch-align;
then count scaffold CAs inside the cylinder (plain) and those NOT within --exclude_dist of a
native non-epitope heavy atom (native-aware carve).

  python3 scripts/dp3_native_cylinder.py \
      --metrics_csv <run>/04_filter/metrics_cylinder_full.csv \
      --dp2_parquet <datasets>/dp2.parquet \
      --native_dir  <abdb cleaned dir> \
      --out_csv     <run>/04_filter/metrics_native_cyl_full.csv \
      --exclude_dist 1.0
"""
import argparse, difflib, gzip, math
from functools import lru_cache
from pathlib import Path
from typing import List

import gemmi
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

RADIUS, HEIGHT, OFFSET = 16.0, 40.0, -4.0

TOKEN_COLS = ["token", "assay_scaffolded_epitope_id"]; AF3_DIR_COLS = ["af3_dir"]
ID_COLS = ["id"]; DP2_EPI_COLS = ["assay_scaffolded_epitope_chunk_resindices"]
DP2_EPI_MPNN = ["scaffolded_epitope_chunk_resindices"]; WS_COLS = ["af3_window_start"]
TRUE_CLASH = ["af3_n_clash_res"]; EPI_RMSD_COLS = ["epitope_chunk_rmsd"]

THREE2ONE = {
 'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E','GLY':'G',
 'HIS':'H','ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S',
 'THR':'T','TRP':'W','TYR':'Y','VAL':'V'}


def parse_index_list(x) -> List[int]:
    if x is None or (isinstance(x, float) and math.isnan(x)): return []
    if isinstance(x, (list, tuple, np.ndarray)): return [int(i) for i in x]
    s = str(x).replace("[", "").replace("]", "").replace(",", " ")
    out = []
    for t in s.split():
        try: out.append(int(t))
        except ValueError: pass
    return out

def first_present(cols, cands):
    for c in cands:
        if c in cols: return c
    return None

def read_gemmi(p: Path):
    s = str(p)
    if s.endswith(".gz"):
        with gzip.open(p, "rt") as f: doc = gemmi.cif.read_string(f.read())
        return gemmi.make_structure_from_block(doc.sole_block())
    if p.suffix.lower() == ".cif":
        return gemmi.make_structure_from_block(gemmi.cif.read(str(p)).sole_block())
    return gemmi.read_structure(str(p))

def chain_by_name(model, name):
    for ch in model:
        if ch.name == name: return ch
    return None

def find_af3_cif(af3_dir: Path):
    for pat in ("*_model.cif", "*_model.cif.gz", "model.cif", "model.cif.gz"):
        h = next(af3_dir.glob(pat), None) or next(af3_dir.rglob(pat), None)
        if h: return h
    return None

def load_af3(af3_dir: Path):
    cif = find_af3_cif(Path(af3_dir))
    if cif is None: return None
    try: st = read_gemmi(cif)
    except Exception: return None
    ch = chain_by_name(st[0], "A") or st[0][0]
    ca, ridx, seq = [], [], []
    for i, res in enumerate(ch):
        a = res.find_atom("CA", altloc="*")
        if a:
            ca.append([a.pos.x, a.pos.y, a.pos.z]); ridx.append(i)
            seq.append(THREE2ONE.get(res.name.upper(), "X"))
    if not ca: return None
    return np.asarray(ca, float), np.asarray(ridx, int), "".join(seq)

@lru_cache(maxsize=128)
def load_native_antigen(pdb: Path):
    """Cached per native pdb (≈59 unique vs ~150k designs)."""
    try: st = read_gemmi(pdb)
    except Exception: return None
    ch = chain_by_name(st[0], "A")
    if ch is None: return None
    ca, seq, heavy, heavy_resi = [], [], [], []
    for ri, res in enumerate(ch):
        a = res.find_atom("CA", altloc="*")
        if not a: continue
        ca.append([a.pos.x, a.pos.y, a.pos.z]); seq.append(THREE2ONE.get(res.name.upper(), "X"))
        ridx = len(ca) - 1
        for at in res:
            if at.element == gemmi.Element("H"): continue
            heavy.append([at.pos.x, at.pos.y, at.pos.z]); heavy_resi.append(ridx)
    if len(ca) < 5: return None
    return (np.asarray(ca, float), "".join(seq),
            np.asarray(heavy, float), np.asarray(heavy_resi, int))


def cylinder_frame(epi_ca, all_ca):
    centroid = epi_ca.mean(0)
    _, _, Vt = np.linalg.svd(epi_ca - centroid); normal = Vt[-1]
    if np.dot(normal, all_ca.mean(0) - centroid) > 0: normal = -normal
    return centroid, normal, centroid + OFFSET * normal

def inside_cylinder(pts, base, normal, r=RADIUS, h=HEIGHT):
    v = pts - base; proj = v @ normal
    dist = np.linalg.norm(v - np.outer(proj, normal), axis=1)
    return (proj >= 0) & (proj <= h) & (dist <= r)

def kabsch(P, Q):
    Pc, Qc = P - P.mean(0), Q - Q.mean(0)
    U, _, Vt = np.linalg.svd(Pc.T @ Qc)
    d = np.sign(np.linalg.det(Vt.T @ U.T)); D = np.diag([1, 1, d])
    R = Vt.T @ D @ U.T
    return R, Q.mean(0) - R @ P.mean(0)

def match_epitope(epi_seq, epi_ca, nseq, nca, min_match=5, max_rmsd=8.0):
    sm = difflib.SequenceMatcher(None, epi_seq, nseq, autojunk=False)
    a, b, size = sm.find_longest_match(0, len(epi_seq), 0, len(nseq))
    if size < min_match: return None
    P = nca[b:b + size]; Q = epi_ca[a:a + size]
    R, t = kabsch(P, Q)
    rmsd = float(np.sqrt((((R @ P.T).T + t - Q) ** 2).sum(1).mean()))
    if rmsd > max_rmsd: return None
    off = b - a
    nepi = {off + k for k in range(len(epi_ca)) if 0 <= off + k < len(nca)}
    return R, t, nepi, rmsd


def process(af3_dir, native_pdb, epi_ris, exclude_dist):
    a = load_af3(af3_dir)
    if a is None: return None
    ca, res_idx, seq = a
    if not epi_ris or max(epi_ris) >= len(ca): return None
    epi_pos = [p for p, ri in enumerate(res_idx) if ri in set(epi_ris)]
    if len(epi_pos) < 3: return None
    epi_ca = ca[epi_pos]; epi_seq = "".join(seq[p] for p in epi_pos)
    centroid, normal, base = cylinder_frame(epi_ca, ca)

    scaf_pos = [p for p in range(len(ca)) if res_idx[p] not in set(epi_ris)]
    scaf_ca = ca[scaf_pos]
    ins_scaf = inside_cylinder(scaf_ca, base, normal)
    n_plain = int(ins_scaf.sum())

    nat = load_native_antigen(Path(native_pdb)) if native_pdb else None
    native_in = 0; n_aware = n_plain
    if nat is not None:
        nca, nseq, nheavy, nresi = nat
        m = match_epitope(epi_seq, epi_ca, nseq, nca)
        if m is not None:
            R, t, nepi_set, _ = m
            nca_al = (R @ nca.T).T + t
            nheavy_al = (R @ nheavy.T).T + t
            nonepi = np.ones(len(nca), bool); nonepi[list(nepi_set)] = False
            nin = inside_cylinder(nca_al, base, normal) & nonepi
            native_in = int(nin.sum())
            heavy_nonepi = nheavy_al[nonepi[nresi]]
            inside_xyz = scaf_ca[ins_scaf]
            if n_plain > 0 and len(heavy_nonepi):
                d, _ = cKDTree(heavy_nonepi).query(inside_xyz, k=1)
                n_aware = int((d > exclude_dist).sum())
    return n_plain, native_in, n_aware


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--metrics_csv", required=True)
    ap.add_argument("--dp2_parquet", required=True)
    ap.add_argument("--native_dir", required=True)
    ap.add_argument("--out_csv", required=True)
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--exclude_dist", type=float, default=1.0)
    ap.add_argument("--gate", type=float, default=0.0,
                    help="only process designs with epitope_chunk_rmsd < gate (0 = ALL; use 0 here)")
    args = ap.parse_args()

    df = pd.read_csv(args.metrics_csv, low_memory=False)
    if args.limit > 0: df = df.head(args.limit).copy()
    tok = first_present(df.columns, TOKEN_COLS); dirc = first_present(df.columns, AF3_DIR_COLS)
    idc = first_present(df.columns, ID_COLS); wsc = first_present(df.columns, WS_COLS)
    erc = first_present(df.columns, EPI_RMSD_COLS)

    dp2 = pd.read_parquet(args.dp2_parquet)
    dp2["assay_scaffolded_epitope_id"] = dp2["assay_scaffolded_epitope_id"].astype(str).str.lower()
    epic = first_present(dp2.columns, DP2_EPI_COLS); epim = first_present(dp2.columns, DP2_EPI_MPNN)
    look = dp2.drop_duplicates("assay_scaffolded_epitope_id").set_index("assay_scaffolded_epitope_id")

    native_index = {p.stem.lower(): p for p in Path(args.native_dir).glob("*.pdb")}

    plain = (pd.to_numeric(df["cylinder_ca_clashes"], errors="coerce").to_numpy()
             if "cylinder_ca_clashes" in df.columns else np.full(len(df), np.nan))
    natin = np.full(len(df), np.nan); aware = np.full(len(df), np.nan)
    n_ok = n_fail = n_skip = 0
    for pos, (_, row) in enumerate(df.iterrows()):
        t = str(row[tok]).lower()
        if t not in look.index or pd.isna(row.get(dirc)): n_fail += 1; continue
        if args.gate > 0:
            erm = pd.to_numeric(row.get(erc), errors="coerce")
            if not (pd.notna(erm) and erm < args.gate): n_skip += 1; continue
        dr = look.loc[t]
        if epic is not None: epi = parse_index_list(dr[epic])
        else:
            ws = int(row[wsc]) if (wsc and pd.notna(row.get(wsc))) else 0
            epi = [ws + i for i in parse_index_list(dr[epim])]
        npdb = native_index.get(str(row[idc]).lower()) if idc else None
        res = process(str(row[dirc]), npdb, epi, args.exclude_dist)
        if res is None: n_fail += 1; continue
        plain[pos], natin[pos], aware[pos] = res
        n_ok += 1
        if (pos + 1) % 5000 == 0:
            print(f"  {pos+1}/{len(df)}  ok={n_ok} fail={n_fail}", flush=True)

    df["cylinder_ca_clashes"] = plain
    df["native_in_cylinder"] = natin
    df["cylinder_native_aware"] = aware
    out = Path(args.out_csv); out.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out, index=False)
    print(f"\nWrote {out}  (ok={n_ok} fail={n_fail} gate-skipped={n_skip})")

    na = pd.to_numeric(df["cylinder_native_aware"], errors="coerce")
    print(f"cylinder_native_aware non-null: {na.notna().sum():,}/{len(df):,}")
    tc = first_present(df.columns, TRUE_CLASH)
    if tc is not None:
        truth = pd.to_numeric(df[tc], errors="coerce")
        for c in ("cylinder_ca_clashes", "cylinder_native_aware"):
            p = pd.to_numeric(df[c], errors="coerce"); m = p.notna() & truth.notna()
            if m.sum() >= 10:
                print(f"  {c:24s} vs {tc}: Pearson {p[m].corr(truth[m]):+.3f}")


if __name__ == "__main__":
    main()
