#!/usr/bin/env python3
"""
07_consolidate.py -- score the 8VDL scaffolding runs and emit the top-10 per run for DP4.

Custom-protocol metrics (8VDL is a chain-C antigen + H/L C7 Fab, not the AbDb chain-A/B-C layout the
episcaf stage05 assumes). Per design:
  - align the predicted epitope CAs onto the native chain-C epitope -> epitope RMSD, and the transform
    that carries the scaffold into the antibody frame;
  - in that frame, count scaffold residues whose heavy atoms clash with the H/L Fab -> af3_n_clash_res
    (the REAL clash, since the antibody is known);
  - epitope PAE from the AF3 confidence; overall RMSD vs the MPNN backbone.
Composite = antibody preset weights (0.35 clash + 0.35 epitope RMSD + 0.15 overall RMSD + 0.15 epitope
PAE), percentile within run, lower-is-better. Top-`topk` per run -> 8-column rows for stage06_assemble.

Runs on the cluster (gemmi + AF3 outputs). Usage (from dp4_8vdl/):
  python scripts/07_consolidate.py --runs epitope20,hotspots,contact --topk 10 \
      --out ../results/dp4_8vdl_top10.csv
"""
from __future__ import annotations
import argparse
import json
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

_HERE = Path(__file__).resolve().parent
_REPO = _HERE.parents[1]
sys.path.insert(0, str(_REPO / "episcaf_analysis"))
sys.path.insert(0, str(_HERE))
import compute_metrics as CM                       # noqa: E402
from contact_epitope import AA3                    # noqa: E402

COMPOSITE = [("af3_n_clash_res", 0.35), ("epitope_chunk_rmsd", 0.35),
             ("overall_rmsd", 0.15), ("epitope_pae", 0.15)]
_CID = re.compile(r"_contig(\d+)")


def island_sizes(field: str):
    return [int(b) - int(a) + 1 for a, b in (s.split("-") for s in str(field).split(";"))]


def construct_positions(scaffold_segs, sizes):
    pos, out = 0, []
    for seg, isz in zip(scaffold_segs, sizes):
        pos += seg
        out.extend(range(pos, pos + isz))
        pos += isz
    return out


def load_native(pdb: Path, resids):
    """native epitope CA coords (in `resids` order) + H/L heavy-atom cloud."""
    st = CM.read_structure(pdb)
    ca = {}
    for res in CM.get_chain(st, "C"):
        c = res.find_atom("CA", altloc="*")
        if c is not None:
            ca[res.seqid.num] = [c.pos.x, c.pos.y, c.pos.z]
    epi_ca = np.array([ca[r] for r in resids], float)
    ab = [CM.heavy_coords_residue(res) for ch in st[0] if ch.name in ("H", "L")
          for res in ch if res.name in AA3]
    ab = [h for h in ab if h.shape[0]]
    return epi_ca, (np.vstack(ab) if ab else np.zeros((0, 3)))


def compute_one(af3_dir, mpnn_pdb, cpos, epi_ca, ab_heavy, cutoff=4.0):
    cif, conf, _ = CM.find_af3_files(Path(af3_dir))
    if cif is None:
        return None
    ch_a = CM.get_chain(CM.read_structure(cif), "A")
    res = list(ch_a)
    seq = CM.chain_seq(ch_a)
    if not cpos or max(cpos) >= len(res):
        return None
    P = np.array([CM.ca_coord(res[i]) for i in cpos], float)
    if P.shape[0] != epi_ca.shape[0] or not np.all(np.isfinite(P)):
        return None
    rec = {"designedSequence": seq, "epitope_chunk_rmsd": CM.rmsd_superpose(P, epi_ca),
           "af3_n_clash_res": None, "epitope_pae": None, "mean_pae": None, "overall_rmsd": None}

    # real H/L clash: predicted scaffold -> native frame via the epitope fit
    if ab_heavy.shape[0]:
        R, t = CM.kabsch_fit(P, epi_ca)
        tree = cKDTree(ab_heavy)
        epi_set = set(cpos)
        n = 0
        for i in range(len(res)):
            if i in epi_set:
                continue
            h = CM.heavy_coords_residue(res[i])
            if h.shape[0] and np.any(tree.query((h @ R) + t, k=1)[0] < cutoff):
                n += 1
        rec["af3_n_clash_res"] = n

    if conf is not None:
        try:
            d = json.loads(Path(conf).read_text())
            pae = next((np.asarray(d[k], float) for k in ("pae", "predicted_aligned_error")
                        if k in d and np.asarray(d[k]).ndim == 2), None)
            if pae is not None:
                rec["mean_pae"] = float(np.nanmean(pae))
                epi = [i for i in cpos if i < pae.shape[0]]
                rec["epitope_pae"] = float(np.nanmean(pae[np.ix_(epi, epi)])) if epi else None
        except Exception:  # noqa: BLE001
            pass

    if mpnn_pdb and Path(mpnn_pdb).exists():
        try:
            ch_m = CM.get_chain(CM.read_structure(Path(mpnn_pdb)), "A")
            mseq = CM.chain_seq(ch_m)
            ws = CM.find_subseq(seq, mseq)
            if ws >= 0:
                rec["overall_rmsd"] = CM.rmsd_superpose(CM.bb_coords_range(ch_m, 0, len(mseq)),
                                                        CM.bb_coords_range(ch_a, ws, ws + len(mseq)))
        except Exception:  # noqa: BLE001
            pass
    return rec


def composite(df):
    sc = pd.Series(0.0, index=df.index)
    for col, w in COMPOSITE:
        sc = sc + (1.0 - df[col].rank(pct=True)).fillna(0.0) * w    # lower metric -> higher score
    return sc


def process_run(run, topk, base):
    contigs = pd.read_csv(base / "01_contigs" / f"{run}.csv")
    resids = [int(x) for x in str(contigs.fixed_resids.iloc[0]).split(",")]
    cpos_by_cid = {int(r.contig_id): construct_positions(
        [int(x) for x in str(r.scaffold_segs).split(",")], island_sizes(r.islands))
        for r in contigs.itertuples(index=False)}
    epi_ca, ab_heavy = load_native(base / "data" / "8VDL.pdb", resids)

    mpnn_idx = {p.stem.lower(): p for p in
                (base / "runs" / run / "03_mpnn" / "mpnn_pdbs").rglob("*_fixed_dldesign_*.pdb")}
    out_root = base / "runs" / run / "04_af3" / "outputs"
    rows, n_seen, n_bad = [], 0, 0
    for d in sorted(p for p in out_root.iterdir() if p.is_dir()):
        if next(d.glob("*_model.cif"), None) is None:
            continue
        n_seen += 1
        mp = mpnn_idx.get(d.name)
        m = _CID.search(d.name)
        if m is None:
            n_bad += 1
            continue
        cpos = cpos_by_cid.get(int(m.group(1)))
        rec = compute_one(d, mp, cpos, epi_ca, ab_heavy)
        if rec is None:
            n_bad += 1
            continue
        rec.update(design_ID=d.name, run=run)
        rows.append(rec)

    df = pd.DataFrame(rows)
    print(f"[{run}] af3 dirs {n_seen}  scored {len(df)}  skipped {n_bad}")
    if df.empty:
        return df
    df["composite"] = composite(df)
    top = df.nlargest(topk, "composite").copy()
    print(f"[{run}] top-{topk}: epiRMSD {top.epitope_chunk_rmsd.min():.2f}-{top.epitope_chunk_rmsd.max():.2f}  "
          f"clash {top.af3_n_clash_res.min()}-{top.af3_n_clash_res.max()}")
    return top


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--runs", default="epitope20,hotspots,contact")
    ap.add_argument("--topk", type=int, default=10)
    ap.add_argument("--base", type=Path, default=_HERE.parent, help="dp4_8vdl/ dir")
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()

    parts = [process_run(r.strip(), args.topk, args.base) for r in args.runs.split(",") if r.strip()]
    top = pd.concat([p for p in parts if not p.empty], ignore_index=True)
    out = pd.DataFrame({
        "sequence": top.designedSequence.str[:103],
        "category": "scaffolded8VDL",
        "model": "RFD",
        "designedSequence": top.designedSequence,
        "designedSequenceLength": top.designedSequence.str.len(),
        "design_ID": top.design_ID,
        "target": "8VDL_" + top.run,
    })
    args.out.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.out, index=False)
    print(f"\nwrote {len(out)} rows ({out.target.value_counts().to_dict()}) -> {args.out}")


if __name__ == "__main__":
    main()
