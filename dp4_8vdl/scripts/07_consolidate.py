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
Composite = the shared `antibody_softgate` scorer -- the SAME soft-gate + global-pass promotion the C1/C2
arms use (John, 2026-07-20: apply the soft-gate clash weighting to 8VDL too, so its clashers are
penalized the way they are everywhere else). 8VDL has the known H/L Fab, so accessibility is the REAL
clash `af3_n_clash_res`, not the cylinder surrogate. Top-`topk` per run -> 8-column rows for
stage06_assemble; `--metrics-out` also dumps every design's scored metrics for cheap re-ranks.

Runs on the cluster (gemmi + AF3 outputs). Usage (from dp4_8vdl/):
  python scripts/07_consolidate.py --runs epitope,hotspots --topk 10 \
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
from score import score                            # noqa: E402  the shared config-driven scorer
from presets import PRESETS                        # noqa: E402

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
    cif, conf, summ = CM.find_af3_files(Path(af3_dir))
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
    # EPITOPEscaffold casing: epitope residues (the fixed contig positions) UPPER, scaffold lower --
    # the same convention every other component ships (John, 2026-07-14).
    epi_set = set(cpos)
    ds = "".join(c.upper() if i in epi_set else c.lower() for i, c in enumerate(seq))
    rec = {"designedSequence": ds, "epitope_chunk_rmsd": CM.rmsd_superpose(P, epi_ca),
           "af3_n_clash_res": None, "epitope_pae": None, "scaffold_pae": None, "mean_pae": None,
           "ptm": None, "overall_rmsd": None}

    # real H/L clash: predicted scaffold -> native frame via the epitope fit
    if ab_heavy.shape[0]:
        R, t = CM.kabsch_fit(P, epi_ca)
        tree = cKDTree(ab_heavy)
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
                scaf = [i for i in range(pae.shape[0]) if i not in epi_set]   # scaffold x scaffold block
                rec["scaffold_pae"] = float(np.nanmean(pae[np.ix_(scaf, scaf)])) if scaf else None
        except Exception:  # noqa: BLE001
            pass

    rec["ptm"] = CM.summary_scalars(Path(summ) if summ else None).get("ptm")

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


def softgate(df):
    """Score with the SAME antibody_softgate preset the C1/C2 arms use: steep soft gates on fold
    quality, a broad heavily-weighted REAL-clash term (af3_n_clash_res -- 8VDL has the known H/L Fab, so
    no cylinder surrogate), and global-pass promotion. Replaces the old scale-blind percentile scorer;
    disable the preset's own top-k select so we rank the full run and cut in process_run."""
    preset = {k: (v.copy() if isinstance(v, dict) else v) for k, v in PRESETS["antibody_softgate"].items()}
    preset["select"] = None
    return score(df, preset)                        # adds 'composite' (+ pass_indicator)


def process_run(run, topk, base):
    contigs = pd.read_csv(base / "01_contigs" / f"{run}.csv")
    resids = [int(x) for x in str(contigs.fixed_resids.iloc[0]).split(",")]
    cpos_by_cid = {int(r.contig_id): construct_positions(
        [int(x) for x in str(r.scaffold_segs).split(",")], island_sizes(r.islands))
        for r in contigs.itertuples(index=False)}
    epi_ca, ab_heavy = load_native(base / "data" / "8VDL.pdb", resids)

    out_root = base / "runs" / run / "04_af3" / "outputs"
    if not out_root.is_dir():
        print(f"[{run}] SKIP -- no AF3 outputs at {out_root} (run its MPNN->AF3 first)")
        return pd.DataFrame()
    mpnn_idx = {p.stem.lower(): p for p in
                (base / "runs" / run / "03_mpnn" / "mpnn_pdbs").rglob("*_fixed_dldesign_*.pdb")}
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
        return df, df
    df = softgate(df)
    df["rank_in_group"] = df["composite"].rank(ascending=False, method="first").astype(int)  # within run
    df["is_global_pass"] = df.get("pass_indicator", pd.Series(0.0, index=df.index)) > 0.5
    top = df.nlargest(topk, "composite").copy()
    npass = int(df["is_global_pass"].sum())
    print(f"[{run}] top-{topk}: epiRMSD {top.epitope_chunk_rmsd.min():.2f}-{top.epitope_chunk_rmsd.max():.2f}  "
          f"clash {top.af3_n_clash_res.min()}-{top.af3_n_clash_res.max()} (was scale-blind percentile); "
          f"{npass} four-filter passers in the run")
    return top, df


def emit(top, out_path):
    """Build the stage06_assemble-shaped selection file from a `top` frame (per-run top-k rows, with the
    softgate `composite`/`rank_in_group`/`is_global_pass` already attached). Shared by the fresh-scoring
    path and the --from-metrics re-rank path so both emit an identical column layout."""
    # 8 standard columns + the 5 scoring columns the library ships (stage06_assemble METRICS names).
    # `sequence` is the plain synthesized 103-mer; `designedSequence` keeps the EPITOPEscaffold casing.
    # No cylinder_clashes column: 8VDL is a known-antibody target, so accessibility is the REAL H/L
    # clash (af3_clashes) and the cylinder surrogate was never computed -- assembly leaves it blank.
    out = pd.DataFrame({
        "sequence": top.designedSequence.str[:103].str.upper(),
        "category": "scaffolded8VDL",
        "model": "RFD",
        "designedSequence": top.designedSequence,
        "designedSequenceLength": top.designedSequence.str.len(),
        "design_ID": top.design_ID,
        "target": "8VDL_" + top.run,
        "epitope_rmsd": top.epitope_chunk_rmsd,
        "overall_rmsd": top.overall_rmsd,
        "epitope_pae": top.epitope_pae,
        "scaffold_pae": top.scaffold_pae,
        "mean_pae": top.mean_pae,
        "ptm": top.ptm,
        "af3_clashes": top.af3_n_clash_res,
        "composite": top.composite,
        "rank_in_group": top.rank_in_group,
        "is_global_pass": top.is_global_pass,
    })
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(out_path, index=False)
    print(f"\nwrote {len(out)} rows ({out.target.value_counts().to_dict()}) -> {out_path}")


def rerank_from_metrics(metrics_csv, topk, runs):
    """Re-slice a deeper (or shallower) top-k from an existing --metrics-out dump -- NO AF3 re-read (the
    dump already holds every design's softgate `composite` and per-run `rank_in_group`). The cut is
    deterministic: `rank_in_group <= topk` is exactly `nlargest(topk, composite)` (method='first'), so
    this reproduces a fresh --topk run bit-for-bit. This is the documented purpose of the metrics dump."""
    df = pd.read_csv(metrics_csv)
    keep = {r.strip() for r in runs.split(",") if r.strip()}
    df = df[df.run.isin(keep)]
    top = df[df.rank_in_group <= topk].sort_values(["run", "rank_in_group"]).reset_index(drop=True)
    for run, g in top.groupby("run"):
        print(f"[{run}] re-rank top-{topk}: epiRMSD {g.epitope_chunk_rmsd.min():.2f}-"
              f"{g.epitope_chunk_rmsd.max():.2f}  clash {int(g.af3_n_clash_res.min())}-"
              f"{int(g.af3_n_clash_res.max())}  ({int(g.is_global_pass.sum())} four-filter passers)")
    return top


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--runs", default="epitope,hotspots",
                    help="the two epitope definitions: `epitope` (contiguous C652-673) and `hotspots` "
                         "(F655/F656/E666)")
    ap.add_argument("--topk", type=int, default=10)
    ap.add_argument("--base", type=Path, default=_HERE.parent, help="dp4_8vdl/ dir")
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--from-metrics", type=Path, default=None,
                    help="re-rank from an existing --metrics-out dump instead of re-reading AF3 -- deepens "
                         "or shrinks the per-run top-k with no cluster round-trip (deterministic re-slice)")
    ap.add_argument("--metrics-out", type=Path, default=None,
                    help="also dump the FULL per-design scored metrics (all designs, not just top-k), so a "
                         "later re-rank needs no AF3 re-read. Default: <out>_allmetrics.csv beside --out")
    args = ap.parse_args()

    if args.from_metrics is not None:
        emit(rerank_from_metrics(args.from_metrics, args.topk, args.runs), args.out)
        return

    parts = [process_run(r.strip(), args.topk, args.base) for r in args.runs.split(",") if r.strip()]
    top = pd.concat([t for t, _ in parts if not t.empty], ignore_index=True)

    fulls = [f for _, f in parts if not f.empty]
    if fulls:
        mout = args.metrics_out or args.out.with_name(args.out.stem + "_allmetrics.csv")
        mout.parent.mkdir(parents=True, exist_ok=True)
        pd.concat(fulls, ignore_index=True).to_csv(mout, index=False)
        print(f"full per-design metrics ({sum(len(f) for f in fulls)} rows) -> {mout}")
    emit(top, args.out)


if __name__ == "__main__":
    main()
