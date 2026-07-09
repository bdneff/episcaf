#!/usr/bin/env python3
"""
stage05_extract_metrics.py -- the SINGLE source-of-truth metrics extractor for the
dual-island per-island design run (RFD3 -> MPNN -> AF3). One pass over the 04_af3 outputs
writes one parquet+csv; every downstream decision (composite scoring, top-5/island) reads
from it and analysis is never rerun. This is our analog of Lawson's dp2.parquet.

WHY a new extractor (neither existing one fits this run):
  - episcaf_analysis/compute_metrics.py is the Lawson-validated antibody extractor, but its
    `run` mode is keyed on dp2's 32-hex tokens and the no-MPNN (RFD3->AF3) layout. We reuse
    its *functions* (RMSD/PAE/clash primitives) but drive them from OUR per-island ledger.
  - episcaf_analysis/build_12mer_metrics.py emits exactly the scorer's column contract
    (epitope_pae, cylinder_*) but parses the 12-mer naming and has no antibody. We mirror its
    PAE-decomposition and cylinder calls.

This run HAS antibodies (the 46 dual-island mAb epitopes), so we compute BOTH:
  - the real DP3 ground-truth clash (af3_n_clash_res, vs native antibody chains B/C), AND
  - the cylinder surrogate (cylinder_ca_clashes / cylinder_native_aware),
giving a bonus validation of the surrogate against truth on a set where both exist.

Definitional choices (match the DP3-validated convention so this set stays comparable):
  - epitope_chunk_rmsd / overall_rmsd : BACKBONE (N,CA,C,O), Kabsch-superposed, exactly as
    compute_metrics.py (NOT build_12mer's CA-only).  [--*_rmsd thresholds: <=1, <=2]
  - cylinder exclude_dist : 1.0 (the "adopted ed=1" of the validated 12-mer preset).
  - epitope set : the WHOLE island span [n_flank, n_flank+island_size), i.e. what stage03
    actually held FIXED during MPNN (load_fixed_lookup), not the contacts-only subset that
    build_dual_island wrote into the ledger's epitope_resindices.

Reference structure for overall_rmsd = the MPNN-designed PDB (the structure AF3 got its
sequence from), matching Lawson's `overall_rmsd` (MPNN chain A backbone vs AF3 window).

Usage (on the cluster):
  python3 scripts/stage05_extract_metrics.py \
      --af3_out_dir  runs/dual_island_rfd3/04_af3/outputs \
      --mpnn_pdb_dir runs/dual_island_rfd3/03_mpnn/mpnn_pdbs \
      --ledger       results/dual_island_designs.csv \
      --native_dir   /tgen_labs/altin/alphafold3/workspace/episcaf-experiments/data/abdb/complex_pdbfiles/cleaned \
      --out          runs/dual_island_rfd3/05_analysis/metrics_dual_island.parquet \
      --workers 16
  # smoke test first:  add  --limit_ids 2   (process two epitopes, ~a few thousand designs)
"""
from __future__ import annotations

import argparse
import json
import re
import sys
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

# Reuse the validated primitives from the antibody extractor + the cylinder core.
_REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(_REPO / "episcaf_analysis"))
import compute_metrics as CM          # noqa: E402  (read_structure, get_chain, chain_seq,
                                      #              bb_coords_*, rmsd_superpose, find_subseq,
                                      #              kabsch_fit, ca_coord, heavy_coords_residue,
                                      #              find_af3_files, summary_scalars)
import native_cylinder_core as CYL    # noqa: E402  (cylinder_frame, native_aware_scaffold_count,
                                      #              count_native_in_cylinder)
from scipy.spatial import cKDTree     # noqa: E402

# predID = <id>__contig<cid>__seed<s>__rep<r>...__fixed_dldesign_<dl>  (lowercased on disk).
PRED_RE = re.compile(
    r"^(?P<id>.+?)__contig(?P<cid>\d+)__seed(?P<seed>\d+)__rep(?P<rep>\d+)"
    r".*_fixed_dldesign_(?P<dl>\d+)$",
    re.IGNORECASE,
)
ISLAND_SEG = re.compile(r"A(\d+)-(\d+)")


def parse_pred(stem: str) -> Optional[Dict[str, Any]]:
    m = PRED_RE.match(stem)
    if not m:
        return None
    return {
        "id": m.group("id"),
        "contig_id": int(m.group("cid")),
        "seed": int(m.group("seed")),
        "rep": int(m.group("rep")),
        "dldesign": int(m.group("dl")),
    }


_CONTIG_TOK = re.compile(r"^(A?)(\d+)-(\d+)$")


def epi_positions_from_contig(contig_string: str) -> Tuple[List[int], List[int]]:
    """Parallel (design_epi, native_epi) 0-based epitope positions from a scaffolding contig.

    Walks the slash-delimited contig (`N-N/Aa-b/gap/Ac-d/C-C`): each A-segment is an epitope island
    whose residues occupy consecutive positions in the scaffolded output (design_epi, cumulative) and
    map to native chain-A positional indices a-1..b-1 (native_epi). Handles ONE island (dual-island
    run) or MANY (whole-epitope C1-103). For a single island this reproduces the old
    range(n_flank, n_flank+island_size) / range(a-1, a-1+island_size) exactly.
    """
    pos = 0
    design_epi: List[int] = []
    native_epi: List[int] = []
    for tok in str(contig_string).split("/"):
        m = _CONTIG_TOK.match(tok)
        if not m:
            raise ValueError(f"unparseable contig token {tok!r} in {contig_string!r}")
        a, b = int(m.group(2)), int(m.group(3))
        if m.group(1) == "A":                         # epitope island
            span = b - a + 1
            design_epi.extend(range(pos, pos + span))
            native_epi.extend(range(a - 1, b))        # native chain-A 0-based positional
            pos += span
        else:                                          # scaffold flank/spacer ('N-N' -> a residues)
            pos += a
    return design_epi, native_epi


def _opt_int(r, key):
    v = r.get(key)
    return int(v) if v is not None and pd.notna(v) else None


def load_ledger(path: Path) -> Dict[Tuple[str, int], Dict[str, Any]]:
    """(id, contig_id) -> {design_epi, native_epi, + optional island_index/size/variant/flanks}.

    Epitope positions come from the contig's A-segments (`epi_positions_from_contig`), so this
    handles BOTH the single-island dual-island ledger and the multi-island whole-epitope ledger
    (which has no n_flank/island_size/island_segment columns). The island_* fields are carried
    through for the dual-island output table and are None for the whole-epitope run.
    """
    df = pd.read_csv(path)
    out: Dict[Tuple[str, int], Dict[str, Any]] = {}
    for _, r in df.iterrows():
        design_epi, native_epi = epi_positions_from_contig(str(r["contig_string"]))
        out[(str(r["id"]), int(r["contig_id"]))] = {
            "design_epi": design_epi, "native_epi": native_epi,
            "island_index": _opt_int(r, "island_index"),
            "island_size": _opt_int(r, "island_size"),
            "variant": _opt_int(r, "variant"),
            "n_flank": _opt_int(r, "n_flank"),
            "c_flank": _opt_int(r, "c_flank"),
        }
    return out


# --------------------------------------------------------------------------- #
# Native complex (loaded ONCE per epitope id, reused across all its designs)
# --------------------------------------------------------------------------- #
def load_native(pdb: Path) -> Dict[str, Any]:
    """Native AbDb complex via gemmi: chain-A CA + heavy (positional), antibody (B/C) heavy."""
    st = CM.read_structure(pdb)
    model = st[0]
    chA = CM.get_chain(st, "A")
    ca_A, heavy_A = [], []
    for res in chA:
        ca = res.find_atom("CA", altloc="*")
        ca_A.append([ca.pos.x, ca.pos.y, ca.pos.z] if ca else [np.nan, np.nan, np.nan])
        h = CM.heavy_coords_residue(res)
        if h.shape[0]:
            heavy_A.append(h)
    ab_heavy = []
    for ch in model:
        if ch.name in ("B", "C"):
            for res in ch:
                h = CM.heavy_coords_residue(res)
                if h.shape[0]:
                    ab_heavy.append(h)
    return {
        "ca_A": np.asarray(ca_A, float),
        "heavy_A": np.vstack(heavy_A) if heavy_A else np.zeros((0, 3)),
        "ab_heavy": np.vstack(ab_heavy) if ab_heavy else np.zeros((0, 3)),
        "n_res_A": len(ca_A),
    }


# --------------------------------------------------------------------------- #
# Per-design metric computation
# --------------------------------------------------------------------------- #
def compute_one(af3_dir: Path, mpnn_pdb: Path, info: Dict[str, Any],
                native: Dict[str, Any], clash_cutoff: float, exclude_dist: float) -> Dict[str, Any]:
    rec: Dict[str, Any] = {
        "overall_rmsd": None, "epitope_chunk_rmsd": None,
        "mean_pae": None, "epitope_pae": None, "scaffold_pae": None,
        "ptm": None, "ranking_score": None, "fraction_disordered": None,
        "has_clash_af3": None, "chain_pair_pae_min": None,
        "af3_n_clash_res": None, "af3_has_clash": None, "af3_clash_status": "not_run",
        "cylinder_ca_clashes": None, "cylinder_native_aware": None,
        "native_in_cylinder": None,
        "af3_window_start": None, "status": "ok",
    }

    cif, conf, summary = CM.find_af3_files(af3_dir)
    if cif is None:
        rec["status"] = "no_af3_cif"
        return rec
    rec.update(CM.summary_scalars(summary))   # ptm, ranking_score, fraction_disordered, ...

    # --- load design (MPNN PDB AF3 saw) and AF3 prediction ---
    try:
        st_m = CM.read_structure(mpnn_pdb)
        st_a = CM.read_structure(cif)
        ch_m = CM.get_chain(st_m, "A")
        ch_a = CM.get_chain(st_a, "A")
    except Exception as e:  # noqa: BLE001
        rec["status"] = f"structure_load_fail: {e}"
        return rec

    mpnn_seq, af3_seq = CM.chain_seq(ch_m), CM.chain_seq(ch_a)
    ws = CM.find_subseq(af3_seq, mpnn_seq)     # AF3 received exactly the design seq -> ws==0
    if ws < 0:
        rec["status"] = "seqmatch_fail"
        return rec
    we = ws + len(mpnn_seq)
    rec["af3_window_start"] = ws

    design_epi = list(info["design_epi"])                  # design/MPNN frame (one or many islands)
    af3_epi = [ws + i for i in design_epi]                  # AF3 frame
    native_epi = list(info["native_epi"])                  # native chain-A positional

    af3_res = list(ch_a)
    if af3_epi and max(af3_epi) >= len(af3_res):
        rec["status"] = f"af3_epi_oob(n={len(af3_res)})"
        return rec

    # --- overall + epitope backbone RMSD (mirror compute_metrics.compute_pair_metrics) ---
    try:
        P = CM.bb_coords_range(ch_m, 0, len(mpnn_seq))
        Q = CM.bb_coords_range(ch_a, ws, we)
        rec["overall_rmsd"] = CM.rmsd_superpose(P, Q)
        P_epi = CM.bb_coords_resindices(ch_m, design_epi)
        Q_epi = CM.bb_coords_resindices(ch_a, af3_epi)
        rec["epitope_chunk_rmsd"] = CM.rmsd_superpose(P_epi, Q_epi)
    except Exception as e:  # noqa: BLE001
        rec["status"] = f"rmsd_fail: {e}"
        return rec

    # --- PAE decomposition (mirror build_12mer_metrics) ---
    if conf is not None:
        try:
            d = json.loads(conf.read_text())
            pae = None
            for k in ("pae", "predicted_aligned_error"):
                if k in d:
                    arr = np.asarray(d[k], float)
                    if arr.ndim == 2:
                        pae = arr
                        break
            if pae is not None:
                n = pae.shape[0]
                epi = [i for i in af3_epi if i < n]
                scaf = [i for i in range(n) if i not in set(epi)]
                rec["mean_pae"] = float(np.nanmean(pae))
                rec["epitope_pae"] = float(np.nanmean(pae[np.ix_(epi, epi)])) if epi else None
                rec["scaffold_pae"] = float(np.nanmean(pae[np.ix_(scaf, scaf)])) if scaf else None
        except Exception:  # noqa: BLE001
            pass

    # --- AF3 CA arrays (for clash transform + cylinder) ---
    af3_ca_all = np.array([CM.ca_coord(r) if CM.ca_coord(r) is not None else [np.nan] * 3
                           for r in af3_res], float)
    # epitope CA pairs (AF3 <-> native), keep only where both exist
    P_pairs, Q_pairs, af3_epi_ca = [], [], []
    for i_af3, i_nat in zip(af3_epi, native_epi):
        p = CM.ca_coord(af3_res[i_af3])
        q = native["ca_A"][i_nat] if i_nat < native["n_res_A"] else None
        af3_epi_ca.append(p if p is not None else [np.nan] * 3)
        if p is not None and q is not None and np.all(np.isfinite(q)):
            P_pairs.append(p)
            Q_pairs.append(q)
    af3_epi_ca = np.asarray(af3_epi_ca, float)

    epi_set = set(af3_epi)
    scaf_mask = np.array([i not in epi_set for i in range(len(af3_res))])
    scaf_ca = af3_ca_all[scaf_mask]

    have_align = len(P_pairs) >= 3 and native["n_res_A"] > 0
    if have_align:
        P_pairs, Q_pairs = np.vstack(P_pairs), np.vstack(Q_pairs)

        # --- real antibody clash (DP3 ground truth): AF3 scaffold -> native frame ---
        if native["ab_heavy"].shape[0] > 0:
            R, t = CM.kabsch_fit(P_pairs, Q_pairs)        # (af3 @ R) + t ~= native
            tree = cKDTree(native["ab_heavy"])
            n_clash = 0
            for i in np.where(scaf_mask)[0]:
                heavy = CM.heavy_coords_residue(af3_res[i])
                if heavy.shape[0] == 0:
                    continue
                heavy_fit = (heavy @ R) + t
                dmin, _ = tree.query(heavy_fit, k=1)
                if np.any(dmin < clash_cutoff):
                    n_clash += 1
            rec["af3_n_clash_res"] = int(n_clash)
            rec["af3_has_clash"] = n_clash > 0
            rec["af3_clash_status"] = "ok"
        else:
            rec["af3_clash_status"] = "no_antibody_atoms"

        # --- cylinder surrogate: native heavy -> AF3 frame, build cylinder in AF3 frame ---
        try:
            R2, t2 = CM.kabsch_fit(Q_pairs, P_pairs)      # (native @ R2) + t2 ~= af3
            native_heavy_af3 = (native["heavy_A"] @ R2) + t2
            base, normal = CYL.cylinder_frame(af3_epi_ca[np.isfinite(af3_epi_ca).all(1)], af3_ca_all[np.isfinite(af3_ca_all).all(1)])
            aware, plain = CYL.native_aware_scaffold_count(
                scaf_ca[np.isfinite(scaf_ca).all(1)], base, normal,
                native_heavy_af3, exclude_dist=exclude_dist)
            rec["cylinder_ca_clashes"] = int(plain)
            rec["cylinder_native_aware"] = int(aware)
            nat_ca_af3 = (native["ca_A"] @ R2) + t2
            nat_is_epi = np.zeros(native["n_res_A"], bool)
            for i in native_epi:
                if i < native["n_res_A"]:
                    nat_is_epi[i] = True
            finite = np.isfinite(nat_ca_af3).all(1)
            rec["native_in_cylinder"] = CYL.count_native_in_cylinder(
                nat_ca_af3[finite], nat_is_epi[finite], base, normal)
        except Exception as e:  # noqa: BLE001
            rec["status"] = f"cylinder_fail: {e}"
    else:
        rec["af3_clash_status"] = "too_few_epitope_pairs"

    return rec


def process_id(idkey: str, native_pdb: str, designs: List[Dict[str, Any]],
               clash_cutoff: float, exclude_dist: float) -> List[Dict[str, Any]]:
    """Process every design of one epitope id, loading the native complex once."""
    try:
        native = load_native(Path(native_pdb))
        native_status = "ok"
    except Exception as e:  # noqa: BLE001
        native, native_status = None, f"native_load_fail: {e}"

    rows = []
    for d in designs:
        base = {
            "predID": d["predID"], "id": d["id"], "contig_id": d["contig_id"],
            "island_index": d["island_index"], "island_size": d["island_size"],
            "variant": d["variant"], "n_flank": d["n_flank"], "c_flank": d["c_flank"],
            "seed": d["seed"], "rep": d["rep"], "dldesign": d["dldesign"],
            "antigen": d["id"],          # per-epitope scope for the scorer
            "af3_dir": d["af3_dir"], "native_status": native_status,
        }
        if native is None:
            base["status"] = native_status
            rows.append(base)
            continue
        info = {"design_epi": d["design_epi"], "native_epi": d["native_epi"]}
        base.update(compute_one(Path(d["af3_dir"]), Path(d["mpnn_pdb"]), info,
                                native, clash_cutoff, exclude_dist))
        rows.append(base)
    return rows


# --------------------------------------------------------------------------- #
def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--af3_out_dir", required=True, help="<run>/04_af3/outputs")
    ap.add_argument("--mpnn_pdb_dir", required=True, help="<run>/03_mpnn/mpnn_pdbs")
    ap.add_argument("--ledger", default="results/dual_island_designs.csv")
    ap.add_argument("--native_dir", required=True, help="AbDb cleaned <id>.pdb dir")
    ap.add_argument("--out", required=True, help="output .parquet (a .csv is written alongside)")
    ap.add_argument("--workers", type=int, default=16)
    ap.add_argument("--clash_cutoff", type=float, default=4.0)
    ap.add_argument("--exclude_dist", type=float, default=1.0)
    ap.add_argument("--limit_ids", type=int, default=0, help="process only N epitope ids (smoke test)")
    args = ap.parse_args()

    ledger = load_ledger(Path(args.ledger))
    print(f"[stage05] ledger contigs: {len(ledger):,}")

    # Index MPNN PDBs by lowercased stem (AF3 dir names are lowercased).
    mpnn_idx = {p.stem.lower(): p for p in Path(args.mpnn_pdb_dir).rglob("*_fixed_dldesign_*.pdb")}
    print(f"[stage05] MPNN PDBs     : {len(mpnn_idx):,}")

    # Walk AF3 output dirs that actually produced a model CIF; group designs by epitope id.
    native_dir = Path(args.native_dir)
    by_id: Dict[str, List[Dict[str, Any]]] = defaultdict(list)
    n_seen = n_noparse = n_nompnn = n_noledger = 0
    for d in Path(args.af3_out_dir).iterdir():
        if not d.is_dir() or next(d.glob("*_model.cif"), None) is None:
            continue
        n_seen += 1
        mp = mpnn_idx.get(d.name)
        if mp is None:
            n_nompnn += 1
            continue
        p = parse_pred(mp.stem)                 # original-case id from MPNN stem
        if p is None:
            n_noparse += 1
            continue
        info = ledger.get((p["id"], p["contig_id"]))
        if info is None:
            n_noledger += 1
            continue
        by_id[p["id"]].append({
            "predID": mp.stem, "af3_dir": str(d), "mpnn_pdb": str(mp),
            **p, **info,
        })

    ids = sorted(by_id)
    if args.limit_ids > 0:
        ids = ids[:args.limit_ids]
    n_designs = sum(len(by_id[i]) for i in ids)
    print(f"[stage05] AF3 dirs w/CIF: {n_seen:,}  (no_mpnn={n_nompnn}, no_parse={n_noparse}, "
          f"no_ledger={n_noledger})")
    print(f"[stage05] epitopes: {len(ids)}  designs to score: {n_designs:,}")

    rows: List[Dict[str, Any]] = []
    with ProcessPoolExecutor(max_workers=args.workers) as ex:
        futs = {ex.submit(process_id, i, str(native_dir / f"{i}.pdb"), by_id[i],
                          args.clash_cutoff, args.exclude_dist): i for i in ids}
        for k, fut in enumerate(as_completed(futs), 1):
            rows.extend(fut.result())
            print(f"[stage05] {k}/{len(ids)} epitopes done ({futs[fut]})", flush=True)

    df = pd.DataFrame(rows)
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    df.to_parquet(out, index=False)
    df.to_csv(out.with_suffix(".csv"), index=False)
    print(f"\n[stage05] wrote {len(df):,} rows -> {out}  (+ .csv)")

    ok = df[df["status"] == "ok"] if "status" in df else df
    print(f"[stage05] status==ok: {len(ok):,}/{len(df):,}")
    if "status" in df and (df["status"] != "ok").any():
        print(df.loc[df["status"] != "ok", "status"].value_counts().head(10).to_string())

    # DP3 four-filter summary (uses the REAL antibody clash) + cylinder agreement.
    if len(ok):
        for c in ("overall_rmsd", "epitope_chunk_rmsd", "mean_pae", "epitope_pae",
                  "ptm", "af3_n_clash_res", "cylinder_ca_clashes", "cylinder_native_aware"):
            if c in ok:
                v = pd.to_numeric(ok[c], errors="coerce")
                print(f"  {c:22s} mean={v.mean():.3f}  median={v.median():.3f}  nonnull={v.notna().sum():,}")
        passing = ok[
            (pd.to_numeric(ok["overall_rmsd"], errors="coerce") <= 2.0)
            & (pd.to_numeric(ok["epitope_chunk_rmsd"], errors="coerce") <= 1.0)
            & (pd.to_numeric(ok["mean_pae"], errors="coerce") < 5.0)
            & (pd.to_numeric(ok["af3_n_clash_res"], errors="coerce") == 0)
        ]
        print(f"\n  DP3 4-filter PASS (real clash): {len(passing):,}  "
              f"({100*len(passing)/max(len(ok),1):.1f}% of ok)")


if __name__ == "__main__":
    main()
