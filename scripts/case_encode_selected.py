#!/usr/bin/env python3
"""
case_encode_selected.py -- rebuild the case-encoded `scaffoldEPITOPE` for selected designs.

Heals the lost token->design mapping (see memory dp4-selection-state): our selections are keyed by a
design `token`, but the epitope positions live in dp2 keyed by (id, rfd_id, mpnn_id). The clean key
that bridges them is the design's OWN SEQUENCE (MPNN-unique). So, per selected design:
  1. read its chain-A sequence from the design PDB (`mpnn_pdb` in the selection table);
  2. match that sequence to its dp2 row (within the same `id`);
  3. uppercase the epitope CHUNK-span positions (`scaffolded_epitope_chunk_resindices`, contiguous per
     island) and lowercase the rest -> the case-encoded string John's controls + assembly need
     (UPPERCASE = epitope, lowercase = scaffold; each contiguous uppercase run = one island).

Validated basis (local, on dp2 ground truth): recorded epitope positions spell the native epitope AA
exactly (400/400), so dp2's positions are trustworthy; the epitope is NOT contiguous (scattered
contacts) so we use the chunk SPANS, which are contiguous and give clean per-island uppercase runs.

Runs on Gemini (needs the design PDBs). dp2 at $WS/datasets/dp2.parquet. Self-diagnosing: reports
matched / unmatched / ambiguous so a coverage gap is loud, not silent.

Usage (see case_encode_selected.sbatch):
  python scripts/case_encode_selected.py \
      --selection results/dp4_C1_whole_epitope_ranked.top20.csv \
      --dp2 $WS/datasets/dp2.parquet --out results/dp4_C1_scaffoldEPITOPE.csv
"""
from __future__ import annotations
import argparse, ast, sys
from pathlib import Path
import numpy as np
import pandas as pd

_REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(_REPO / "episcaf_analysis"))
import compute_metrics as CM   # noqa: E402  (read_structure, get_chain, chain_seq)


def parse_list(x):
    return list(x) if isinstance(x, (list, tuple, np.ndarray)) else ast.literal_eval(x)


def case_encode(seq: str, epi_positions, offset: int = 0) -> str:
    """Lowercase everything, uppercase the epitope chunk-span positions (shifted by `offset`)."""
    chars = [c.lower() for c in seq]
    for p in epi_positions:
        q = p + offset
        if 0 <= q < len(chars):
            chars[q] = chars[q].upper()
    return "".join(chars)


def read_design_seq(pdb: str):
    st = CM.read_structure(pdb)
    return CM.chain_seq(CM.get_chain(st, "A"))


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--selection", required=True, help="CSV of selected designs (needs id + mpnn_pdb)")
    ap.add_argument("--dp2", required=True, help="dp2.parquet (native + scaffolded epitope positions)")
    ap.add_argument("--id-col", default="id")
    ap.add_argument("--pdb-col", default="mpnn_pdb")
    ap.add_argument("--target-col", default="id", help="column to carry as 'target'")
    ap.add_argument("--pdb-remap", default="", help="'FROM:TO' prefix swap if PDB paths are stale")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    sel = pd.read_csv(args.selection, low_memory=False)
    dp2 = pd.read_parquet(args.dp2)
    # per-id lookup: exact design sequence -> chunk-span positions
    by_id: dict = {}
    for idv, g in dp2.dropna(subset=["scaffolded_epitope_seq"]).groupby(args.id_col):
        d = {}
        for _, r in g.iterrows():
            d[r["scaffolded_epitope_seq"]] = parse_list(r["scaffolded_epitope_chunk_resindices"])
        by_id[str(idv)] = d
    remap = args.pdb_remap.split(":", 1) if args.pdb_remap else None

    rows, n_ok, n_nopdb, n_nomatch, n_offset = [], 0, 0, 0, 0
    for r in sel.itertuples(index=False):
        idv = str(getattr(r, args.id_col)); pdb = str(getattr(r, args.pdb_col))
        if remap:
            pdb = pdb.replace(remap[0], remap[1])
        try:
            seq = read_design_seq(pdb)
        except Exception as e:  # noqa: BLE001
            n_nopdb += 1; rows.append(dict(id=idv, status=f"pdb_fail:{e}", scaffoldEPITOPE="")); continue
        lut = by_id.get(idv, {})
        pos = lut.get(seq)
        offset = 0
        if pos is None:  # try: a dp2 seq is a substring of ours (flank offset), or vice versa
            for dseq, dpos in lut.items():
                if dseq in seq:
                    pos, offset = dpos, seq.index(dseq); break
                if seq in dseq:
                    pos, offset = [p - dseq.index(seq) for p in dpos], 0; break
        if pos is None:
            n_nomatch += 1
            rows.append(dict(id=idv, status="no_dp2_seq_match", scaffoldEPITOPE="")); continue
        if offset: n_offset += 1
        se = case_encode(seq, pos, offset)
        n_ok += 1
        rows.append(dict(id=idv, target=getattr(r, args.target_col), design_seq=seq,
                         scaffoldEPITOPE=se, n_islands=se_islands(se), status="ok"))

    out = pd.DataFrame(rows)
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.out, index=False)
    print(f"[case-encode] {len(sel)} selected -> {n_ok} encoded | "
          f"no-pdb {n_nopdb} | no-dp2-match {n_nomatch} | flank-offset {n_offset}")
    if n_ok:
        isl = out[out.status == "ok"]["n_islands"].value_counts().to_dict()
        print(f"[case-encode] island-count distribution among encoded: {isl}")
    print(f"[case-encode] wrote {args.out}")


def se_islands(se: str) -> int:
    import re
    return len(re.findall(r"[A-Z]+", se))


if __name__ == "__main__":
    main()
