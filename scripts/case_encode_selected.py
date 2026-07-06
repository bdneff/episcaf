#!/usr/bin/env python3
"""
case_encode_selected.py -- rebuild the case-encoded `scaffoldEPITOPE` for selected designs.

Full mechanism + why-it-works in docs/CASE_ENCODING.md. In short: we reran RFdiffusion3 on Lawson's
contigs, so each design's epitope span is a contig property. The metrics driver
(`compute_metrics.py::run_metrics`, ~line 700-712) joins each design to `dp2` by
`token == dp2.assay_scaffolded_epitope_id` and takes the epitope positions from
`dp2.scaffolded_epitope_chunk_resindices` (0-based indices into the design chain, contiguous per
island). We do exactly that here, then uppercase those positions in the design's own chain-A sequence
(read from its PDB) -> the case-encoded string (UPPERCASE = epitope, lowercase = scaffold; each
contiguous uppercase run = one island) that C6 (`build_c6_mutants.py`) and assembly consume.

IMPORTANT: join by TOKEN, not by sequence -- `dp2.scaffolded_epitope_seq` is LAWSON's sequence, not
ours (different MPNN run); only the contig-determined POSITIONS transfer. And use the dp2 that carries
OUR tokens (repo_refactored/datasets or $WS/datasets), not the local known_antigen dp2 (only ~149 match).

Runs on Gemini (needs the design PDBs). Self-diagnosing: prints token coverage + failure buckets.
"""
from __future__ import annotations
import argparse, ast, re, sys
from pathlib import Path
import numpy as np
import pandas as pd

_REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(_REPO / "episcaf_analysis"))
import compute_metrics as CM   # noqa: E402  (read_structure, get_chain, chain_seq)


def parse_list(x):
    if isinstance(x, (list, tuple, np.ndarray)):
        return [int(v) for v in x]
    if x is None or (isinstance(x, float) and np.isnan(x)):
        return []
    return [int(v) for v in ast.literal_eval(x)] if isinstance(x, str) else list(x)


def case_encode(seq: str, epi_positions) -> str:
    chars = [c.lower() for c in seq]
    for p in epi_positions:
        if 0 <= p < len(chars):
            chars[p] = chars[p].upper()
    return "".join(chars)


def n_islands(se: str) -> int:
    return len(re.findall(r"[A-Z]+", se))


def read_design_seq(pdb: str) -> str:
    st = CM.read_structure(Path(pdb))          # CM.read_structure needs a Path (uses .suffix)
    return CM.chain_seq(CM.get_chain(st, "A"))


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--selection", required=True, help="CSV of selected designs (needs token + mpnn_pdb)")
    ap.add_argument("--dp2", required=True, help="dp2.parquet that carries OUR tokens (assay_scaffolded_epitope_id)")
    ap.add_argument("--token-col", default="token")
    ap.add_argument("--pdb-col", default="mpnn_pdb")
    ap.add_argument("--target-col", default="id")
    ap.add_argument("--pdb-remap", default="", help="'FROM:TO' prefix swap if PDB paths are stale")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    sel = pd.read_csv(args.selection, low_memory=False)
    dp2 = pd.read_parquet(args.dp2)
    dp2["assay_scaffolded_epitope_id"] = dp2["assay_scaffolded_epitope_id"].astype(str).str.lower()
    # token -> epitope chunk-span positions (contig-determined; same for our rerun as Lawson's)
    by_token = {}
    for tok, chunks in zip(dp2["assay_scaffolded_epitope_id"], dp2["scaffolded_epitope_chunk_resindices"]):
        if tok not in by_token:
            by_token[tok] = parse_list(chunks)

    sel_tokens = sel[args.token_col].astype(str).str.lower()
    cover = sel_tokens.isin(by_token).sum()
    print(f"[case-encode] dp2 tokens cover {cover}/{len(sel_tokens)} selected tokens "
          f"({100*cover/max(len(sel_tokens),1):.0f}%)")
    if cover == 0:
        print("[case-encode] WARNING: 0 token overlap -- wrong dp2 (Lawson-token version?). "
              "Point --dp2 at the run's dp2 (repo_refactored/datasets or $WS/datasets).")

    remap = args.pdb_remap.split(":", 1) if args.pdb_remap else None
    rows, n_ok, n_notok, n_nopdb, n_oob = [], 0, 0, 0, 0
    for r in sel.itertuples(index=False):
        tok = str(getattr(r, args.token_col)).lower()
        pos = by_token.get(tok)
        if pos is None:
            n_notok += 1; rows.append(dict(token=tok, status="no_dp2_token", scaffoldEPITOPE="")); continue
        pdb = str(getattr(r, args.pdb_col))
        if remap:
            pdb = pdb.replace(remap[0], remap[1])
        try:
            seq = read_design_seq(pdb)
        except Exception as e:  # noqa: BLE001
            n_nopdb += 1; rows.append(dict(token=tok, status=f"pdb_fail:{e}", scaffoldEPITOPE="")); continue
        if pos and max(pos) >= len(seq):
            n_oob += 1
            rows.append(dict(token=tok, status=f"pos_oob(max={max(pos)},len={len(seq)})", scaffoldEPITOPE=""))
            continue
        se = case_encode(seq, pos)
        n_ok += 1
        rows.append(dict(token=tok, target=getattr(r, args.target_col), design_seq=seq,
                         scaffoldEPITOPE=se, n_islands=n_islands(se), status="ok"))

    out = pd.DataFrame(rows)
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.out, index=False)
    print(f"[case-encode] {len(sel)} selected -> {n_ok} encoded | "
          f"no-token {n_notok} | no-pdb {n_nopdb} | pos-oob {n_oob}")
    if n_ok:
        print(f"[case-encode] island-count distribution among encoded: "
              f"{out[out.status=='ok']['n_islands'].value_counts().to_dict()}")
    print(f"[case-encode] wrote {args.out}")


if __name__ == "__main__":
    main()
