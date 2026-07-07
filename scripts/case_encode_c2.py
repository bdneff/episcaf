#!/usr/bin/env python3
"""
case_encode_c2.py -- case-encoded scaffoldEPITOPE for C2 (single-island scaffolds). GEMINI.

C2 designs are our dual-island run (not in dp2), so we localize the epitope from the LEDGER fields
carried in the selection table: a single island occupies design positions [n_flank, n_flank+island_size)
(manuscript sec:methods; same formula stage05 uses). The design SEQUENCE is read from the AF3 output
(`af3_dir` chain A). Uppercase the island span, lowercase the rest -> the case-encoded sequence
(epitope UPPER, scaffold lower; one island -> one uppercase run). C2 is natively 103 (no trim).

Runs on Gemini (needs the af3 outputs). Self-diagnosing.

Usage (see case_encode_c2.sbatch):
  python scripts/case_encode_c2.py \
      --selection results/dp4_C2_single_island_ranked.top20.csv \
      --out results/dp4_C2_scaffoldEPITOPE.csv
"""
from __future__ import annotations
import argparse, re, sys
from pathlib import Path
import pandas as pd

_REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(_REPO / "episcaf_analysis"))
import compute_metrics as CM   # noqa: E402  (find_af3_files, read_structure, get_chain, chain_seq)


def case_encode(seq: str, start: int, length: int) -> str:
    chars = [c.lower() for c in seq]
    for i in range(start, start + length):
        if 0 <= i < len(chars):
            chars[i] = chars[i].upper()
    return "".join(chars)


def read_af3_seq(af3_dir: str) -> str:
    cif, _, _ = CM.find_af3_files(Path(af3_dir))
    if cif is None:
        raise FileNotFoundError(f"no af3 cif under {af3_dir}")
    return CM.chain_seq(CM.get_chain(CM.read_structure(cif), "A"))


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--selection", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--af3-remap", default="", help="'FROM:TO' prefix swap if af3_dir paths are stale")
    args = ap.parse_args()

    d = pd.read_csv(args.selection, low_memory=False)
    remap = args.af3_remap.split(":", 1) if args.af3_remap else None
    rows, n_ok, n_noaf3, n_oob = [], 0, 0, 0
    for r in d.itertuples(index=False):
        af3 = str(getattr(r, "af3_dir"))
        if remap:
            af3 = af3.replace(remap[0], remap[1])
        try:
            seq = read_af3_seq(af3)
        except Exception as e:  # noqa: BLE001
            n_noaf3 += 1
            rows.append(dict(predID=getattr(r, "predID", ""), status=f"af3_fail:{e}", scaffoldEPITOPE="")); continue
        ws = int(getattr(r, "af3_window_start", 0) or 0)
        start = ws + int(getattr(r, "n_flank")); length = int(getattr(r, "island_size"))
        if start + length > len(seq):
            n_oob += 1
            rows.append(dict(predID=getattr(r, "predID", ""),
                             status=f"island_oob(start={start},len={length},seq={len(seq)})",
                             scaffoldEPITOPE="")); continue
        se = case_encode(seq, start, length)
        n_ok += 1
        rows.append(dict(predID=getattr(r, "predID", ""), id=getattr(r, "id", ""),
                         island_index=getattr(r, "island_index", ""), target=getattr(r, "id", ""),
                         design_seq=seq, scaffoldEPITOPE=se,
                         n_islands=len(re.findall(r"[A-Z]+", se)), status="ok"))

    out = pd.DataFrame(rows)
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.out, index=False)
    print(f"[c2] {len(d)} designs -> {n_ok} encoded | af3-fail {n_noaf3} | island-oob {n_oob}")
    if n_ok:
        ok = out[out.status == "ok"]
        print(f"[c2] island-count: {ok.n_islands.value_counts().to_dict()} (expect all 1); "
              f"design_seq len: {ok.design_seq.str.len().value_counts().to_dict()}")
    print(f"[c2] wrote {args.out}")


if __name__ == "__main__":
    main()
