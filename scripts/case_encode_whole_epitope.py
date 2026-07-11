#!/usr/bin/env python3
"""
case_encode_whole_epitope.py -- case-encode the native-103 C1 designs (whole-epitope run).

Unlike the old C1 (which recovered epitope positions via a token->dp2 join), these are FRESH
designs whose epitope positions are known directly from their contig: the A-segments give the
0-based construct positions (`epi_positions_from_contig`, the same primitive stage05 uses). For
each selected design we read its AlphaFold3 chain-A sequence and uppercase those positions, so
each contiguous uppercase run is one epitope island. Feeds C6 + assembly (the `scaffoldEPITOPE`
column), same as the other components.

Runs on the cluster (needs gemmi + the AF3 outputs). Usage:
  python scripts/case_encode_whole_epitope.py \
      --selected results/dp4_C1_whole_epitope_ranked.top20.csv \
      --ledger   results/whole_epitope_designs.csv \
      --out      results/dp4_C1_scaffoldEPITOPE.csv
"""
from __future__ import annotations
import argparse
import sys
from pathlib import Path

import pandas as pd

_REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(_REPO / "scripts"))
sys.path.insert(0, str(_REPO / "episcaf_analysis"))
from stage05_extract_metrics import epi_positions_from_contig   # noqa: E402  (design_epi, native_epi)
import compute_metrics as CM                                    # noqa: E402  (read_structure, get_chain, chain_seq, find_af3_files)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--selected", required=True, help="C1 top-n CSV from stage06_select")
    ap.add_argument("--ledger", default="results/whole_epitope_designs.csv")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    sel = pd.read_csv(args.selected)
    led = pd.read_csv(args.ledger)
    # (id, contig_id) -> design_epi (0-based construct positions of the epitope islands)
    epi = {(str(r.id), int(r.contig_id)): epi_positions_from_contig(str(r.contig_string))[0]
           for r in led.itertuples(index=False)}

    rows, n_bad = [], 0
    for r in sel.itertuples(index=False):
        key = (str(r.id), int(r.contig_id))
        de = epi.get(key)
        if de is None:
            n_bad += 1
            continue
        cif, _, _ = CM.find_af3_files(Path(r.af3_dir))
        if cif is None:
            n_bad += 1
            continue
        seq = CM.chain_seq(CM.get_chain(CM.read_structure(cif), "A"))
        if de and max(de) >= len(seq):
            n_bad += 1
            continue
        chars = [c.lower() for c in seq]
        for i in de:
            chars[i] = chars[i].upper()
        se = "".join(chars)
        rows.append({
            "id": r.id, "contig_id": int(r.contig_id),
            "predID": getattr(r, "predID", ""),
            "token": getattr(r, "predID", ""),        # id column C6/assembly key on
            "target": r.id,
            "scaffoldEPITOPE": se,
            "epitope_seq": "".join(se[i] for i in de),  # audit: the uppercased epitope
        })

    out = pd.DataFrame(rows)
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.out, index=False)
    lens = out.scaffoldEPITOPE.str.len()
    print(f"[case-encode] {len(sel)} selected -> {len(out)} encoded  (skipped {n_bad})")
    print(f"[case-encode] length range {lens.min()}-{lens.max()} (expect 103); "
          f"islands/design: {out.scaffoldEPITOPE.str.findall(r'[A-Z]+').str.len().describe()[['min','50%','max']].to_dict()}")
    print(f"[case-encode] wrote {args.out}")


if __name__ == "__main__":
    main()
