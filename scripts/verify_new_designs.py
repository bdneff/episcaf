#!/usr/bin/env python3
"""
verify_new_designs.py -- confirm DP4 C1/C5/C6 ship OUR new RFD3 designs, not Lawson's dp2 designs.

We reran RFD3+MPNN on the (whole-epitope) contigs; the epitope POSITIONS are shared with Lawson's dp2
(same contigs) but the SCAFFOLD SEQUENCES should be our own MPNN output, i.e. DIFFERENT from Lawson's
stored `scaffolded_epitope_seq`. This checks exactly that, per selected design:

  our_seq   = design_seq (read from our run's PDB, in dp4_C1_scaffoldEPITOPE.csv)
  lawson_seq = dp2.scaffolded_epitope_seq for the same token (assay_scaffolded_epitope_id)

If we were mistakenly using Lawson's designs, the two would be IDENTICAL. Ours should DIFFER (while the
epitope residues themselves match, since the epitope is held fixed). Reports identical vs different, and
as a sanity check that the epitope AAs still agree (so the difference is scaffold, not a mismatch).

Usage (Gemini, full set):
  python scripts/verify_new_designs.py \
      --scaffoldepitope results/dp4_C1_scaffoldEPITOPE.csv \
      --dp2 /tgen_labs/altin/alphafold3/workspace/episcaf_v2_bneff/datasets/dp2.parquet
"""
from __future__ import annotations
import argparse, ast, sys
from pathlib import Path
import numpy as np
import pandas as pd


def parse_list(x):
    if isinstance(x, (list, tuple, np.ndarray)): return [int(v) for v in x]
    if x is None or (isinstance(x, float) and np.isnan(x)): return []
    return [int(v) for v in ast.literal_eval(x)] if isinstance(x, str) else list(x)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--scaffoldepitope", required=True, help="dp4_C1_scaffoldEPITOPE.csv (token + design_seq)")
    ap.add_argument("--dp2", required=True)
    ap.add_argument("--token-col", default="token")
    args = ap.parse_args()

    ce = pd.read_csv(args.scaffoldepitope, low_memory=False)
    dp2 = pd.read_parquet(args.dp2)
    dp2["assay_scaffolded_epitope_id"] = dp2["assay_scaffolded_epitope_id"].astype(str).str.lower()
    seq_lut, chunk_lut = {}, {}
    for tok, s, ch in zip(dp2["assay_scaffolded_epitope_id"], dp2["scaffolded_epitope_seq"],
                          dp2["scaffolded_epitope_chunk_resindices"]):
        if tok not in seq_lut:
            seq_lut[tok] = s; chunk_lut[tok] = parse_list(ch)

    n = n_cov = n_identical = n_different = n_epi_ok = 0
    for r in ce.itertuples(index=False):
        tok = str(getattr(r, args.token_col)).lower()
        ours = str(getattr(r, "design_seq"))
        n += 1
        law = seq_lut.get(tok)
        if law is None:
            continue
        n_cov += 1
        if ours == law:
            n_identical += 1
        else:
            n_different += 1
            # sanity: do the epitope residues (chunk positions) still agree between ours and Lawson's?
            pos = [p for p in chunk_lut.get(tok, []) if p < len(ours) and p < len(law)]
            if pos and all(ours[p].upper() == law[p].upper() for p in pos):
                n_epi_ok += 1

    print(f"[verify] selected designs: {n:,} | dp2 token coverage: {n_cov:,}")
    print(f"[verify]   IDENTICAL to Lawson (would mean using dp2 designs): {n_identical:,}")
    print(f"[verify]   DIFFERENT from Lawson (our own MPNN scaffold):      {n_different:,}")
    print(f"[verify]   of the different, epitope residues still match Lawson: {n_epi_ok:,}/{n_different:,}")
    verdict = ("OUR OWN designs (all differ from Lawson)" if n_identical == 0 and n_different > 0
               else "MIXED / SOME are Lawson's -- investigate" if n_identical else "no coverage")
    print(f"[verify] VERDICT: {verdict}")


if __name__ == "__main__":
    main()
