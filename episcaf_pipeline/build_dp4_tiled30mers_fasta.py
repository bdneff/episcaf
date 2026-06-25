#!/usr/bin/env python3
"""
build_dp4_tiled30mers_fasta.py -- DP4 tiled-30mer linear-control library from the PDB FASTA
sequences (gap-free), the corrected source.

Why a new builder: build_dp4_tiled30mers.py tiled dp2 `antigen_seq` (the crystal-RESOLVED
chain), which is truncated at the termini and, for one antigen (6okm), splices across an
internal unresolved gap (a fake junction) -- see scripts/verify_antigen_seq_vs_pdb_fasta.py and
results/antigen_seq_vs_fasta.csv. The PDB FASTA chain is the full deposited sequence with NO
gaps; we tile that instead.

Sequence rule (per antigen):
  - fetch the RCSB FASTA, pick the antigen chain by matching dp2 `antigen_seq` (the chain that
    contains it; episcaf_pipeline.pdb_fasta.pick_antigen_chain);
  - native window = FASTA_chain[start:] where `start` is where `antigen_seq` begins. This drops
    any N-terminal expression tag/leader (matching the DP2 native start, which began at the
    resolved residue, not the tag) while RECOVERING the gap-free C-terminus the resolved
    sequence had lost. The window is a contiguous slice of a gap-free FASTA -> gap-free.
  - antigen 6okm has no exact substring match (internal splice in antigen_seq); it is tiled from
    the full best-match chain and flagged status=needs_review for a manual look.

The 3 tiled controls (1d2k/4wat/6m0j) are John-provided full sequences, kept as-is.

Output (same 8-column DP2 annotated format as the approved file; target = dp2 patch id):
  data/libraries/dp4_tiled30mers_fasta.csv
  data/libraries/dp4_tiled30mers_fasta_summary.csv   (per-antigen: old/new length, residues
                                                       recovered, tiles, protein, organism, status)

Reproducible: RCSB responses are cached under data/sequences/pdb_fasta_cache/ (committed-ignored
or kept), so re-runs are deterministic and need no network. Each tile is asserted to be a
substring of its source FASTA chain (no fake junctions).

Usage:
  python -m episcaf_pipeline.build_dp4_tiled30mers_fasta
"""
from __future__ import annotations
import argparse
from pathlib import Path
import pandas as pd

from episcaf_pipeline.tile_fasta import read_fasta, tile_sequence, CONSTRUCT_PREFIX
from episcaf_pipeline.pdb_fasta import fetch_fasta, parse_chains, pick_antigen_chain

import sys
ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
from configs.paths import DP2_PARQUET_LOCAL  # noqa: E402

MER, STEP, CATEGORY, ID_PREFIX = 30, 6, "tiled30mer", "DP4"
TILED_ANTIGENS = {"1d2k", "4wat", "6m0j"}


def antigen_fasta_records(dp2_path):
    """For each dp2 mAb antigen, yield (patch_id, native_window_seq, meta)."""
    dp2 = pd.read_parquet(dp2_path).drop_duplicates("id")
    for _, r in dp2.iterrows():
        pid, aseq = str(r["id"]), str(r["antigen_seq"])
        chains = parse_chains(fetch_fasta(pid[:4]))
        pick = pick_antigen_chain(aseq, chains)
        full = pick["seq"]
        native = full[pick["start"]:] if pick["start"] is not None else full
        meta = dict(old_len=len(aseq), new_len=len(native), fasta_len=len(full),
                    recovered=len(native) - len(aseq), protein=pick["protein"],
                    organism=pick["organism"], status=pick["status"], chain_seq=full)
        yield pid, native, meta


def build(dp2_path, control_fasta):
    rows, summ = [], []
    # mAb antigens from PDB FASTA (gap-free), then the 3 controls
    sources = list(antigen_fasta_records(dp2_path))
    for rid, seq in read_fasta(control_fasta):
        sources.append((rid, seq, dict(old_len="", new_len=len(seq), fasta_len="",
                                       recovered="", protein="(control)", organism="",
                                       status="control", chain_seq=seq)))
    for rid, seq, meta in sources:
        n = 0
        for start, tile in tile_sequence(seq, MER, STEP, include_cterm=True):
            assert tile in meta["chain_seq"], f"{rid}: tile not a substring of source (gap!)"
            rows.append(dict(sequence=CONSTRUCT_PREFIX + tile, category=CATEGORY, model="(none)",
                             designedSequence=tile, designedSequenceLength=len(tile),
                             design_ID=start, target=rid))
            n += 1
        summ.append(dict(target=rid, **{k: meta[k] for k in
                    ("old_len", "new_len", "recovered", "protein", "organism", "status")},
                    tiles=n))

    df = pd.DataFrame(rows)
    df.insert(0, "library_member", [f"{ID_PREFIX}_{i}" for i in range(1, len(df) + 1)])
    df = df[["library_member", "sequence", "category", "model", "designedSequence",
             "designedSequenceLength", "design_ID", "target"]]
    summary = pd.DataFrame(summ)
    return df, summary


def main():
    data = Path(__file__).resolve().parent.parent / "data"
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dp2", type=Path, default=DP2_PARQUET_LOCAL)
    ap.add_argument("--controls", type=Path, default=data / "sequences/control_antigens.fasta")
    ap.add_argument("--out", type=Path, default=data / "libraries/dp4_tiled30mers_fasta.csv")
    args = ap.parse_args()

    df, summary = build(args.dp2, args.controls)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.out, index=False)
    summ_path = args.out.with_name(args.out.stem + "_summary.csv")
    summary.to_csv(summ_path, index=False)

    mab = summary[~summary.target.isin(TILED_ANTIGENS) & (summary.status != "control")]
    rec = pd.to_numeric(mab["recovered"], errors="coerce")
    print(f"{len(summary)} antigens ({len(mab)} mAb + {len(summary)-len(mab)} control) "
          f"-> {len(df)} tiled-30mers")
    print(f"  residues recovered vs resolved antigen_seq: total {int(rec.sum())}, "
          f"median {rec.median():.0f}, max {int(rec.max())}")
    nr = summary[summary.status == "needs_review"]
    if len(nr):
        print(f"  NEEDS REVIEW ({len(nr)}): {', '.join(nr.target)} "
              f"(no exact FASTA substring -- internal gap in resolved seq; tiled from full chain)")
    print(f"  library -> {args.out}")
    print(f"  summary -> {summ_path}")


if __name__ == "__main__":
    main()
