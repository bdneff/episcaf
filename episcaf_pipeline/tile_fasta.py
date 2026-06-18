#!/usr/bin/env python3
"""
tile_fasta.py  --  tile protein sequences into overlapping k-mers and place each in the
assay construct. Throw any FASTA at it; specify the mer length and the step size.

Construct (constant filler + TEV-like protease site + the tile), e.g. for a 30mer:
  GSGAGSGAGSGAGSGAGSGAGSGAGSGAGSGAGSGAGSGAGSGAGSGAGSGAGSGAGSGAGSGAGENLYFQGA[30mer]
The 73-residue prefix is constant, so every member is synthesized as a constant-length
construct (103 for 30mers) and the protease site lets it cleave to the bare k-mer in the
final step -- which is what matches the linear-epitope assay.

Output columns match the DP2 annotated library:
  library_member, sequence, category, model, designedSequence, designedSequenceLength,
  design_ID, target
  - sequence            : prefix + tile (the full construct)
  - category            : e.g. "tiled30mer"
  - model               : "(none)" for unscaffolded tiles (controls)
  - designedSequence    : the bare tile (what gets displayed after cleavage)
  - design_ID           : tile start position in the antigen, 1-indexed
  - target              : FASTA record id (the antigen)
  - library_member      : <id_prefix>_<running counter>; final numbering is set when the
                          full library is assembled across categories, so treat this as
                          provisional unless --start-id is coordinated.

Tiling matches the existing tiled_library_*: windows at start = 1, 1+step, 1+2*step, ...
while start+mer-1 <= L (the C-terminal remainder shorter than one step is not padded).

Usage:
  python -m episcaf_pipeline.tile_fasta --fasta antigens.fasta --mer 30 --step 6 \
      --category tiled30mer --id-prefix DP4 --out runs/.../tiled30mers.csv
"""
from __future__ import annotations
import argparse
from pathlib import Path
import pandas as pd

CONSTRUCT_PREFIX = ("GSGAGSGAGSGAGSGAGSGAGSGAGSGAGSGAGSGAGSGAGSGAGSGAGSGAGSGAGSGAGSGAG"
                    "ENLYFQGA")   # 73 aa: GS filler + ENLYFQGA protease site (verified vs DP2)


def read_fasta(path: Path):
    """Yield (id, sequence). id = first whitespace-delimited token of the header."""
    rid, seq = None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if rid is not None:
                    yield rid, "".join(seq)
                rid = line[1:].strip().split()[0]
                seq = []
            elif line.strip():
                seq.append(line.strip())
    if rid is not None:
        yield rid, "".join(seq)


def read_parquet_records(path: Path, id_col: str, seq_col: str, dedupe: bool):
    """Yield (id, sequence) from a parquet/ledger (e.g. dp2.parquet antigen_seq).
    dedupe collapses to unique sequences (tile each antigen protein once)."""
    df = pd.read_parquet(path)
    for c in (id_col, seq_col):
        if c not in df.columns:
            raise SystemExit(f"column {c!r} not in {path} (have: {list(df.columns)})")
    sub = df[[id_col, seq_col]].dropna()
    if dedupe:
        sub = sub.drop_duplicates(subset=[seq_col])     # one row per distinct antigen
    else:
        sub = sub.drop_duplicates(subset=[id_col])
    for _, r in sub.iterrows():
        yield str(r[id_col]), str(r[seq_col])


def tile_sequence(seq: str, mer: int, step: int, include_cterm: bool = True):
    """Yield (start_1indexed, tile) overlapping windows.
    If include_cterm and the last regular window does not reach the C-terminus, append a
    final window ending at the C-terminus (covers the trailing residues), matching the
    convention in the existing 6M0J 12-mer library."""
    L = len(seq)
    if L < mer:
        return
    last = None
    for i in range(0, L - mer + 1, step):
        yield i + 1, seq[i:i + mer]
        last = i
    if include_cterm and last is not None and last + mer < L:
        i = L - mer
        yield i + 1, seq[i:i + mer]


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    src = ap.add_mutually_exclusive_group(required=True)
    src.add_argument("--fasta", type=Path, help="input FASTA (e.g. the 3 tiled antigens)")
    src.add_argument("--parquet", type=Path, help="ledger with antigen sequences (e.g. dp2.parquet)")
    ap.add_argument("--id-col", default="id", help="parquet: id column")
    ap.add_argument("--seq-col", default="antigen_seq", help="parquet: sequence column")
    ap.add_argument("--dedupe", action="store_true",
                    help="parquet: tile each distinct antigen sequence once (unique proteins)")
    ap.add_argument("--mer", required=True, type=int)
    ap.add_argument("--step", required=True, type=int)
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--category", default=None, help="default: tiled<mer>mer")
    ap.add_argument("--model", default="(none)", help="(none) for unscaffolded controls")
    ap.add_argument("--id-prefix", default="DP4")
    ap.add_argument("--start-id", type=int, default=1)
    ap.add_argument("--construct-prefix", default=CONSTRUCT_PREFIX,
                    help="set to empty string to emit bare tiles (no construct)")
    ap.add_argument("--no-cterm", action="store_true",
                    help="do not add a final C-terminal coverage tile")
    args = ap.parse_args()

    records = (read_parquet_records(args.parquet, args.id_col, args.seq_col, args.dedupe)
               if args.parquet else read_fasta(args.fasta))
    category = args.category or f"tiled{args.mer}mer"
    rows, n = [], args.start_id
    n_prot = 0
    for rid, seq in records:
        n_prot += 1
        if len(seq) < args.mer:
            print(f"[warn] {rid}: length {len(seq)} < mer {args.mer}, skipped")
            continue
        for start, tile in tile_sequence(seq, args.mer, args.step, include_cterm=not args.no_cterm):
            rows.append(dict(
                library_member=f"{args.id_prefix}_{n}",
                sequence=args.construct_prefix + tile,
                category=category,
                model=args.model,
                designedSequence=tile,
                designedSequenceLength=len(tile),
                design_ID=start,
                target=rid,
            ))
            n += 1

    df = pd.DataFrame(rows, columns=["library_member", "sequence", "category", "model",
                                     "designedSequence", "designedSequenceLength",
                                     "design_ID", "target"])
    args.out.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.out, index=False)
    print(f"{n_prot} proteins -> {len(df):,} {category} members  ->  {args.out}")


if __name__ == "__main__":
    main()
