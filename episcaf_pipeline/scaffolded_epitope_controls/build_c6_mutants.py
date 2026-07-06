#!/usr/bin/env python3
"""
build_c6_mutants.py -- DP4 scaffolded-epitope controls (component C6).

Python port of John Altin's DP3 R scripts (reference_dp3/), with his two DP4 modifications:
  1. island->Ala substitutes EVERY epitope residue of the focal island (DP3 did every other);
  2. scaffold-disruption keeps only _scaffoldMutX4 (drop _scaffoldMutX1).
Plus: the scaffold 6-mer placement is SEEDED (reproducible libraries) and FAIL-SOFT (a design that
can't fit 4 disruptions is logged and skipped, not a hard error like the R default).

This is NOT new scaffolding -- it substitutes residues on an EXISTING design's case-encoded sequence.
Input `scaffoldEPITOPE` string: UPPERCASE = epitope residue, lowercase = scaffold; each contiguous
uppercase run = one epitope island. That casing IS the island annotation (no residue-index map needed).

Flavors emitted per design (the unmutated original is NOT emitted -- it is already shipped as C1):
  <id>_epitopeIsland1Mut>A   island 1 -> all Ala
  <id>_epitopeIsland2Mut>A   island 2 -> all Ala   (only if a 2nd island exists; else skipped)
  <id>_scaffoldMutX4         4 lowercase 6-mers (each >=4 from any epitope residue) -> "PPDDGG"

Usage:
  # reproduce/sanity-check against John's DP3 input:
  python episcaf_pipeline/scaffolded_epitope_controls/build_c6_mutants.py \
      --input episcaf_pipeline/scaffolded_epitope_controls/reference_dp3/scaffolded403.csv \
      --id-col Design_ID --target-col Target --seq-col scaffoldEPITOPE \
      --out results/dp4_C6_controls.csv
  # DP3-compat (every-other Ala + include X1) to validate the port matches John exactly:
  ... --every-other --include-x1
"""
from __future__ import annotations
import argparse, random, re, sys
from pathlib import Path
import pandas as pd

REPLACEMENT = "PPDDGG"       # John's structure-disrupting hexamer
SYNTH_LEN = 103              # synthesized peptide length (John: substring(seq, 1, 103))


def islands(seq: str):
    """Contiguous uppercase runs = epitope islands; returns list of (start, end_exclusive)."""
    return [(m.start(), m.end()) for m in re.finditer(r"[A-Z]+", seq)]


def mutate_to_A(seq: str, n: int, every_other: bool = False) -> str:
    """Substitute the n-th island (1-based) to alanine. '' if island n doesn't exist."""
    runs = islands(seq)
    if n > len(runs):
        return ""
    s, e = runs[n - 1]
    chars = list(seq)
    for k, i in enumerate(range(s, e)):
        if every_other and k % 2 != 0:      # DP3: only odd (1-based) positions within the chunk
            continue
        chars[i] = "A"
    return "".join(chars)


def mutate6(seq: str, rng: random.Random, n: int = 4, max_tries: int = 1000):
    """Replace n lowercase 6-mers (each position >=4 away from any uppercase/epitope residue) with
    REPLACEMENT. Returns the mutated string, or None if it can't place all n (fail-soft)."""
    chars = list(seq)
    for _ in range(n):
        placed = False
        for _attempt in range(max_tries):
            L = len(chars)
            is_lower = [c.isalpha() and c.islower() for c in chars]
            caps = [i for i, c in enumerate(chars) if c.isalpha() and c.isupper()]
            valid = []
            for i in range(0, L - 5):
                six = range(i, i + 6)
                if not all(is_lower[j] for j in six):
                    continue
                if caps and any(abs(j - c) <= 4 for j in six for c in caps):
                    continue
                valid.append(i)
            if valid:
                start = rng.choice(valid)
                chars[start:start + 6] = list(REPLACEMENT)
                placed = True
                break
        if not placed:
            return None
    return "".join(chars)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--input", required=True, help="CSV/parquet with a case-encoded scaffoldEPITOPE column")
    ap.add_argument("--id-col", default="Design_ID")
    ap.add_argument("--target-col", default="Target")
    ap.add_argument("--seq-col", default="scaffoldEPITOPE")
    ap.add_argument("--out", required=True)
    ap.add_argument("--drop-targets", default="", help="comma prefixes of target ids to exclude (e.g. 4xwo,7a3t)")
    ap.add_argument("--seed", type=int, default=0, help="RNG seed for scaffold 6-mer placement")
    ap.add_argument("--every-other", action="store_true", help="DP3-compat: Ala every OTHER island residue")
    ap.add_argument("--include-x1", action="store_true", help="also emit _scaffoldMutX1 (DP3-compat)")
    args = ap.parse_args()

    p = Path(args.input)
    df = pd.read_parquet(p) if p.suffix == ".parquet" else pd.read_csv(p, low_memory=False)
    for c in (args.id_col, args.seq_col):
        if c not in df.columns:
            sys.exit(f"[c6] column {c!r} not in input. have: {list(df.columns)}")
    tgt = args.target_col if args.target_col in df.columns else None
    if args.drop_targets and tgt:
        pref = tuple(p.strip().lower() for p in args.drop_targets.split(",") if p.strip())
        n0 = len(df)
        df = df[~df[tgt].astype(str).str.lower().str.startswith(pref)].copy()
        print(f"[c6] dropped {n0-len(df)} rows with target in {pref} -> {len(df)} base designs")
    rng = random.Random(args.seed)

    rows, n_isl2, n_x4_fail = [], 0, 0
    scaff_ns = ([1, 4] if args.include_x1 else [4])
    for r in df.itertuples(index=False):
        did = getattr(r, args.id_col); seq = str(getattr(r, args.seq_col))
        target = getattr(r, tgt) if tgt else ""
        if not seq or seq == "nan":
            continue
        flavors = [(f"{did}_epitopeIsland1Mut>A", mutate_to_A(seq, 1, args.every_other)),
                   (f"{did}_epitopeIsland2Mut>A", mutate_to_A(seq, 2, args.every_other))]
        for nX in scaff_ns:
            flavors.append((f"{did}_scaffoldMutX{nX}", mutate6(seq, rng, n=nX)))
        for name, mut in flavors:
            if not mut:                       # '' (no such island) or None (X4 couldn't place)
                if name.endswith("Island2Mut>A"): pass
                elif "scaffoldMut" in name: n_x4_fail += 1
                continue
            if name.endswith("Island2Mut>A"): n_isl2 += 1
            up = mut.upper()
            rows.append(dict(sequence=up[:SYNTH_LEN], category="scaffoldedAbEpitope", model="RFD",
                             designedSequence=up, designedSequenceLength=len(up),
                             design_ID=name, target=target, scaffoldEPITOPE=mut))

    out = pd.DataFrame(rows)
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.out, index=False)
    n_in = len(df)
    print(f"[c6] {n_in} input designs -> {len(out)} control constructs")
    print(f"[c6]   island1->A: {(out.design_ID.str.endswith('Island1Mut>A')).sum()}  "
          f"island2->A: {n_isl2} (dual-island only)  "
          f"scaffoldMutX4: {(out.design_ID.str.contains('scaffoldMutX4')).sum()}")
    if n_x4_fail:
        print(f"[c6]   WARNING: {n_x4_fail} designs could not place 4 disruption hexamers (skipped X4)")
    print(f"[c6] wrote {args.out}")


if __name__ == "__main__":
    main()
