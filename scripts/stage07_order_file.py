#!/usr/bin/env python3
"""Stage 07c -- emit (and verify) the Twist synthesis order file.

The NN selector (step 2) writes `<LIB>_best_encodings`, which already carries the
adapter-flanked oligo in its "Nucleotide Encoding w/ Adapters" column. The synthesis
order file is therefore a two-column slice of it, not a further encoding step:

    "Seq ID"                          -> Seq ID
    "Nucleotide Encoding w/ Adapters" -> nucleotide_encoding_with_twist_adapters

What this script adds is the checking. An oligo order is expensive and irreversible, and
the one thing that can go silently wrong is the adapter: the two encoder installs ship
DIFFERENT --adapter defaults (master's are 19-mers, DP3's synthesized oligos carry
20-mers), so an unpinned run yields oligos two bases short with no error. Every row is
therefore checked to (a) start/end with the expected adapters, (b) be the expected total
length, and (c) have its core translate back to exactly the peptide it claims to encode.

Usage
-----
Emit the order file from a best-encodings table (checks every row):

    python scripts/stage07_order_file.py \
        --best-encodings DP4_best_encodings \
        --peptides data/libraries/dp4_named_peptides.csv \
        --out data/libraries/dp4_order_file.csv

Verify an existing order file (e.g. the DP3 one, as a self-test of this checker):

    python scripts/stage07_order_file.py --verify \
        --order-file episcaf_pipeline/oligo_encoding/examples/DP3_order_file.csv \
        --peptides episcaf_pipeline/oligo_encoding/examples/DP3_named_peptides.csv
"""

from __future__ import annotations

import argparse
import csv
import sys

# The Twist adapters. Default = the 20-mers John CONFIRMED (2026-07-16) as the intended DP4 adapters
# (ACCTATACTTCCAAGGCGCA / GGTGACTCTCTGTCTTGGCT), the same ones the DP3 order file carried -> 349-nt
# oligos, one base under Twist's next price tier at 350. NOT the GitHub-master 19-mer default (347 nt).
STD_PREFIX = "ACCTATACTTCCAAGGCGCA"
STD_SUFFIX = "GGTGACTCTCTGTCTTGGCT"

ID_COL = "Seq ID"
OLIGO_COL = "Nucleotide Encoding w/ Adapters"
OUT_OLIGO_COL = "nucleotide_encoding_with_twist_adapters"

_BASES = "TCAG"
_AAS = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODON_TABLE = {
    b1 + b2 + b3: _AAS[i]
    for i, (b1, b2, b3) in enumerate(
        (b1, b2, b3) for b1 in _BASES for b2 in _BASES for b3 in _BASES
    )
}


def translate(dna: str) -> str:
    if len(dna) % 3:
        raise ValueError(f"coding sequence is not a multiple of 3 ({len(dna)} nt)")
    return "".join(CODON_TABLE[dna[i : i + 3]] for i in range(0, len(dna), 3))


def load_peptides(path: str) -> dict[str, str]:
    """The encoder input: `name,seq`, no header."""
    peptides: dict[str, str] = {}
    with open(path, newline="") as fh:
        for row in csv.reader(fh):
            if not row:
                continue
            if len(row) != 2:
                raise SystemExit(f"{path}: expected 2 columns (name,seq), got {len(row)}: {row!r}")
            peptides[row[0]] = row[1]
    return peptides


def peptide_name(seq_id: str) -> str:
    """`DP4_1234_00050` -> `DP4_1234`. The selector appends `_<encoding index>`."""
    return seq_id.rsplit("_", 1)[0]


def check_rows(rows, peptides, prefix, suffix, label):
    """Check every oligo's adapters, length, and that its core translates to its peptide."""
    expected_len = None
    errors: list[str] = []
    seen: set[str] = set()

    for row in rows:
        seq_id, oligo = row["id"], row["oligo"]
        name = peptide_name(seq_id)
        pep = peptides.get(name)

        if pep is None:
            errors.append(f"{seq_id}: no peptide named {name!r} in the peptide file")
            continue
        if name in seen:
            errors.append(f"{name}: more than one encoding (expected exactly one per peptide)")
            continue
        seen.add(name)

        if not oligo.startswith(prefix):
            errors.append(f"{seq_id}: 5' adapter mismatch (starts {oligo[:len(prefix)]!r})")
            continue
        if not oligo.endswith(suffix):
            errors.append(f"{seq_id}: 3' adapter mismatch (ends {oligo[-len(suffix):]!r})")
            continue

        want_len = len(prefix) + 3 * len(pep) + len(suffix)
        if len(oligo) != want_len:
            errors.append(f"{seq_id}: oligo is {len(oligo)} nt, expected {want_len}")
            continue
        if expected_len is None:
            expected_len = len(oligo)
        elif len(oligo) != expected_len:
            errors.append(f"{seq_id}: oligo is {len(oligo)} nt, but others are {expected_len}")
            continue

        core = oligo[len(prefix) : len(oligo) - len(suffix)]
        try:
            prot = translate(core)
        except (ValueError, KeyError) as exc:
            errors.append(f"{seq_id}: cannot translate core ({exc})")
            continue
        if prot != pep:
            errors.append(f"{seq_id}: core translates to a different peptide than {name}")

    missing = set(peptides) - seen
    if missing:
        sample = ", ".join(sorted(missing)[:5])
        errors.append(f"{len(missing)} peptide(s) have no encoding (e.g. {sample})")

    print(f"{label}: {len(rows)} oligos, {len(seen)} distinct peptides, {expected_len} nt each")
    print(f"  5' adapter {prefix} ({len(prefix)} nt)")
    print(f"  3' adapter {suffix} ({len(suffix)} nt)")

    if errors:
        print(f"\nFAILED -- {len(errors)} problem(s):", file=sys.stderr)
        for err in errors[:20]:
            print(f"  {err}", file=sys.stderr)
        if len(errors) > 20:
            print(f"  ... and {len(errors) - 20} more", file=sys.stderr)
        raise SystemExit(1)

    print(f"  OK: every core translates to its peptide; one encoding per peptide.")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--best-encodings", help="the step-2 output table")
    ap.add_argument("--order-file", help="an existing order file (with --verify)")
    ap.add_argument("--peptides", required=True, help="the encoder input: name,seq, no header")
    ap.add_argument("--out", help="order file to write")
    ap.add_argument("--verify", action="store_true", help="check only, write nothing")
    ap.add_argument("--prefix", default=STD_PREFIX, help="5' Twist adapter (default: DP4 20-mer)")
    ap.add_argument("--suffix", default=STD_SUFFIX, help="3' Twist adapter (default: DP4 20-mer)")
    args = ap.parse_args()

    peptides = load_peptides(args.peptides)

    if args.verify and args.order_file:
        with open(args.order_file, newline="") as fh:
            rows = [
                {"id": r[ID_COL], "oligo": r[OUT_OLIGO_COL]} for r in csv.DictReader(fh)
            ]
        check_rows(rows, peptides, args.prefix, args.suffix, args.order_file)
        return

    if not args.best_encodings:
        raise SystemExit("need --best-encodings (or --verify with --order-file)")

    with open(args.best_encodings, newline="") as fh:
        reader = csv.DictReader(fh)
        for col in (ID_COL, OLIGO_COL):
            if col not in (reader.fieldnames or []):
                raise SystemExit(
                    f"{args.best_encodings}: missing column {col!r}. Found: {reader.fieldnames}"
                )
        rows = [{"id": r[ID_COL], "oligo": r[OLIGO_COL]} for r in reader]

    check_rows(rows, peptides, args.prefix, args.suffix, args.best_encodings)

    if args.verify:
        return
    if not args.out:
        raise SystemExit("need --out (or pass --verify to check only)")

    with open(args.out, "w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow([ID_COL, OUT_OLIGO_COL])
        for row in rows:
            writer.writerow([row["id"], row["oligo"]])
    print(f"\nwrote {args.out}  ({len(rows)} oligos)")


if __name__ == "__main__":
    main()
