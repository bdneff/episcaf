#!/usr/bin/env python3
"""
04_make_af3_jsons.py

Robust AF3 JSON generator from RFD3 *.cif.gz outputs.

- Reads manifest.csv with columns: pred_id,rfd3_cif_gz
- Extracts sequence from mmCIF _atom_site table (robust to nonstandard CIFs)
- Writes AF3 JSON files matching your working schema

No external deps required.
"""

import argparse
import csv
import gzip
import json
from pathlib import Path

AA3_TO_1 = {
    "ALA":"A","ARG":"R","ASN":"N","ASP":"D","CYS":"C",
    "GLN":"Q","GLU":"E","GLY":"G","HIS":"H","ILE":"I",
    "LEU":"L","LYS":"K","MET":"M","PHE":"F","PRO":"P",
    "SER":"S","THR":"T","TRP":"W","TYR":"Y","VAL":"V",
    # Common alternatives / unknowns
    "MSE":"M",  # selenomethionine
    "SEC":"U",  # selenocysteine
    "PYL":"O",  # pyrrolysine
    "ASX":"B",  # Asp/Asn ambiguous
    "GLX":"Z",  # Glu/Gln ambiguous
    "UNK":"X",
}

def _split_cif_tokens(line: str):
    # mmCIF tokenization for simple atom_site loops (no quoted strings expected here)
    # Works for typical RFD3 outputs.
    return line.strip().split()

def extract_sequence_from_cif_gz_atomsite(cif_gz: Path) -> str:
    """
    Parse mmCIF _atom_site loop to reconstruct polymer sequence.
    Prefers auth_* fields if present; otherwise label_*.
    """
    with gzip.open(cif_gz, "rt") as f:
        lines = f.readlines()

    # Find the _atom_site loop header
    atom_site_cols = []
    in_loop = False
    collecting_cols = False
    data_start_idx = None

    for i, line in enumerate(lines):
        s = line.strip()
        if not s:
            continue

        if s == "loop_":
            # start of a loop; we don't know which yet
            in_loop = True
            collecting_cols = True
            atom_site_cols = []
            data_start_idx = None
            continue

        if in_loop and collecting_cols:
            if s.startswith("_atom_site."):
                atom_site_cols.append(s)
                continue
            # once we hit a non-_atom_site line, either this loop isn't atom_site,
            # or we've finished collecting cols and are at data.
            if atom_site_cols:
                data_start_idx = i
                collecting_cols = False
                break
            else:
                # loop_, but not atom_site; reset
                in_loop = False
                collecting_cols = False
                atom_site_cols = []
                continue

    if not atom_site_cols or data_start_idx is None:
        raise RuntimeError(f"No _atom_site loop found in {cif_gz}")

    col_index = {c: idx for idx, c in enumerate(atom_site_cols)}

    # Choose best available columns for chain/residue identification
    chain_col = None
    for c in ("_atom_site.auth_asym_id", "_atom_site.label_asym_id"):
        if c in col_index:
            chain_col = c
            break

    seqid_col = None
    for c in ("_atom_site.auth_seq_id", "_atom_site.label_seq_id"):
        if c in col_index:
            seqid_col = c
            break

    comp_col = None
    for c in ("_atom_site.auth_comp_id", "_atom_site.label_comp_id"):
        if c in col_index:
            comp_col = c
            break

    if chain_col is None or seqid_col is None or comp_col is None:
        raise RuntimeError(
            f"Missing required atom_site cols in {cif_gz} "
            f"(chain={chain_col}, seqid={seqid_col}, comp={comp_col})"
        )

    # Collect residues in encounter order (atom_site contains many atoms per residue)
    residues_by_chain = {}  # chain -> ordered dict-like list of (seqid, comp)
    seen = set()            # (chain, seqid)

    for line in lines[data_start_idx:]:
        s = line.strip()
        if not s:
            continue
        if s.startswith("_") or s == "loop_" or s.startswith("data_"):
            # End of loop data
            break
        if s.startswith("#"):
            break

        toks = _split_cif_tokens(line)
        if len(toks) < len(atom_site_cols):
            # mmCIF can wrap long lines; RFD3 outputs typically don't.
            # Skip safely.
            continue

        chain = toks[col_index[chain_col]]
        seqid = toks[col_index[seqid_col]]
        comp = toks[col_index[comp_col]]

        if seqid == "." or seqid == "?":
            continue

        key = (chain, seqid)
        if key in seen:
            continue
        seen.add(key)

        residues_by_chain.setdefault(chain, []).append((seqid, comp))

    if not residues_by_chain:
        raise RuntimeError(f"Failed to extract residues from {cif_gz}")

    # Prefer chain "A" if present; else first chain encountered (sorted for determinism)
    chain = "A" if "A" in residues_by_chain else sorted(residues_by_chain.keys())[0]
    residue_list = residues_by_chain[chain]

    seq = []
    for _, comp in residue_list:
        comp = comp.upper()
        seq.append(AA3_TO_1.get(comp, "X"))

    s = "".join(seq).strip()
    if not s:
        raise RuntimeError(f"Failed to extract sequence from {cif_gz}")
    return s

def make_af3_payload(name: str, seq: str, seed: int) -> dict:
    msa = f">query\n{seq}\n"
    return {
        "name": name,
        "modelSeeds": [seed],
        "sequences": [
            {
                "protein": {
                    "id": "A",
                    "sequence": seq,
                    "unpairedMsa": msa,
                    "pairedMsa": msa,
                    "templates": []
                }
            }
        ],
        "dialect": "alphafold3",
        "version": 1
    }

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--manifest", required=True, help="CSV with columns: pred_id,rfd3_cif_gz")
    ap.add_argument("--out_dir", required=True, help="Output dir for AF3 JSONs")
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--limit", type=int, default=0)
    args = ap.parse_args()

    manifest = Path(args.manifest).resolve()
    out_dir = Path(args.out_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    n = 0
    with open(manifest, newline="") as f:
        r = csv.DictReader(f)
        for row_i, row in enumerate(r, start=1):
            if args.limit and row_i > args.limit:
                break

            pred_id = row["pred_id"].strip()
            cif_gz = Path(row["rfd3_cif_gz"]).resolve()

            seq = extract_sequence_from_cif_gz_atomsite(cif_gz)
            payload = make_af3_payload(pred_id, seq, args.seed)

            out_path = out_dir / f"{pred_id}.json"
            with open(out_path, "w") as out:
                json.dump(payload, out, indent=2)

            n += 1

    print(f"Wrote {n} AF3 JSONs to: {out_dir}")

if __name__ == "__main__":
    main()

