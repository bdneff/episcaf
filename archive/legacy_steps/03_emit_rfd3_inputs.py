#!/usr/bin/env python3
"""
03_emit_rfd3_inputs.py

Reads the expanded contigs table (01_contigs/contigs.parquet) and emits one
RFD3 input JSON per design into 02_rfd3/inputs/.

Key behavior you asked for:
- Uses contig_string from parquet (no regeneration)
- Converts contig_string "15-15/A1-16/..." -> "15-15,A1-16,..." (RFD3 style)
- Adds select_fixed_atoms automatically (BKBN) from:
    1) assay_scaffolded_epitope_resindices (preferred)
    2) scaffolded_epitope_resindices
    3) epitope_resindices
  Assumes indices are 0-based positions on antigen chain A and maps to A(i+1).
- Writes absolute input paths (SLURM-proof)
- Also writes a manifest.csv with design_id,json_path

Example output JSON (one per design):
{
  "DESIGN_ID": {
    "input": "/abs/path/to/cleaned/2h32_0P.pdb",
    "input_structure": "/abs/path/to/cleaned/2h32_0P.pdb",
    "contig": "15-15,A1-16,57-57,A81-81,15-15",
    "length": "104-104",
    "seed": 0,
    "dump_trajectory": false,
    "prevalidate_inputs": true,
    "select_fixed_atoms": { "A1": "BKBN", ... },
    "select_unfixed_sequence": false
  }
}
"""

import argparse
import csv
import json
import os
import re
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Tuple

import pandas as pd


# -----------------------------
# helpers
# -----------------------------
_INT_RE = re.compile(r"-?\d+")


def _parse_index_list(x: Any) -> List[int]:
    """
    Parse list-like columns that may be:
      - python list
      - numpy array
      - string like "[ 0  1  2 15 80]"
      - string like "[0, 1, 2, 15, 80]"
    Returns a list[int].
    """
    if x is None or (isinstance(x, float) and pd.isna(x)):
        return []
    if isinstance(x, (list, tuple)):
        out = []
        for v in x:
            try:
                out.append(int(v))
            except Exception:
                pass
        return out
    # pandas may give numpy scalars/arrays; treat as iterable if possible
    try:
        # numpy arrays have .tolist()
        if hasattr(x, "tolist"):
            return _parse_index_list(x.tolist())
    except Exception:
        pass

    s = str(x).strip()
    if not s:
        return []
    nums = _INT_RE.findall(s)
    return [int(n) for n in nums]


def _contig_to_rfd3(contig_string: str) -> str:
    """
    Convert parquet contig_string format (slash-separated) to RFD3 format (comma-separated).
    Example: "15-15/A1-16/57-57/A81-81/15-15" -> "15-15,A1-16,57-57,A81-81,15-15"
    """
    if contig_string is None:
        return ""
    s = str(contig_string).strip()
    if not s:
        return ""
    # Accept either "/" or "," already
    if "/" in s and "," not in s:
        parts = [p.strip() for p in s.split("/") if p.strip()]
        return ",".join(parts)
    # If it already has commas, keep it (just normalize spaces)
    parts = [p.strip() for p in s.replace("/", ",").split(",") if p.strip()]
    return ",".join(parts)


def _length_to_rfd3(contig_length: Any) -> str:
    """
    contig_length sometimes already looks like "104-104". Keep as string.
    """
    if contig_length is None or (isinstance(contig_length, float) and pd.isna(contig_length)):
        return ""
    return str(contig_length).strip()


def _abs_cleaned_pdb(cleaned_pdb_dir: Path, abdb_id: str) -> str:
    """
    Build absolute path to cleaned ABDB PDB file.
    """
    fname = f"{abdb_id}.pdb"
    p = (cleaned_pdb_dir / fname).resolve()
    return str(p)


def _fixed_atoms_from_indices(indices_0based: Iterable[int], chain: str = "A") -> Dict[str, str]:
    """
    Map 0-based positions -> PDB residue numbers starting at 1: i -> (i+1)
    Returns {"A149": "BKBN", ...}
    """
    fixed: Dict[str, str] = {}
    for i in sorted(set(int(v) for v in indices_0based)):
        resid = i + 1
        fixed[f"{chain}{resid}"] = "BKBN"
    return fixed


def _pick_epitope_indices(row: Dict[str, Any]) -> Tuple[List[int], str]:
    """
    Choose the best available epitope index column.
    Returns (indices, source_colname)
    """
    for col in ("assay_scaffolded_epitope_resindices", "scaffolded_epitope_resindices", "epitope_resindices"):
        if col in row:
            vals = _parse_index_list(row.get(col))
            if vals:
                return vals, col
    return [], "none"


# -----------------------------
# main
# -----------------------------
def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--contigs_parquet", required=True, help="01_contigs/contigs.parquet")
    ap.add_argument("--out_dir", required=True, help="02_rfd3/inputs (directory for JSON inputs)")
    ap.add_argument("--cleaned_pdb_dir", required=True, help="Dir containing cleaned/*.pdb (ABDB cleaned)")
    ap.add_argument("--limit", type=int, default=0, help="If >0, emit only first N designs (quick test)")
    ap.add_argument("--dump_trajectory", action="store_true", help="Set dump_trajectory=true in JSON")
    ap.add_argument("--prevalidate_inputs", action="store_true", help="Set prevalidate_inputs=true in JSON (recommended)")
    args = ap.parse_args()

    contigs_path = Path(args.contigs_parquet).resolve()
    out_dir = Path(args.out_dir).resolve()
    cleaned_pdb_dir = Path(args.cleaned_pdb_dir).resolve()

    out_dir.mkdir(parents=True, exist_ok=True)

    df = pd.read_parquet(contigs_path)

    required = ["design_id", "id", "contig_string", "contig_length", "seed"]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns in contigs parquet: {missing}")

    if args.limit and args.limit > 0:
        df = df.head(args.limit).copy()

    manifest_path = out_dir.parent / "inputs_manifest.csv"

    n_written = 0
    with open(manifest_path, "w", newline="") as mf:
        w = csv.writer(mf)
        w.writerow(["design_id", "json_path"])

        for _, r in df.iterrows():
            row = r.to_dict()
            design_id = str(row["design_id"])
            abdb_id = str(row["id"])

            contig = _contig_to_rfd3(row.get("contig_string"))
            length = _length_to_rfd3(row.get("contig_length"))

            # absolute input paths
            input_pdb = _abs_cleaned_pdb(cleaned_pdb_dir, abdb_id)
            if not os.path.exists(input_pdb):
                raise FileNotFoundError(f"Missing cleaned PDB for id={abdb_id}: {input_pdb}")

            # fixed atoms from epitope indices
            ep_inds, src = _pick_epitope_indices(row)
            select_fixed_atoms = _fixed_atoms_from_indices(ep_inds, chain="A")

            # Build spec

            spec = {
                "input": input_pdb,
                "contig": contig,
                "length": length,
                "select_fixed_atoms": select_fixed_atoms,
                "select_unfixed_sequence": False,
                }
            payload = {design_id: spec}

            json_path = out_dir / f"{design_id}.json"
            with open(json_path, "w") as f:
                json.dump(payload, f, indent=2, sort_keys=False)

            w.writerow([design_id, str(json_path)])
            n_written += 1

    print(f"Wrote {n_written} JSON inputs into: {out_dir}")
    print(f"Wrote manifest: {manifest_path}")


if __name__ == "__main__":
    main()

