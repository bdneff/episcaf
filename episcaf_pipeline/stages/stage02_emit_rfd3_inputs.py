#!/usr/bin/env python3
"""Stage 03: Emit one RFD3 input JSON per expanded contig row.

Inputs:
- contigs.parquet (output of stage02)
- a cleaned input PDB (typically antigen-only) to provide coordinates to RFD3

Outputs:
- 02_rfd3/inputs/<design_id>.json
- 02_rfd3/inputs_manifest.csv  (design_id,json_path)

Behavior:
- uses contig_string from parquet (no regeneration)
- converts contig_string '.../.../...' -> '..., ..., ...' (RFD3 style)
- builds select_fixed_atoms (BKBN) from residue indices columns:
    1) assay_scaffolded_epitope_resindices
    2) scaffolded_epitope_resindices
    3) epitope_resindices
  indices are assumed 0-based on chain A; mapped to A(i+1).
"""

from __future__ import annotations

import csv
import json
import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List


def parse_contig_length(x):
    """Parse contig_length which may be int-like or a range string like "104-104"."""
    if x is None:
        raise ValueError("contig_length is None")
    # polars may give python int, numpy int, or str
    if isinstance(x, (int,)):
        return int(x)
    s = str(x).strip()
    if "-" in s:
        a, b = s.split("-", 1)
        a = int(a.strip())
        b = int(b.strip())
        return b  # choose max of range; for 104-104 returns 104
    return int(s)

import pandas as pd

from episcaf_pipeline.utils import abs_path, parse_index_list


def _resolve_input_pdb(args, row_dict):
    """Return absolute path to input PDB for this row.

    If args.pdb_dir is set, expect <id>.pdb in that directory.
    Otherwise fall back to args.input_pdb.
    """
    if getattr(args, "pdb_dir", None):
        base = str(row_dict.get("id"))
        cand = Path(args.pdb_dir) / f"{base}.pdb"
        if not cand.exists():
            # allow already-suffixed or alternate naming
            cand2 = Path(args.pdb_dir) / base
            if cand2.exists():
                cand = cand2
            else:
                raise FileNotFoundError(
                    f"Could not find PDB for id={base}: tried {cand} and {cand2}"
                )
        return abs_path(cand)
    if getattr(args, "input_pdb", None):
        return abs_path(Path(args.input_pdb))
    raise ValueError("Must provide either pdb_dir or input_pdb")

log = logging.getLogger(__name__)

@dataclass
class EmitRFD3Args:
    contigs_parquet: Path
    input_pdb: Path | None
    pdb_dir: Path | None
    out_dir: Path
    manifest_csv: Path
    dump_trajectory: bool = False
    prevalidate_inputs: bool = True


def _fixed_atoms_from_row(r: Dict[str, Any]) -> Dict[str, str]:
    cols = [
        "assay_scaffolded_epitope_resindices",
        "scaffolded_epitope_resindices",
        "epitope_resindices",
    ]
    idx: List[int] = []
    for c in cols:
        if c in r and r[c] is not None:
            parsed = parse_index_list(r[c])
            if parsed:
                idx = parsed
                break
    # map 0-based -> A(i+1)
    fixed = {}
    for i in idx:
        if i < 0:
            continue
        fixed[f"A{i+1}"] = "BKBN"
    return fixed


def emit_rfd3_inputs(args: EmitRFD3Args) -> None:
    df = pd.read_parquet(args.contigs_parquet)
    required = ["design_id", "id", "contig_rfd3", "contig_length", "seed"]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns in contigs parquet: {missing}")

    args.out_dir.mkdir(parents=True, exist_ok=True)
    args.manifest_csv.parent.mkdir(parents=True, exist_ok=True)


    with open(args.manifest_csv, "w", newline="") as mf:
        w = csv.writer(mf)
        w.writerow(["design_id", "json_path"])
        n = 0
        for _, row in df.iterrows():
            r = row.to_dict()
            input_pdb_abs = _resolve_input_pdb(args, r)
            design_id = str(r["design_id"])
            contig = str(r["contig_rfd3"])
            length = parse_contig_length(r["contig_length"])
            seed = int(r["seed"])

            payload = {
                design_id: {
                    "input": input_pdb_abs,
                    "contig": contig,
                    "length": f"{length}-{length}",
                    "select_fixed_atoms": _fixed_atoms_from_row(r),
                    "select_unfixed_sequence": False,
                }
            }

            out_json = args.out_dir / f"{design_id}.json"
            out_json.write_text(json.dumps(payload, indent=2))
            w.writerow([design_id, abs_path(out_json)])
            n += 1

    log.info("Emitted %d RFD3 JSONs -> %s", n, args.out_dir)
    log.info("Manifest -> %s", args.manifest_csv)
