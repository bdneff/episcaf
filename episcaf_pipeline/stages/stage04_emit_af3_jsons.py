#!/usr/bin/env python3
"""Stage 04: Create AlphaFold3 JSON inputs from RFD3 output structures.

By default this stage scans:
  <run_dir>/02_rfd3/outputs/**/*.cif.gz
and writes:
  <run_dir>/03_af3/inputs/<pred_id>.json
plus:
  <run_dir>/03_af3/inputs_manifest.csv  (pred_id,json_path,rfd3_cif_gz)

pred_id defaults to the cif.gz stem (filename without .cif.gz).

This is intentionally lightweight and reproducible: given the same RFD3 outputs,
the AF3 JSONs are deterministic.
"""

from __future__ import annotations

import csv
import gzip
import json
import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Tuple

from episcaf_pipeline.utils import abs_path

log = logging.getLogger(__name__)

AA3_TO_1 = {
    "ALA":"A","ARG":"R","ASN":"N","ASP":"D","CYS":"C",
    "GLN":"Q","GLU":"E","GLY":"G","HIS":"H","ILE":"I",
    "LEU":"L","LYS":"K","MET":"M","PHE":"F","PRO":"P",
    "SER":"S","THR":"T","TRP":"W","TYR":"Y","VAL":"V",
    "MSE":"M","SEC":"U","PYL":"O","ASX":"B","GLX":"Z","UNK":"X",
}

@dataclass
class EmitAF3Args:
    rfd3_outputs_dir: Path
    out_dir: Path
    manifest_csv: Path
    seed: int = 42
    limit: int = 0  # 0 = no limit


def _split_cif_tokens(line: str) -> List[str]:
    return line.strip().split()

def extract_sequence_from_cif_gz_atomsite(cif_gz: Path) -> str:
    with gzip.open(cif_gz, "rt") as f:
        lines = f.readlines()

    atom_site_cols = []
    in_loop = False
    collecting_cols = False
    data_start_idx = None

    for i, line in enumerate(lines):
        s = line.strip()
        if not s:
            continue
        if s == "loop_":
            in_loop = True
            collecting_cols = True
            atom_site_cols = []
            data_start_idx = None
            continue
        if in_loop and collecting_cols:
            if s.startswith("_atom_site."):
                atom_site_cols.append(s)
                continue
            if atom_site_cols:
                data_start_idx = i
                collecting_cols = False
                break
            in_loop = False
            collecting_cols = False
            atom_site_cols = []

    if not atom_site_cols or data_start_idx is None:
        raise RuntimeError(f"No _atom_site loop found in {cif_gz}")

    col_index = {c: j for j, c in enumerate(atom_site_cols)}

    # choose columns
    chain_col = None
    for c in ("_atom_site.auth_asym_id", "_atom_site.label_asym_id"):
        if c in col_index:
            chain_col = c; break
    seqid_col = None
    for c in ("_atom_site.auth_seq_id", "_atom_site.label_seq_id"):
        if c in col_index:
            seqid_col = c; break
    comp_col = None
    for c in ("_atom_site.auth_comp_id", "_atom_site.label_comp_id"):
        if c in col_index:
            comp_col = c; break

    if chain_col is None or seqid_col is None or comp_col is None:
        raise RuntimeError(f"Missing required atom_site cols in {cif_gz}")

    residues = []
    seen = set()

    for line in lines[data_start_idx:]:
        s = line.strip()
        if not s:
            continue
        if s.startswith("_") or s == "loop_" or s.startswith("data_"):
            break
        toks = _split_cif_tokens(line)
        if len(toks) < len(atom_site_cols):
            continue
        chain = toks[col_index[chain_col]]
        seqid = toks[col_index[seqid_col]]
        comp = toks[col_index[comp_col]]
        key = (chain, seqid)
        if key in seen:
            continue
        seen.add(key)
        residues.append((chain, seqid, comp))

    # prefer chain A if present, else most-populated chain
    by_chain: Dict[str, List[str]] = {}
    for chain, _, comp in residues:
        by_chain.setdefault(chain, []).append(comp)

    if "A" in by_chain:
        comps = by_chain["A"]
    else:
        comps = max(by_chain.values(), key=len)

    seq = "".join(AA3_TO_1.get(c.upper(), "X") for c in comps)
    if not seq:
        raise RuntimeError(f"Failed to extract sequence from {cif_gz}")
    return seq

def make_af3_payload(name: str, seq: str, seed: int) -> dict:
    msa = f">query\n{seq}\n"
    return {
        "name": name,
        "modelSeeds": [seed],
        "sequences": [
            {"protein": {"id": "A", "sequence": seq, "unpairedMsa": msa, "pairedMsa": msa, "templates": []}}
        ],
        "dialect": "alphafold3",
        "version": 1
    }

def emit_af3_jsons(args: EmitAF3Args) -> None:
    args.out_dir.mkdir(parents=True, exist_ok=True)
    args.manifest_csv.parent.mkdir(parents=True, exist_ok=True)

    cif_files = sorted(args.rfd3_outputs_dir.rglob("*.cif.gz"))
    if args.limit and args.limit > 0:
        cif_files = cif_files[: args.limit]

    with open(args.manifest_csv, "w", newline="") as mf:
        w = csv.writer(mf)
        w.writerow(["pred_id", "json_path", "rfd3_cif_gz"])
        n=0
        for cif_gz in cif_files:
            pred_id = cif_gz.name.replace(".cif.gz","")
            seq = extract_sequence_from_cif_gz_atomsite(cif_gz)
            payload = make_af3_payload(pred_id, seq, args.seed)
            out_path = args.out_dir / f"{pred_id}.json"
            out_path.write_text(json.dumps(payload, indent=2))
            w.writerow([pred_id, abs_path(out_path), abs_path(cif_gz)])
            n += 1

    log.info("Wrote %d AF3 JSONs -> %s", n, args.out_dir)
    log.info("Manifest -> %s", args.manifest_csv)
