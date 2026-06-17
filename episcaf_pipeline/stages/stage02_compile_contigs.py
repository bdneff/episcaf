#!/usr/bin/env python3
"""Stage 02: Compile an expanded contig table from the design ledger.

Inputs:
- designs parquet (source-of-truth) with at least: id, contig_string, contig_length.
Outputs:
- contigs.parquet: one row per (design, seed, rep), with a stable design_id.

Notes:
- This stage deliberately does **not** depend on downstream result columns.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from pathlib import Path
from typing import List

import pandas as pd

from episcaf_pipeline.schema import SAFE_COLS
from episcaf_pipeline.utils import contig_to_rfd3

log = logging.getLogger(__name__)

@dataclass
class CompileContigsArgs:
    in_parquet: Path
    out_parquet: Path
    seeds: List[int]
    reps: int
    max_rows: int = 0


def compile_contigs(args: CompileContigsArgs) -> None:
    df = pd.read_parquet(args.in_parquet)

    required = ["id", "contig_id", "contig_string", "contig_length"]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns in parquet: {missing}")

    # Keep only SAFE columns that exist
    keep = [c for c in SAFE_COLS if c in df.columns]
    df = df[keep].copy()

    if args.max_rows and args.max_rows > 0:
        df = df.head(args.max_rows).copy()

    # Collapse to ONE contig per design id (Lawson parquet contains many contigs per id).
    # Note: contig_id is only locally unique within each id, so do NOT dedupe on contig_id alone.
    before = len(df)
    df = df.sort_values(['id', 'contig_id']).drop_duplicates(subset=['id','contig_id']).copy()
    log.info('[stage02] Collapsed %d -> %d rows by unique (id, contig_id)', before, len(df))

    # Collapse Lawson-style expanded parquet to ONE row per unique contig.
    # Lawson parquet is expanded across rfd/mpnn fanout; we only want unique contigs.
    # Use a safe key: (assay_scaffolded_epitope_id or id) + contig_id to avoid collisions.
    key0 = "assay_scaffolded_epitope_id" if "assay_scaffolded_epitope_id" in df.columns else "id"
    before = len(df)
    df = df.drop_duplicates(subset=[key0, "contig_id"]).copy()
    log.info("[stage02] Collapsed %d -> %d unique contigs (key=%s+contig_id)", before, len(df), key0)

    rows = []
    for _, r in df.iterrows():
        base_id = r.get("assay_scaffolded_epitope_id", None) or r["id"]
        contig_id = int(r["contig_id"])
        for seed in args.seeds:
            for rep in range(args.reps):
                design_id = f"{base_id}__contig{contig_id}__seed{seed}__rep{rep}"
                out = r.to_dict()
                out["design_id"] = design_id
                out["seed"] = int(seed)
                out["rep"] = int(rep)
                out["contig_rfd3"] = contig_to_rfd3(out["contig_string"])
                rows.append(out)

    out_df = pd.DataFrame(rows)
    args.out_parquet.parent.mkdir(parents=True, exist_ok=True)
    out_df.to_parquet(args.out_parquet, index=False)
    log.info("Wrote %d expanded contig rows -> %s", len(out_df), args.out_parquet)
