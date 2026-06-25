#!/usr/bin/env python3
"""
verify_antigen_seq_vs_pdb_fasta.py -- does dp2 `antigen_seq` match the PDB FASTA, and does it
splice across internal unresolved gaps?

John's question (2026-06-25): the tiled-30mer antigen sequences should come from the PDB FASTA
files (full deposited sequence, NO gaps) -- unlike the .pdb ATOM records, which drop unresolved
residues. Our current dp4_tiled30mers.csv was tiled from dp2 `antigen_seq`. This script checks,
for every mAb antigen, whether that `antigen_seq` is a CLEAN CONTIGUOUS WINDOW of the RCSB FASTA
chain (resolved-but-honest: only termini trimmed) or whether it SPLICES across an internal gap
(a fake junction that would create a non-native tile).

Method (no cluster; RCSB web only): fetch https://www.rcsb.org/fasta/entry/<PDB>, cache it, and
for each antigen pick the FASTA chain that best contains `antigen_seq`. Classify:
  clean_substring : antigen_seq is an exact substring of the FASTA chain  -> no internal gap,
                    only N/C-terminal trim (report how many residues trimmed each end)
  internal_gap    : antigen_seq aligns to the FASTA chain with >=1 internal deletion (fake join)
  mismatch        : no good alignment (chain-mapping or sequence problem -- inspect by hand)

Outputs results/antigen_seq_vs_fasta.csv and a console summary. Reproducible: re-runs read the
FASTA cache under data/sequences/pdb_fasta_cache/ so no network is needed the second time.

Usage:
  python scripts/verify_antigen_seq_vs_pdb_fasta.py
"""
from __future__ import annotations
import sys, time, difflib, urllib.request, urllib.error
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
from configs.paths import DP2_PARQUET_LOCAL  # noqa: E402

CACHE = ROOT / "data/sequences/pdb_fasta_cache"
OUT = ROOT / "results/antigen_seq_vs_fasta.csv"
FASTA_URL = "https://www.rcsb.org/fasta/entry/{pdb}"


def fetch_fasta(pdb: str) -> str:
    """RCSB FASTA for one entry, cached to disk for deterministic re-runs."""
    CACHE.mkdir(parents=True, exist_ok=True)
    cf = CACHE / f"{pdb.lower()}.fasta"
    if cf.exists():
        return cf.read_text()
    url = FASTA_URL.format(pdb=pdb.upper())
    for attempt in range(3):
        try:
            with urllib.request.urlopen(url, timeout=30) as r:
                txt = r.read().decode()
            cf.write_text(txt)
            time.sleep(0.2)
            return txt
        except urllib.error.URLError:
            if attempt == 2:
                raise
            time.sleep(1.0)
    return ""


def parse_chains(fasta_txt: str):
    """Return list of (header, seq) records."""
    recs, hdr, seq = [], None, []
    for line in fasta_txt.splitlines():
        if line.startswith(">"):
            if hdr is not None:
                recs.append((hdr, "".join(seq)))
            hdr, seq = line[1:], []
        elif line.strip():
            seq.append(line.strip())
    if hdr is not None:
        recs.append((hdr, "".join(seq)))
    return recs


def name_from_header(hdr: str) -> tuple[str, str]:
    """'7OX3_3|Chain C|Interleukin-9|Homo sapiens (9606)' -> ('Interleukin-9','Homo sapiens (9606)')."""
    parts = hdr.split("|")
    prot = parts[2].strip() if len(parts) > 2 else ""
    org = parts[3].strip() if len(parts) > 3 else ""
    return prot, org


def classify(antigen_seq: str, chains):
    """Pick the chain best matching antigen_seq; classify the relationship."""
    # 1) exact substring anywhere -> clean window
    for hdr, seq in chains:
        idx = seq.find(antigen_seq)
        if idx >= 0:
            prot, org = name_from_header(hdr)
            return dict(status="clean_substring", chain=hdr.split("|")[1] if "|" in hdr else "",
                        len_fasta=len(seq), trim_n=idx, trim_c=len(seq) - idx - len(antigen_seq),
                        n_internal_gaps=0, protein=prot, organism=org)
    # 2) best difflib alignment -> internal gap or mismatch
    best = None
    for hdr, seq in chains:
        sm = difflib.SequenceMatcher(None, antigen_seq, seq, autojunk=False)
        ratio = sm.ratio()
        if best is None or ratio > best[0]:
            best = (ratio, hdr, seq, sm)
    ratio, hdr, seq, sm = best
    blocks = [b for b in sm.get_matching_blocks() if b.size > 0]
    # count internal gaps in antigen_seq coverage (jumps in the b-index between consecutive blocks)
    gaps = 0
    for p, q in zip(blocks, blocks[1:]):
        if (q.a - (p.a + p.size)) > 0 or (q.b - (p.b + p.size)) > 0:
            gaps += 1
    prot, org = name_from_header(hdr)
    status = "internal_gap" if ratio >= 0.95 else "mismatch"
    return dict(status=status, chain=hdr.split("|")[1] if "|" in hdr else "",
                len_fasta=len(seq), trim_n="", trim_c="", n_internal_gaps=gaps,
                protein=prot, organism=org, ratio=round(ratio, 3))


def main():
    dp2 = pd.read_parquet(DP2_PARQUET_LOCAL).drop_duplicates("id")
    rows = []
    for _, r in dp2.iterrows():
        pid = str(r["id"])
        pdb = pid[:4]
        aseq = str(r["antigen_seq"])
        try:
            chains = parse_chains(fetch_fasta(pdb))
        except Exception as e:  # noqa: BLE001
            rows.append(dict(id=pid, pdb=pdb, len_antigen_seq=len(aseq), status=f"fetch_error:{e}"))
            continue
        info = classify(aseq, chains)
        rows.append(dict(id=pid, pdb=pdb, len_antigen_seq=len(aseq), **info))
        print(f"{pid:9s} {info['status']:15s} "
              f"aseq={len(aseq):4d} fasta={info.get('len_fasta','?')} "
              f"trimN={info.get('trim_n','')} trimC={info.get('trim_c','')} "
              f"gaps={info.get('n_internal_gaps','')}  {info.get('protein','')}")
    df = pd.DataFrame(rows)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUT, index=False)
    print("\n=== summary ===")
    print(df["status"].value_counts().to_string())
    print(f"\nwrote {OUT}  ({len(df)} antigens)")


if __name__ == "__main__":
    main()
