"""Reusable RCSB PDB-FASTA helpers for the tiled-antigen library.

Antigen sequences for the tiled-30mer controls come from the PDB FASTA (the full deposited
per-chain sequence, gap-free) -- not the crystal-resolved ATOM records, which drop unresolved
residues. These helpers fetch+cache the FASTA, parse chains, read the protein name/organism
from the header, and pick the antigen chain by matching a reference (resolved) sequence.

Deterministic re-runs: fetched FASTAs are cached under data/sequences/pdb_fasta_cache/, so a
second run needs no network.
"""
from __future__ import annotations
import time, difflib, urllib.request, urllib.error
from pathlib import Path

CACHE = Path(__file__).resolve().parent.parent / "data/sequences/pdb_fasta_cache"
FASTA_URL = "https://www.rcsb.org/fasta/entry/{pdb}"


def fetch_fasta(pdb: str) -> str:
    """RCSB FASTA text for one entry, cached to disk."""
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
    """[(header, seq), ...] from FASTA text."""
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
    """'7OX3_3|Chain C|Interleukin-9|Homo sapiens (9606)' -> ('Interleukin-9', 'Homo sapiens (9606)')."""
    parts = hdr.split("|")
    prot = parts[2].strip() if len(parts) > 2 else ""
    org = parts[3].strip() if len(parts) > 3 else ""
    return prot, org


def pick_antigen_chain(ref_seq: str, chains):
    """Pick the FASTA chain that is the antigen, by matching the resolved reference `ref_seq`.

    Returns dict: header, seq (full FASTA chain), start (index where ref_seq begins, or None),
    protein, organism, status ('clean' if ref_seq is an exact substring -> no internal gap;
    'needs_review' if only a fuzzy alignment, e.g. ref_seq spliced an internal gap), ratio.
    """
    for hdr, seq in chains:
        idx = seq.find(ref_seq)
        if idx >= 0:
            prot, org = name_from_header(hdr)
            return dict(header=hdr, seq=seq, start=idx, protein=prot, organism=org,
                        status="clean", ratio=1.0)
    # fuzzy fallback (e.g. ref_seq has an internal splice that isn't in any chain verbatim)
    best = max(chains, key=lambda hs: difflib.SequenceMatcher(None, ref_seq, hs[1],
                                                              autojunk=False).ratio())
    hdr, seq = best
    ratio = difflib.SequenceMatcher(None, ref_seq, seq, autojunk=False).ratio()
    prot, org = name_from_header(hdr)
    return dict(header=hdr, seq=seq, start=None, protein=prot, organism=org,
                status="needs_review", ratio=round(ratio, 3))


ORGANISM_SHORT = {  # extend as needed; falls back to the raw organism string
    "Homo sapiens": "human", "Mus musculus": "mouse", "Rattus norvegicus": "rat",
}


def organism_short(org: str) -> str:
    base = org.split("(")[0].strip()
    return ORGANISM_SHORT.get(base, base.replace(" ", "_"))
