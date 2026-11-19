"""DiscoTope-3.0 per-residue epitope-likelihood provider for the live tool backend.

Scoring an uploaded antigen is a CAPABILITY, not a guarantee. The local/static tool ships with no
scorer configured, so `score_cif` returns None and the viewer shows the full growPatch candidate set
unranked ("DiscoTope scores unavailable"). A hosted backend can turn scoring on without any code
change, via environment variables:

  EPISCAF_DISCOTOPE=biolib            Run DiscoTope-3.0 on BioLib's cloud for each uploaded antigen.
                                      Needs `pip install pybiolib` (imports as `biolib`); the ESM-IF1
                                      compute is remote, no local GPU. Off by default.
  EPISCAF_DISCOTOPE_SCORED_PDB=PATH   Read per-residue scores from a pre-baked DiscoTope PDB's CA
                                      B-factor column (0-100, divided by EPISCAF_DISCOTOPE_SCALE,
                                      default 100). For a host that bakes scores ahead of time, or a
                                      local fixture when exercising the ranked path.

Scores are returned as {pdb_residue_number: likelihood in 0..1}, min-max normalized within the chain,
so they rank regions within one antigen. They are NEVER derived from an uploaded PDB's own B-factor
column -- only from a DiscoTope run or an explicitly pre-baked DiscoTope score file.

Mirrors scripts/webtool_run_discotope.py (the offline baker) for the BioLib call and CSV columns.
"""
from __future__ import annotations

import os
import re


def available() -> bool:
    """True when a scorer is configured (so callers can message the capability without running it)."""
    return bool(os.environ.get("EPISCAF_DISCOTOPE_SCORED_PDB")
                or os.environ.get("EPISCAF_DISCOTOPE"))


def _minmax(raw: dict) -> dict | None:
    """Min-max normalize raw per-residue scores to 0..1 within the chain (as the bake does)."""
    if not raw:
        return None
    vals = list(raw.values())
    lo, hi = min(vals), max(vals)
    rng = hi - lo
    if rng <= 0:
        return {k: 0.0 for k in raw}
    return {k: (v - lo) / rng for k, v in raw.items()}


def _from_scored_pdb(path: str, chain: str) -> dict | None:
    """Read normalized scores from a pre-baked DiscoTope PDB's CA B-factor column (0-100 / scale)."""
    try:
        from Bio.PDB import PDBParser
    except Exception:
        return None
    scale = float(os.environ.get("EPISCAF_DISCOTOPE_SCALE", "100"))
    try:
        model = next(iter(PDBParser(QUIET=True).get_structure("scored", path)))
    except Exception:
        return None
    if chain not in model:
        return None
    out = {res.id[1]: res["CA"].get_bfactor() / scale for res in model[chain] if res.has_id("CA")}
    return out or None


def _from_biolib(struct_path: str, chain: str) -> dict | None:
    """Run DiscoTope-3.0 on BioLib and return min-max-normalized per-residue scores, or None.

    Any failure (no client, network error, unexpected output) returns None so the live tool simply
    falls back to the unscored "Candidates" state rather than erroring.
    """
    try:
        import csv as _csv
        import glob
        import tempfile
        import biolib  # pybiolib; optional
    except Exception:
        return None
    try:
        app = biolib.load("DTU/DiscoTope-3")
        job = app.cli(args=f"-f {struct_path} --struc_type solved --out_dir output")
        out = tempfile.mkdtemp(prefix="episcaf_disco_")
        job.save_files(out)
        hits = glob.glob(out + "/**/*discotope3.csv", recursive=True)
        if not hits:
            return None
        raw: dict = {}
        with open(hits[0]) as fh:
            for row in _csv.DictReader(fh):
                ch = (row.get("chain") or row.get("chain_id") or "").strip()
                if ch and chain and ch != chain:
                    continue
                rid, sc = row.get("res_id"), row.get("DiscoTope-3.0_score")
                if rid is None or sc in (None, ""):
                    continue
                m = re.match(r"-?\d+", str(rid))
                if m:
                    raw[int(m.group())] = float(sc)
        return _minmax(raw)
    except Exception:
        return None


def score_cif(cif_path: str, chain: str) -> dict | None:
    """Per-residue DiscoTope-3.0 likelihood (0..1, within-chain normalized) for an uploaded antigen.

    Returns None when no scorer is configured or scoring fails -- the viewer then shows the full
    growPatch candidate set unranked.
    """
    scored = os.environ.get("EPISCAF_DISCOTOPE_SCORED_PDB")
    if scored:
        s = _from_scored_pdb(scored, chain)
        if s is not None:
            return s
    if os.environ.get("EPISCAF_DISCOTOPE") == "biolib":
        return _from_biolib(cif_path, chain)
    return None
