"""Per-residue DiscoTope scores for the selected chain, keyed by author residue label."""
import csv
import logging
import math
import os
from pathlib import Path
import tempfile

from Bio.PDB import PDBIO
from Bio.PDB.Chain import Chain
from Bio.SeqUtils import seq1

from episcaf_analysis.growpatch import read as reader


def available():
    return bool(os.getenv("EPISCAF_DISCOTOPE_SCORED_PDB") or os.getenv("EPISCAF_DISCOTOPE") == "biolib")


def _residues(path, chain):
    model = next(iter(reader._parse_structure(str(path))))
    return [r for r in model[chain] if r.id[0] == " " and r.resname in reader.MAXASA and r.has_id("CA")]


def _label(residue):
    return str(residue.id[1]) + residue.id[2].strip()


def _from_scored_pdb(path, uploaded, chain):
    # A configured score file belongs to one structure. Chain name and numbering
    # alone cannot identify it when another researcher uploads a different PDB.
    def signature(residues):
        return [(_label(r), r.resname, [(a.name, tuple(a.coord)) for a in r]) for r in residues]

    scored = _residues(path, chain)
    if signature(scored) != signature(_residues(uploaded, chain)):
        return None
    scale = float(os.getenv("EPISCAF_DISCOTOPE_SCALE", "100"))
    if not math.isfinite(scale) or scale <= 0:
        return None
    return {_label(r): r["CA"].bfactor / scale for r in scored if math.isfinite(r["CA"].bfactor)} or None


def _from_biolib(struct_path, chain, structure_type):
    import biolib

    # The provider accepts PDB and returns a per-chain CSV. Send one chain with
    # simple IDs so relative and PDB-numbered outputs map to the same residues,
    # including insertions and author IDs from mmCIF that PDB cannot represent.
    selected = Chain("A")
    mapping = {}
    for residue in _residues(struct_path, chain):
        index = len(mapping) + 1
        mapping[index] = (_label(residue), seq1(residue.resname))
        copy = residue.copy()
        copy.id = (" ", index, " ")
        selected.add(copy)
    if not mapping:
        return None
    with tempfile.TemporaryDirectory(prefix="episcaf_disco_") as tmp:
        pdb = Path(tmp) / "antigen.pdb"
        writer = PDBIO()
        writer.set_structure(selected)
        writer.save(str(pdb))
        job = biolib.load("DTU/DiscoTope-3").cli(args=["-f", str(pdb), "--struc_type", structure_type])
        if job.get_exit_code() != 0:
            raise ValueError("DiscoTope exited unsuccessfully")
        job.save_files(tmp)
        result, = Path(tmp).rglob("*discotope3.csv")
        raw = {}
        with result.open(newline="") as fh:
            for row in csv.DictReader(fh):
                label, amino = mapping[int(row["res_id"])]
                if row["residue"].strip().upper() != amino or label in raw:
                    raise ValueError("DiscoTope residue mapping mismatch")
                score = float(row["DiscoTope-3.0_score"])
                if math.isfinite(score):
                    raw[label] = score
    if not raw:
        return None
    lo, hi = min(raw.values()), max(raw.values())
    return {key: (v - lo) / (hi - lo) if hi > lo else 0.0 for key, v in raw.items()}


def score_cif(cif_path, chain, structure_type="solved"):
    if structure_type not in {"solved", "alphafold"}:
        raise ValueError("Select solved or AlphaFold for DiscoTope.")
    scored = os.getenv("EPISCAF_DISCOTOPE_SCORED_PDB")
    if scored:
        try:
            values = _from_scored_pdb(scored, cif_path, chain)
            if values:
                return values
        except (OSError, ValueError, KeyError, StopIteration):
            logging.getLogger(__name__).warning("Configured DiscoTope score file unavailable")
    if os.getenv("EPISCAF_DISCOTOPE") == "biolib":
        try:
            return _from_biolib(cif_path, chain, structure_type)
        except Exception as exc:
            # Never log uploaded coordinates, provider response bodies or credentials.
            logging.getLogger(__name__).warning("DiscoTope ranking unavailable (%s)", type(exc).__name__)
    return None
