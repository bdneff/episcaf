#!/usr/bin/env python3
"""
stage04_af3_emit_jsons.py

Generate AF3 JSON inputs (04_af3) from MPNN-designed PDB files (03_mpnn).

Usage:
    python scripts/stage04_af3_emit_jsons.py \
        --mpnn_pdb_dir runs/<run>/03_mpnn/mpnn_pdbs \
        --out_dir      runs/<run>/04_af3/inputs \
        --seed         1
"""

import argparse
import json
import logging
from pathlib import Path

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
log = logging.getLogger(__name__)

AA3_TO_1 = {
    "ALA":"A","ARG":"R","ASN":"N","ASP":"D","CYS":"C",
    "GLN":"Q","GLU":"E","GLY":"G","HIS":"H","ILE":"I",
    "LEU":"L","LYS":"K","MET":"M","PHE":"F","PRO":"P",
    "SER":"S","THR":"T","TRP":"W","TYR":"Y","VAL":"V",
    "MSE":"M","SEC":"U","PYL":"O","ASX":"B","GLX":"Z","UNK":"X",
}


def extract_sequence_from_pdb(pdb_path: Path) -> str:
    """Extract one-letter sequence from chain A of an all-atom PDB."""
    seen = []
    seen_keys = set()
    with open(pdb_path) as f:
        for line in f:
            if not line.startswith("ATOM"):
                continue
            if line[21] != "A":
                continue
            res_name = line[17:20].strip()
            res_seq  = line[22:26].strip()
            key = res_seq
            if key in seen_keys:
                continue
            seen_keys.add(key)
            seen.append(res_name)
    if not seen:
        raise RuntimeError(f"No chain A ATOM records in {pdb_path}")
    return "".join(AA3_TO_1.get(r.upper(), "X") for r in seen)


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
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mpnn_pdb_dir", required=True)
    parser.add_argument("--out_dir",      required=True)
    parser.add_argument("--seed",         type=int, default=1,
                        help="AF3 model seed (default: 1, matching Lawson)")
    args = parser.parse_args()

    mpnn_pdb_dir = Path(args.mpnn_pdb_dir)
    out_dir      = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    all_pdbs = sorted(mpnn_pdb_dir.rglob("*_fixed_dldesign_*.pdb"))
    log.info(f"Found {len(all_pdbs)} MPNN PDBs")

    n_written = n_skipped = n_failed = 0
    for pdb_path in all_pdbs:
        pred_id  = pdb_path.stem
        out_json = out_dir / f"{pred_id}.json"

        if out_json.exists():
            n_skipped += 1
            continue

        try:
            seq = extract_sequence_from_pdb(pdb_path)
        except Exception as e:
            log.warning(f"FAIL {pdb_path.name}: {e}")
            n_failed += 1
            continue

        out_json.write_text(json.dumps(make_af3_payload(pred_id, seq, args.seed), indent=2))
        n_written += 1

        if n_written % 10000 == 0:
            log.info(f"Progress: {n_written} written, {n_skipped} skipped, {n_failed} failed")

    log.info(f"Done: {n_written} written, {n_skipped} skipped, {n_failed} failed")
    log.info(f"Output: {out_dir}")

    total = n_written + n_skipped
    chunk_size = 100
    n_tasks = (total + chunk_size - 1) // chunk_size
    log.info(f"Total JSONs: {total}")
    log.info(f"To submit: sbatch --array=1-{n_tasks} scripts/stage04_af3_array.sbatch {out_dir.parent.parent}")


if __name__ == "__main__":
    main()
