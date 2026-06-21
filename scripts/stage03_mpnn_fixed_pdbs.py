#!/usr/bin/env python3
"""
stage03_mpnn_fixed_pdbs.py -- RFD3 CIF outputs (02_rfd3) -> backbone PDBs with FIXED epitope
remarks (03_mpnn/fixed_pdbs), sourcing the fixed residues from OUR ledger. Part 1 of stage 03.

This is the dual-island analog of the generic strip-to-fixed-PDB step. The naive version reads
the epitope's FIXED positions from dp2's `scaffolded_epitope_chunk_resindices` (the original
dual-island scaffolded coordinates), which are wrong here: each island is re-scaffolded into a
NEW 103-mer contig `N/Aa-b/C`, so in the RFD3 output the island lands at residues N+1 .. N+span
(1-based). We compute exactly that from `results/dual_island_designs.csv` (the n_flank and
island_size of each contig), keyed by design_id, so ProteinMPNN holds the epitope fixed and
redesigns only the scaffold.

Pipeline position (stage 03_mpnn, part 1): 02_rfd3 (done) -> [THIS] strip-to-backbone+FIXED ->
stage03_mpnn_submit (ProteinMPNN) -> stage04_af3_emit_jsons -> stage04_af3_array. Output: one
`<model_stem>_fixed.pdb` per RFD3 model CIF, which stage03_mpnn_submit then globs.

Usage:
  python3 scripts/stage03_mpnn_fixed_pdbs.py \
      --rfd3_outputs_dir runs/dual_island_rfd3/02_rfd3/outputs \
      --ledger           results/dual_island_designs.csv \
      --outdir           runs/dual_island_rfd3/03_mpnn/fixed_pdbs \
      --n_workers        8
"""
import argparse
import gzip
import logging
import re
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import pandas as pd

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
log = logging.getLogger(__name__)

BB_ATOMS = ("N", "CA", "C", "O")
# design_id looks like  2h32_0P__contig3__seed0__rep0  (id before __contig, contig number after)
DESIGN_ID = re.compile(r"^(?P<id>.+?)__contig(?P<cid>\d+)__seed")


def parse_design_id(stem: str):
    """(id, contig_id) from an RFD3 output cif stem; None if it doesn't match."""
    m = DESIGN_ID.match(stem)
    if not m:
        return None
    return m.group("id"), int(m.group("cid"))


def load_fixed_lookup(ledger_csv: Path) -> dict:
    """(id, contig_id) -> 0-based fixed residue indices in the scaffolded 103-mer.

    The island occupies output residues [n_flank .. n_flank+island_size-1] (0-based); we fix
    the whole island motif (matching Lawson's chunk-level FIXED), so MPNN designs only scaffold.
    """
    df = pd.read_csv(ledger_csv)
    lookup = {}
    for _, r in df.iterrows():
        n_flank, span = int(r["n_flank"]), int(r["island_size"])
        lookup[(str(r["id"]), int(r["contig_id"]))] = list(range(n_flank, n_flank + span))
    return lookup


def cif_to_fixed_pdb(cif_path: Path, fixed_ris: list, out_pdb: Path) -> str:
    """RFD3 CIF -> backbone PDB with PDBinfo-LABEL FIXED remarks (1-based). 'ok' or error."""
    import gemmi  # cluster-only; imported here so the pure logic above is testable without it
    try:
        with gzip.open(cif_path, "rt") as f:
            doc = gemmi.cif.read_string(f.read())
        st = gemmi.make_structure_from_block(doc.sole_block())
    except Exception as e:  # noqa: BLE001
        return f"load_fail: {e}"

    chain = st[0][0]
    lines, serial = [], 1
    for res_idx, res in enumerate(chain):
        for atom_name in BB_ATOMS:
            atom = res.find_atom(atom_name, altloc="*")
            if atom is None:
                continue
            p = atom.pos
            lines.append(
                f"ATOM  {serial:5d}  {atom_name:<3s} {res.name:3s} A{res_idx+1:4d}    "
                f"{p.x:8.3f}{p.y:8.3f}{p.z:8.3f}  1.00  0.00           {atom_name[0]:1s}"
            )
            serial += 1
    for ri in fixed_ris:                       # FIXED remarks are 1-based for PyRosetta
        lines.append(f"REMARK PDBinfo-LABEL:{ri+1:5d} FIXED")

    out_pdb.parent.mkdir(parents=True, exist_ok=True)
    out_pdb.write_text("\n".join(lines) + "\n")
    return "ok"


def _process(args):
    cif_str, fixed_ris, out_str = args
    out_pdb = Path(out_str)
    if out_pdb.exists():
        return out_pdb.name, "skipped"
    return out_pdb.name, cif_to_fixed_pdb(Path(cif_str), fixed_ris, out_pdb)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--rfd3_outputs_dir", required=True)
    ap.add_argument("--ledger", default="results/dual_island_designs.csv")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--n_workers", type=int, default=8)
    args = ap.parse_args()

    lookup = load_fixed_lookup(Path(args.ledger))
    log.info(f"Loaded fixed-residue lookup for {len(lookup)} contigs from {args.ledger}")

    cifs = sorted(Path(args.rfd3_outputs_dir).glob("*.cif.gz"))
    log.info(f"Found {len(cifs)} RFD3 model CIFs in {args.rfd3_outputs_dir}")

    work, n_nokey = [], 0
    for cif in cifs:
        stem = cif.name[:-len(".cif.gz")]
        key = parse_design_id(stem)
        if key is None or key not in lookup:
            n_nokey += 1
            continue
        out_pdb = Path(args.outdir) / f"{stem}_fixed.pdb"
        work.append((str(cif), lookup[key], str(out_pdb)))
    if n_nokey:
        log.warning(f"{n_nokey} CIFs had no parseable/known design_id and were skipped")
    log.info(f"Converting {len(work)} CIFs with {args.n_workers} workers")

    n_ok = n_skip = n_fail = 0
    with ProcessPoolExecutor(max_workers=args.n_workers) as ex:
        futs = [ex.submit(_process, w) for w in work]
        for i, fut in enumerate(as_completed(futs), 1):
            name, status = fut.result()
            if status == "ok":
                n_ok += 1
            elif status == "skipped":
                n_skip += 1
            else:
                n_fail += 1
                log.warning(f"FAIL {name}: {status}")
            if i % 5000 == 0:
                log.info(f"  {i}/{len(work)}  ok={n_ok} skip={n_skip} fail={n_fail}")

    log.info(f"Done: {n_ok} ok, {n_skip} skipped, {n_fail} failed -> {args.outdir}")


if __name__ == "__main__":
    main()
