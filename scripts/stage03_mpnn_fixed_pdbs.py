#!/usr/bin/env python3
"""
stage03_mpnn_fixed_pdbs.py -- RFD3 CIF outputs (02_rfd3) -> backbone PDBs with FIXED epitope
remarks (03_mpnn/fixed_pdbs), sourcing the fixed residues from OUR ledger. Part 1 of stage 03.

The epitope's FIXED positions are NOT dp2's `scaffolded_epitope_chunk_resindices` (the original
coordinates): each epitope is re-scaffolded into a NEW 103-mer contig, so in the RFD3 output the
islands land at positions determined by the contig's flanks/gaps. We compute those by walking the
contig's A-segments (`fixed_from_contig`), keyed by design_id, so ProteinMPNN holds the epitope
fixed and redesigns only the scaffold. This handles BOTH the single-island dual-island ledger
(`results/dual_island_designs.csv`) and the multi-island whole-epitope C1-103 ledger
(`results/whole_epitope_designs.csv`, which has no n_flank/island_size columns).

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


_CONTIG_TOK = re.compile(r"^(A?)(\d+)-(\d+)$")


def fixed_from_contig(contig_string: str) -> list:
    """0-based positions of every epitope (A-segment) residue in the scaffolded output construct.

    Walks the slash-delimited contig (``N-N/Aa-b/gap/Ac-d/C-C``), accumulating the output position;
    each A-segment's whole span is fixed, so ProteinMPNN preserves the epitope and redesigns only
    scaffold. Handles ONE island (dual-island run) or MANY (whole-epitope run) uniformly. For a
    single island this equals the old ``n_flank .. n_flank+island_size`` -- so dual-island is
    unchanged -- while multi-island contigs simply contribute one span per island.
    """
    pos, fixed = 0, []
    for tok in str(contig_string).split("/"):
        m = _CONTIG_TOK.match(tok)
        if not m:
            raise ValueError(f"unparseable contig token {tok!r} in {contig_string!r}")
        a, b = int(m.group(2)), int(m.group(3))
        if m.group(1) == "A":                       # epitope island: fix its whole span
            fixed.extend(range(pos, pos + (b - a + 1)))
            pos += b - a + 1
        else:                                        # scaffold flank/spacer ('N-N' -> a residues)
            pos += a
    return fixed


def load_fixed_lookup(ledger_csv: Path) -> dict:
    """(id, contig_id) -> 0-based fixed epitope positions in the scaffolded 103-mer.

    Parsed from each contig's A-segments (`fixed_from_contig`), so it works for the single-island
    dual-island ledger AND the multi-island whole-epitope ledger (which has no n_flank/island_size).
    """
    df = pd.read_csv(ledger_csv)
    lookup = {}
    for _, r in df.iterrows():
        lookup[(str(r["id"]), int(r["contig_id"]))] = fixed_from_contig(r["contig_string"])
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
