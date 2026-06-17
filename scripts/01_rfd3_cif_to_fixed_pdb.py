#!/usr/bin/env python3
"""
01_rfd3_cif_to_fixed_pdb.py

Converts all RFD3 all-atom CIF outputs to backbone PDB files with
FIXED residue remarks for ProteinMPNN, using epitope indices from dp2.

Usage:
    python 01_rfd3_cif_to_fixed_pdb.py \
        --metrics_csv  runs/run_test_rfd3_nompmn/04_filter/metrics_full.csv \
        --dp2_parquet  datasets/dp2.parquet \
        --outdir       runs/run_rfd3_mpnn/01_fixed_pdbs \
        --n_workers    8

Outputs one _fixed.pdb per RFD3 CIF in outdir, ready for dl_interface_design.py.
"""

import argparse
import gzip
import logging
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import gemmi
import pandas as pd

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
log = logging.getLogger(__name__)

BB_ATOMS = ("N", "CA", "C", "O")


def cif_to_fixed_pdb(cif_path: Path, fixed_ris: list, out_pdb: Path) -> str:
    """
    Convert a RFD3 CIF to a backbone PDB with FIXED remarks.
    Returns 'ok' or an error string.
    """
    try:
        with gzip.open(cif_path, "rt") as f:
            doc = gemmi.cif.read_string(f.read())
        st = gemmi.make_structure_from_block(doc.sole_block())
    except Exception as e:
        return f"load_fail: {e}"

    chain = st[0][0]
    residues = list(chain)

    lines = []
    atom_serial = 1
    for res_idx, res in enumerate(residues):
        for atom_name in BB_ATOMS:
            atom = res.find_atom(atom_name, altloc="*")
            if atom is None:
                continue
            p = atom.pos
            lines.append(
                f"ATOM  {atom_serial:5d}  {atom_name:<3s} {res.name:3s} A{res_idx+1:4d}    "
                f"{p.x:8.3f}{p.y:8.3f}{p.z:8.3f}  1.00  0.00           {atom_name[0]:1s}"
            )
            atom_serial += 1

    # FIXED remarks are 1-based for PyRosetta
    for ri in fixed_ris:
        lines.append(f"REMARK PDBinfo-LABEL:{ri+1:5d} FIXED")

    out_pdb.parent.mkdir(parents=True, exist_ok=True)
    out_pdb.write_text("\n".join(lines) + "\n")
    return "ok"


def process_row(args):
    tok, cif_path_str, fixed_ris, out_pdb_str = args
    cif_path = Path(cif_path_str)
    out_pdb = Path(out_pdb_str)

    if out_pdb.exists():
        return tok, "skipped"

    if not cif_path.exists():
        return tok, "missing_cif"

    status = cif_to_fixed_pdb(cif_path, fixed_ris, out_pdb)
    return tok, status


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--metrics_csv",  required=True, help="metrics_full.csv from RFD3 run")
    parser.add_argument("--dp2_parquet",  required=True, help="dp2.parquet (Lawson's dataset)")
    parser.add_argument("--outdir",       required=True, help="Output directory for fixed PDBs")
    parser.add_argument("--n_workers",    type=int, default=4, help="Parallel workers (default: 4)")
    parser.add_argument("--id_filter",    default=None, help="Optional: only process this id (e.g. 7ox3_0P)")
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # --- load metrics ---
    log.info(f"Loading metrics from {args.metrics_csv}")
    df = pd.read_csv(args.metrics_csv)
    df["assay_scaffolded_epitope_id"] = df["assay_scaffolded_epitope_id"].str.lower()
    df = df[df["status"] == "ok"].copy()

    if args.id_filter:
        df = df[df["id"] == args.id_filter]
        log.info(f"Filtered to id={args.id_filter}: {len(df)} rows")

    log.info(f"Total RFD3 designs to convert: {len(df)}")

    # --- load dp2 epitope indices ---
    log.info(f"Loading dp2 from {args.dp2_parquet}")
    dp2 = pd.read_parquet(args.dp2_parquet)
    dp2["assay_scaffolded_epitope_id"] = dp2["assay_scaffolded_epitope_id"].astype(str).str.lower()

    # one row per token is sufficient — indices are same for all rows sharing a token
    dp2_lookup = (
        dp2.drop_duplicates("assay_scaffolded_epitope_id")
        .set_index("assay_scaffolded_epitope_id")["scaffolded_epitope_chunk_resindices"]
    )

    # --- build work items ---
    work = []
    n_no_dp2 = 0
    for _, row in df.iterrows():
        tok = row["assay_scaffolded_epitope_id"]

        if tok not in dp2_lookup.index:
            n_no_dp2 += 1
            continue

        fixed_ris = list(dp2_lookup[tok])
        cif_path = row["rfd3_path"]
        pred = int(row["pred"])
        out_pdb = outdir / f"{tok}_pred{pred}_fixed.pdb"

        work.append((tok, cif_path, fixed_ris, str(out_pdb)))

    if n_no_dp2 > 0:
        log.warning(f"{n_no_dp2} rows had no matching dp2 entry and will be skipped")

    log.info(f"Processing {len(work)} designs with {args.n_workers} workers")

    # --- run in parallel ---
    n_ok = n_skip = n_fail = 0
    with ProcessPoolExecutor(max_workers=args.n_workers) as ex:
        futures = {ex.submit(process_row, w): w[0] for w in work}
        for i, fut in enumerate(as_completed(futures), 1):
            tok, status = fut.result()
            if status == "ok":
                n_ok += 1
            elif status == "skipped":
                n_skip += 1
            else:
                n_fail += 1
                log.warning(f"FAIL {tok}: {status}")

            if i % 500 == 0:
                log.info(f"Progress: {i}/{len(work)}  ok={n_ok}  skip={n_skip}  fail={n_fail}")

    log.info(f"Done: {n_ok} ok, {n_skip} skipped, {n_fail} failed")
    log.info(f"Output PDBs: {outdir}")


if __name__ == "__main__":
    main()
