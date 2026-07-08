#!/usr/bin/env python3
"""
04_make_fixed_pdbs.py -- 8VDL RFD3 CIF outputs -> backbone PDBs with the epitope FIXED, ready for
ProteinMPNN. This is the ONLY 8VDL-specific step in the MPNN->AF3 leg: it computes the multi-island
FIXED positions (the prior-art minibinder script hardcoded a single 12-mer at the start, wrong for
both the 20-mer epitope and the two-island hotspots), then reuses episcaf's proven `cif_to_fixed_pdb`
primitive. Downstream is episcaf's existing, just-run tooling:

    python scripts/stage03_mpnn_submit.py   --fixed_pdb_dir <out_dir> --outdir <mpnn_pdbs> --tag 8vdl_<run>
    python scripts/stage04_af3_emit_jsons.py --mpnn_pdb_dir <mpnn_pdbs> --out_dir <run_dir>/04_af3/inputs
    sbatch --array=1-N scripts/stage04_af3_array.sbatch <run_dir>

Fixed positions: an island in the RFD3 output construct occupies the residues that follow the scaffold
placed before it, cumulatively. For scaffold segments [s0, s1, ...] (N-flank, inter-island gaps,
C-flank) and island sizes [m1, m2, ...] in native order, island i starts (0-based) at
s0+m1+...+s_{i-1}+m_{i-1}+s_{i-1}. We fix every motif residue, so MPNN redesigns only the scaffold.

Output PDBs are named `<cif_stem>_fixed.pdb` so `stage03_mpnn_submit` (globs `*_fixed.pdb`) and
`stage04_af3_emit_jsons` (globs `*_fixed_dldesign_*.pdb`) pick them up unchanged.

Usage (from repo root, per run):
  python dp4_8vdl/scripts/04_make_fixed_pdbs.py \
      --contigs_csv      dp4_8vdl/01_contigs/epitope20.csv \
      --rfd3_outputs_dir dp4_8vdl/02_rfd3/epitope20/outputs \
      --out_dir          dp4_8vdl/runs/epitope20/03_mpnn/fixed_pdbs
"""
from __future__ import annotations
import argparse
import sys
from multiprocessing import Pool
from pathlib import Path

import pandas as pd

# Reuse episcaf's proven CIF->FIXED-backbone-PDB primitive (gemmi imported lazily inside it, so this
# import works off-cluster; only the actual conversion needs gemmi).
_REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(_REPO / "scripts"))
from stage03_mpnn_fixed_pdbs import cif_to_fixed_pdb  # noqa: E402


def island_sizes(islands_field: str) -> list[int]:
    """'655-656;666-666' -> [2, 1]."""
    sizes = []
    for seg in str(islands_field).split(";"):
        a, b = (int(x) for x in seg.split("-"))
        sizes.append(b - a + 1)
    return sizes


def fixed_positions(scaffold_segs: list[int], sizes: list[int]) -> list[int]:
    """0-based construct positions of every motif residue (see module docstring)."""
    pos, fixed = 0, []
    for seg, isz in zip(scaffold_segs, sizes):     # last (trailing C-flank) seg is unused
        pos += seg
        fixed.extend(range(pos, pos + isz))
        pos += isz
    return fixed


def build_lookup(contigs_csv: Path) -> dict[str, list[int]]:
    df = pd.read_csv(contigs_csv)
    lookup = {}
    for _, r in df.iterrows():
        design_id = f"{r['target']}_contig{int(r['contig_id']):04d}"
        segs = [int(x) for x in str(r["scaffold_segs"]).split(",")]
        lookup[design_id] = fixed_positions(segs, island_sizes(r["islands"]))
    return lookup


def _match_design_id(cif: Path, design_ids) -> str | None:
    """RFD3 writes outputs/<design_id>/<...>.cif.gz; fall back to stem-prefix match."""
    if cif.parent.name in design_ids:
        return cif.parent.name
    stem = cif.name[: -len(".cif.gz")]
    hits = [d for d in design_ids if stem.startswith(d)]
    return max(hits, key=len) if hits else None


def _process(job):
    cif_str, fixed_ris, out_str = job
    out_pdb = Path(out_str)
    if out_pdb.exists():
        return out_pdb.name, "skipped"
    return out_pdb.name, cif_to_fixed_pdb(Path(cif_str), fixed_ris, out_pdb)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--contigs_csv", required=True, type=Path)
    ap.add_argument("--rfd3_outputs_dir", required=True, type=Path)
    ap.add_argument("--out_dir", required=True, type=Path)
    ap.add_argument("--n_workers", type=int, default=8)
    args = ap.parse_args()

    lookup = build_lookup(args.contigs_csv)
    args.out_dir.mkdir(parents=True, exist_ok=True)

    cifs = sorted(args.rfd3_outputs_dir.rglob("*.cif.gz"))
    jobs, unmatched = [], 0
    for cif in cifs:
        did = _match_design_id(cif, lookup)
        if did is None:
            unmatched += 1
            continue
        out_pdb = args.out_dir / f"{cif.name[:-len('.cif.gz')]}_fixed.pdb"
        jobs.append((str(cif), lookup[did], str(out_pdb)))

    print(f"contigs: {len(lookup)} design_ids  |  CIFs found: {len(cifs)}  "
          f"matched: {len(jobs)}  unmatched: {unmatched}")
    if unmatched:
        print(f"  WARNING: {unmatched} CIFs did not match any design_id in the contigs CSV")

    with Pool(args.n_workers) as pool:
        results = pool.map(_process, jobs)
    ok = sum(1 for _, s in results if s == "ok")
    skip = sum(1 for _, s in results if s == "skipped")
    fail = [(n, s) for n, s in results if s not in ("ok", "skipped")]
    print(f"fixed PDBs: ok={ok} skipped={skip} failed={len(fail)} -> {args.out_dir}")
    for n, s in fail[:10]:
        print(f"  FAIL {n}: {s}")
    print(f"\nnext:  python scripts/stage03_mpnn_submit.py --fixed_pdb_dir {args.out_dir} "
          f"--outdir {args.out_dir.parent}/mpnn_pdbs --tag 8vdl")


if __name__ == "__main__":
    main()
