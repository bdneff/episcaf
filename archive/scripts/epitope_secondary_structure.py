#!/usr/bin/env python3
"""
epitope_secondary_structure.py

Characterize the secondary structure composition of epitope chunk residues
for each unique epitope target, then compare between passing and failing designs.

Tests the hypothesis that passing designs are enriched for helical epitopes,
and that the pipeline is biased against non-helical epitope structures.

Requires DSSP to be installed:
    conda install -c conda-forge dssp
    or: apt-get install dssp

Usage:
    python scripts/epitope_secondary_structure.py \
        --metrics_csv  runs/run_rfd3_mpnn/04_filter/metrics_partial.csv \
        --dp2_parquet  datasets/dp2.parquet \
        --true_dir     /tgen_labs/altin/alphafold3/workspace/episcaf-experiments/data/abdb/complex_pdbfiles/cleaned \
        --out_csv      runs/run_rfd3_mpnn/epitope_ss_analysis.csv
"""

import argparse
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

# try MDAnalysis DSSP first, fall back to Biopython
try:
    import MDAnalysis as mda
    from MDAnalysis.analysis.dssp import DSSP as MDA_DSSP
    USE_MDA = True
except ImportError:
    USE_MDA = False

try:
    from Bio.PDB import PDBParser, DSSP as Bio_DSSP
    USE_BIO = True
except ImportError:
    USE_BIO = False

if not USE_MDA and not USE_BIO:
    raise ImportError("Either MDAnalysis (with DSSP) or Biopython must be installed")


# DSSP codes grouped into broad categories
HELIX  = {'H', 'G', 'I'}   # alpha, 3-10, pi helix
STRAND = {'E', 'B'}         # beta strand, isolated bridge
LOOP   = {'T', 'S', 'C', '-', ' '}  # turn, bend, coil


def ss_to_category(code: str) -> str:
    if code in HELIX:
        return 'helix'
    if code in STRAND:
        return 'strand'
    return 'loop'


def get_epitope_ss_mda(pdb_path: Path, epitope_ris: list) -> dict:
    """Get secondary structure of epitope residues using MDAnalysis DSSP."""
    u = mda.Universe(str(pdb_path))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        dssp = MDA_DSSP(u).run()

    # dssp.results.dssp is shape (n_frames, n_residues)
    ss_per_res = dssp.results.dssp[0]  # first (only) frame

    codes = []
    for ri in epitope_ris:
        if ri < len(ss_per_res):
            codes.append(str(ss_per_res[ri]))

    return summarize_ss(codes)


def get_epitope_ss_bio(pdb_path: Path, epitope_ris: list, pdb_id: str) -> dict:
    """Get secondary structure of epitope residues using Biopython DSSP."""
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure(pdb_id, str(pdb_path))
    model = structure[0]

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        dssp = Bio_DSSP(model, str(pdb_path))

    # get all residues in chain A in order
    chainA_residues = list(model['A'].get_residues())

    codes = []
    for ri in epitope_ris:
        if ri < len(chainA_residues):
            res = chainA_residues[ri]
            key = (res.get_full_id()[2], res.get_id())
            try:
                ss = dssp[key][2]
                codes.append(ss)
            except KeyError:
                codes.append('-')

    return summarize_ss(codes)


def summarize_ss(codes: list) -> dict:
    """Summarize a list of DSSP codes into fraction helix/strand/loop."""
    if not codes:
        return {'n_residues': 0, 'n_helix': 0, 'n_strand': 0, 'n_loop': 0,
                'frac_helix': np.nan, 'frac_strand': np.nan, 'frac_loop': np.nan,
                'dominant_ss': 'unknown', 'ss_codes': ''}

    n = len(codes)
    n_helix  = sum(1 for c in codes if c in HELIX)
    n_strand = sum(1 for c in codes if c in STRAND)
    n_loop   = n - n_helix - n_strand

    frac_helix  = n_helix  / n
    frac_strand = n_strand / n
    frac_loop   = n_loop   / n

    dominant = max(['helix', 'strand', 'loop'],
                   key=lambda x: {'helix': frac_helix, 'strand': frac_strand, 'loop': frac_loop}[x])

    return {
        'n_residues':  n,
        'n_helix':     n_helix,
        'n_strand':    n_strand,
        'n_loop':      n_loop,
        'frac_helix':  frac_helix,
        'frac_strand': frac_strand,
        'frac_loop':   frac_loop,
        'dominant_ss': dominant,
        'ss_codes':    ''.join(codes),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--metrics_csv", required=True)
    parser.add_argument("--dp2_parquet", required=True)
    parser.add_argument("--true_dir",    required=True)
    parser.add_argument("--out_csv",     required=True)
    args = parser.parse_args()

    true_dir = Path(args.true_dir)
    out_csv  = Path(args.out_csv)

    # --- Load data ---
    print("Loading metrics CSV ...")
    df = pd.read_csv(args.metrics_csv, low_memory=False)

    passing_mask = (
        (df['mean_pae'] < 5.0) &
        (df['overall_rmsd'] <= 2.0) &
        (df['epitope_chunk_rmsd'] <= 1.0) &
        (df['af3_n_clash_res'].fillna(1) == 0)
    )
    passing = df[passing_mask]
    print(f"  Total designs: {len(df):,}  Passing: {len(passing):,}")

    print("Loading dp2 ...")
    dp2 = pd.read_parquet(args.dp2_parquet)
    dp2['assay_scaffolded_epitope_id'] = dp2['assay_scaffolded_epitope_id'].astype(str).str.lower()
    dp2_tok = dp2.drop_duplicates('assay_scaffolded_epitope_id').set_index('assay_scaffolded_epitope_id')

    # --- Per-epitope pass counts ---
    pass_counts  = passing.groupby('id').size().rename('n_passing')
    total_counts = df.groupby('id').size().rename('n_total')
    epi_stats = pd.concat([total_counts, pass_counts], axis=1).fillna(0)
    epi_stats['pass_rate'] = epi_stats['n_passing'] / epi_stats['n_total']

    # --- Run DSSP for each unique epitope ---
    print(f"\nRunning DSSP on {len(epi_stats)} epitope targets ...")
    print(f"Using {'MDAnalysis' if USE_MDA else 'Biopython'} for DSSP")

    ss_rows = []
    for pid, row in epi_stats.iterrows():
        pdb_path = true_dir / f"{pid}.pdb"
        if not pdb_path.exists():
            print(f"  SKIP {pid}: PDB not found")
            continue

        # get epitope residue indices from dp2 — look up by id
        dp2_matches = dp2[dp2['id'] == pid]
        if len(dp2_matches) == 0:
            print(f"  SKIP {pid}: not in dp2")
            continue

        epi_ris = [int(i) for i in dp2_matches.iloc[0]['epitope_chunk_resindices']]

        try:
            if USE_MDA:
                ss = get_epitope_ss_mda(pdb_path, epi_ris)
            else:
                ss = get_epitope_ss_bio(pdb_path, epi_ris, pid)

            ss_row = {
                'id':         pid,
                'n_total':    int(row['n_total']),
                'n_passing':  int(row['n_passing']),
                'pass_rate':  row['pass_rate'],
                **ss
            }
            ss_rows.append(ss_row)
            print(f"  {pid:15s}  passes={int(row['n_passing']):4d}  "
                  f"helix={ss['frac_helix']:.2f}  strand={ss['frac_strand']:.2f}  "
                  f"loop={ss['frac_loop']:.2f}  dominant={ss['dominant_ss']}")

        except Exception as e:
            print(f"  FAIL {pid}: {e}")

    out_df = pd.DataFrame(ss_rows)
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    out_df.to_csv(out_csv, index=False)
    print(f"\nWrote: {out_csv}")

    # --- Summary ---
    print(f"\n=== SECONDARY STRUCTURE VS PASS RATE ===")
    print(out_df[['id', 'n_passing', 'pass_rate', 'frac_helix',
                  'frac_strand', 'frac_loop', 'dominant_ss']].sort_values(
        'pass_rate', ascending=False).to_string(index=False))

    print(f"\n=== MEAN PASS RATE BY DOMINANT SS ===")
    print(out_df.groupby('dominant_ss').agg(
        n_epitopes=('id', 'count'),
        mean_pass_rate=('pass_rate', 'mean'),
        total_passing=('n_passing', 'sum')
    ).to_string())

    # correlation
    corr_helix  = out_df['frac_helix'].corr(out_df['pass_rate'])
    corr_strand = out_df['frac_strand'].corr(out_df['pass_rate'])
    corr_loop   = out_df['frac_loop'].corr(out_df['pass_rate'])
    print(f"\nCorrelation with pass rate:")
    print(f"  frac_helix  : {corr_helix:.3f}")
    print(f"  frac_strand : {corr_strand:.3f}")
    print(f"  frac_loop   : {corr_loop:.3f}")


if __name__ == "__main__":
    main()
