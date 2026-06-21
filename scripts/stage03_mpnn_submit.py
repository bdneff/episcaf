#!/usr/bin/env python3
"""
stage03_mpnn_submit.py

Splits the FIXED backbone PDBs (03_mpnn/fixed_pdbs) into batches and submits a SLURM array
job to run dl_interface_design.py (ProteinMPNN) on each batch, writing designed PDBs to
03_mpnn/mpnn_pdbs.

Usage (run after stage03_mpnn_fixed_pdbs.py completes):
    python scripts/stage03_mpnn_submit.py \
        --fixed_pdb_dir runs/<run>/03_mpnn/fixed_pdbs \
        --outdir        runs/<run>/03_mpnn/mpnn_pdbs \
        --batch_size    500 \
        --dry_run       # omit to actually submit
"""

import argparse
import os
import subprocess
from pathlib import Path

MPNN_SCRIPT = "/tgen_labs/altin/alphafold3/workspace/dl_binder_design/mpnn_fr/dl_interface_design.py"
MPNN_ENV    = "/tgen_labs/altin/alphafold3/miniconda3/envs/proteinmpnn_binder_design"

SBATCH_TEMPLATE = """\
#!/bin/bash
#SBATCH --job-name=mpnn_batch_{batch_id:04d}
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --gres=gpu:1
#SBATCH --time=2:00:00
#SBATCH --output=mpnn_batch_{batch_id:04d}_%A.out

mkdir -p {outdir}

{python} {mpnn_script} \\
    -pdbdir  {pdbdir} \\
    -outpdbdir {outdir} \\
    -seqs_per_struct 8 \\
    -temperature 0.1 \\
    -relax_cycles 0
"""


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fixed_pdb_dir", required=True)
    parser.add_argument("--outdir",        required=True)
    parser.add_argument("--batch_size",    type=int, default=500)
    parser.add_argument("--dry_run",       action="store_true",
                        help="Print sbatch scripts without submitting")
    args = parser.parse_args()

    fixed_pdb_dir = Path(args.fixed_pdb_dir).resolve()
    outdir        = Path(args.outdir).resolve()
    # batch staging lives under the same run dir as outdir (no hardcoded run name)
    batch_dir     = outdir.parent / "02_mpnn_batches"

    batch_dir.mkdir(parents=True, exist_ok=True)

    # collect all fixed PDBs
    all_pdbs = sorted(fixed_pdb_dir.glob("*_fixed.pdb"))
    print(f"Found {len(all_pdbs)} fixed PDBs")

    if len(all_pdbs) == 0:
        print("ERROR: no fixed PDBs found — has 01_rfd3_cif_to_fixed_pdb.sh finished?")
        return

    # split into batches — each batch gets its own input subdir
    batches = [all_pdbs[i:i+args.batch_size]
               for i in range(0, len(all_pdbs), args.batch_size)]
    print(f"Splitting into {len(batches)} batches of ~{args.batch_size}")

    python = f"{MPNN_ENV}/bin/python"

    submitted = []
    for batch_id, batch_pdbs in enumerate(batches):
        # create a subdirectory with symlinks to this batch's PDBs
        pdbdir  = batch_dir / f"batch_{batch_id:04d}"
        pdbdir.mkdir(exist_ok=True)

        # symlink each PDB into the batch dir
        for pdb in batch_pdbs:
            link = pdbdir / pdb.name
            if not link.exists():
                link.symlink_to(pdb.resolve())

        batch_outdir = outdir / f"batch_{batch_id:04d}"

        script = SBATCH_TEMPLATE.format(
            batch_id    = batch_id,
            pdbdir      = pdbdir,
            outdir      = batch_outdir,
            mpnn_script = MPNN_SCRIPT,
            python      = python,
        )

        script_path = batch_dir / f"batch_{batch_id:04d}.sh"
        script_path.write_text(script)

        if args.dry_run:
            print(f"[dry_run] Would submit: {script_path}")
        else:
            result = subprocess.run(
                ["sbatch", str(script_path)],
                capture_output=True, text=True
            )
            job_id = result.stdout.strip()
            print(f"Submitted batch {batch_id:04d}: {job_id}")
            submitted.append(job_id)

    if not args.dry_run:
        print(f"\nSubmitted {len(submitted)} jobs")
        print("Monitor with: squeue -u $USER")
        print(f"Output will be in: {outdir}/batch_*/")


if __name__ == "__main__":
    main()
