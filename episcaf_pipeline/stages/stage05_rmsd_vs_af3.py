#!/usr/bin/env python3
"""Stage 05: RMSD comparison between RFD3 and AF3 outputs.

This stage currently wraps the proven legacy implementation (gemmi + MDAnalysis),
but exposes a stable CLI API and standard output locations.

Outputs:
- <run_dir>/04_analysis/rmsd_vs_af3_all.csv
- <run_dir>/04_analysis/rmsd_vs_af3_best_per_run.csv

If you later want to fully inline the implementation, replace the wrapper call
with the function code from legacy_steps/05_rmsd_vs_af3.py.
"""

from __future__ import annotations

import logging
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path

log = logging.getLogger(__name__)

@dataclass
class RMSDArgs:
    legacy_script: Path
    rfd3_outputs_root: Path
    af3_outputs_root: Path
    out_all: Path
    out_best: Path
    tmp_pdb_dir: Path
    verbose: bool = False


def run_rmsd_vs_af3(args: RMSDArgs) -> None:
    args.out_all.parent.mkdir(parents=True, exist_ok=True)
    args.tmp_pdb_dir.mkdir(parents=True, exist_ok=True)

    cmd = [
        sys.executable,
        str(args.legacy_script),
        "--rfd3_outputs_root", str(args.rfd3_outputs_root),
        "--af3_outputs_root", str(args.af3_outputs_root),
        "--out_all", str(args.out_all),
        "--out_best", str(args.out_best),
        "--tmp_pdb_dir", str(args.tmp_pdb_dir),
    ]
    if args.verbose:
        cmd.append("--verbose")

    log.info("Running: %s", " ".join(cmd))
    subprocess.run(cmd, check=True)
