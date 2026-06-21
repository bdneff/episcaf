#!/usr/bin/env python3
"""episcaf_pipeline command-line interface.

Design principles:
- One canonical dataset parquet (design ledger).
- Each run is a *snapshot* + generated artifacts in a standard layout.
- Stages are deterministic and parameterized by (dataset, run_dir).

Stage numbers match the run-dir layout they produce:
  00_input -> 01_design -> 02_rfd3 -> 03_mpnn -> 04_af3 -> 05_analysis
(03_mpnn is driven by scripts/stage03_mpnn_*.py using cluster ProteinMPNN tooling.)

Typical usage:

  # init: copy dataset into run/00_input/
  python -m episcaf_pipeline init --dataset datasets/dp2.parquet --run_dir runs/run_20260220_120000

  # stage01: compile expanded contigs into 01_design (adds seeds/reps + contig_rfd3)
  python -m episcaf_pipeline stage01 --run_dir runs/run_... --seeds 0,1,2,3 --reps 1

  # stage02: emit RFD3 JSONs into 02_rfd3
  python -m episcaf_pipeline stage02 --run_dir runs/run_... --input_pdb /abs/path/to/antigen_clean.pdb

  # run RFD3, then 03_mpnn (scripts/), then:
  # stage04: emit AF3 JSONs into 04_af3
  python -m episcaf_pipeline stage04 --run_dir runs/run_...

  # stage05: analyze RMSD into 05_analysis (requires gemmi + MDAnalysis)
  python -m episcaf_pipeline stage05 --run_dir runs/run_...

  # or init + stage01 + stage02 in one shot:
  python -m episcaf_pipeline prep --dataset datasets/dp2.parquet --run_dir runs/run_... --input_pdb ... --seeds 0,1,2,3 --reps 1

"""

from __future__ import annotations

import argparse
import shutil
from pathlib import Path

from episcaf_pipeline.utils import setup_logging
from episcaf_pipeline.paths import ensure_run_layout, RunPaths
from episcaf_pipeline.stages.stage01_compile_contigs import CompileContigsArgs, compile_contigs
from episcaf_pipeline.stages.stage02_emit_rfd3_inputs import EmitRFD3Args, emit_rfd3_inputs
from episcaf_pipeline.stages.stage04_emit_af3_jsons import EmitAF3Args, emit_af3_jsons
from episcaf_pipeline.stages.stage05_rmsd_vs_af3 import RMSDArgs, run_rmsd_vs_af3


def _parse_int_list(s: str) -> list[int]:
    return [int(x.strip()) for x in s.split(",") if x.strip() != ""]


def cmd_init(args: argparse.Namespace) -> None:
    rp = ensure_run_layout(Path(args.run_dir))
    dataset = Path(args.dataset).expanduser().resolve()
    if not dataset.exists():
        raise FileNotFoundError(dataset)
    out = rp.input_dir / "designs.parquet"
    if out.exists() and not args.force:
        raise SystemExit(f"{out} already exists (use --force to overwrite)")
    shutil.copy2(dataset, out)
    print(f"Initialized run: {rp.run_dir}")
    print(f"  dataset snapshot: {out}")


def cmd_stage01(args: argparse.Namespace) -> None:
    rp = ensure_run_layout(Path(args.run_dir))
    in_parquet = rp.input_dir / "designs.parquet"
    if getattr(args, "in_parquet", ""):
        in_parquet = Path(args.in_parquet).expanduser().resolve()
    out_parquet = rp.contigs_parquet
    if getattr(args, "out_parquet", ""):
        out_parquet = Path(args.out_parquet).expanduser().resolve()

    compile_contigs(CompileContigsArgs(
        in_parquet=in_parquet,
        out_parquet=out_parquet,
        seeds=_parse_int_list(args.seeds),
        reps=int(args.reps),
        max_rows=int(args.max_rows),
    ))
    print(f"Wrote contigs: {out_parquet}")


def cmd_stage02(args: argparse.Namespace) -> None:
    rp = ensure_run_layout(Path(args.run_dir))
    contigs = rp.contigs_parquet
    if getattr(args, "contigs_parquet", ""):
        contigs = Path(args.contigs_parquet).expanduser().resolve()

    emit_rfd3_inputs(EmitRFD3Args(
        contigs_parquet=contigs,
        input_pdb=(Path(args.input_pdb).expanduser().resolve() if getattr(args, "input_pdb", "") else None),
        pdb_dir=(Path(args.pdb_dir).expanduser().resolve() if getattr(args, "pdb_dir", "") else None),
        out_dir=rp.rfd3_inputs_dir,
        manifest_csv=rp.rfd3_manifest_csv,
        dump_trajectory=args.dump_trajectory,
        prevalidate_inputs=not args.no_prevalidate,
    ))
    print(f"RFD3 inputs: {rp.rfd3_inputs_dir}")
    print(f"Manifest   : {rp.rfd3_manifest_csv}")


def cmd_stage04(args: argparse.Namespace) -> None:
    rp = ensure_run_layout(Path(args.run_dir))
    emit_af3_jsons(EmitAF3Args(
        rfd3_outputs_dir=rp.rfd3_outputs_dir,
        out_dir=rp.af3_inputs_dir,
        manifest_csv=rp.af3_dir / "inputs_manifest.csv",
        seed=int(args.seed),
        limit=int(args.limit),
    ))
    print(f"AF3 inputs: {rp.af3_inputs_dir}")
    print(f"Manifest : {rp.af3_dir / 'inputs_manifest.csv'}")


def cmd_stage05(args: argparse.Namespace) -> None:
    rp = ensure_run_layout(Path(args.run_dir))
    repo_root = Path(args.repo_root).expanduser().resolve() if args.repo_root else Path(__file__).resolve().parents[1]
    legacy_script = repo_root / "legacy_steps" / "05_rmsd_vs_af3.py"
    if not legacy_script.exists():
        raise SystemExit(f"Legacy RMSD script not found at: {legacy_script}")

    run_rmsd_vs_af3(RMSDArgs(
        legacy_script=legacy_script,
        rfd3_outputs_root=rp.rfd3_outputs_dir,
        af3_outputs_root=rp.af3_outputs_dir,
        out_all=rp.rmsd_all_csv,
        out_best=rp.rmsd_best_csv,
        tmp_pdb_dir=repo_root / "tmp_pdb_for_rmsd",
        verbose=bool(args.verbose),
    ))
    print(f"Wrote: {rp.rmsd_all_csv}")
    print(f"Wrote: {rp.rmsd_best_csv}")


def cmd_prep(args: argparse.Namespace) -> None:
    # init + stage01 (contigs) + stage02 (rfd3 inputs)
    cmd_init(args)
    cmd_stage01(args)
    cmd_stage02(args)


def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(prog="episcaf_pipeline")
    ap.add_argument("-v", "--verbose", action="count", default=0, help="Increase verbosity (-v, -vv)")
    sub = ap.add_subparsers(dest="cmd", required=True)

    ap_init = sub.add_parser("init", help="Create a new run directory and snapshot the dataset parquet into 00_input/")
    ap_init.add_argument("--dataset", required=True, help="Path to canonical dataset parquet")
    ap_init.add_argument("--run_dir", required=True, help="Where to create the run directory")
    ap_init.add_argument("--force", action="store_true", help="Overwrite existing 00_input/designs.parquet")
    ap_init.set_defaults(func=cmd_init)

    ap_s1 = sub.add_parser("stage01", help="Compile expanded contigs table into 01_design (adds design_id, seeds, reps)")
    ap_s1.add_argument("--run_dir", required=True)
    ap_s1.add_argument("--in_parquet", default="", help="Optional override of input parquet (default: run/00_input/designs.parquet)")
    ap_s1.add_argument("--out_parquet", default="", help="Optional override of output parquet")
    ap_s1.add_argument("--seeds", default="0")
    ap_s1.add_argument("--reps", default=1, type=int)
    ap_s1.add_argument("--max_rows", default=0, type=int)
    ap_s1.set_defaults(func=cmd_stage01)

    ap_s2 = sub.add_parser("stage02", help="Emit RFD3 JSON inputs + manifest into 02_rfd3")
    ap_s2.add_argument("--run_dir", required=True)
    ap_s2.add_argument("--contigs_parquet", default="", help="Optional override of contigs parquet")
    mx = ap_s2.add_mutually_exclusive_group(required=True)
    mx.add_argument("--input_pdb", default="", help="Single cleaned antigen-only PDB used by RFD3 (toy runs)")
    mx.add_argument("--pdb_dir", default="", help="Directory of cleaned antigen-only PDBs; per-row PDB is <id>.pdb")
    ap_s2.add_argument("--dump_trajectory", action="store_true")
    ap_s2.add_argument("--no_prevalidate", action="store_true")
    ap_s2.set_defaults(func=cmd_stage02)

    # stage03 = MPNN lives in scripts/ (stage03_mpnn_*.py; uses cluster ProteinMPNN tooling)

    ap_s4 = sub.add_parser("stage04", help="Emit AF3 JSON inputs from RFD3 outputs into 04_af3 (RFD3-direct; no MPNN)")
    ap_s4.add_argument("--run_dir", required=True)
    ap_s4.add_argument("--seed", type=int, default=42)
    ap_s4.add_argument("--limit", type=int, default=0)
    ap_s4.set_defaults(func=cmd_stage04)

    ap_s5 = sub.add_parser("stage05", help="Analyze RMSD between RFD3 and AF3 outputs")
    ap_s5.add_argument("--run_dir", required=True)
    ap_s5.add_argument("--repo_root", default="", help="Optional override if running from elsewhere")
    ap_s5.add_argument("--verbose", action="store_true")
    ap_s5.set_defaults(func=cmd_stage05)

    ap_prep = sub.add_parser("prep", help="init + stage01 + stage02 (prepare run for RFD3 generation)")
    ap_prep.add_argument("--dataset", required=True)
    ap_prep.add_argument("--run_dir", required=True)
    ap_prep.add_argument("--force", action="store_true")
    ap_prep.add_argument("--seeds", default="0")
    ap_prep.add_argument("--reps", default=1, type=int)
    ap_prep.add_argument("--max_rows", default=0, type=int)
    mxp = ap_prep.add_mutually_exclusive_group(required=True)
    mxp.add_argument("--input_pdb", default="")
    mxp.add_argument("--pdb_dir", default="")
    ap_prep.add_argument("--dump_trajectory", action="store_true")
    ap_prep.add_argument("--no_prevalidate", action="store_true")
    ap_prep.set_defaults(func=cmd_prep)

    return ap


def main() -> None:
    ap = build_parser()
    args = ap.parse_args()
    setup_logging(args.verbose if isinstance(args.verbose, int) else 0)
    args.func(args)
