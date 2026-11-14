#!/usr/bin/env python3
"""Thin CLI for the IPD-loop round driver (episcaf_pipeline/ipd_driver.py).

Drives ONE round transition per call -- the GPU steps (RFd3/MPNN/AF3) run on the cluster between calls.

  # round 0: stage de-novo RFd3 inputs from an antigen + a growPatch designs table
  python scripts/ipd_round.py init --base-dir runs/ipd_8pww --antigen antigen_clean.pdb \
      --designs-table results/growpatch/designs_8pww.csv

  # (or let it enumerate epitopes with growPatch first)
  python scripts/ipd_round.py init --base-dir runs/ipd_8pww --antigen antigen_clean.pdb \
      --cif antigen.cif --chain A --pdb-id 8pww --total-len 103 --variants 4

  # ...run the cluster steps for round 0, score them, then:
  python scripts/ipd_round.py advance --base-dir runs/ipd_8pww --scored-csv round00_scored.csv

Config dials map 1:1 onto IPDConfig (see episcaf_pipeline/ipd.py).
"""
from __future__ import annotations

import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from episcaf_pipeline import ipd_driver          # noqa: E402
from episcaf_pipeline.ipd import IPDConfig        # noqa: E402


def _cfg_from_args(a) -> IPDConfig:
    noise = tuple(int(x) for x in str(a.noise_levels).split(",") if str(x).strip())
    return IPDConfig(
        score_key=a.score_key, better=a.better, group_key=a.group_key, pass_key=a.pass_key,
        topn_winners=a.topn_winners, budget_per_epitope=a.budget, mode=a.mode,
        children_fraction=a.children_fraction, children_step=a.children_step,
        children_cap=a.children_cap, noise_levels=noise, stop_eps=a.stop_eps, max_rounds=a.max_rounds,
    )


def _add_cfg_args(ap: argparse.ArgumentParser) -> None:
    ap.add_argument("--score-key", default="composite")
    ap.add_argument("--better", default="high", choices=["high", "low"])
    ap.add_argument("--group-key", default="id", help="the epitope/target column")
    ap.add_argument("--pass-key", default="is_pass", help="truthy success column ('' to ignore)")
    ap.add_argument("--topn-winners", type=int, default=3)
    ap.add_argument("--budget", type=int, default=64, help="designs per still-open epitope per round")
    ap.add_argument("--mode", default="adaptive", choices=["adaptive", "fixed"])
    ap.add_argument("--children-fraction", type=float, default=0.5, help="fixed mode: flat child share")
    ap.add_argument("--children-step", type=float, default=0.25, help="adaptive: child share per streak")
    ap.add_argument("--children-cap", type=float, default=0.6)
    ap.add_argument("--noise-levels", default="5,10", help="partial-diffusion re-noise steps (RFd3 partial_T)")
    ap.add_argument("--stop-eps", type=float, default=0.0)
    ap.add_argument("--max-rounds", type=int, default=4)


def cmd_init(a) -> None:
    designs_table = a.designs_table
    if not designs_table:
        if not a.cif:
            raise SystemExit("init needs --designs-table or --cif (to enumerate epitopes with growPatch)")
        # build the growPatch designs table first, then init from it
        from episcaf_pipeline import build_growpatch_designs as bgp
        designs_table = os.path.join(a.base_dir, "growpatch_designs.csv")
        os.makedirs(a.base_dir, exist_ok=True)
        argv = ["--cif", a.cif, "--chain", a.chain, "--total-len", str(a.total_len),
                "--variants", str(a.variants), "--seed", str(a.seed), "--out", designs_table]
        if a.pdb_id:
            argv += ["--pdb-id", a.pdb_id]
        if a.top_k:
            argv += ["--top-k", str(a.top_k)]
        bgp.main(argv)

    state = ipd_driver.init_round0(a.base_dir, a.antigen, designs_table, _cfg_from_args(a))
    run0 = ipd_driver._round_dir(a.base_dir, 0)
    print(f"IPD round 0 staged under {run0}")
    print(f"  epitopes      : {len(state['epitopes'])}  ({', '.join(sorted(state['epitopes']))})")
    print(f"  de-novo inputs: {state['rounds'][0]['n_denovo']} RFd3 JSONs -> {run0}/02_rfd3/inputs/")
    print(f"  state         : {os.path.join(a.base_dir, ipd_driver.STATE_FILE)}")
    print(f"  summary       : {run0}/ROUND_SUMMARY.md")
    print("NEXT (cluster): run RFd3 on the inputs, AOMclash pre-filter, MPNN+AF3, score; then `advance`.")


def cmd_advance(a) -> None:
    res = ipd_driver.advance(a.base_dir, a.scored_csv, parent_outputs_dir=a.parent_outputs_dir)
    run = res["run_dir"]
    cont = [g for g, s in res["summary"].items() if s["status"] == "continue"]
    stop = [g for g, s in res["summary"].items() if s["status"] == "stop"]
    print(f"IPD advanced to round {res['round_idx']:02d} -> {run}")
    print(f"  continue: {len(cont)}  {sorted(cont)}")
    print(f"  stop    : {len(stop)}  {sorted(stop)}")
    print(f"  summary : {run}/ROUND_SUMMARY.md")
    if res["complete"]:
        print("ALL epitopes stopped -- IPD loop complete.")
    else:
        print("NEXT (cluster): RFd3 on the new inputs, AOMclash pre-filter, MPNN+AF3, score; then `advance`.")


def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(prog="ipd_round", description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)

    ai = sub.add_parser("init", help="round 0: stage de-novo RFd3 inputs + write state")
    ai.add_argument("--base-dir", required=True, help="IPD run root (holds round00/, round01/, state)")
    ai.add_argument("--antigen", required=True, help="cleaned antigen PDB (RFd3 de-novo input + epitope frame)")
    ai.add_argument("--designs-table", default="", help="prebuilt growPatch designs table (CSV)")
    ai.add_argument("--cif", default="", help="antigen mmCIF to enumerate epitopes with growPatch")
    ai.add_argument("--chain", default="A")
    ai.add_argument("--pdb-id", default="")
    ai.add_argument("--total-len", type=int, default=103)
    ai.add_argument("--variants", type=int, default=4)
    ai.add_argument("--top-k", type=int, default=0)
    ai.add_argument("--seed", type=int, default=0)
    _add_cfg_args(ai)
    ai.set_defaults(func=cmd_init)

    aa = sub.add_parser("advance", help="one round transition: scored CSV -> next round's RFd3 inputs")
    aa.add_argument("--base-dir", required=True)
    aa.add_argument("--scored-csv", required=True, help="finished round's scored metrics CSV")
    aa.add_argument("--parent-outputs-dir", default="",
                    help="where the finished round's RFd3 backbones live (default: <base>/round<prev>/02_rfd3/outputs)")
    aa.set_defaults(func=cmd_advance)
    return ap


def main(argv=None) -> None:
    a = build_parser().parse_args(argv)
    if getattr(a, "parent_outputs_dir", "") == "":
        a.parent_outputs_dir = None
    a.func(a)


if __name__ == "__main__":
    main()
