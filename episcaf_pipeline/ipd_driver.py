#!/usr/bin/env python3
"""IPD-loop ROUND DRIVER: orchestrate one round transition at a time.

The IPD ("iterative protein design") loop runs rounds of
    generate (RFd3) -> AOMclash backbone pre-filter -> MPNN -> AF3 -> score -> seed next round.
The pure CONTROLLER that decides *what* the next round should be lives in `episcaf_pipeline/ipd.py`
(`plan_round`); THIS module is the local DRIVER that turns those decisions into concrete RFd3 input
artifacts and keeps the round-state on disk. It does **not** touch the controller's logic -- it drives it.

Why "one round transition at a time": the generate/MPNN/AF3 steps are GPU work that runs on the
cluster *between* driver calls. So the driver has exactly two entry points:

  * `init_round0`  -- stage de-novo RFd3 inputs for round 0 from an antigen + chosen epitopes
                      (or a prebuilt growPatch designs table), and write the round-state file.
  * `advance`      -- given the just-finished round's SCORED metrics CSV, call `ipd.plan_round`,
                      then for each still-open epitope emit the next round's RFd3 inputs
                      (n_denovo fresh de-novo contigs + partial-diffusion children) and rewrite state.

WHAT RUNS WHERE
  LOCAL  (this driver): build contigs, emit RFd3 input JSONs, bookkeeping, scoring-CSV -> plan.
  CLUSTER (user runs, between driver calls): RFd3 (02_rfd3/outputs), AOMclash pre-filter, MPNN, AF3,
          metrics + score -> the scored CSV fed back into `advance`.

The driver reuses the existing stages where they already apply:
  * de-novo inputs  : `build_growpatch_designs` (optional, --cif) -> `stage01_compile_contigs`
                      -> `stage02_emit_rfd3_inputs` (standard motif-scaffolding JSONs vs. the antigen).
  * children inputs : a dedicated partial-diffusion emitter (stage02 cannot do partial diffusion);
                      see `partial_diffusion_rfd3_payload` -- its RFd3 field names are FLAGGED unverified.
"""
from __future__ import annotations

import csv
import json
import logging
import re
from dataclasses import asdict
from pathlib import Path
from typing import Any

from episcaf_pipeline.ipd import IPDConfig, plan_round
from episcaf_pipeline.paths import ensure_run_layout
from episcaf_pipeline.stages.stage01_compile_contigs import CompileContigsArgs, compile_contigs
from episcaf_pipeline.stages.stage02_emit_rfd3_inputs import EmitRFD3Args, emit_rfd3_inputs
from episcaf_pipeline.utils import abs_path

# design_epi_positions: where a design's fixed epitope lands in the scaffold (chain A, 1-based),
# derived from the contig alone -- the same helper the AOMclash pre-filter/scorer use.
from episcaf_analysis.aomclash.score import design_epi_positions

log = logging.getLogger(__name__)

STATE_FILE = "ipd_state.json"
# design_id shape produced by stage01 and consumed by the AOMclash scripts:
#   <id>__contig<cid>__seed<s>__rep<r>[...trailing RFd3/MPNN suffix...]
_DESIGN_ID = re.compile(r"^(?P<id>.+?)__contig(?P<cid>\d+)__seed(?P<seed>\d+)__rep(?P<rep>\d+)")

# Columns of the per-round designs table (what the AOMclash pre-filter / scorer consume):
# id, contig_id, contig_string, contig_length, growpatch_epitope map a design_id back to its
# epitope + contig so the scorer can place the epitope. The extra bookkeeping columns
# (kind/parent/noise) record how each design was spawned.
DESIGNS_COLS = ["id", "contig_id", "contig_string", "contig_length", "growpatch_epitope",
                "scaffolded_epitope_resindices", "epitope_chunks", "kind", "parent", "noise"]


# --------------------------------------------------------------------------------------------------
# PARTIAL-DIFFUSION CHILD INPUT  (the one new RFd3 mapping; isolated + FLAGGED)
# --------------------------------------------------------------------------------------------------
def partial_diffusion_rfd3_payload(design_id: str, parent_pdb: str, epi_pos_1based: list[int],
                                   length: int, noise: int) -> dict[str, Any]:
    """Build ONE RFd3 partial-diffusion input-JSON entry for a child of a winning backbone.

    A child re-noises a winner's backbone part-way (`noise` steps) and denoises again with the
    epitope motif held fixed -- exploring the neighbourhood of something that already worked. The
    shape mirrors stage02's motif-scaffolding JSON (input / contig / length / select_fixed_atoms)
    plus the partial-diffusion controls.

    !!! FLAG -- UNVERIFIED RFd3 SCHEMA !!!
    The partial-diffusion FIELD NAMES below (`partial_T`, `noise_steps`, `input_pdb`,
    `partial_diffusion`) are a *sensible guess*; they MUST be confirmed against the actual RFd3(v3)
    partial-diffusion config before a real run. Everything partial-diffusion-specific is confined to
    THIS function on purpose, so confirming the schema is a one-function edit. The `input`/`contig`/
    `length`/`select_fixed_atoms` keys match stage02 and are not in doubt.

    Args:
        design_id      : the child's design_id (<id>__contig<cid>__seed<s>__rep<r>).
        parent_pdb     : absolute path to the winner's RFd3 backbone PDB (produced on the cluster).
        epi_pos_1based : the epitope's residue ids in the backbone (chain A, 1-based) to hold fixed.
        length         : the backbone length (partial diffusion preserves it).
        noise          : partial-diffusion re-noise step count (IPDConfig.noise_levels / RFd3 partial_T).
    """
    select_fixed_atoms = {f"A{p}": "BKBN" for p in epi_pos_1based}
    return {
        design_id: {
            # --- same meaning as stage02 (confirmed) ---
            "input": parent_pdb,                 # the backbone to re-noise (NOT the native antigen)
            "contig": f"A1-{length}",            # keep the whole parent chain A; partial diff fills gaps
            "length": f"{length}-{length}",
            "select_fixed_atoms": select_fixed_atoms,
            "select_unfixed_sequence": False,
            # --- partial-diffusion controls (FLAGGED: names unverified, confirm vs RFd3 config) ---
            "partial_diffusion": True,
            "partial_T": int(noise),             # assumed RFd3 name for the re-noise step count ...
            "noise_steps": int(noise),           # ... alias in case the config uses this name instead
            "input_pdb": parent_pdb,             # assumed duplicate key some RFd3 configs expect
        }
    }


# --------------------------------------------------------------------------------------------------
# state I/O
# --------------------------------------------------------------------------------------------------
def _round_dir(base_dir: Path, idx: int) -> Path:
    return Path(base_dir) / f"round{idx:02d}"


def load_state(base_dir: str | Path) -> dict:
    path = Path(base_dir) / STATE_FILE
    if not path.exists():
        raise FileNotFoundError(f"no IPD state at {path}; run `init` first")
    return json.loads(path.read_text())


def save_state(base_dir: str | Path, state: dict) -> None:
    path = Path(base_dir) / STATE_FILE
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(state, indent=2))


def _cfg_from_state(state: dict) -> IPDConfig:
    d = dict(state["config"])
    if "noise_levels" in d:
        d["noise_levels"] = tuple(d["noise_levels"])   # JSON list -> tuple
    return IPDConfig(**d)


# --------------------------------------------------------------------------------------------------
# designs-table helpers
# --------------------------------------------------------------------------------------------------
def _epi_resindices_0based(contig_string: str) -> str:
    """0-based epitope positions in the scaffold, as a comma string (stage02 maps i -> A(i+1))."""
    return ",".join(str(p - 1) for p in design_epi_positions(contig_string, resid_offset=1))


def _length_int(contig_length: Any) -> int:
    s = str(contig_length).strip()
    return int(s.split("-")[0]) if "-" in s else int(s)


def _epi_col(row: dict) -> str:
    for c in ("growpatch_epitope", "epitope"):
        if c in row and row[c]:
            return c
    raise ValueError("designs table row has no epitope column (growpatch_epitope|epitope)")


def _load_designs_pool(designs_table: str | Path) -> dict[str, list[dict]]:
    """Group a build_growpatch-style designs table into a per-epitope de-novo contig pool.

    Returns {id: [ {contig_id, contig_string, contig_length, epitope, resind, epitope_chunks}, ... ]}.
    """
    rows = list(csv.DictReader(open(designs_table)))
    if not rows:
        raise SystemExit(f"designs table {designs_table} is empty")
    epi_col = _epi_col(rows[0])
    pool: dict[str, list[dict]] = {}
    for r in rows:
        pool.setdefault(r["id"], []).append(dict(
            contig_id=int(r["contig_id"]),
            contig_string=r["contig_string"],
            contig_length=r["contig_length"],
            epitope=r[epi_col],
            resind=_epi_resindices_0based(r["contig_string"]),
            epitope_chunks=r.get("epitope_chunks", ""),
        ))
    return pool


def _write_designs_table(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=DESIGNS_COLS)
        w.writeheader()
        for r in rows:
            w.writerow({c: r.get(c, "") for c in DESIGNS_COLS})


# --------------------------------------------------------------------------------------------------
# emission
# --------------------------------------------------------------------------------------------------
def _emit_denovo(run_dir: Path, antigen_pdb: str, denovo_rows: list[dict], seed: int) -> list[str]:
    """Emit de-novo RFd3 inputs via the real stage01 -> stage02 path. Returns design_ids emitted.

    Each de-novo row carries a UNIQUE contig_id within its epitope, so stage01 (seeds=[seed], reps=1)
    collapses to exactly one design per row: design_id = <id>__contig<cid>__seed<seed>__rep0.
    `scaffolded_epitope_resindices` (added here from the contig) is what stage02 turns into the BKBN
    select_fixed_atoms, so growPatch contigs get correct fixed atoms without native resindex columns.
    """
    import pandas as pd
    rp = ensure_run_layout(run_dir)
    df = pd.DataFrame([{
        "id": r["id"],
        "contig_id": int(r["contig_id"]),
        "contig_string": r["contig_string"],
        "contig_length": r["contig_length"],
        "scaffolded_epitope_resindices": r["scaffolded_epitope_resindices"],
        "epitope_chunks": r.get("epitope_chunks", ""),
    } for r in denovo_rows])
    in_parquet = rp.input_dir / "designs.parquet"
    df.to_parquet(in_parquet, index=False)

    compile_contigs(CompileContigsArgs(in_parquet=in_parquet, out_parquet=rp.contigs_parquet,
                                       seeds=[seed], reps=1))
    emit_rfd3_inputs(EmitRFD3Args(contigs_parquet=rp.contigs_parquet, input_pdb=Path(antigen_pdb),
                                  pdb_dir=None, out_dir=rp.rfd3_inputs_dir,
                                  manifest_csv=rp.rfd3_manifest_csv))
    with open(rp.rfd3_manifest_csv) as fh:
        return [r["design_id"] for r in csv.DictReader(fh)]


def _emit_children(run_dir: Path, child_rows: list[dict]) -> list[tuple[str, str]]:
    """Emit one partial-diffusion RFd3 input JSON per child row. Returns [(design_id, json_path)].

    Each child row carries: design_id, parent_pdb (abs path to the winner backbone on the cluster),
    contig_string (inherited from the parent -> same epitope positions), length, noise.
    """
    rp = ensure_run_layout(run_dir)
    out = []
    for r in child_rows:
        epi_pos = design_epi_positions(r["contig_string"], resid_offset=1)
        payload = partial_diffusion_rfd3_payload(r["design_id"], r["parent_pdb"], epi_pos,
                                                  _length_int(r["contig_length"]), int(r["noise"]))
        out_json = rp.rfd3_inputs_dir / f"{r['design_id']}.json"
        out_json.write_text(json.dumps(payload, indent=2))
        out.append((r["design_id"], abs_path(out_json)))
    return out


def _append_manifest(run_dir: Path, entries: list[tuple[str, str]]) -> None:
    """Append child (design_id, json_path) rows to the RFd3 manifest, creating it if needed."""
    rp = ensure_run_layout(run_dir)
    mf = rp.rfd3_manifest_csv
    new = not mf.exists()
    with open(mf, "a", newline="") as fh:
        w = csv.writer(fh)
        if new:
            w.writerow(["design_id", "json_path"])
        for did, jp in entries:
            w.writerow([did, jp])


# --------------------------------------------------------------------------------------------------
# round 0: init
# --------------------------------------------------------------------------------------------------
def init_round0(base_dir: str | Path, antigen_pdb: str, designs_table: str | Path,
                cfg: IPDConfig) -> dict:
    """Stage de-novo RFd3 inputs for round 0 from a (growPatch-style) designs table + the antigen,
    and write the round-state file. Returns the state dict."""
    base_dir = Path(base_dir)
    antigen_abs = abs_path(antigen_pdb)
    pool = _load_designs_pool(designs_table)

    # per-epitope controller state (plan_round's prev_state) + driver bookkeeping
    epitopes: dict[str, dict] = {}
    registry: dict[str, dict] = {}       # "<id>__contig<cid>" -> {contig_string, epitope, length}
    denovo_pool: dict[str, list[dict]] = {}
    round0_rows: list[dict] = []
    for g, contigs in pool.items():
        denovo_pool[g] = contigs
        next_cid = max(c["contig_id"] for c in contigs) + 1
        epitopes[g] = {"best": None, "streak": 0, "status": "open", "reason": "", "next_cid": next_cid}
        for c in contigs:
            key = f"{g}__contig{c['contig_id']}"
            registry[key] = {"contig_string": c["contig_string"], "epitope": c["epitope"],
                             "length": _length_int(c["contig_length"])}
            round0_rows.append({
                "id": g, "contig_id": c["contig_id"], "contig_string": c["contig_string"],
                "contig_length": c["contig_length"], "growpatch_epitope": c["epitope"],
                "scaffolded_epitope_resindices": c["resind"], "epitope_chunks": c["epitope_chunks"],
                "kind": "denovo", "parent": "", "noise": "",
            })

    run_dir = _round_dir(base_dir, 0)
    _write_designs_table(run_dir / "designs_table.csv", round0_rows)
    design_ids = _emit_denovo(run_dir, antigen_abs, round0_rows, seed=0)

    state = {
        "antigen_pdb": antigen_abs,
        "config": {**asdict(cfg), "noise_levels": list(cfg.noise_levels)},
        "round_idx": 0,
        "epitopes": epitopes,
        "denovo_pool": denovo_pool,
        "design_registry": registry,
        "rounds": [{"round_idx": 0, "run_dir": str(run_dir),
                    "designs_table": str(run_dir / "designs_table.csv"),
                    "n_denovo": len(design_ids), "n_child": 0,
                    "epitopes_open": sorted(epitopes), "epitopes_stopped": []}],
    }
    save_state(base_dir, state)
    _write_summary(run_dir, 0, cfg, {g: {"status": "init", "n_denovo": len(denovo_pool[g]),
                                         "n_child": 0, "children_fraction": 0.0,
                                         "reason": "round-0 de-novo exploration"} for g in epitopes},
                   total_denovo=len(design_ids), total_child=0)
    return state


# --------------------------------------------------------------------------------------------------
# advance: one round transition
# --------------------------------------------------------------------------------------------------
def _parse_group_cid(design_id: str) -> tuple[str, str] | None:
    m = _DESIGN_ID.match(design_id)
    return (m.group("id"), m.group("cid")) if m else None


def _lookup_parent(state: dict, design_id: str) -> dict | None:
    gc = _parse_group_cid(design_id)
    if not gc:
        return None
    return state["design_registry"].get(f"{gc[0]}__contig{gc[1]}")


def _take_denovo(state: dict, g: str, n: int, seed: int) -> list[dict]:
    """n fresh de-novo design rows for epitope g: cycle its round-0 contig pool with new unique cids."""
    pool = state["denovo_pool"][g]
    rows = []
    for i in range(n):
        base = pool[i % len(pool)]
        cid = state["epitopes"][g]["next_cid"]
        state["epitopes"][g]["next_cid"] = cid + 1
        key = f"{g}__contig{cid}"
        state["design_registry"][key] = {"contig_string": base["contig_string"],
                                          "epitope": base["epitope"],
                                          "length": _length_int(base["contig_length"])}
        rows.append({
            "id": g, "contig_id": cid, "contig_string": base["contig_string"],
            "contig_length": base["contig_length"], "growpatch_epitope": base["epitope"],
            "scaffolded_epitope_resindices": _epi_resindices_0based(base["contig_string"]),
            "epitope_chunks": base.get("epitope_chunks", ""),
            "kind": "denovo", "parent": "", "noise": "",
        })
    return rows


def _make_children(state: dict, g: str, children_spec: list[dict], seed: int,
                   parent_outputs: Path) -> list[dict]:
    """Expand plan_round's per-(winner,noise) child specs into concrete child design rows."""
    rows = []
    for spec in children_spec:
        parent_did, noise, count = spec["parent"], spec["noise"], spec["count"]
        meta = _lookup_parent(state, parent_did)
        if meta is None:
            log.warning("child parent %s not in registry; skipping", parent_did)
            continue
        parent_pdb = abs_path(parent_outputs / f"{parent_did}.pdb")
        for _ in range(count):
            cid = state["epitopes"][g]["next_cid"]
            state["epitopes"][g]["next_cid"] = cid + 1
            did = f"{g}__contig{cid}__seed{seed}__rep0"
            # child inherits the parent's contig_string -> identical epitope positions (partial
            # diffusion preserves the backbone/length), so the scorer maps it the same way.
            state["design_registry"][f"{g}__contig{cid}"] = {
                "contig_string": meta["contig_string"], "epitope": meta["epitope"],
                "length": meta["length"]}
            rows.append({
                "id": g, "contig_id": cid, "contig_string": meta["contig_string"],
                "contig_length": f"{meta['length']}-{meta['length']}",
                "growpatch_epitope": meta["epitope"],
                "scaffolded_epitope_resindices": _epi_resindices_0based(meta["contig_string"]),
                "epitope_chunks": "",
                "kind": "child", "parent": parent_did, "noise": noise,
                "design_id": did, "parent_pdb": parent_pdb,
            })
    return rows


def advance(base_dir: str | Path, scored_csv: str | Path,
            parent_outputs_dir: str | Path | None = None) -> dict:
    """One round transition: read the finished round's scored CSV, plan + emit the next round.

    Returns {"round_idx", "run_dir", "plan", "summary", "complete"}.
    """
    base_dir = Path(base_dir)
    state = load_state(base_dir)
    cfg = _cfg_from_state(state)
    last_round = state["round_idx"]
    new_round = last_round + 1

    scored = list(csv.DictReader(open(scored_csv)))
    if not scored:
        raise SystemExit(f"scored CSV {scored_csv} is empty")

    # prev_state for the controller is its own {best, streak}; our extra keys are dropped.
    prev_state = {g: {"best": e.get("best"), "streak": e.get("streak", 0)}
                  for g, e in state["epitopes"].items()}
    plan, next_state = plan_round(scored, prev_state, last_round, cfg)

    prev_outputs = (Path(parent_outputs_dir) if parent_outputs_dir
                    else _round_dir(base_dir, last_round) / "02_rfd3" / "outputs")
    run_dir = _round_dir(base_dir, new_round)

    denovo_rows: list[dict] = []
    child_rows: list[dict] = []
    summary: dict[str, dict] = {}
    for g, p in plan.items():
        # carry the controller's updated best/streak back into persistent state
        if g in next_state:
            state["epitopes"].setdefault(g, {}).update(best=next_state[g]["best"],
                                                        streak=next_state[g]["streak"])
        if p["status"] == "stop":
            state["epitopes"][g]["status"] = "stopped"
            state["epitopes"][g]["reason"] = p["reason"]
            summary[g] = {"status": "stop", "reason": p["reason"], "n_denovo": 0, "n_child": 0,
                          "children_fraction": 0.0}
            continue
        state["epitopes"][g]["status"] = "open"
        dn = _take_denovo(state, g, p["n_denovo"], seed=new_round)
        ch = _make_children(state, g, p["children"], seed=new_round, parent_outputs=prev_outputs)
        denovo_rows += dn
        child_rows += ch
        summary[g] = {"status": "continue", "reason": p["reason"], "n_denovo": len(dn),
                      "n_child": len(ch), "children_fraction": p["children_fraction"]}

    all_rows = denovo_rows + child_rows
    _write_designs_table(run_dir / "designs_table.csv", all_rows)

    # emit: de-novo via stage01/stage02, children via the partial-diffusion emitter (append manifest)
    if denovo_rows:
        _emit_denovo(run_dir, state["antigen_pdb"], denovo_rows, seed=new_round)
    else:
        ensure_run_layout(run_dir)  # still create the layout + an empty manifest
        _append_manifest(run_dir, [])
    if child_rows:
        _append_manifest(run_dir, _emit_children(run_dir, child_rows))

    complete = len(all_rows) == 0
    state["round_idx"] = new_round
    state["rounds"].append({
        "round_idx": new_round, "run_dir": str(run_dir),
        "designs_table": str(run_dir / "designs_table.csv"),
        "n_denovo": len(denovo_rows), "n_child": len(child_rows),
        "epitopes_open": sorted(g for g, s in summary.items() if s["status"] == "continue"),
        "epitopes_stopped": sorted(g for g, s in summary.items() if s["status"] == "stop"),
    })
    save_state(base_dir, state)
    _write_summary(run_dir, new_round, cfg, summary, total_denovo=len(denovo_rows),
                   total_child=len(child_rows), complete=complete)
    return {"round_idx": new_round, "run_dir": str(run_dir), "plan": plan,
            "summary": summary, "complete": complete}


# --------------------------------------------------------------------------------------------------
# human-readable per-round summary
# --------------------------------------------------------------------------------------------------
def _write_summary(run_dir: Path, round_idx: int, cfg: IPDConfig, summary: dict,
                   total_denovo: int, total_child: int, complete: bool = False) -> None:
    run_dir = Path(run_dir)
    run_dir.mkdir(parents=True, exist_ok=True)
    lines = [f"# IPD round {round_idx:02d} summary",
             f"mode={cfg.mode}  budget/epitope={cfg.budget_per_epitope}  "
             f"topn_winners={cfg.topn_winners}  noise_levels={list(cfg.noise_levels)}  "
             f"max_rounds={cfg.max_rounds}",
             f"emitted this round: {total_denovo} de-novo + {total_child} children "
             f"= {total_denovo + total_child} RFd3 inputs",
             ""]
    for g in sorted(summary):
        s = summary[g]
        if s["status"] == "stop":
            lines.append(f"  [STOP]     {g}: {s['reason']}  (emits nothing)")
        elif s["status"] == "init":
            lines.append(f"  [INIT]     {g}: {s['n_denovo']} de-novo contigs  ({s['reason']})")
        else:
            lines.append(f"  [CONTINUE] {g}: {s['n_denovo']} de-novo + {s['n_child']} children "
                         f"(children_fraction={s['children_fraction']}; {s['reason']})")
    lines.append("")
    if complete:
        lines.append("ALL epitopes stopped -- IPD loop complete; nothing emitted for this round.")
    else:
        lines += [
            "NEXT (cluster, user-run, in order):",
            f"  1. RFd3        : generate backbones from {run_dir}/02_rfd3/inputs/ "
            f"-> {run_dir}/02_rfd3/outputs/",
            "  2. AOMclash    : scripts/aomclash_prefilter.py --backbones <outputs> "
            f"--contigs {run_dir}/designs_table.csv --antigen <antigen>",
            "  3. MPNN + AF3  : on the kept backbones (stage03_mpnn_* -> stage04_af3_*)",
            "  4. score       : metrics + composite -> a scored CSV (design_id,id,composite,is_pass)",
            f"  5. advance     : scripts/ipd_round.py advance --base-dir {run_dir.parent} "
            "--scored-csv <scored.csv>",
        ]
    (run_dir / "ROUND_SUMMARY.md").write_text("\n".join(lines) + "\n")
