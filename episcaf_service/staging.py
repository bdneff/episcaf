"""Shared staging logic for the episcaf tool (local stdlib server + cluster FastAPI service).

This is the single source of truth for the two CPU helpers both servers use -- enumerating candidate
epitopes for the viewer (`viewer_data`) and turning chosen epitopes into RFd3 contigs (`stage_run`) --
plus the two things the cluster service adds on top: persisting a staged run under $WS
(`persist_run`) and emitting the exact cluster command sequence (`build_submit_sh`).

Nothing here runs a GPU step. It stages the job; the user submits it on the cluster.

Every command `build_submit_sh` emits was verified against the real script interfaces
(`scripts/ipd_round.py init`, `aomclash_prefilter.py`, `aomclash_score_designs.py`,
`stage06_select.py`, `stage04_af3_array.sbatch`) and the documented NN_stage run-dir convention in
`docs/PIPELINE.md`. The ProteinMPNN->AF3 middle is the sequence the repo itself marks
"cluster-specific, not verified off-cluster" (`run_whole_epitope_mpnn_af3.sh`); it is referenced,
not reinvented.
"""
from __future__ import annotations

import csv as _csv
import datetime
import io
import json
import os
import re

from episcaf_analysis.growpatch import read as R
from episcaf_pipeline import build_growpatch_designs as BGD
from episcaf_service import disco as _disco


def _resi(label):
    m = re.match(r"-?\d+", str(label))
    return int(m.group()) if m else None


def viewer_data(cif_path, chain, sasa_cut, footprint_radius, methods, top_k, mode="growpatch"):
    """Residues + candidate patches for the viewer (like docs/data/growpatch_*.js).

    mode 'growpatch' -> enumerate candidate epitopes with growPatch (lone antigen), ranked by mean
                        DiscoTope score when a scorer is configured (see episcaf_service.disco).
    mode 'antibody'  -> the epitope is already defined by the bound partner; return it as one patch,
                        no enumeration (growPatch unnecessary).

    Each growPatch candidate carries `resis` (the full gap-filled residue set sent to scaffold design)
    and `resis_rank` (the growth-selected surface residues used for ranking). The full candidate set
    is generated and deduplicated, then ranked, and only then truncated to `top_k` -- so the top-k are
    the highest-ranked candidates, never an early geometric slice.
    """
    ag = R.read_antigen(cif_path, chain, relsasa=True)
    residues = [dict(resi=_resi(ag.label(c, s)), resn=r["resname"],
                     rel=round(float(ag.rs(r["key"])), 3),
                     surf=bool(ag.rs(r["key"]) >= sasa_cut),
                     p=[round(float(x), 2) for x in r["rep"]])
                for r in ag.residues for (c, s) in [r["key"]]]
    if mode == "antibody":
        det = R.detect_complex(cif_path, chain)
        epi = det["epitope"]
        resis = sorted({_resi(ag.label(c, s)) for (c, s) in ag.parse_epitope(epi)}) if epi else []
        patches = ([dict(i=0, methods=["antibody-contact"], epitope=epi, n_res=len(resis),
                         n_isl=len([t for t in epi.split(";") if t]), resis=resis)] if resis else [])
        extra = dict(mode="antibody", partner_chains=det["partner_chains"], has_disco=False)
        return dict(chain=chain, sasa_cut=sasa_cut, n_residues=len(residues),
                    n_surface=sum(r["surf"] for r in residues), n_patches=len(patches),
                    residues=residues, patches=patches, **extra)

    # lone antigen: full growPatch set (no early truncation), then rank, then truncate.
    cands, _ = R.enumerate_candidates(cif_path, chain, antigen=ag, methods=methods,
                                      sasa_cut=sasa_cut, footprint_radius=footprint_radius)
    patches = [dict(i=p["patch"], epitope=p["epitope"], n_res=p["n_res"], n_isl=p["n_islands"],
                    resis=sorted({_resi(ag.label(c, s)) for (c, s) in p["keys"]}),
                    resis_rank=sorted({_resi(ag.label(c, s)) for (c, s) in p["keys_raw"]}))
               for p in cands]

    scores = _disco.score_cif(cif_path, chain)
    has_disco = scores is not None
    if has_disco:
        for p in patches:
            vals = [scores[r] for r in p["resis_rank"] if r in scores]
            p["disco_mean"] = round(sum(vals) / len(vals), 4) if vals else None
            p["disco_n"] = len(vals)
            p["disco_total"] = len(p["resis_rank"])
        # rank within the antigen by descending mean; No-score candidates last; stable id tie-break.
        patches.sort(key=lambda p: (p["disco_mean"] is None, -(p["disco_mean"] or 0.0), p["i"]))
    if top_k:
        patches = patches[:top_k]
    return dict(chain=chain, sasa_cut=sasa_cut, n_residues=len(residues),
                n_surface=sum(r["surf"] for r in residues), n_patches=len(patches),
                has_disco=has_disco, residues=residues, patches=patches, mode="growpatch")


def stage_run(cif_path, chain, epitopes, cfg):
    """Build the RFd3 contigs for the chosen epitopes + a run spec the cluster consumes.

    Returns (run_spec dict, contigs_csv str, n_contigs, skipped). The heavy run (RFd3->MPNN->AF3) is
    cluster/GPU work; this stages the job, it does not run it.
    """
    chains = R.parse_chains(chain)
    cands = [dict(patch=i, epitope=e, n_islands=len([t for t in e.split(";") if t]),
                  n_res=0, methods=[]) for i, e in enumerate(epitopes)]
    rows, skipped = BGD.build_rows(cfg["name"], chains, cands, cfg["total_len"],
                                   cfg["variants"], cfg["seed"])
    cols = ["id", "contig_id", "contig_string", "contig_length", "epitope_chunks", "source",
            "pdb_id", "chains", "growpatch_epitope", "growpatch_methods", "growpatch_n_res", "variant"]
    buf = io.StringIO()
    w = _csv.DictWriter(buf, fieldnames=cols)
    w.writeheader()
    w.writerows({k: r.get(k) for k in cols} for r in rows)
    n_epi = len(epitopes)
    run_spec = {
        "created": datetime.datetime.now().isoformat(timespec="seconds"),
        "antigen": cfg["name"], "chain": chain, "mode": cfg.get("mode", "growpatch"),
        "n_epitopes": n_epi, "epitopes": epitopes,
        "generation": {"flank_variants": cfg["variants"], "total_len": cfg["total_len"],
                       "rfd3_designs_per_contig": cfg["rfd3"], "mpnn_seqs_per_backbone": cfg["mpnn"],
                       "helix_bias": cfg.get("helix_bias", 0.0)},   # RFd3 secondary-structure bias (0 = neutral)
        "budget": cfg["budget"],
        "n_contigs": len(rows),
        "designs_at_full_fanout": len(rows) * cfg["rfd3"] * cfg["mpnn"],
        "note": "Stage only. Run RFd3->MPNN->AF3 on the cluster; score + view the metrics CSV in "
                "the Results tab. accessibility metric = AOMclash (aomclash_pred).",
    }
    return run_spec, buf.getvalue(), len(rows), skipped


# --------------------------------------------------------------------------------------------------
# Cluster service additions: persist a staged run under $WS and emit the real submit script.
# --------------------------------------------------------------------------------------------------

def build_submit_sh(run_dir, input_name, preset, settings):
    """The exact cluster command sequence for a staged run, as a runnable bash script.

    `settings` keys: ipd(bool), ipdrounds(int), ipd_mode(str), prefilter(bool), keepfrac(float),
    topn(int), group(str). All flags verified against the real script `--help` interfaces.
    """
    rd = run_dir                      # absolute, under $WS
    inp = f"$RUN_DIR/{input_name}"
    ipd = bool(settings.get("ipd"))
    rounds = int(settings.get("ipdrounds", 1))
    ipd_mode = settings.get("ipd_mode", "adaptive")
    topn = int(settings.get("topn", 10))
    group = settings.get("group", "id")

    L = ["#!/usr/bin/env bash",
         f"# episcaf staged run -- {os.path.basename(rd)}",
         "# Generated by episcaf_service /submit. Run on the cluster (gemini). This stages & scores;",
         "# the RFd3 and AF3 steps print the sbatch array you submit (you stay in control of the GPU).",
         "set -euo pipefail",
         "conda activate ~/rfd3/env/rfd3_py312",
         f'export RUN_DIR="{rd}"',
         "# input.pdb, designs.csv, run_spec.json are already staged in $RUN_DIR by the service.",
         ""]

    # 1. RFd3 generation, driven by the staged designs table (the one entry that consumes designs.csv).
    L += ["# 1. RFd3 backbones from the staged contigs"
          + (f" (IPD loop: {ipd_mode} children, partial diffusion, <= {rounds} rounds)" if ipd
             else " (single round)"),
          f'python scripts/ipd_round.py init --base-dir "$RUN_DIR" --antigen "{inp}" \\',
          f'    --designs-table "$RUN_DIR/designs.csv" --mode {ipd_mode} --max-rounds {rounds if ipd else 1}',
          '#    submit the RFd3 array the init prints, e.g.:',
          '#    sbatch --array=1-N episcaf_pipeline/hpc/sbatch/rfd3_array.sbatch "$RUN_DIR"',
          ""]
    step = 2

    # 2. optional AOMclash backbone pre-filter (verified flags).
    if settings.get("prefilter"):
        keep = float(settings.get("keepfrac", 0.5))
        L += [f"# {step}. AOMclash backbone pre-filter: keep best {round(keep*100)}% per epitope before MPNN/AF3",
              'python scripts/aomclash_prefilter.py --backbones "$RUN_DIR/02_rfd3/outputs" \\',
              f'    --contigs "$RUN_DIR/designs.csv" --antigen "{inp}" --keep-frac {keep} \\',
              '    --out "$RUN_DIR/aomclash_prefilter.csv"',
              ""]
        step += 1

    # 3. MPNN -> AF3. Repo marks this sequence cluster-specific / unverified off-cluster.
    L += [f"# {step}. kept backbones -> ProteinMPNN -> AlphaFold3",
          "#    cluster-specific; follow scripts/run_whole_epitope_mpnn_af3.sh for the exact MPNN/AF3",
          "#    env on this cluster. The canonical sequence (paths follow the NN_stage convention):",
          'python scripts/stage03_mpnn_fixed_pdbs.py --run "$RUN_DIR"'
          + ('   # if prefiltered, point at the .pass backbones' if settings.get("prefilter") else ""),
          'python scripts/stage03_mpnn_submit.py --fixed_pdb_dir "$RUN_DIR/03_mpnn/fixed_pdbs" \\',
          '    --outdir "$RUN_DIR/03_mpnn/out"',
          'python scripts/stage04_af3_emit_jsons.py --mpnn_pdb_dir "$RUN_DIR/03_mpnn/out" \\',
          '    --out_dir "$RUN_DIR/04_af3/in" --seed 1',
          'sbatch --array=1-N scripts/stage04_af3_array.sbatch "$RUN_DIR"',
          ""]
    step += 1

    # 4. AOMclash on the AF3 designs, then score + select (verified flags).
    L += [f"# {step}. AOMclash accessibility on the AF3 designs, then composite score + select (top {topn})",
          'python scripts/aomclash_score_designs.py --designs "$RUN_DIR/04_af3/outputs" \\',
          f'    --contigs "$RUN_DIR/designs.csv" --antigen "{inp}" \\',
          '    --merge "$RUN_DIR/05_analysis/metrics.csv" --out "$RUN_DIR/05_analysis/metrics.csv"',
          f'python scripts/stage06_select.py --preset {preset} \\',
          '    --metrics-csv "$RUN_DIR/05_analysis/metrics.csv" \\',
          f'    --group {group} --topk {topn} --out "$RUN_DIR/05_analysis/ranked.csv"',
          "#    or fetch the scored table from the service: GET /results/<run_id>",
          ""]
    return "\n".join(L)


def persist_run(ws_root, structure_text, fmt, chain, mode, epitopes, cfg, preset, settings):
    """Write a staged run under $WS and return its manifest. Durable artifacts, never /scratch.

    Creates $WS/runs/<run_id>/ with input.(pdb|cif), designs.csv, run_spec.json, submit.sh.
    Returns {run_id, run_dir, submit_sh, n_contigs, designs_at_full_fanout, skipped, n_epitopes}.
    """
    stamp = datetime.datetime.utcnow().strftime("%Y%m%dT%H%M%SZ")
    safe = re.sub(r"[^A-Za-z0-9_.-]", "_", cfg.get("name", "antigen")) or "antigen"
    run_id = f"{safe}_{stamp}"
    run_dir = os.path.join(ws_root, "runs", run_id)
    os.makedirs(run_dir, exist_ok=True)

    input_name = "input.cif" if fmt == "cif" else "input.pdb"
    input_path = os.path.join(run_dir, input_name)
    with open(input_path, "w") as fh:
        fh.write(structure_text)

    # Resolve epitopes: known-antibody -> the contact epitope; lone antigen -> the chosen patches.
    if mode == "antibody" and not epitopes:
        det = R.detect_complex(input_path, chain)
        epitopes = [det["epitope"]] if det.get("epitope") else []
    if not epitopes:
        raise ValueError("no epitopes to stage (lone antigen needs chosen patches; "
                         "complex needs a detectable contact epitope)")

    run_spec, csv_text, n_contigs, skipped = stage_run(input_path, chain or "A", epitopes, cfg)
    with open(os.path.join(run_dir, "designs.csv"), "w") as fh:
        fh.write(csv_text)

    submit_sh = build_submit_sh(run_dir, input_name, preset, settings)
    sh_path = os.path.join(run_dir, "submit.sh")
    with open(sh_path, "w") as fh:
        fh.write(submit_sh + "\n")
    os.chmod(sh_path, 0o755)

    run_spec["run_id"] = run_id
    run_spec["preset"] = preset
    run_spec["settings"] = settings
    with open(os.path.join(run_dir, "run_spec.json"), "w") as fh:
        json.dump(run_spec, fh, indent=2)

    return dict(run_id=run_id, run_dir=run_dir, submit_sh=submit_sh, n_contigs=n_contigs,
                designs_at_full_fanout=run_spec["designs_at_full_fanout"],
                skipped=len(skipped), n_epitopes=len(epitopes))
