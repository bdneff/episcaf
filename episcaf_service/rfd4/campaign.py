"""Single-pass campaign: antigen + epitope -> RFD4 scaffolds -> clash filter -> AF2 -> composite
-> ranked results. Chains the tested primitives (generate, read, clash, fold, score). The seeded
partial-diffusion round loop (phase 3) layers on top of evaluate_design + the generate/fold calls.

evaluate_design() is pure (structures in, metrics out) and unit-tested. run_campaign() adds the
subprocess calls to rfd4-generate and the structure-prediction step, which only run on Gemini.

Residue-numbering note (validate on first real run): the design CIF is 0-indexed (label_seq_id);
ColabFold PDBs are 1-indexed, so af2_resnum = design_seq_id + af2_offset (default 1). Native residues
come from the design's own src->native map.
"""
from __future__ import annotations

import csv
import json
import math
import subprocess
import threading
import time
from dataclasses import dataclass
from pathlib import Path
from zipfile import ZIP_DEFLATED, ZipFile

import numpy as np

from episcaf_service.rfd4 import clash, generate, read, score, timing, topology, viewer
from episcaf_service.rfd4.config import DEFAULT, EngineConfig


def epitope_native_residues(epitope) -> list:
    """Ordered [(chain, resnum)] the epitope string names, matching RFD4's src label order."""
    out = []
    for isl in topology.parse_islands(epitope) if isinstance(epitope, str) else epitope:
        out.extend((isl.chain, r) for r in range(isl.start, isl.end + 1))
    return out


def _src_to_native(design: "read.RFD4Design", epitope_native: list) -> dict:
    """Map RFD4's (possibly renumbered) src residue numbers -> the real native (chain, resnum).
    atomworks may renumber the epitope on PDB->CIF; the src residues sorted still correspond to the
    epitope residues in input order, so zip them."""
    src_sorted = sorted({rn for _ch, rn in design.native_of.values()})
    return {src: epitope_native[i] for i, src in enumerate(src_sorted) if i < len(epitope_native)}


@dataclass
class DesignResult:
    design_id: str
    sequence: str
    epitope_rmsd: float
    overall_rmsd: float
    overall_pae: float
    af_clashes: int
    af_clashes_strict: int
    composite: float
    backbone_clashes: int
    model_pdb: str
    aom_s_p1_free: float = None     # no-antibody accessibility (AOMclash); None when an antibody is used


def _paired_epitope_backbone(design: read.RFD4Design, native: read.Structure, epitope_native: list):
    """design & native epitope backbones paired atom-for-atom (N,CA,C,O per residue), + the seq_ids
    used. `epitope_native` is the ordered [(chain, resnum)] from the epitope string, used to remap
    RFD4's (possibly renumbered) src residues onto the real native structure."""
    remap = _src_to_native(design, epitope_native)
    dbb, nbb, used = [], [], []
    for s in sorted(design.epitope_seq_ids):
        if s not in design.native_of:
            continue
        _ch, src_rn = design.native_of[s]
        if src_rn not in remap:
            continue
        nchain, nresnum = remap[src_rn]
        d = design.backbone_xyz([s])
        n = native.backbone_xyz(nchain, [nresnum])
        if d.shape == (4, 3) and n.shape == (4, 3):
            dbb.append(d)
            nbb.append(n)
            used.append(s)
    if not used:
        return np.empty((0, 3)), np.empty((0, 3)), []
    return np.vstack(dbb), np.vstack(nbb), used


def backbone_prefilter_clashes(design: read.RFD4Design, native: read.Structure,
                               antibody_chains, epitope_native: list, cfg: EngineConfig = DEFAULT) -> int:
    """Clash count on the RFD4 backbone (pre-AF2), aligning the design epitope onto the native."""
    dbb, nbb, _ = _paired_epitope_backbone(design, native, epitope_native)
    if len(dbb) == 0:
        return 0
    scaffold_xyz, scaffold_resid = design.heavy_xyz_with_resid(design.scaffold_seq_ids)
    antibody = native.heavy_xyz(antibody_chains) if antibody_chains else np.empty((0, 3))
    antigen = native.heavy_xyz([c for c in {a.chain for a in native.atoms} if c not in set(antibody_chains or ())])
    shadow, _strict = clash.count_clashing_residues(dbb, nbb, scaffold_xyz, scaffold_resid,
                                                    antibody, antigen, cfg.clash)
    return shadow


def evaluate_design(design: read.RFD4Design, native: read.Structure, af2_struct: read.Structure,
                    mean_pae: float, antibody_chains, epitope_native: list, cfg: EngineConfig = DEFAULT,
                    af2_chain: str = "A") -> dict:
    """The four metrics + composite for one AF2-folded design. Pure (no I/O).

    The AF2 model and the design fold the SAME residue sequence, so pair them by ordinal position
    (i-th residue <-> i-th residue) rather than residue-number arithmetic -- correct whether the
    design CIF is 0- or 1-indexed and whatever the AF2 numbering is."""
    _dbb, nbb, used = _paired_epitope_backbone(design, native, epitope_native)
    af2_resnums = sorted({a.seq_id for a in af2_struct.atoms if a.chain == af2_chain})
    ordinal = {sid: i for i, sid in enumerate(design.seq_ids)}

    def af2_res(sid):
        i = ordinal.get(sid)
        return af2_resnums[i] if i is not None and i < len(af2_resnums) else None

    af2_epi_bb = af2_struct.backbone_xyz(af2_chain, [r for r in (af2_res(s) for s in used) if r is not None])
    epitope_rmsd = score.epitope_backbone_rmsd(af2_epi_bb, nbb) if len(nbb) and af2_epi_bb.shape == nbb.shape else float("nan")
    overall_rmsd = score.overall_ca_rmsd(af2_struct.ca_xyz([af2_chain]), design.ca_xyz())

    # clash on the AF2 model: align its epitope backbone onto native, test scaffold heavy vs antibody
    scaffold_af2 = {r for r in (af2_res(s) for s in design.scaffold_seq_ids) if r is not None}
    sc_atoms = [a for a in af2_struct.atoms if a.chain == af2_chain and a.seq_id in scaffold_af2
                and a.element.upper() != "H"]
    sc_xyz = np.array([a.xyz for a in sc_atoms], float) if sc_atoms else np.empty((0, 3))
    sc_resid = np.array([a.seq_id for a in sc_atoms])
    antibody = native.heavy_xyz(antibody_chains) if antibody_chains else np.empty((0, 3))
    antigen = native.heavy_xyz([c for c in {a.chain for a in native.atoms} if c not in set(antibody_chains or ())])
    af_clashes, af_clashes_strict = clash.count_clashing_residues(af2_epi_bb, nbb, sc_xyz, sc_resid,
                                                                  antibody, antigen, cfg.clash)

    comp = score.composite(af_clashes, epitope_rmsd, overall_rmsd, mean_pae, cfg.composite)
    return {"epitope_rmsd": epitope_rmsd, "overall_rmsd": overall_rmsd, "overall_pae": mean_pae,
            "af_clashes": af_clashes, "af_clashes_strict": af_clashes_strict, "composite": comp}


_CSV_COLS = ["design_id", "sequence", "epitope_rmsd", "overall_rmsd", "overall_pae",
             "af_clashes", "af_clashes_strict", "composite", "backbone_clashes", "aom_s_p1_free"]


def write_results(out_dir: Path, results: list[DesignResult]):
    ranked = sorted(results, key=lambda r: (r.composite is None, -(r.composite or -1e9)))
    csv_path = out_dir / "final_results.csv"
    with csv_path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=_CSV_COLS)
        w.writeheader()
        for r in ranked:
            w.writerow({c: getattr(r, c) for c in _CSV_COLS})
    with ZipFile(out_dir / "designs.zip", "w", compression=ZIP_DEFLATED) as z:
        z.write(csv_path, "final_results.csv")
        for extra in ("timings.json", "campaign.json"):   # ship the benchmark breakdown in the bundle
            if (out_dir / extra).is_file():
                z.write(out_dir / extra, extra)
        for r in ranked:
            if r.model_pdb and Path(r.model_pdb).is_file():
                z.write(r.model_pdb, f"designs/{r.design_id}.pdb")
            for ext in ("pdb", "json"):                # per-design viewer assets (if emitted)
                vf = out_dir / "viewer" / f"{r.design_id}.{ext}"
                if vf.is_file():
                    z.write(vf, f"viewer/{r.design_id}.{ext}")
    return csv_path, ranked


@dataclass
class Candidate:
    """One RFD4 backbone before the structure-prediction step: its id, the parsed design, its
    backbone clash count, and whether it cleared the pre-filter."""
    design_id: str
    design: "read.RFD4Design"
    backbone_clashes: int
    passes: bool


def select_scaffolds(pool: list[Candidate], target_n: int) -> list[Candidate]:
    """Pick exactly `target_n` scaffolds to carry into the structure-prediction step (AF2 — the
    expensive stage). Clash-free backbones come first (fewest clashes first); if fewer than N pass,
    fill with the lowest-clash of the rest so AF2 always runs on N. Returns fewer than N only when
    fewer than N were generated."""
    ordered = sorted(pool, key=lambda c: (not c.passes, c.backbone_clashes))
    return ordered[:target_n]


def wave_size(target_n: int, n_passed: int, n_generated: int,
              designs_per_job: int, max_tries: int) -> int:
    """How many backbones to request in the next generation wave. The passer shortfall inflated by
    the survival rate seen so far (x1.2 margin), at least one job, clamped to the remaining try
    budget. Returns 0 when N passers already exist or the try budget is spent -- the loop's stop
    condition. Backbone-clashing designs don't count against N, so a low survival rate just sizes a
    bigger next wave, up to `max_tries` total (then `select_scaffolds` fills the rest)."""
    remaining = max_tries - n_generated
    want = target_n - n_passed
    if want <= 0 or remaining <= 0:
        return 0
    if n_generated == 0:
        size = target_n
    else:
        survival = max(n_passed, 1) / n_generated
        size = math.ceil(want / survival * 1.2)
    return min(max(size, designs_per_job), remaining)


def fold_selected(selected, fold_one, workers: int, gpus: list) -> list:
    """Run `fold_one(cand, gpu)` over `selected`, serially when workers<=1 (gpu=None) or concurrently
    across `workers` GPUs otherwise. Returns results in input order. Each concurrently-running task
    holds a distinct GPU id from `gpus` (or its worker index when fewer gpus than workers are known),
    returned to the pool when it finishes — so no two live folds share a GPU."""
    if workers <= 1 or len(selected) <= 1:
        return [fold_one(c, None) for c in selected]
    import concurrent.futures
    import queue as _queue
    gpu_q = _queue.Queue()
    for i in range(workers):
        gpu_q.put(gpus[i] if i < len(gpus) else str(i))
    def worker(c):
        g = gpu_q.get()
        try:
            return fold_one(c, g)
        finally:
            gpu_q.put(g)
    with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as ex:
        return list(ex.map(worker, selected))


def run_campaign(antigen_path: str, epitope: str, out_dir: str, *, predictor, native_path: str = None,
                 antibody_chains=(), rfd4_cmd=("rfd4-generate",), num_replicates: int = None,
                 cfg: EngineConfig = DEFAULT, recorder: "timing.Recorder" = None,
                 progress=None) -> list[DesignResult]:
    """Full single-pass run on Gemini. Writes final_results.csv + designs.zip under out_dir.

    `predictor` is a fold.StructurePredictor (AF2/ColabFold by default). A timing.Recorder times every
    stage (generation, per-design backbone prefilter, structure prediction, scoring) -> timings.json,
    so a benchmark run reports GPU-hours per task type for the parallelization work."""
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    rec = recorder or timing.Recorder()

    def _p(label):                                          # step-aware status for the My-runs pill
        if progress:
            try:
                progress(label)
            except Exception:
                pass

    native = read.read_pdb(native_path or antigen_path)
    gen = cfg.generation
    # The user asks for `target_n` scaffolds (clamped to the UI cap). Backbone-clashing designs are
    # cheap to throw away, so they don't count against N: over-generate up to `max_tries` RFD4 designs.
    target_n = max(1, min(num_replicates or gen.target_scaffolds, gen.max_scaffolds))
    max_tries = max(target_n, gen.max_generation_tries)

    epitope_native = epitope_native_residues(epitope)      # ordered native residues for src remapping
    islands = topology.parse_islands(epitope)

    # Antibody auto-detection (John §1.4 contact fallback): no antibody chains given, but the uploaded
    # complex may contain a Fab -> detect it and score clashes against it (method=antibody). A true lone
    # antigen detects nothing -> the no-antibody path (AOMclash/zero). Keeps scored clashes consistent
    # with the viewer, and lets the growPatch avenue score accessibility when a complex is uploaded.
    if not antibody_chains:
        nat_atoms = viewer.parse_pdb(native_path or antigen_path)
        antigen_chains = list(dict.fromkeys(i.chain for i in islands))
        epi_set = set(epitope_native)
        epi_heavy = np.array([[a["x"], a["y"], a["z"]] for a in nat_atoms
                              if (a["chain"], a["resi"]) in epi_set and a["elem"].upper() != "H"], float)
        antibody_chains = tuple(viewer.detect_antibody_chains(nat_atoms, antigen_chains, epi_heavy))

    # No-antibody accessibility (AOMclash): no Fab -> read the occupancy map so each design still gets an
    # accessibility score (aom_s_p1_free) instead of nothing. Scorer is numpy/scipy + the vendored
    # aom_map.pkl (fast, local). How this feeds the composite is a separate call (discuss w/ John); here
    # we only compute + report it. Graceful: any failure (deps/map/frame) just leaves it None.
    aom_mod = aom_target = None
    if not antibody_chains:
        try:
            from episcaf_analysis.aomclash import score as aom_mod  # noqa: PLC0415
            aom_target = aom_mod.build_target(native_path or antigen_path, epitope)
        except Exception:
            aom_mod = aom_target = None

    compact = gen.compact_input if gen.compact_input is not None else topology.needs_compact_input(islands)
    if compact:
        # multi-chain / high-resnum epitope -> feed RFD4 a single-chain, renumbered motif. Scoring still
        # uses the original native structure + epitope string (the motif map remaps by order).
        compact_islands = generate.write_compact_antigen(
            antigen_path, islands, out / "compact_antigen.pdb", gen.compact_gap)
        antigen_cif = generate.ensure_cif(str(out / "compact_antigen.pdb"), out)
        spec = generate.condition_spec(antigen_cif, compact_islands, gen)
    else:
        antigen_cif = generate.ensure_cif(antigen_path, out)   # RFD4 needs CIF; converts a PDB upload
        spec = generate.condition_spec(antigen_cif, epitope, gen)
    spec_path = out / "spec.yaml"
    generate.write_spec(spec, spec_path)
    gen_out = out / "rfd4"

    # --- generation budget: waves of RFD4 + backbone pre-filter until N clash-free backbones (or the
    # try cap). Clashers don't count against N; we only pay the structure-prediction step on N. ---
    _p("generating scaffolds")
    pool: list[Candidate] = []
    n_generated = 0
    wave = 0
    while True:
        n_passed = sum(1 for c in pool if c.passes)
        batch = wave_size(target_n, n_passed, n_generated, gen.designs_per_job, max_tries)
        if batch <= 0:
            break
        wave_out = gen_out / f"wave_{wave}"
        argv = list(rfd4_cmd) + generate.generate_argv(str(spec_path), str(wave_out), batch)[1:]
        with rec.stage("generation", wave=wave, n_requested=batch):
            subprocess.run(argv, check=True)
        for cif in sorted(wave_out.glob("inference_samples/*/*_generated.cif")):
            design = read.read_rfd4_cif(str(cif))
            did = f"scaffold_{n_generated}"                 # global id across waves (unique for AF2)
            n_generated += 1
            with rec.stage("backbone_prefilter", did):
                bb_clash = backbone_prefilter_clashes(design, native, antibody_chains, epitope_native, cfg)
            pool.append(Candidate(did, design, bb_clash,
                                  clash.passes_backbone_prefilter(bb_clash, cfg.clash)))
        wave += 1

    selected = select_scaffolds(pool, target_n)
    n_passed = sum(1 for c in pool if c.passes)

    # --- structure-prediction step on the selected N, then score. The per-design fold is independent,
    # so it runs serially (fold_workers=1, default) or concurrently across the job's GPUs (one colabfold
    # subprocess per GPU). This is the stage the future cross-node SLURM array would shard instead. ---
    view_lock = threading.Lock()
    n_total = len(selected)
    done = [0]
    done_lock = threading.Lock()
    _p(f"predicting structures 0/{n_total}")

    def fold_and_score(cand: Candidate, gpu=None) -> DesignResult:
        design, did, bb_clash = cand.design, cand.design_id, cand.backbone_clashes
        t0 = time.perf_counter()                            # timed by hand so the fold event carries
        fold = predictor.predict(design.sequence, did, str(out / "af2"), gpu=gpu)   # the MSA/model split
        elapsed = time.perf_counter() - t0
        detail = timing.parse_colabfold_log(Path(fold.out_dir) / "log.txt") if fold.out_dir else {"log_found": False}
        if "model_seconds" in detail:                       # MSA = stage wall - model inference
            detail["msa_seconds"] = round(max(0.0, elapsed - detail["model_seconds"]), 3)
        rec.record("structure_prediction", elapsed, did, gpu=gpu, **detail)
        with rec.stage("scoring", did):
            af2_struct = read.read_pdb(fold.model_pdb)
            m = evaluate_design(design, native, af2_struct, fold.mean_pae, antibody_chains, epitope_native, cfg)
        # viewer assets: AF2 model aligned onto the native epitope + native antigen/antibody (detected
        # if antibody_chains is empty). Serialised: all designs share one viewer/native.pdb.
        with view_lock:
            cfg_v = viewer.write_viewer_files(out, did, native_path or antigen_path, fold.model_pdb,
                                              epitope, antibody_chains or None)
        aom = None
        if aom_target is not None and cfg_v is not None:   # no-antibody accessibility for this design
            try:
                epi_pos = cfg_v.get("epitopeScaffoldOrdered") or cfg_v["epitopeScaffold"]
                s = aom_mod.aomclash_score(aom_target, fold.model_pdb, epi_pos)
                aom = round(float(s), 4) if s is not None else None
            except Exception:
                aom = None
        with done_lock:
            done[0] += 1
            _p(f"predicting structures {done[0]}/{n_total}")
        return DesignResult(design_id=did, sequence=design.sequence, backbone_clashes=bb_clash,
                            model_pdb=fold.model_pdb, aom_s_p1_free=aom, **m)

    workers = max(1, cfg.predictor.fold_workers)
    results = fold_selected(selected, fold_and_score, workers, timing.visible_gpus())
    rec.note(target_scaffolds=target_n, n_generated=n_generated, n_passed_prefilter=n_passed,
             n_selected=len(selected), n_folded=len(results), n_scored=len(results), fold_workers=workers)
    (out / "campaign.json").write_text(json.dumps(
        {"antigen": antigen_path, "epitope": epitope, "target_scaffolds": target_n,
         "antibody_chains": list(antibody_chains), "n_generated": n_generated,
         "n_passed_prefilter": n_passed, "n_scored": len(results)}, indent=2))
    rec.write(out / "timings.json")                        # before write_results so it lands in the zip
    write_results(out, results)
    return results


def main():
    import argparse
    from episcaf_service.rfd4 import fold
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--antigen", required=True, help="antigen structure fed to RFD4 (epitope chain)")
    p.add_argument("--epitope", required=True, help="islands, e.g. 'A:205-216;A:278-292'")
    p.add_argument("--native", help="full native complex for clash/RMSD (default: --antigen)")
    p.add_argument("--antibody-chains", default="", help="comma-separated antibody chain ids")
    p.add_argument("--out", required=True)
    p.add_argument("--num-replicates", "--num-scaffolds", dest="num_replicates", type=int,
                   help="N scaffolds to fold (default 32, UI cap 128); RFD4 over-generates to find N "
                        "clash-free backbones, so this is not the number of RFD4 designs made")
    p.add_argument("--rfd4-cmd", default="rfd4-generate", help="generator command (shlex-split)")
    p.add_argument("--colabfold", default="colabfold_batch")
    p.add_argument("--af2-data", required=True, help="ColabFold params dir (--data)")
    args = p.parse_args()

    import shlex
    predictor = fold.ColabFoldAF2(data_dir=args.af2_data, colabfold_bin=args.colabfold, cfg=DEFAULT.predictor)
    ab = tuple(c for c in args.antibody_chains.split(",") if c)
    results = run_campaign(args.antigen, args.epitope, args.out, predictor=predictor, native_path=args.native,
                           antibody_chains=ab, rfd4_cmd=tuple(shlex.split(args.rfd4_cmd)),
                           num_replicates=args.num_replicates, cfg=DEFAULT)
    ranked = sorted((r for r in results if r.composite is not None), key=lambda r: -r.composite)
    print(f"{len(results)} designs scored; top composite "
          f"{ranked[0].composite:.3f} ({ranked[0].design_id})" if ranked else "no designs scored")


if __name__ == "__main__":
    main()
