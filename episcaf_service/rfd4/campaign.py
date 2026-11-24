"""Single-pass campaign: antigen + epitope -> RFD4 scaffolds -> clash filter -> AF2 -> composite
-> ranked results. Chains the tested primitives (generate, read, clash, fold, score). The seeded
partial-diffusion round loop (phase 3) layers on top of evaluate_design + the generate/fold calls.

evaluate_design() is pure (structures in, metrics out) and unit-tested. run_campaign() adds the
subprocess calls to rfd4-generate and the oracle, which only run on Gemini.

Residue-numbering note (validate on first real run): the design CIF is 0-indexed (label_seq_id);
ColabFold PDBs are 1-indexed, so af2_resnum = design_seq_id + af2_offset (default 1). Native residues
come from the design's own src->native map.
"""
from __future__ import annotations

import csv
import json
import subprocess
from dataclasses import dataclass
from pathlib import Path
from zipfile import ZIP_DEFLATED, ZipFile

import numpy as np

from episcaf_service.rfd4 import clash, generate, read, score
from episcaf_service.rfd4.config import DEFAULT, EngineConfig


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


def _paired_epitope_backbone(design: read.RFD4Design, native: read.Structure):
    """design & native epitope backbones paired atom-for-atom (N,CA,C,O per residue), + the seq_ids used."""
    dbb, nbb, used = [], [], []
    for s in sorted(design.epitope_seq_ids):
        if s not in design.native_of:
            continue
        ch, rn = design.native_of[s]
        d = design.backbone_xyz([s])
        n = native.backbone_xyz(ch, [rn])
        if d.shape == (4, 3) and n.shape == (4, 3):
            dbb.append(d)
            nbb.append(n)
            used.append(s)
    if not used:
        return np.empty((0, 3)), np.empty((0, 3)), []
    return np.vstack(dbb), np.vstack(nbb), used


def backbone_prefilter_clashes(design: read.RFD4Design, native: read.Structure,
                               antibody_chains, cfg: EngineConfig = DEFAULT) -> int:
    """Clash count on the RFD4 backbone (pre-AF2), aligning the design epitope onto the native."""
    dbb, nbb, _ = _paired_epitope_backbone(design, native)
    if len(dbb) == 0:
        return 0
    scaffold_xyz, scaffold_resid = design.heavy_xyz_with_resid(design.scaffold_seq_ids)
    antibody = native.heavy_xyz(antibody_chains) if antibody_chains else np.empty((0, 3))
    antigen = native.heavy_xyz([c for c in {a.chain for a in native.atoms} if c not in set(antibody_chains or ())])
    shadow, _strict = clash.count_clashing_residues(dbb, nbb, scaffold_xyz, scaffold_resid,
                                                    antibody, antigen, cfg.clash)
    return shadow


def evaluate_design(design: read.RFD4Design, native: read.Structure, af2_struct: read.Structure,
                    mean_pae: float, antibody_chains, cfg: EngineConfig = DEFAULT,
                    af2_chain: str = "A", af2_offset: int = 1) -> dict:
    """The four metrics + composite for one AF2-folded design. Pure (no I/O)."""
    _dbb, nbb, used = _paired_epitope_backbone(design, native)
    af2_epi_bb = af2_struct.backbone_xyz(af2_chain, [s + af2_offset for s in used])

    epitope_rmsd = score.epitope_backbone_rmsd(af2_epi_bb, nbb) if len(nbb) and af2_epi_bb.shape == nbb.shape else float("nan")
    overall_rmsd = score.overall_ca_rmsd(af2_struct.ca_xyz([af2_chain]), design.ca_xyz())

    # clash on the AF2 model: align its epitope backbone onto native, test scaffold heavy vs antibody
    scaffold_ids = [s + af2_offset for s in design.scaffold_seq_ids]
    sc_atoms = [a for a in af2_struct.atoms if a.chain == af2_chain and a.seq_id in set(scaffold_ids)
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
             "af_clashes", "af_clashes_strict", "composite", "backbone_clashes"]


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
        for r in ranked:
            if r.model_pdb and Path(r.model_pdb).is_file():
                z.write(r.model_pdb, f"designs/{r.design_id}.pdb")
    return csv_path, ranked


def run_campaign(antigen_path: str, epitope: str, out_dir: str, *, oracle, native_path: str = None,
                 antibody_chains=(), rfd4_cmd=("rfd4-generate",), num_replicates: int = None,
                 cfg: EngineConfig = DEFAULT) -> list[DesignResult]:
    """Full single-pass run on Gemini. Writes final_results.csv + designs.zip under out_dir."""
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    native = read.read_pdb(native_path or antigen_path)
    n = num_replicates or cfg.generation.pilot_denovo_jobs * cfg.generation.designs_per_job

    spec = generate.condition_spec(antigen_path, epitope, cfg.generation)
    spec_path = out / "spec.yaml"
    generate.write_spec(spec, spec_path)
    gen_out = out / "rfd4"
    argv = list(rfd4_cmd) + generate.generate_argv(str(spec_path), str(gen_out), n)[1:]
    subprocess.run(argv, check=True)

    results: list[DesignResult] = []
    for cif in sorted(gen_out.glob("inference_samples/*/*_generated.cif")):
        design = read.read_rfd4_cif(str(cif))
        did = cif.parent.name
        bb_clash = backbone_prefilter_clashes(design, native, antibody_chains, cfg)
        if not clash.passes_backbone_prefilter(bb_clash, cfg.clash):
            continue
        fold = oracle.predict(design.sequence, did, str(out / "af2"))
        af2_struct = read.read_pdb(fold.model_pdb)
        m = evaluate_design(design, native, af2_struct, fold.mean_pae, antibody_chains, cfg)
        results.append(DesignResult(design_id=did, sequence=design.sequence, backbone_clashes=bb_clash,
                                    model_pdb=fold.model_pdb, **m))
    write_results(out, results)
    (out / "campaign.json").write_text(json.dumps(
        {"antigen": antigen_path, "epitope": epitope, "n_requested": n, "n_scored": len(results)}, indent=2))
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
    p.add_argument("--num-replicates", type=int)
    p.add_argument("--rfd4-cmd", default="rfd4-generate", help="generator command (shlex-split)")
    p.add_argument("--colabfold", default="colabfold_batch")
    p.add_argument("--af2-data", required=True, help="ColabFold params dir (--data)")
    args = p.parse_args()

    import shlex
    oracle = fold.ColabFoldAF2(data_dir=args.af2_data, colabfold_bin=args.colabfold, cfg=DEFAULT.oracle)
    ab = tuple(c for c in args.antibody_chains.split(",") if c)
    results = run_campaign(args.antigen, args.epitope, args.out, oracle=oracle, native_path=args.native,
                           antibody_chains=ab, rfd4_cmd=tuple(shlex.split(args.rfd4_cmd)),
                           num_replicates=args.num_replicates, cfg=DEFAULT)
    ranked = sorted((r for r in results if r.composite is not None), key=lambda r: -r.composite)
    print(f"{len(results)} designs scored; top composite "
          f"{ranked[0].composite:.3f} ({ranked[0].design_id})" if ranked else "no designs scored")


if __name__ == "__main__":
    main()
