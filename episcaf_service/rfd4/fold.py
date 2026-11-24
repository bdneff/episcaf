"""The structure-prediction oracle — a pluggable foundation-model seam (Brandon's call).

One interface, swappable backends: AF2/ColabFold (default; the composite is calibrated to it, so we
validate the port here), then AF3 (already on Gemini) or Boltz-2 (Apache, commercial-safe). The
backend is a `config.OracleConfig` dial; swapping it needs the composite dials re-validated for that
oracle's PAE/confidence scale.

AF2 backend = localcolabfold `colabfold_batch` (5 models, 3 recycles, no relax, MSA on, no templates).
Output per query (confirmed on Gemini 2026-11-24):
  <name>_unrelaxed_rank_001_*.pdb              the kept model
  <name>_scores_rank_001_*.json               plddt (list), ptm, pae (NxN) / max_pae
  <name>_predicted_aligned_error_v1.json      the PAE matrix (fallback for mean PAE)
We score on `mean_pae` = mean of the PAE matrix (John's overall_pAE), plus the model PDB for RMSD.
"""
from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from episcaf_service.rfd4.config import DEFAULT, OracleConfig


@dataclass(frozen=True)
class FoldResult:
    model_pdb: str          # path to the rank_001 model
    mean_pae: float         # overall_pAE for the composite
    plddt: float            # mean pLDDT (diagnostic)
    ptm: float              # diagnostic


def af2_argv(fasta_path: str, out_dir: str, data_dir: str, colabfold_bin: str = "colabfold_batch",
             cfg: OracleConfig = DEFAULT.oracle) -> list[str]:
    """The `colabfold_batch` command. MSA on by default; num_relax=0 => no --amber; templates off."""
    argv = [colabfold_bin,
            "--num-models", str(cfg.num_models),
            "--num-recycle", str(cfg.num_recycles),
            "--data", data_dir]
    if cfg.num_relax > 0:
        argv += ["--amber", "--num-relax", str(cfg.num_relax)]
    argv += [fasta_path, out_dir]
    return argv


def _one(paths):
    hits = sorted(paths)
    if not hits:
        return None
    return hits[0]


def _mean_pae(out: Path, name: str, scores: dict) -> float:
    """Mean of the PAE matrix: prefer the scores json's pae, else the predicted_aligned_error file."""
    pae = scores.get("pae") or scores.get("predicted_aligned_error")
    if pae is None:
        paej = _one(out.glob(f"{name}_predicted_aligned_error*.json"))
        if paej is not None:
            obj = json.loads(Path(paej).read_text())
            # colabfold writes a dict, or sometimes a list wrapping one
            if isinstance(obj, list) and obj:
                obj = obj[0]
            pae = obj.get("predicted_aligned_error") or obj.get("pae")
    if pae is None:
        return float("nan")
    return float(np.asarray(pae, dtype=float).mean())


def ingest_af2(out_dir: str, name: str) -> FoldResult:
    """Read the rank_001 model + metrics for one query from a colabfold_batch output dir."""
    out = Path(out_dir)
    pdb = _one(out.glob(f"{name}_unrelaxed_rank_001_*.pdb")) or _one(out.glob(f"{name}_*rank_001*.pdb"))
    scores_path = _one(out.glob(f"{name}_scores_rank_001_*.json"))
    if pdb is None or scores_path is None:
        raise FileNotFoundError(f"no rank_001 model/scores for {name!r} under {out_dir}")
    scores = json.loads(Path(scores_path).read_text())
    plddt_raw = scores.get("plddt")
    plddt = float(np.mean(plddt_raw)) if isinstance(plddt_raw, (list, tuple)) else float(plddt_raw or "nan")
    ptm = float(scores.get("ptm", float("nan")))
    return FoldResult(model_pdb=str(pdb), mean_pae=_mean_pae(out, name, scores), plddt=plddt, ptm=ptm)


class Oracle:
    """Interface: predict(sequence, name, workdir) -> FoldResult. Backends implement _run()."""

    def predict(self, sequence: str, name: str, workdir: str) -> FoldResult:
        raise NotImplementedError


@dataclass
class ColabFoldAF2(Oracle):
    data_dir: str
    colabfold_bin: str = "colabfold_batch"
    cfg: OracleConfig = DEFAULT.oracle

    def predict(self, sequence: str, name: str, workdir: str) -> FoldResult:
        import subprocess
        work = Path(workdir)
        work.mkdir(parents=True, exist_ok=True)
        fasta = work / f"{name}.fasta"
        fasta.write_text(f">{name}\n{sequence}\n")
        out = work / f"{name}_af2"
        argv = af2_argv(str(fasta), str(out), self.data_dir, self.colabfold_bin, self.cfg)
        subprocess.run(argv, check=True)
        return ingest_af2(str(out), name)


def get_oracle(cfg: OracleConfig = DEFAULT.oracle, **kwargs) -> Oracle:
    """Factory: the oracle backend named by cfg.backend. af3/boltz are future plug-ins."""
    if cfg.backend == "af2":
        return ColabFoldAF2(cfg=cfg, **kwargs)
    raise NotImplementedError(f"oracle backend {cfg.backend!r} not implemented yet (af3/boltz are planned)")
