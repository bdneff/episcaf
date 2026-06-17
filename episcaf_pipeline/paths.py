"""Run directory layout helpers."""

from __future__ import annotations
from dataclasses import dataclass
from pathlib import Path

@dataclass(frozen=True)
class RunPaths:
    run_dir: Path

    @property
    def input_dir(self) -> Path: return self.run_dir / "00_input"
    @property
    def design_dir(self) -> Path: return self.run_dir / "01_design"
    @property
    def rfd3_dir(self) -> Path: return self.run_dir / "02_rfd3"
    @property
    def af3_dir(self) -> Path: return self.run_dir / "03_af3"
    @property
    def analysis_dir(self) -> Path: return self.run_dir / "04_analysis"

    # stage-specific subdirs/files
    @property
    def contigs_parquet(self) -> Path: return self.design_dir / "contigs.parquet"

    @property
    def rfd3_inputs_dir(self) -> Path: return self.rfd3_dir / "inputs"
    @property
    def rfd3_manifest_csv(self) -> Path: return self.rfd3_dir / "inputs_manifest.csv"
    @property
    def rfd3_outputs_dir(self) -> Path: return self.rfd3_dir / "outputs"
    @property
    def rfd3_logs_dir(self) -> Path: return self.rfd3_dir / "logs"

    @property
    def af3_inputs_dir(self) -> Path: return self.af3_dir / "inputs"
    @property
    def af3_outputs_dir(self) -> Path: return self.af3_dir / "outputs"
    @property
    def af3_logs_dir(self) -> Path: return self.af3_dir / "logs"

    @property
    def results_parquet(self) -> Path: return self.analysis_dir / "results.parquet"
    @property
    def rmsd_all_csv(self) -> Path: return self.analysis_dir / "rmsd_vs_af3_all.csv"
    @property
    def rmsd_best_csv(self) -> Path: return self.analysis_dir / "rmsd_vs_af3_best_per_run.csv"


def ensure_run_layout(run_dir: Path) -> RunPaths:
    run_dir = Path(run_dir).resolve()
    rp = RunPaths(run_dir)
    for d in [rp.input_dir, rp.design_dir, rp.rfd3_dir, rp.af3_dir, rp.analysis_dir,
              rp.rfd3_inputs_dir, rp.rfd3_outputs_dir, rp.rfd3_logs_dir,
              rp.af3_inputs_dir, rp.af3_outputs_dir, rp.af3_logs_dir]:
        d.mkdir(parents=True, exist_ok=True)
    return rp
