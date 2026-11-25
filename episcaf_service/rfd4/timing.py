"""Per-stage wall-clock instrumentation for a campaign run -> GPU-hours per task type.

Why: the single-pass campaign runs generation once and then folds each surviving design serially
(campaign.run_campaign), so we don't know where a run's wall-time goes -- the thing we need before
designing the AF2 fan-out (see memory af2-parallelization-plan). This records a timed event per
stage, captures the GPU/SLURM context once, and writes a machine-readable `timings.json` that ships
back in the result bundle so a benchmark run reports its own breakdown (no SLURM-log digging).

GPU-hours note: the job holds one GPU for its whole wall and runs serially, so for a single-GPU run
`gpu_hours == wall_hours` per task type. If a run ever uses N GPUs at once, multiply that type's
wall by its concurrency -- recorded as `gpus` so the reader can.
"""
from __future__ import annotations

import json
import os
import platform
import re
import statistics
import subprocess
import time
from contextlib import contextmanager
from pathlib import Path


def gpu_context() -> dict:
    """SLURM + GPU facts captured once, so a timing record is attributable to a machine/GPU."""
    ctx = {
        "slurm_job_id": os.environ.get("SLURM_JOB_ID"),
        "node": platform.node(),
        "cuda_visible_devices": os.environ.get("CUDA_VISIBLE_DEVICES"),
    }
    try:
        out = subprocess.run(
            ["nvidia-smi", "--query-gpu=name,memory.total", "--format=csv,noheader"],
            capture_output=True, text=True, timeout=15)
        if out.returncode == 0:
            ctx["gpus"] = [ln.strip() for ln in out.stdout.strip().splitlines() if ln.strip()]
    except Exception:
        pass
    return ctx


# ColabFold writes a per-run log.txt; these pull the MSA cost and per-model inference cost out of it.
# Best-effort and version-tolerant: a miss just leaves the field absent (the stage total is always
# reliable from the wall-clock around the subprocess).
_MODEL_TOOK = re.compile(r"\b(model_\d+(?:_[\w.]+)?)\b.*?\btook\s+([\d.]+)s", re.I)
_MSA_TOOK = re.compile(r"\b(?:getting msa|running msa|msa generation|pairing)\b.*?([\d.]+)s", re.I)


def parse_colabfold_log(log_path: str | Path) -> dict:
    """Pull {msa_seconds?, models:[{name,seconds}], model_seconds} from a ColabFold log.txt."""
    try:
        text = Path(log_path).read_text(errors="replace")
    except OSError:
        return {}
    models = [{"name": m.group(1), "seconds": float(m.group(2))} for m in _MODEL_TOOK.finditer(text)]
    out: dict = {}
    if models:
        out["models"] = models
        out["model_seconds"] = round(sum(m["seconds"] for m in models), 3)
    msa = _MSA_TOOK.search(text)
    if msa:
        out["msa_seconds"] = float(msa.group(1))
    return out


class Recorder:
    """Collects timed stage events for one campaign; `stage()` wraps a block, `write()` dumps JSON."""

    def __init__(self):
        self.events: list[dict] = []
        self.gpu = gpu_context()
        self._start = time.time()
        self.counts: dict = {}   # free-form run facts (n_generated, n_folded, ...)

    def record(self, name: str, seconds: float, design_id: str | None = None, **extra):
        """Append one timed event directly (use when the duration is measured by the caller, e.g. to
        attach a parsed MSA/model split to a fold)."""
        ev = {"stage": name, "seconds": round(seconds, 3)}
        if design_id is not None:
            ev["design_id"] = design_id
        ev.update(extra)
        self.events.append(ev)

    @contextmanager
    def stage(self, name: str, design_id: str | None = None, **extra):
        t0 = time.perf_counter()
        try:
            yield
        finally:
            self.record(name, time.perf_counter() - t0, design_id, **extra)

    def note(self, **facts):
        self.counts.update(facts)

    def _per_type(self) -> dict:
        by: dict[str, list[float]] = {}
        for ev in self.events:
            by.setdefault(ev["stage"], []).append(ev["seconds"])
        out = {}
        for name, secs in by.items():
            total = sum(secs)
            out[name] = {
                "count": len(secs),
                "total_seconds": round(total, 3),
                "gpu_hours": round(total / 3600.0, 4),
                "mean_seconds": round(statistics.mean(secs), 3),
                "median_seconds": round(statistics.median(secs), 3),
                "min_seconds": round(min(secs), 3),
                "max_seconds": round(max(secs), 3),
            }
        return out

    def summary(self) -> dict:
        wall = time.time() - self._start
        per = self._per_type()
        accounted = sum(t["total_seconds"] for t in per.values())
        return {
            "gpu": self.gpu,
            "counts": self.counts,
            "wall_seconds": round(wall, 3),
            "wall_hours": round(wall / 3600.0, 4),
            "accounted_seconds": round(accounted, 3),
            "unaccounted_seconds": round(wall - accounted, 3),  # startup, weight loads, I/O, overhead
            "per_task_type": per,
            "note": "single-GPU serial run: gpu_hours == wall_hours per type; "
                    "unaccounted = process startup / model-weight loads / I/O outside timed stages.",
        }

    def write(self, path: str | Path):
        out = {"summary": self.summary(), "events": self.events}
        Path(path).write_text(json.dumps(out, indent=2))
        return out
