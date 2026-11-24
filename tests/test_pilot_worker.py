#!/usr/bin/env python3
"""Offline tests for the Gemini data-mover worker (episcaf_service/gcloud/pilot_worker.py).

GCS and Slurm are faked, so no cluster, no gcloud and no network are touched. Covers: run
classification (new / in-flight / done), auto-submit writing a live status.json, squeue/sacct ->
stage mapping, that a published run is never resubmitted, that --dry-run writes nothing, and that a
finished run's backbones get published to the result bucket.

Run:  /Users/bneff/micromamba/envs/bcell-repair/bin/python tests/test_pilot_worker.py
"""
import importlib.util
import io
import json
import os
import shutil
import sys
import tempfile
from collections import namedtuple
from pathlib import Path
from zipfile import ZipFile

ROOT = Path(__file__).resolve().parents[1]
WORKER_PATH = ROOT / "episcaf_service" / "gcloud" / "pilot_worker.py"

REQ = "episcaf-requests-test"
RES = "episcaf-results-test"
_INBOX = tempfile.mkdtemp(prefix="episcaf_worker_inbox_")
os.environ.update(EPISCAF_REQUEST_BUCKET=REQ, EPISCAF_RESULT_BUCKET=RES, EPISCAF_INBOX=_INBOX)

spec = importlib.util.spec_from_file_location("pilot_worker", WORKER_PATH)
pw = importlib.util.module_from_spec(spec)
spec.loader.exec_module(pw)

RID = "a" * 32
OWNER = "alice"

CP = namedtuple("CP", ["returncode", "stdout", "stderr"])


class FakeGcloud:
    """In-memory gs:// store that answers the worker's `storage ls|cat|cp` calls."""

    def __init__(self):
        self.store = {}  # uri -> bytes

    def __call__(self, *args, check=False, capture=True):
        a = list(args)
        if a[:2] == ["storage", "ls"]:
            pattern = a[2]
            prefix, _, suffix = pattern.partition("**")
            hits = [u for u in self.store
                    if u.startswith(prefix) and u.endswith(suffix) and "**" not in (prefix + suffix)]
            return CP(0 if hits else 1, "\n".join(sorted(hits)), "")
        if a[:2] == ["storage", "cat"]:
            data = self.store.get(a[2])
            return CP(0, data.decode() if data is not None else "", "") if data is not None else CP(1, "", "")
        if a[:2] == ["storage", "cp"]:
            rest = a[2:]
            once = "--if-generation-match=0" in rest
            paths = [x for x in rest if not x.startswith("--")]
            src, dst = paths[0], paths[1]
            if src.startswith("gs://"):  # download
                data = self.store.get(src)
                if data is None:
                    return CP(1, "", "")
                Path(dst).write_bytes(data)
                return CP(0, "", "")
            if once and dst in self.store:  # create-once lost the race
                return CP(1, "", "exists")
            self.store[dst] = Path(src).read_bytes()
            return CP(0, "", "")
        return CP(1, "", "")


def fresh(fake=None):
    """Reset worker module state for a test: clean inbox, fresh fake gcloud, flags off."""
    pw.DRY_RUN = False
    shutil.rmtree(_INBOX, ignore_errors=True)
    os.makedirs(_INBOX, exist_ok=True)
    pw.gcloud = fake if fake is not None else FakeGcloud()
    return pw.gcloud


def seed_request(fake, owner=OWNER, rid=RID):
    fake.store[f"gs://{REQ}/inbox/{owner}/{rid}/reference.zip"] = b"PK-fake-zip"


def designs_zip(design_ids=("d1", "d2")):
    retained = "design_id,score\n" + "\n".join(f"{d},0.1" for d in design_ids) + "\n"
    stream = io.BytesIO()
    with ZipFile(stream, "w") as archive:
        archive.writestr("retained.csv", retained)
        for d in design_ids:
            archive.writestr(f"backbones/{d}.pdb", "ATOM\n")
    return stream.getvalue(), retained.encode()


def test_discovery_classifies():
    fake = fresh()
    seed_request(fake, rid="a" * 32)
    seed_request(fake, rid="b" * 32)
    seed_request(fake, rid="c" * 32)
    fake.store[f"gs://{RES}/{'b' * 32}/status.json"] = b"{}"
    fake.store[f"gs://{RES}/{'c' * 32}/terminal.json"] = b"{}"
    phases = {rid: phase for phase, _owner, rid in pw.discover()}
    assert phases["a" * 32] == "new", phases
    assert phases["b" * 32] == "inflight", phases
    assert phases["c" * 32] == "done", phases


def test_submit_new_writes_status():
    fake = fresh()
    seed_request(fake)
    pw.slurm_submit = lambda rid: "12345"
    pw.tick()
    status = json.loads(fake.store[f"gs://{RES}/{RID}/status.json"])
    assert status["stage"] == "submitted" and status["job_id"] == "12345", status
    assert (Path(_INBOX) / RID / "request.zip").exists()


def test_poll_running_updates_status():
    fake = fresh()
    seed_request(fake)
    fake.store[f"gs://{RES}/{RID}/status.json"] = json.dumps(
        {"stage": "submitted", "job_id": "12345", "updated": 0}).encode()
    pw.slurm_state = lambda job: ("running", "RUNNING")
    pw.tick()
    status = json.loads(fake.store[f"gs://{RES}/{RID}/status.json"])
    assert status["stage"] == "running" and status["job_id"] == "12345", status


def test_done_is_never_resubmitted():
    fake = fresh()
    seed_request(fake)
    fake.store[f"gs://{RES}/{RID}/terminal.json"] = json.dumps({"stage": "completed"}).encode()

    def boom(rid):
        raise AssertionError("a completed run must not be resubmitted")

    pw.slurm_submit = boom
    pw.tick()  # must not raise


def test_dry_run_writes_nothing():
    fake = fresh()
    seed_request(fake)
    pw.DRY_RUN = True
    pw.slurm_submit = lambda rid: (_ for _ in ()).throw(AssertionError("no sbatch in dry-run"))
    pw.tick()
    assert f"gs://{RES}/{RID}/status.json" not in fake.store


def test_completion_publishes_backbones():
    fake = fresh()
    seed_request(fake)
    fake.store[f"gs://{RES}/{RID}/status.json"] = json.dumps(
        {"stage": "running", "job_id": "12345", "updated": 0}).encode()
    run_dir = Path(_INBOX) / RID
    run_dir.mkdir(parents=True, exist_ok=True)
    zip_bytes, retained = designs_zip()
    (run_dir / "designs.zip").write_bytes(zip_bytes)
    (run_dir / "retained.csv").write_bytes(retained)
    (run_dir / "terminal.json").write_text(json.dumps(
        {"id": RID, "stage": "completed", "job_id": "12345"}))
    pw.tick()
    for name in ("retained.csv", "designs.zip", "terminal.json"):
        assert f"gs://{RES}/{RID}/{name}" in fake.store, name
    published = json.loads(fake.store[f"gs://{RES}/{RID}/terminal.json"])
    assert published["stage"] == "completed", published


if __name__ == "__main__":
    for name, fn in sorted((n, f) for n, f in globals().items() if n.startswith("test_")):
        fn()
        print("ok", name)
    print("pilot worker OK")
