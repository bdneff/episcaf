#!/usr/bin/env python3
"""episcaf Gemini data-mover worker: one poll tick.

Runs on the Gemini ``data-mover`` node, fired every ~5 min by ``scrontab`` (Slurm cron). Each tick
is self-contained and idempotent; scrontab skips a new tick while the previous one is still running.
A tick:

  1. **Discovers** pilot runs by listing the request bucket (``inbox/<owner>/<rid>/reference.zip``)
     and the result bucket's per-run markers, and classifies each run NEW / IN-FLIGHT / DONE.
  2. **Submits** each NEW reference run to Slurm: fetches ``reference.zip`` from GCS into
     ``~/.episcaf-pilot/inbox/<rid>/request.zip`` and ``sbatch``es ``handoff.sbatch`` to ``gpu-a100``.
  3. **Polls** each IN-FLIGHT run with ``sacct``/``squeue`` and keeps a live ``status.json`` in the
     result bucket; on completion it publishes the backbones, on Slurm failure it writes a terminal
     failure. ``terminal.json`` present = done, never resubmit.

This is the Gemini-native replacement for the laptop ``pilot_transfer.py`` hops. Only this node
reaches GCS; ``gcloud`` runs through the ``google/cloud-sdk`` singularity image with a
buckets-scoped service-account key (see ``run_worker.sh`` and ``episcaf_service/gcloud/README.md``).

Standalone + stdlib-only on purpose: the data mover has no project virtualenv. Run it as a script,
never ``-m``. Idempotency is belt-and-suspenders: scrontab serialises ticks, a create-once
``status.json`` claims a run before ``sbatch``, and the compute job itself refuses a duplicate run id.

Config (environment, set by ``run_worker.sh``):
  EPISCAF_REQUEST_BUCKET   request bucket name, e.g. ``episcaf-requests-123456789`` (required)
  EPISCAF_RESULT_BUCKET    result bucket name,  e.g. ``episcaf-results-123456789``  (required)
  EPISCAF_GCLOUD           gcloud invocation, default ``gcloud`` (run_worker.sh sets the
                           singularity-wrapped form); shlex-split
  EPISCAF_INBOX            local inbox root, default ``~/.episcaf-pilot/inbox``
  EPISCAF_SBATCH           path to handoff.sbatch, default ``~/episcaf-advisor/episcaf_service/handoff.sbatch``
  EPISCAF_GRACE_SECONDS    seconds a COMPLETED-but-resultless job is tolerated before it is failed,
                           default 900

Usage:  python3 pilot_worker.py tick [--dry-run]
"""
import argparse
import csv
import io
import json
import os
import re
import shlex
import subprocess
import sys
import tempfile
import time
from pathlib import Path
from zipfile import ZipFile

RID_RE = re.compile(r"[a-f0-9]{32}")
OWNER_RE = re.compile(r"[a-zA-Z0-9_-]{1,48}")
JOB_RE = re.compile(r"[0-9]+")

# Slurm states that mean the job is finished and failed (no usable result).
SLURM_FAILED = {"FAILED", "TIMEOUT", "OUT_OF_MEMORY", "NODE_FAIL", "CANCELLED",
                "BOOT_FAIL", "DEADLINE", "PREEMPTED", "REVOKED", "SPECIAL_EXIT"}
SLURM_PENDING = {"PENDING", "CONFIGURING", "RESV_DEL_HOLD", "REQUEUE_HOLD",
                 "REQUEUED", "RESIZING", "SUSPENDED"}
SLURM_RUNNING = {"RUNNING", "COMPLETING", "STAGE_OUT", "SIGNALING"}

REQUEST_BUCKET = os.environ.get("EPISCAF_REQUEST_BUCKET", "")
RESULT_BUCKET = os.environ.get("EPISCAF_RESULT_BUCKET", "")
GCLOUD = shlex.split(os.environ.get("EPISCAF_GCLOUD", "gcloud"))
INBOX = Path(os.environ.get("EPISCAF_INBOX", str(Path.home() / ".episcaf-pilot" / "inbox")))
SBATCH = Path(os.environ.get(
    "EPISCAF_SBATCH", str(Path.home() / "episcaf-advisor" / "episcaf_service" / "handoff.sbatch")))
GRACE_SECONDS = float(os.environ.get("EPISCAF_GRACE_SECONDS", "900"))

DRY_RUN = False


def log(message):
    print(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {message}", flush=True)


# --- gcloud storage (via the configured, possibly singularity-wrapped, invocation) --------------

def gcloud(*args, check=False, capture=True):
    kw = {"text": True}
    if capture:
        kw["stdout"] = subprocess.PIPE
        kw["stderr"] = subprocess.PIPE
    return subprocess.run([*GCLOUD, *args], check=check, **kw)


def gs_ls(bucket, pattern):
    """Return the object URIs under gs://bucket/pattern; empty list when nothing matches."""
    result = gcloud("storage", "ls", f"gs://{bucket}/{pattern}")
    if result.returncode != 0:
        return []  # "matched no objects" (nothing staged yet) or a transient error; next tick retries.
    return [line.strip() for line in (result.stdout or "").splitlines() if line.strip()]


def gs_read_json(uri):
    """Download and parse a JSON object, or None when it is absent / unreadable."""
    result = gcloud("storage", "cat", uri)
    if result.returncode != 0 or not (result.stdout or "").strip():
        return None
    try:
        return json.loads(result.stdout)
    except json.JSONDecodeError:
        return None


def gs_cp_down(uri, dest):
    gcloud("storage", "cp", uri, str(dest), check=True)


def gs_write_json(uri, obj, once=False):
    """Write a small JSON object to GCS. once=True creates immutably (claim / terminal); else overwrites."""
    with tempfile.NamedTemporaryFile("w", suffix=".json", delete=False) as handle:
        json.dump(obj, handle)
        temp = handle.name
    try:
        args = ["storage", "cp"]
        if once:
            args.append("--if-generation-match=0")
        args += [temp, uri]
        result = gcloud(*args)
        return result.returncode == 0
    finally:
        os.unlink(temp)


def gs_upload_once(path, uri):
    """Immutable upload used for result artifacts; a retry must match the existing bytes."""
    result = gcloud("storage", "cp", "--if-generation-match=0", str(path), uri, capture=False)
    if result.returncode != 0:
        # The object already exists: a retry is only valid if the bytes are identical.
        check = subprocess.run([*GCLOUD, "storage", "cat", uri], stdout=subprocess.PIPE)
        if check.returncode != 0 or check.stdout != path.read_bytes():
            raise RuntimeError(f"A different result already exists at {uri}")


# --- result-artifact validation (mirrors pilot_transfer.publish; keep the two in sync) ----------

def validate_designs(run_dir, rid, job):
    """Assert terminal.json matches this run+job and the retained table agrees with designs.zip."""
    record = json.loads((run_dir / "terminal.json").read_text())
    if record.get("id") != rid or str(record.get("job_id")) != str(job) or record.get("stage") != "completed":
        raise RuntimeError("terminal.json does not match this run and Slurm receipt")
    with ZipFile(run_dir / "designs.zip") as archive:
        retained = archive.read("retained.csv")
        models = list(csv.DictReader(io.StringIO(retained.decode())))
        if (retained != (run_dir / "retained.csv").read_bytes() or len(models) > 8 or
                sorted("backbones/" + row["design_id"] + ".pdb" for row in models) !=
                sorted(name for name in archive.namelist() if name.startswith("backbones/"))):
            raise RuntimeError("the retained table and structure archive disagree")


# --- Slurm ---------------------------------------------------------------------------------------

def slurm_submit(rid):
    """sbatch handoff.sbatch from the run's inbox dir so SLURM_SUBMIT_DIR points at request.zip."""
    run_dir = INBOX / rid
    result = subprocess.run(["sbatch", "--parsable", str(SBATCH)], cwd=str(run_dir),
                            stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if result.returncode != 0:
        raise RuntimeError("sbatch failed: " + (result.stderr or result.stdout).strip())
    job = (result.stdout or "").strip().split(";")[0]  # --parsable emits "jobid[;cluster]"
    if not JOB_RE.fullmatch(job):
        raise RuntimeError("unexpected sbatch output: " + (result.stdout or "").strip())
    return job


def slurm_state(job):
    """Return a coarse Slurm state for the job: 'pending' | 'running' | 'done' | 'failed' | 'unknown'.

    sacct is authoritative once the job is known to the accounting db; squeue covers the brief window
    before it lands there. The raw state can carry a suffix (e.g. 'CANCELLED by 1234'); take the head.
    """
    acct = subprocess.run(["sacct", "-j", job, "-n", "-X", "-o", "State"],
                          stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    raw = ""
    if acct.returncode == 0:
        for line in (acct.stdout or "").splitlines():
            token = line.strip().split()[0] if line.strip() else ""
            if token:
                raw = token
                break
    if not raw:
        queue = subprocess.run(["squeue", "-j", job, "-h", "-o", "%T"],
                               stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        raw = (queue.stdout or "").strip().split("\n")[0].strip()
    if not raw:
        return "unknown", ""
    state = raw.upper()
    if state == "COMPLETED":
        return "done", state
    if state in SLURM_FAILED:
        return "failed", state
    if state in SLURM_RUNNING:
        return "running", state
    if state in SLURM_PENDING:
        return "pending", state
    return "unknown", state


# --- discovery -----------------------------------------------------------------------------------

def discover():
    """Classify every staged reference run as ('new'|'inflight'|'done', owner, rid)."""
    runs = {}
    for uri in gs_ls(REQUEST_BUCKET, "inbox/**/reference.zip"):
        # gs://<req>/inbox/<owner>/<rid>/reference.zip
        parts = uri.rstrip("/").split("/")
        if len(parts) < 4 or parts[-1] != "reference.zip":
            continue
        owner, rid = parts[-3], parts[-2]
        if OWNER_RE.fullmatch(owner) and RID_RE.fullmatch(rid):
            runs[rid] = owner
    # One listing of the result bucket tells which runs already have markers.
    markers = {}
    for uri in gs_ls(RESULT_BUCKET, "**"):
        parts = uri.rstrip("/").split("/")
        if len(parts) >= 2 and RID_RE.fullmatch(parts[-2]):
            markers.setdefault(parts[-2], set()).add(parts[-1])
    classified = []
    for rid, owner in runs.items():
        names = markers.get(rid, set())
        if "terminal.json" in names:
            phase = "done"
        elif "status.json" in names:
            phase = "inflight"
        else:
            phase = "new"
        classified.append((phase, owner, rid))
    return classified


# --- per-run actions -----------------------------------------------------------------------------

def submit(owner, rid):
    results = f"gs://{RESULT_BUCKET}/{rid}"
    now = time.time()
    if DRY_RUN:
        log(f"DRY-RUN would submit run {rid} (owner {owner})")
        return
    # Claim the run with a create-once status.json before sbatch so an overlapping tick cannot
    # double-submit. If the claim loses the race, another tick owns it; poll it next time.
    if not gs_write_json(results + "/status.json",
                         {"stage": "submitted", "job_id": None, "updated": now}, once=True):
        log(f"run {rid} already claimed; will poll")
        return
    run_dir = INBOX / rid
    run_dir.mkdir(parents=True, exist_ok=True)
    gs_cp_down(f"gs://{REQUEST_BUCKET}/inbox/{owner}/{rid}/reference.zip", run_dir / "request.zip")
    try:
        job = slurm_submit(rid)
    except Exception as exc:  # leave the claim's status.json so the dashboard shows the stall.
        log(f"run {rid} sbatch error: {exc}")
        gs_write_json(results + "/status.json",
                      {"stage": "submit error", "job_id": None, "updated": time.time(),
                       "reason": str(exc)})
        return
    gs_write_json(results + "/status.json",
                  {"stage": "submitted", "job_id": job, "updated": time.time()})
    log(f"run {rid} submitted as Slurm job {job}")


def poll(owner, rid):
    results = f"gs://{RESULT_BUCKET}/{rid}"
    status = gs_read_json(results + "/status.json") or {}
    job = status.get("job_id")
    run_dir = INBOX / rid

    # A crashed or lost submit (claimed, no job id): the compute job refuses a duplicate run id, so a
    # resubmit is safe. Retry it.
    if not job:
        log(f"run {rid} claimed without a job id; resubmitting")
        submit(owner, rid)
        return

    # The compute job writes terminal.json into the shared-home inbox on success: publish it.
    if (run_dir / "terminal.json").exists():
        try:
            publish(rid, job)
            log(f"run {rid} published (Slurm job {job})")
        except Exception as exc:
            log(f"run {rid} publish error: {exc}")
        return

    phase, raw = slurm_state(job)
    if phase == "failed":
        if DRY_RUN:
            log(f"DRY-RUN would fail run {rid}: Slurm {raw}")
            return
        mark_failed(rid, job, f"The Slurm job ended in state {raw}.")
        log(f"run {rid} failed (Slurm {raw})")
        return
    if phase == "done":
        # Slurm says COMPLETED but no terminal.json landed. Allow a grace window for shared-home
        # propagation, then treat it as a failure so the run does not hang forever.
        age = time.time() - float(status.get("updated", time.time()))
        if age > GRACE_SECONDS:
            if DRY_RUN:
                log(f"DRY-RUN would fail run {rid}: completed with no result bundle")
            else:
                mark_failed(rid, job, "The Slurm job completed but produced no result bundle.")
                log(f"run {rid} completed with no result after grace; failed")
        else:
            log(f"run {rid} completed; awaiting result bundle ({int(age)}s)")
        return
    stage = {"pending": "queued on cluster", "running": "running"}.get(phase, "submitted")
    if DRY_RUN:
        log(f"DRY-RUN run {rid} -> {stage} (Slurm {raw or phase})")
        return
    gs_write_json(results + "/status.json",
                  {"stage": stage, "job_id": job, "updated": time.time()})
    log(f"run {rid} -> {stage} (Slurm {raw or phase})")


def publish(rid, job):
    run_dir = INBOX / rid
    validate_designs(run_dir, rid, job)
    results = f"gs://{RESULT_BUCKET}/{rid}"
    if DRY_RUN:
        log(f"DRY-RUN would publish run {rid}")
        return
    gs_upload_once(run_dir / "retained.csv", results + "/retained.csv")
    gs_upload_once(run_dir / "designs.zip", results + "/designs.zip")
    gs_upload_once(run_dir / "terminal.json", results + "/terminal.json")  # terminal last = fully published


def mark_failed(rid, job, reason):
    results = f"gs://{RESULT_BUCKET}/{rid}"
    now = time.time()
    terminal = {"id": rid, "stage": "failed", "job_id": job, "updated": now,
                "completed": now, "reason": reason}
    # Preserve a prior identical failure on retry; a different terminal is a conflict to surface.
    existing = gs_read_json(results + "/terminal.json")
    if existing:
        if existing.get("stage") == "failed" and str(existing.get("job_id")) == str(job):
            return
        log(f"run {rid} already has a terminal result; leaving it")
        return
    gs_write_json(results + "/terminal.json", terminal, once=True)


# --- tick ----------------------------------------------------------------------------------------

def tick():
    if not REQUEST_BUCKET or not RESULT_BUCKET:
        sys.exit("Set EPISCAF_REQUEST_BUCKET and EPISCAF_RESULT_BUCKET.")
    runs = discover()
    counts = {"new": 0, "inflight": 0, "done": 0}
    for phase, owner, rid in runs:
        counts[phase] = counts.get(phase, 0) + 1
    log(f"tick: {counts['new']} new, {counts['inflight']} in-flight, {counts['done']} done"
        + (" (dry-run)" if DRY_RUN else ""))
    for phase, owner, rid in runs:
        try:
            if phase == "new":
                submit(owner, rid)
            elif phase == "inflight":
                poll(owner, rid)
        except Exception as exc:  # one bad run must not abort the tick.
            log(f"run {rid} tick error: {exc}")


def main(argv=None):
    global DRY_RUN
    parser = argparse.ArgumentParser(description="episcaf Gemini data-mover worker")
    parser.add_argument("command", choices=["tick"])
    parser.add_argument("--dry-run", action="store_true",
                        help="classify and print planned actions without submitting or writing")
    args = parser.parse_args(argv)
    DRY_RUN = args.dry_run
    os.umask(0o077)
    if args.command == "tick":
        tick()


if __name__ == "__main__":
    main()
