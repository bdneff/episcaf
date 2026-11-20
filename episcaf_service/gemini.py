"""Operator-started dispatcher for the private 3HFM reference pilot."""
import argparse
import csv
from datetime import datetime
import json
import logging
import os
from pathlib import Path
import re
import socket
import subprocess
import time

from episcaf_service import pilot
from configs.paths import WORKSPACE
from episcaf_service.lysozyme_job import ROOT

BACKEND = "gemini-lysozyme"
POLL_SECONDS = 1800


def run_root():
    configured = Path(os.environ["EPISCAF_GEMINI_RUNS"])
    root = configured.resolve()
    legacy = (Path(os.getenv("EPISCAF_WS", WORKSPACE)) / "runs").resolve()
    if not configured.is_absolute() or any(root.is_relative_to(p) or p.is_relative_to(root)
                                            for p in (ROOT, legacy)):
        raise ValueError("Choose a separate absolute private directory outside the checkout and legacy runs.")
    if not root.is_dir() or root.stat().st_mode & 0o077:
        raise ValueError("The shared pilot run directory must be private to the operator (mode 0700).")
    return root


def run_dir(rid):
    if not re.fullmatch(r"[0-9a-f]{32}", rid):
        raise ValueError("Invalid run ID")
    return run_root() / rid


def submit(row):
    rid = row["id"]
    directory = run_dir(rid)
    # Commit the claim before contacting Slurm. A second dispatcher or a web
    # retry cannot submit the same run. Uncertain submissions are reconciled by name.
    with pilot.database() as db:
        claimed = db.execute("UPDATE runs SET stage='submitting',updated=? WHERE id=? AND stage='requested'",
                             (time.time(), rid)).rowcount
    if not claimed:
        return
    try:
        directory.mkdir(mode=0o700)
        (directory / "selection.json").write_text(json.dumps({"regions": json.loads(row["plan"])["regions"]}))
        args = ["sbatch", "--parsable", "--job-name", "ep-" + rid, "--chdir", str(ROOT),
                "--output", str(directory / "slurm-%j.log"), str(ROOT / "episcaf_service/lysozyme_job.sbatch"),
                str(directory)]
        result = subprocess.run(args, capture_output=True, text=True, timeout=90)
    except subprocess.TimeoutExpired:
        return  # accounting may already have accepted this job
    except OSError as exc:
        finish(rid, "failed", f"Could not stage or submit the job: {exc.strerror}")
        return
    if result.returncode:
        finish(rid, "failed", "Slurm rejected submission: " + result.stderr.strip()[:400])
        return
    jid = result.stdout.strip().partition(";")[0]
    if jid.isdecimal():
        record_job(rid, jid)


def record_job(rid, jid):
    with pilot.database() as db:
        db.execute("UPDATE runs SET stage='queued',reason=NULL,job_id=?,updated=? WHERE id=? AND stage='submitting'",
                   (jid, time.time(), rid))


def finish(rid, stage, reason=None, result=None, completed=None):
    with pilot.database() as db:
        db.execute("""UPDATE runs SET stage=?,reason=?,result=?,updated=?,completed=?,estimated_start=NULL
                      WHERE id=? AND stage NOT IN ('completed','failed')""",
                   (stage, reason, result, time.time(), completed, rid))


def timestamp(value):
    try:
        return datetime.fromisoformat(value).timestamp()
    except ValueError:
        return None


def slurm(*args):
    output = subprocess.run(args, text=True, capture_output=True, check=True, timeout=45).stdout
    return [line.split("|") for line in output.splitlines() if line]


def reconcile(row):
    rid, jid = row["id"], row["job_id"]
    directory = run_dir(rid)
    if not jid:
        date = datetime.fromtimestamp(row["submitted"] - 86400).strftime("%Y-%m-%d")
        records = slurm("sacct", "-X", "-n", "-P", "-S", date, "--name", "ep-"+rid,
                        "-o", "JobIDRaw,JobName%64")
        matches = {i for i, name in records if name == "ep-"+rid and i.isdecimal()}
        if len(matches) == 1:
            record_job(rid, matches.pop())
        elif time.time() - row["submitted"] >= POLL_SECONDS:
            with pilot.database() as db:
                db.execute("""UPDATE runs SET reason=?,updated=? WHERE id=? AND stage='submitting'""",
                           ("The lab is checking the Slurm submission receipt before retrying.",
                            time.time(), rid))
        return
    records = slurm("sacct", "-X", "-n", "-P", "-j", jid, "-o", "JobIDRaw,State,ExitCode,End")
    for number, state, code, end in records:
        if number != jid:
            continue
        state = state.split()[0].rstrip("+")
        if state == "COMPLETED" and code == "0:0":
            ranked = directory / "05_analysis/ranked.csv"
            if not (directory / "designs.zip").is_file() or not ranked.is_file():
                finish(rid, "failed", "Job ended before the scored designs were published.")
                return
            text = ranked.read_text()
            if not list(csv.DictReader(text.splitlines())):
                finish(rid, "failed", "The ranked metrics table is empty.")
                return
            finish(rid, "completed", result=text, completed=timestamp(end))
            return
        if state in {"CANCELLED","FAILED","TIMEOUT","OUT_OF_MEMORY","NODE_FAIL","BOOT_FAIL",
                     "PREEMPTED","DEADLINE","REVOKED","COMPLETED"}:
            stage_path = directory / "stage.txt"
            step = stage_path.read_text().strip() if stage_path.exists() else "setup"
            finish(rid, "failed", f"Job {jid} ended during {step} ({state}, exit {code}).")
            return
    # A job may disappear between the accounting and queue queries. Query by
    # user and exact name, then match its ID, avoiding an invalid-job-id error.
    entries = slurm("squeue", "--me", "--name", "ep-"+rid, "-h", "-o", "%i|%T|%S")
    for number, state, start in entries:
        if number == jid:
            detail = directory / "stage.txt"
            stage = (detail.read_text().strip() if state == "RUNNING" and detail.exists() else state.lower())
            estimate = timestamp(start) if state == "PENDING" else None
            with pilot.database() as db:
                db.execute("""UPDATE runs SET stage=?,updated=?,estimated_start=? WHERE id=?
                              AND stage NOT IN ('failed','completed')""", (stage, time.time(), estimate, rid))


def dispatch(poll):
    with pilot.database() as db:
        rows = db.execute("""SELECT id,stage,job_id,submitted,plan FROM runs WHERE backend=?
                             AND submitted IS NOT NULL AND stage NOT IN ('completed','failed')
                             ORDER BY submitted""", (BACKEND,)).fetchall()
    for row in rows:
        try:
            if row["stage"] == "requested":
                submit(row)
            elif poll:
                reconcile(row)
        except (OSError, subprocess.SubprocessError, ValueError):
            logging.exception("Dispatcher could not refresh run %s", row["id"])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--once", action="store_true", help="dispatch and refresh once, then exit")
    args = parser.parse_args()
    if pilot.mode() != BACKEND or os.environ.get("EPISCAF_GEMINI_HOST") not in {socket.gethostname(),socket.getfqdn()}:
        parser.error("Set the private live mode and EPISCAF_GEMINI_HOST to the approved local submission hostname.")
    run_root()
    pilot.initialize()
    last_poll = -POLL_SECONDS
    while True:
        poll = time.monotonic() - last_poll >= POLL_SECONDS
        dispatch(poll)
        if args.once:
            break
        if poll:
            last_poll = time.monotonic()
        time.sleep(10)
