"""Owner-scoped run review, submission, and result retrieval for the private pilot."""
import argparse
import csv
import hashlib
import io
import json
import os
import secrets
import sqlite3
import time
import uuid
from dataclasses import dataclass
from contextlib import contextmanager
from pathlib import Path
from typing import Literal, Protocol

from fastapi import APIRouter, Depends, HTTPException, Request, Response
from fastapi.responses import FileResponse
from pydantic import BaseModel, Field

from episcaf_service import staging

router = APIRouter(prefix="/pilot")
COOKIE = "episcaf_tester"
SUMMARY = "backend,id,filename,chain,stage,submitted,updated,job_id,estimated_start,completed,reason"


HOSTED = {"cloud-lysozyme", "manual-gemini"}


def mode():
    configured = os.getenv("EPISCAF_PILOT_MODE")
    if configured in HOSTED and os.getenv("K_SERVICE", "episcaf-pilot") == "episcaf-pilot":
        return configured
    if configured in {"simulated", "gemini-lysozyme"} and not os.getenv("K_SERVICE"):
        return configured
    return "staging"


def enabled():
    return mode() != "staging"


@contextmanager
def database():
    path = Path(os.environ["EPISCAF_PILOT_DB"])
    db = sqlite3.connect(path, timeout=15)
    db.row_factory = sqlite3.Row
    try:
        with db:
            yield db
    finally:
        db.close()


def initialize():
    # SQLite journals inherit the database permissions. Keep uploaded structures private on disk.
    path = Path(os.environ["EPISCAF_PILOT_DB"])
    path.touch(mode=0o600, exist_ok=True)
    with database() as db:
        db.executescript("""
            CREATE TABLE IF NOT EXISTS testers(owner TEXT PRIMARY KEY, digest TEXT UNIQUE);
            CREATE TABLE IF NOT EXISTS runs(
                id TEXT PRIMARY KEY, owner TEXT NOT NULL, filename TEXT, chain TEXT,
                structure TEXT, plan TEXT, stage TEXT, submitted REAL, updated REAL,
                job_id TEXT, estimated_start REAL, completed REAL, reason TEXT, result TEXT);
            CREATE INDEX IF NOT EXISTS by_owner ON runs(owner,submitted);
        """)
        db.execute("BEGIN IMMEDIATE")
        if "backend" not in {r["name"] for r in db.execute("PRAGMA table_info(runs)")}:
            db.execute("ALTER TABLE runs ADD COLUMN backend TEXT NOT NULL DEFAULT 'simulated'")


def issue(owner):
    token = secrets.token_urlsafe(32)
    with database() as db:
        db.execute("INSERT INTO testers VALUES(?,?) ON CONFLICT(owner) DO UPDATE SET digest=excluded.digest", (owner, digest(token)))
    return token


def digest(token):
    return hashlib.sha256(token.encode()).hexdigest()


def tester(token):
    if mode() in HOSTED:
        from episcaf_service import cloud
        return cloud.tester(token)
    if not enabled():
        raise HTTPException(404)
    with database() as db:
        row = db.execute("SELECT owner FROM testers WHERE digest=?", (digest(token),)).fetchone()
    if not row:
        raise HTTPException(401, "Sign in with your current tester token.")
    return row["owner"]


def owner(request: Request):
    return tester(request.cookies.get(COOKIE, ""))


@router.post("/session")
def signin(request: Request, response: Response):
    scheme, _, token = request.headers.get("Authorization", "").partition(" ")
    name = tester(token if scheme.lower() == "bearer" else "")
    response.set_cookie(COOKIE, token, httponly=True, secure=bool(os.getenv("K_SERVICE")) or request.url.scheme == "https",
                        samesite="strict", path="/pilot")
    return {"owner": name}


@router.delete("/session")
def signout(response: Response):
    response.delete_cookie(COOKIE, path="/pilot")
    return {"ok": True}


class RegionSelection(BaseModel):
    filename: str = Field(min_length=1, max_length=255)
    chain: str = Field(min_length=1, max_length=32)
    regions: list[str] = Field(min_length=1, max_length=1000)


class Selection(RegionSelection):
    route: Literal["growpatch", "antibody"] = "growpatch"
    structure: str = Field(min_length=1)


def prepare(s: RegionSelection, trials=8):
    try:
        return staging.candidate_plan(s.filename, s.chain, s.regions, trials)
    except ValueError as exc:
        raise HTTPException(422, str(exc)) from exc


@router.get("/reference.pdb")
def reference(route: Literal["antibody","growpatch"] = "growpatch"):
    if mode() in HOSTED and route == "antibody":
        from episcaf_service.lysozyme_job import REFERENCE
        return FileResponse(REFERENCE, media_type="chemical/x-pdb", filename="3HFM-complex.pdb")
    if mode() not in HOSTED | {"gemini-lysozyme"}:
        raise HTTPException(404)
    from episcaf_service.lysozyme_job import reference_pdb
    return Response(reference_pdb(), media_type="chemical/x-pdb",
                    headers={"Content-Disposition": 'attachment; filename="3HFM-lysozyme.pdb"'})


@router.post("/reviews")
def review(selection: Selection, name=Depends(owner)):
    backend = mode()
    if backend in HOSTED:
        from episcaf_service import cloud
        return cloud.review(selection, name)
    if backend == "gemini-lysozyme":
        from episcaf_service.lysozyme_job import reference_pdb, validate_regions
        selected_atoms = [line.rstrip() for line in selection.structure.splitlines()
                          if line.startswith("ATOM  ") and line[21:22] == "Y"]
        reference_atoms = [line.rstrip() for line in reference_pdb().splitlines() if line.startswith("ATOM  ")]
        if selection.chain != "Y" or selected_atoms != reference_atoms:
            backend = "gemini-reviewed"
        else:
            try:
                validate_regions(selection.regions)
            except ValueError as exc:
                raise HTTPException(422, str(exc)) from exc
    plan = prepare(selection)
    if backend == "gemini-lysozyme":
        from scripts.stage03_mpnn_submit import SEQS_PER_STRUCT
        plan["sequence_predictions"] = plan["total"] * SEQS_PER_STRUCT
    rid = uuid.uuid4().hex
    with database() as db:
        db.execute("""INSERT INTO runs(id,owner,filename,chain,structure,plan,stage,backend)
                      VALUES(?,?,?,?,?,?,'review',?)""",
                   (rid, name, selection.filename, selection.chain, selection.structure, json.dumps(plan), backend))
    return dict(id=rid, backend=backend, **plan)


@dataclass(frozen=True)
class JobStatus:
    stage: str
    estimated_start: float | None = None
    completed: float | None = None
    reason: str | None = None


class Scheduler(Protocol):
    poll_seconds: int

    def submit(self, run_id: str) -> str: ...
    def status(self, job_id: str, submitted: float, now: float) -> JobStatus: ...


class SimulatedScheduler:
    poll_seconds = 1

    def submit(self, run_id):
        outcome = "fail" if os.getenv("EPISCAF_SIM_FAIL") == "1" else "ok"
        return f"sim-{outcome}-{run_id}"

    def status(self, job_id, submitted, now):
        age = now - submitted
        if age < 2:
            return JobStatus("queued", estimated_start=submitted + 2)
        if age < 6:
            return JobStatus("running")
        if job_id.startswith("sim-fail-"):
            return JobStatus("failed", reason="The simulated worker stopped before producing results (test failure).")
        return JobStatus("completed", completed=submitted + 6)


scheduler: Scheduler = SimulatedScheduler()


def simulated_csv(plan):
    buf = io.StringIO()
    writer = csv.writer(buf)
    writer.writerow(["source", "epitope", "backbone_trials"])
    for region in plan["regions"]:
        writer.writerow(["simulated", region, plan["trials_per_candidate"]])
    return buf.getvalue()


def refresh():
    with database() as db:
        active = db.execute("""SELECT id,job_id,submitted,stage,estimated_start
                               FROM runs WHERE backend='simulated' AND stage IN ('queued','running')""").fetchall()
        for row in active:
            now = time.time()
            status = scheduler.status(row["job_id"], row["submitted"], now)
            if status.stage == row["stage"] and status.estimated_start == row["estimated_start"]:
                continue
            result = None
            if status.stage == "completed":
                plan = db.execute("SELECT plan FROM runs WHERE id=?", (row["id"],)).fetchone()
                result = simulated_csv(json.loads(plan[0]))
            db.execute("""UPDATE runs SET stage=?,reason=?,updated=?,estimated_start=?,completed=?,result=?
                          WHERE id=? AND stage=?""",
                       (status.stage, status.reason, now, status.estimated_start, status.completed,
                        result, row["id"], row["stage"]))


def owned(db, rid, name, columns=SUMMARY, include_review=False):
    row = db.execute(f"SELECT {columns} FROM runs WHERE id=? AND owner=?" +
                     ("" if include_review else " AND submitted IS NOT NULL"), (rid, name)).fetchone()
    if not row:
        raise HTTPException(404, "Run unavailable.")
    return row


@router.post("/runs/{rid}")
def launch(rid: str, name=Depends(owner)):
    if mode() in HOSTED:
        from episcaf_service import cloud
        return cloud.launch(name, rid)
    with database() as db:
        db.execute("BEGIN IMMEDIATE")
        row = owned(db, rid, name, include_review=True)
        if row["backend"] != mode() and not (row["backend"] == "gemini-reviewed" and mode() == "gemini-lysozyme"):
            raise HTTPException(409, "Prepare a new review on this backend.")
        if row["stage"] == "review" and mode() == "gemini-lysozyme":
            now = time.time()
            stage = "awaiting lab review" if row["backend"] == "gemini-reviewed" else "requested"
            db.execute("UPDATE runs SET stage=?,submitted=?,updated=? WHERE id=?", (stage, now, now, rid))
        elif row["stage"] == "review":
            now = time.time()
            db.execute("""UPDATE runs SET stage='queued',submitted=?,updated=?,job_id=? WHERE id=?""",
                       (now, now, scheduler.submit(rid), rid))
        return dict(owned(db, rid, name))


@router.get("/runs")
def runs(name=Depends(owner)):
    if mode() in HOSTED:
        from episcaf_service import cloud
        return cloud.runs(name)
    with database() as db:
        return [dict(r) for r in db.execute(f"SELECT {SUMMARY} FROM runs WHERE owner=? AND "
                                            "submitted IS NOT NULL ORDER BY submitted DESC", (name,))]


@router.get("/runs/{rid}")
def status(rid: str, name=Depends(owner)):
    if mode() in HOSTED:
        from episcaf_service import cloud
        return cloud.status(name, cloud.owned(name, rid))
    with database() as db:
        return dict(owned(db, rid, name))


@router.get("/runs/{rid}/results")
def results(rid: str, download: bool = False, name=Depends(owner)):
    if mode() in HOSTED:
        from episcaf_service import cloud
        return cloud.result(name, rid, "retained.csv", download)
    with database() as db:
        row = owned(db, rid, name, "backend,stage,result,reason")
    if row["stage"] != "completed":
        raise HTTPException(409, row["reason"] or "Results are pending.")
    backend = row["backend"]
    headers = {"Content-Disposition": f'attachment; filename="{backend}-{rid}.csv"'} if download else {}
    return Response(row["result"], media_type="text/csv", headers=headers)


@router.get("/runs/{rid}/request.zip")
def request_files(rid: str, name=Depends(owner)):
    if mode() != "manual-gemini":
        raise HTTPException(404)
    from episcaf_service import handoff
    return handoff.download(name, rid)


@router.get("/runs/{rid}/designs.zip")
def designs(rid: str, name=Depends(owner)):
    if mode() in HOSTED:
        from episcaf_service import cloud
        return cloud.result(name, rid, "designs.zip", True)
    from episcaf_service.gemini import run_dir
    with database() as db:
        row = owned(db, rid, name, "backend,stage")
    if row["backend"] not in {"gemini-lysozyme","gemini-reviewed"} or row["stage"] != "completed":
        raise HTTPException(404)
    path = run_dir(rid) / "designs.zip"
    if not path.is_file():
        raise HTTPException(409, "The scored design archive is unavailable.")
    return FileResponse(path, media_type="application/zip", filename=f"designs-{rid}.zip")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Issue (or rotate) a tester token, or revoke it.")
    parser.add_argument("action", choices=["issue", "revoke"])
    parser.add_argument("owner")
    args = parser.parse_args()
    initialize()
    if args.action == "issue":
        print(issue(args.owner))
    else:
        with database() as conn:
            conn.execute("DELETE FROM testers WHERE owner=?", (args.owner,))
