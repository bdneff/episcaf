"""Execute one staged 3HFM reference request inside Brandon's Slurm allocation."""
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import time
from zipfile import ZipFile

from episcaf_service.reference import matches_reference
from configs.paths import WORKSPACE
from episcaf_service import cloud_backbones as worker


def execute(inbox, workspace, job):
    rid = inbox.name
    if not re.fullmatch("[0-9a-f]{32}", rid) or not job.isdecimal():
        raise ValueError("Invalid run directory or Slurm ID.")
    # mkdir without exist_ok prevents a second submission from repeating inference.
    parent = workspace / "pilot" / "requests"
    parent.mkdir(parents=True, exist_ok=True)
    directory = parent / rid
    directory.mkdir(mode=0o700)
    with ZipFile(inbox / "request.zip") as archive:
        manifest = json.loads(archive.read("request.json"))
        upload = archive.read("upload")
        plan = json.loads(archive.read("plan.json"))
    if (manifest["id"] != rid or manifest["backend"] != "manual-gemini" or manifest["chain"] != "Y"
            or manifest["attempts"] != worker.ATTEMPTS or manifest["retained_limit"] != worker.LIMIT
            or hashlib.sha256(upload).hexdigest() != manifest["sha256"]):
        raise ValueError("The staged request failed the reference run checks.")
    route = manifest["route"]
    region, = manifest["regions"]
    worker.region_for(route, region)
    if plan != worker.reference_plan([region]):
        raise ValueError("The staged trial plan differs from this checkout.")
    if not matches_reference(upload.decode(), route, "Y"):
        raise ValueError("This job accepts the checked-in 3HFM lysozyme reference only.")
    shutil.copy2(inbox / "request.zip", directory / "request.zip")
    count = worker.generate(directory, route, region, source="gemini-backbone")
    for name in ("retained.csv", "designs.zip"):
        shutil.copyfile(directory / name, inbox / name)
    now = time.time()
    result = dict(id=rid, stage="completed", job_id=job, updated=now, completed=now,
                  reason=f"{count} of {worker.ATTEMPTS} backbones retained.")
    (directory / "terminal.json").write_text(json.dumps(result))
    temp = inbox / "terminal.json.tmp"
    temp.write_text(json.dumps(result))
    temp.replace(inbox / "terminal.json")


if __name__ == "__main__":
    os.umask(0o077)
    execute(Path(os.environ["SLURM_SUBMIT_DIR"]), WORKSPACE, os.environ["SLURM_JOB_ID"])
