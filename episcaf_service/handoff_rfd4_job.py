"""Execute one staged RFD4-scaffold request inside Brandon's Slurm allocation.

Parallels handoff_job.py (the RFd3 reference path), but runs the RFD4 campaign engine on a general
antigen + epitope instead of the locked 3HFM reference. Reads request.zip from the worker inbox
(SLURM_SUBMIT_DIR), runs episcaf_service.rfd4.campaign, and copies the ranked results back as
retained.csv + designs.zip + terminal.json — the same filenames the data-mover worker publishes.

Tool paths come from the sbatch env (handoff_rfd4.sbatch): EPISCAF_COLABFOLD, EPISCAF_AF2_DATA,
EPISCAF_RFD4_CMD. Weights/tools live on tgen_labs (home is quota-bound).
"""
import hashlib
import json
import os
import re
import shutil
import shlex
import time
from pathlib import Path
from zipfile import ZipFile

from configs.paths import WORKSPACE
from episcaf_service.rfd4 import campaign, fold
from episcaf_service.rfd4.config import DEFAULT

WORKFLOW = "rfd4-scaffold-v1"


def execute(inbox, workspace, job):
    rid = inbox.name
    if not re.fullmatch("[0-9a-f]{32}", rid) or not job.isdecimal():
        raise ValueError("Invalid run directory or Slurm ID.")
    # mkdir without exist_ok: a second submission for the same rid refuses to repeat inference.
    parent = workspace / "pilot" / "rfd4"
    parent.mkdir(parents=True, exist_ok=True)
    directory = parent / rid
    directory.mkdir(mode=0o700)

    with ZipFile(inbox / "request.zip") as archive:
        manifest = json.loads(archive.read("request.json"))
        upload = archive.read("upload")
    if (manifest["id"] != rid or manifest.get("workflow") != WORKFLOW
            or hashlib.sha256(upload).hexdigest() != manifest["sha256"]):
        raise ValueError("The staged request failed the RFD4 run checks.")
    regions = manifest["regions"]
    if not regions:
        raise ValueError("No epitope regions in the staged request.")
    epitope = ";".join(regions)
    antibody_chains = tuple(manifest.get("antibody_chains") or ())

    antigen = directory / "antigen.pdb"
    antigen.write_bytes(upload)

    predictor = fold.ColabFoldAF2(data_dir=os.environ["EPISCAF_AF2_DATA"],
                                  colabfold_bin=os.environ["EPISCAF_COLABFOLD"], cfg=DEFAULT.predictor)
    results = campaign.run_campaign(
        str(antigen), epitope, str(directory), predictor=predictor, native_path=str(antigen),
        antibody_chains=antibody_chains,
        rfd4_cmd=tuple(shlex.split(os.environ.get("EPISCAF_RFD4_CMD", "rfd4-generate"))),
        num_replicates=manifest.get("num_replicates"), cfg=DEFAULT)

    # publish under the names the data-mover worker already handles
    shutil.copyfile(directory / "final_results.csv", inbox / "retained.csv")
    shutil.copyfile(directory / "designs.zip", inbox / "designs.zip")

    scored = sum(1 for r in results if r.composite is not None)
    now = time.time()
    result = dict(id=rid, stage="completed", job_id=job, updated=now, completed=now,
                  reason=f"{scored} of {len(results)} designs scored.")
    (directory / "terminal.json").write_text(json.dumps(result))
    temp = inbox / "terminal.json.tmp"
    temp.write_text(json.dumps(result))
    temp.replace(inbox / "terminal.json")


if __name__ == "__main__":
    os.umask(0o077)
    execute(Path(os.environ["SLURM_SUBMIT_DIR"]), WORKSPACE, os.environ["SLURM_JOB_ID"])
