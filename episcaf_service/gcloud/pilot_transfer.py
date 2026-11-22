"""Transfer a staged request and publish its Slurm receipt or results from Brandon's laptop."""
import argparse
import csv
import io
import json
import os
from pathlib import Path
import re
import subprocess
import tempfile
import time
from zipfile import ZipFile


def call(*args, **kwargs):
    return subprocess.run(args, check=True, **kwargs)


def text_output(*args):
    return call(*args, stdout=subprocess.PIPE, text=True).stdout.strip()


def upload_once(path, uri):
    # A repeated publication must match the existing immutable object.
    if subprocess.run(["gcloud","storage","cp","--if-generation-match=0",str(path),uri]).returncode:
        if call("gcloud","storage","cat",uri,stdout=subprocess.PIPE).stdout != path.read_bytes():
            raise ValueError("A different result already exists at " + uri)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=["fetch","review","receipt","publish","fail"])
    parser.add_argument("owner")
    parser.add_argument("run")
    parser.add_argument("job", nargs="?")
    parser.add_argument("--reason", help="Failure explanation visible to the tester")
    args = parser.parse_args()
    rid, job = args.run, args.job
    if not re.fullmatch("[a-zA-Z0-9_-]{1,48}", args.owner) or not re.fullmatch("[a-f0-9]{32}",rid):
        parser.error("Invalid owner or run ID.")
    if args.action in ("receipt","publish","fail") and not re.fullmatch("[0-9]+",job or ""):
        parser.error("Supply the numeric Slurm job ID.")
    if args.action == "fail" and not args.reason:
        parser.error("Supply --reason with the stage and Slurm exit status.")
    project = text_output("gcloud","config","get-value","project")
    number = text_output("gcloud","projects","describe",project,"--format=value(projectNumber)")
    inbox = ".episcaf-pilot/inbox/" + rid
    object_dir = f"gs://episcaf-requests-{number}/inbox/{args.owner}/{rid}"
    results = f"gs://episcaf-results-{number}/{rid}"
    os.umask(0o077)
    with tempfile.TemporaryDirectory() as temp:
        directory = Path(temp)
        if args.action in ("fetch","review"):
            name = "review.zip" if args.action == "review" else "reference.zip"
            path = directory / "request.zip"
            call("gcloud","storage","cp",object_dir+"/"+name,str(path))
            # The remote shell only sees the validated hex ID. Publish the file atomically.
            with path.open("rb") as source:
                call("ssh","gemini",f'configs={inbox}; umask 077; mkdir -p "$configs" && cd "$configs" && '
                     'tmp=$(mktemp .request.XXXXXX) && cat > "$tmp" && mv "$tmp" request.zip',stdin=source)
            print(f"Gemini inbox: ~/{inbox}/request.zip")
            if args.action == "fetch":
                print(f"On Gemini: cd ~/{inbox} && sbatch --parsable ~/episcaf-advisor/episcaf_service/handoff.sbatch")
            return
        # Check that the operator selected a submitted reference run owned by this tester.
        text_output("gcloud","storage","objects","describe",object_dir+"/reference.zip","--format=value(name)")
        receipt = directory / "receipt.json"
        receipt.write_text(json.dumps(dict(stage="queued", job_id=job)))
        if args.action == "receipt":
            upload_once(receipt, results+"/receipt.json")
            return
        terminal = directory / "terminal.json"
        if args.action == "publish":
            for name in ("terminal.json","retained.csv","designs.zip"):
                call("scp",f"gemini:{inbox}/{name}",str(directory/name))
            record = json.loads(terminal.read_text())
            if record["id"] != rid or record["job_id"] != job or record["stage"] != "completed":
                raise ValueError("The result does not match this run and Slurm receipt.")
            with ZipFile(directory/"designs.zip") as archive:
                retained = archive.read("retained.csv")
                models = list(csv.DictReader(io.StringIO(retained.decode())))
                if (retained != (directory/"retained.csv").read_bytes() or len(models)>8 or
                    sorted("backbones/"+r["design_id"]+".pdb" for r in models) !=
                    sorted(n for n in archive.namelist() if n.startswith("backbones/"))):
                    raise ValueError("The retained table and structure archive disagree.")
            upload_once(receipt, results+"/receipt.json")
            for name in ("retained.csv","designs.zip"):
                upload_once(directory/name, results+"/"+name)
        else:
            # Preserve the first failure timestamp when the same command is retried.
            existing = subprocess.run(["gcloud","storage","cat",results+"/terminal.json"],capture_output=True)
            if existing.returncode == 0:
                old = json.loads(existing.stdout)
                if old["stage"] == "failed" and old["job_id"] == job and old["reason"] == args.reason:
                    return
                raise ValueError("The run already has a terminal result.")
            upload_once(receipt, results+"/receipt.json")
            now = time.time()
            terminal.write_text(json.dumps(dict(stage="failed",job_id=job,updated=now,
                                                completed=now,reason=args.reason)))
        upload_once(terminal, results+"/terminal.json")


if __name__ == "__main__":
    main()
