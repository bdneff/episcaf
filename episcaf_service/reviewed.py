"""Operator tools for inspecting private requests and returning existing lab results."""
import argparse
import csv
import hashlib
import io
import json
import math
import re
import time
from pathlib import Path
from zipfile import ZIP_DEFLATED, ZipFile

from episcaf_service import pilot
from episcaf_service.gemini import run_dir


def pending(db, rid):
    row = db.execute("""SELECT * FROM runs WHERE id=? AND backend='gemini-reviewed'
                        AND stage='awaiting lab review'""", (rid,)).fetchone()
    if row is None:
        raise ValueError("Select a pending lab review.")
    return row


def export(rid):
    directory = run_dir(rid)
    with pilot.database() as db:
        row = pending(db, rid)
        plan = json.loads(row["plan"])
        data = row["structure"].encode()
        directory.mkdir(mode=0o700, exist_ok=True)
        suffix = ".cif" if row["filename"].lower().endswith((".cif", ".mmcif")) else ".pdb"
        (directory / ("upload" + suffix)).write_bytes(data)
        manifest = dict(run_id=rid, owner=row["owner"], filename=row["filename"], chain=row["chain"],
                        regions=plan["regions"], total=plan["total"], sha256=hashlib.sha256(data).hexdigest())
        (directory / "request.json").write_text(json.dumps(manifest, indent=2))
    return directory


def publish(rid, ranked, models, job_id):
    # A submitted run owns exactly one archive. Serialize concurrent operator
    # publications with the terminal DB update and publish the file atomically.
    directory = run_dir(rid)
    with pilot.database() as db:
        db.execute("BEGIN IMMEDIATE")
        row = pending(db, rid)
        regions = set(json.loads(row["plan"])["regions"])
        with ranked.open(newline="") as fh:
            reader = csv.DictReader(fh)
            columns = reader.fieldnames or []
            required = {"design_id", "sequence", "composite", "epitope", "run_id"}
            if not required.issubset(columns) or not job_id.isdecimal():
                raise ValueError("Supply a numeric Slurm ID and a CSV with " + ", ".join(sorted(required)))
            rows = list(reader)
        if not rows:
            raise ValueError("The lab results CSV is empty.")
        files = {}
        for record in rows:
            design = record["design_id"]
            if not re.fullmatch(r"[a-zA-Z0-9_-]+", design) or design in files:
                raise ValueError("Each design needs a unique model filename.")
            if (record["run_id"] != rid or record["epitope"] not in regions or
                    not re.fullmatch("[ACDEFGHIKLMNPQRSTVWY]+", record["sequence"]) or
                    not math.isfinite(float(record["composite"]))):
                raise ValueError("Check the run ID, selected epitope, sequence and composite for each design.")
            model = models / (design + ".cif")
            if not model.is_file() or not model.stat().st_size:
                raise ValueError(f"Missing predicted structure for {design}.")
            files[design] = model
            record["source"] = "lab-reviewed"
        text = io.StringIO()
        writer = csv.DictWriter(text, fieldnames=columns if "source" in columns else [*columns, "source"])
        writer.writeheader()
        writer.writerows(rows)
        directory.mkdir(mode=0o700, exist_ok=True)
        temp = directory / "designs.zip.tmp"
        with ZipFile(temp, "w", compression=ZIP_DEFLATED) as archive:
            archive.writestr("ranked.csv", text.getvalue())
            for design, model in files.items():
                archive.write(model, f"models/{design}.cif")
        temp.replace(directory / "designs.zip")
        now = time.time()
        db.execute("""UPDATE runs SET stage='completed',completed=?,updated=?,job_id=?,result=?
                      WHERE id=?""", (now, now, job_id, text.getvalue(), rid))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    commands.add_parser("list")
    commands.add_parser("export").add_argument("run_id")
    p = commands.add_parser("publish")
    p.add_argument("run_id")
    p.add_argument("--ranked", type=Path, required=True)
    p.add_argument("--models", type=Path, required=True)
    p.add_argument("--job-id", required=True)
    p = commands.add_parser("decline")
    p.add_argument("run_id")
    p.add_argument("reason")
    args = parser.parse_args()
    if pilot.mode() != "gemini-lysozyme":
        parser.error("Use the private Gemini pilot environment.")
    if args.command == "list":
        with pilot.database() as db:
            for row in db.execute("""SELECT id, owner, filename, chain FROM runs
                                     WHERE stage='awaiting lab review' ORDER BY submitted"""):
                print(json.dumps(dict(row)))
    elif args.command == "export":
        print(export(args.run_id))
    elif args.command == "publish":
        publish(args.run_id, args.ranked, args.models, args.job_id)
    else:
        with pilot.database() as db:
            db.execute("BEGIN IMMEDIATE")
            pending(db, args.run_id)
            db.execute("UPDATE runs SET stage='failed',reason=?,updated=? WHERE id=?",
                       (args.reason, time.time(), args.run_id))


if __name__ == "__main__":
    main()
