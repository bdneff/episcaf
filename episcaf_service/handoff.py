"""Stage private request bundles for operator submission on Gemini."""
import hashlib
import io
import json
from pathlib import Path
import time
from zipfile import ZIP_DEFLATED, ZipFile

from fastapi import HTTPException, Response
from google.api_core.exceptions import NotFound, PreconditionFailed

from episcaf_service import cloud

READY = "awaiting Brandon’s submission"


def bundle_blob(owner, record):
    name = {"manual-gemini": "reference.zip", "manual-review": "review.zip"}.get(record["backend"])
    if name is None:
        raise HTTPException(404, "No staged request for this run.")
    return cloud.store().blob(f"inbox/{owner}/{record['id']}/{name}")


def download(owner, rid):
    record = cloud.owned(owner, rid)
    try:
        data = bundle_blob(owner, record).download_as_bytes()
    except NotFound as exc:
        raise HTTPException(409, "The request files are not available yet.") from exc
    return Response(data, media_type="application/zip",
                    headers={"Content-Disposition": f'attachment; filename="request-{rid}.zip"'})


def bundle(record, owner, upload):
    ref = record["backend"] == "manual-gemini"
    manifest = {key: record[key] for key in ("id", "filename", "chain", "route", "regions", "backend")}
    manifest.update(owner=owner, sha256=hashlib.sha256(upload).hexdigest(),
                    attempts=cloud.worker.ATTEMPTS if ref else None,
                    retained_limit=cloud.worker.LIMIT if ref else None)
    stream = io.BytesIO()
    with ZipFile(stream, "w", compression=ZIP_DEFLATED) as archive:
        archive.writestr("request.json", json.dumps(manifest, indent=2))
        archive.writestr("upload", upload)
        if ref:
            archive.writestr("plan.json", json.dumps(record["plan"], indent=2))
            archive.write(Path(__file__).with_name("handoff.sbatch"), "handoff.sbatch")
    return stream.getvalue()


def launch(owner, r, blob):
    if r["stage"] == "review":
        now = time.time()
        r = dict(r, stage="preparing handoff", submitted=now, updated=now)
        try:
            cloud.write(blob, r, if_generation_match=blob.generation)
        except PreconditionFailed:
            r = cloud.owned(owner, r["id"])
    return status(owner, r)


def status(owner, r):
    rid = r["id"]
    if r["stage"] == "preparing handoff":
        blob = cloud.record_blob(owner, rid)
        r = cloud.read(blob)
        if r["stage"] == "preparing handoff":
            reference = r["backend"] == "manual-gemini"
            upload = cloud.store().blob(f"uploads/{owner}/{rid}").download_as_bytes()
            out = bundle_blob(owner, r)
            try:
                out.upload_from_string(bundle(r, owner, upload), content_type="application/zip",
                                       if_generation_match=0)
            except PreconditionFailed:
                pass  # Another request already staged the immutable bundle.
            r.update(stage=READY if reference else "awaiting lab review", updated=time.time())
            try:
                cloud.write(blob, r, if_generation_match=blob.generation)
            except PreconditionFailed:
                r = cloud.owned(owner, rid)
    if r["backend"] == "manual-review":
        return r
    # Only the operator can write result objects. Completion is published last.
    for marker in ("terminal", "receipt"):
        blob = cloud.outputs().blob(f"{rid}/{marker}.json")
        update = cloud.read(blob)
        if update:
            if marker == "receipt" and blob.time_created:
                update["updated"] = blob.time_created.timestamp()
            return dict(r, **{k: update[k] for k in ("stage","job_id","updated","reason","completed")
                              if k in update})
    return r
