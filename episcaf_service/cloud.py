"""Private Cloud Run requests, on-demand Google Batch jobs and private bucket results."""
from datetime import datetime
from functools import lru_cache
import hmac
import json
import os
import re
import time
import uuid

from fastapi import HTTPException, Response

from episcaf_service import cloud_backbones as worker
from episcaf_service import pilot
from episcaf_service.reference import matches_reference

MODE = "cloud-lysozyme"
PENDING = "cloud-reviewed"


@lru_cache(maxsize=1)
def client():
    from google.cloud import storage
    return storage.Client()


def store():
    return client().bucket(os.environ["EPISCAF_REQUEST_BUCKET"])


def outputs():
    return client().bucket(os.environ["EPISCAF_RESULT_BUCKET"])


def write(blob, value, **kwargs):
    blob.upload_from_string(json.dumps(value), content_type="application/json", **kwargs)


def read(blob):
    from google.api_core.exceptions import NotFound
    try:
        return json.loads(blob.download_as_bytes())
    except NotFound:
        return None


def tester(token):
    owner = token.partition(".")[0]
    if re.fullmatch("[a-zA-Z0-9_-]{1,48}", owner):
        record = read(store().blob("testers/"+owner))
        if record and hmac.compare_digest(record.get("digest") or "", pilot.digest(token)):
            return owner
    raise HTTPException(401, "Sign in with your current tester token.")


def record_blob(owner, rid):
    if not re.fullmatch("[0-9a-f]{32}", rid):
        raise HTTPException(404, "Run unavailable.")
    return store().blob(f"runs/{owner}/{rid}")


def review(s, owner):
    route = s.route
    enabled = matches_reference(s.structure, route, s.chain)
    regions = s.regions
    if enabled:
        if len(regions) != 1:
            raise HTTPException(422, "Choose one region per reference run.")
        try:
            worker.region_for(route, regions[0])
        except ValueError as exc:
            raise HTTPException(422, str(exc)) from exc
    plan = pilot.prepare(s, trials=worker.ATTEMPTS if enabled else 8)
    if enabled:
        plan["retained_limit"] = worker.LIMIT
    rid = uuid.uuid4().hex
    backend = pilot.mode() if enabled else PENDING
    if not enabled and pilot.mode() == "manual-gemini":
        backend = "manual-review"
    record = dict(id=rid, filename=s.filename, chain=s.chain, route=route, regions=regions, backend=backend,
                  stage="review", submitted=None, updated=None, job_id=None, reason=None,
                  completed=None, estimated_start=None)
    if backend == "manual-gemini":
        record["plan"] = worker.reference_plan(regions)
    # Keep the uploaded bytes private. The worker can only read its checked-in reference.
    store().blob(f"uploads/{owner}/{rid}").upload_from_string(s.structure, if_generation_match=0)
    write(record_blob(owner,rid), record, if_generation_match=0)
    return dict(**plan, id=rid, backend=backend)


def batch_parent():
    return ("https://batch.googleapis.com/v1/projects/"+os.environ["EPISCAF_PROJECT"]+
            "/locations/"+os.environ["EPISCAF_REGION"]+"/jobs")


def batch_request(method, suffix="", **kwargs):
    import google.auth
    from google.auth.transport.requests import AuthorizedSession
    from requests.exceptions import RequestException
    credentials, _ = google.auth.default(scopes=["https://www.googleapis.com/auth/cloud-platform"])
    try:
        with AuthorizedSession(credentials) as session:
            res = session.request(method, batch_parent()+suffix, timeout=30, **kwargs)
    except RequestException as exc:
        raise HTTPException(503, "Google Batch could not be reached. Retry this run.") from exc
    if res.status_code not in (400,403,404,409) and not res.ok:
        raise HTTPException(503, f"Google Batch returned HTTP {res.status_code}. Check the pilot service permissions and GPU quota.")
    return res


def job_spec(r):
    spec = {"runnables":[{"container":{"imageUri":os.environ["EPISCAF_GPU_IMAGE"],
                "commands":["--route",r["route"],"--region",r["regions"][0],"--run",r["id"],
                            "--bucket",os.environ["EPISCAF_RESULT_BUCKET"]], "options":"--gpus all"}}],
            "maxRunDuration":"3600s","maxRetryCount":0,
            "computeResource":{"cpuMilli":4000,"memoryMib":14000}}
    instance = {"installGpuDrivers":True,"policy":{"machineType":"g2-standard-8",
                "bootDisk":{"sizeGb":150,"type":"pd-balanced"}}}
    return {"taskGroups":[{"taskCount":1, "parallelism":1, "taskSpec":spec}],
            "allocationPolicy":{"instances":[instance],
                                "serviceAccount":{"email":os.environ["EPISCAF_WORKER_SA"]}},
            "logsPolicy":{"destination":"CLOUD_LOGGING"}}


def owned(owner, rid):
    r = read(record_blob(owner, rid))
    if not r or not r["submitted"]:
        raise HTTPException(404, "Run unavailable.")
    return r


def launch(owner, rid):
    from google.api_core.exceptions import PreconditionFailed
    blob = record_blob(owner, rid)
    # Compare-and-swap review -> submitting. A REST retry uses the same Batch job ID.
    r = read(blob)
    if r is None:
        raise HTTPException(404, "Run unavailable.")
    generation = blob.generation
    if r["backend"] in {"manual-gemini", "manual-review"}:
        from episcaf_service import handoff
        return handoff.launch(owner, r, blob)
    if r["stage"] == "review" and pilot.mode() != MODE:
        raise HTTPException(409, "Prepare a new review for manual Gemini submission.")
    if r["stage"] == "review":
        now = time.time()
        r.update(submitted=now, updated=now, stage="submitting" if r["backend"]==MODE else "awaiting lab review")
        try:
            write(blob,r,if_generation_match=generation)
        except PreconditionFailed:
            r = owned(owner,rid)
    return status(owner, r)


def stamp(value, fallback):
    return datetime.fromisoformat(value.replace("Z","+00:00")).timestamp() if value else fallback


def status(owner, r):
    from google.api_core.exceptions import PreconditionFailed
    if r["backend"] in {"manual-gemini", "manual-review"}:
        from episcaf_service import handoff
        return handoff.status(owner, r)
    if pilot.mode() != MODE and r["stage"] == "submitting":
        return dict(r, reason="Check the earlier Google Batch submission with the operator.")
    if r["stage"] in ("completed","failed","awaiting lab review"):
        return r
    blob = record_blob(owner, r["id"])
    if r["stage"] != "submitting" and time.time() - r["updated"] < 30:
        return r
    # download_as_bytes populates the generation from the response. Guard the
    # write so an older poll cannot overwrite completion recorded by a peer.
    r = read(blob)
    generation = blob.generation
    if r["stage"] in ("completed", "failed"):
        return r
    job_id = "ep-"+r["id"]
    now = time.time()
    done = read(outputs().blob(r["id"]+"/done.json"))
    if done is not None:
        r.update(stage="completed", completed=done["completed"],
                 reason=f"{done['retained']} of {done['attempted']} backbones retained.")
    elif r["stage"] == "submitting" and now-r["submitted"] < 86400:
        res = batch_request("POST", params={"jobId":job_id, "requestId":str(uuid.UUID(r["id"]))},
                            json=job_spec(r))
        if res.status_code in (400,403,404):
            r.update(stage="failed", completed=now,
                     reason=f"Google Batch rejected submission (HTTP {res.status_code}). Check the project, region, job configuration and service account.")
        else:
            r.update(stage="queued",reason=None)
    else:
        res = batch_request("GET","/"+job_id)
        if res.status_code in (400,403):
            return dict(r, reason="The service could not read the Batch job status. Check its permissions.")
        if res.status_code == 404:
            r.update(stage="failed", completed=now, reason="Google Batch has no record of this job. Check its ID in Cloud Logging.")
        else:
            job = res.json()
            info = job["status"]
            state = info["state"]
            if state in {"FAILED","CANCELLED","DELETION_IN_PROGRESS","SUCCEEDED"}:
                events = info.get("statusEvents",[])
                reason = events[-1].get("description") if events else None
                if state == "SUCCEEDED":
                    reason = "The job ended before publishing its CSV and structure archive."
                r.update(stage="failed", reason=reason or f"Google Batch reported {state}. Check the job log.",
                         completed=stamp(job.get("updateTime"),now))
            else:
                r["stage"] = "generating and filtering backbones" if state == "RUNNING" else "queued"
    r.update(job_id=job_id,updated=now)
    try:
        write(blob, r, if_generation_match=generation)
    except PreconditionFailed:
        return owned(owner,r["id"])
    return r


def runs(owner):
    records = (read(blob) for blob in store().list_blobs(prefix=f"runs/{owner}/"))
    return sorted((status(owner,r) for r in records if r and r["submitted"]),
                  key=lambda r:r["submitted"], reverse=True)


def result(owner, rid, name, download):
    r = status(owner,owned(owner,rid))
    if r["stage"] != "completed":
        raise HTTPException(409,r["reason"] or "Results are pending.")
    content = outputs().blob(rid+"/"+name).download_as_bytes()
    headers = {"Content-Disposition":f'attachment; filename="{rid}-{name}"'} if download else {}
    return Response(content, media_type="text/csv" if name.endswith(".csv") else "application/zip", headers=headers)
