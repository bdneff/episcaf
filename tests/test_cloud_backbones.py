"""Run the private web and worker contracts with local GCS, Batch and GPU stand-ins."""
import csv
import gzip
import io
import json
import os
import secrets
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch
from zipfile import ZipFile

from Bio.PDB import MMCIFIO, PDBParser
from fastapi.testclient import TestClient
from google.api_core.exceptions import NotFound, PreconditionFailed

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from episcaf_service import cloud, cloud_backbones as worker, lysozyme_job as reference
from episcaf_service.app import app
from scripts.stage05_extract_metrics import epi_positions_from_contig


class Bucket:
    def __init__(self):
        self.objects = {}

    def blob(self, key):
        return Blob(self.objects, key)

    def list_blobs(self, prefix):
        return [self.blob(k) for k in self.objects if k.startswith(prefix)]


class Blob:
    def __init__(self, objects, key):
        self.objects = objects
        self.key = key
        self.generation = None
        self.time_created = None

    def exists(self):
        return self.key in self.objects

    def reload(self):
        if not self.exists():
            raise NotFound("missing")
        self.generation = self.objects[self.key][0]

    def download_as_bytes(self):
        self.reload()
        return self.objects[self.key][1]

    def upload_from_string(self, value, if_generation_match=None, **kwargs):
        generation = self.objects.get(self.key, (0,None))[0]
        if if_generation_match is not None and if_generation_match != generation:
            raise PreconditionFailed("object changed")
        self.objects[self.key] = generation+1, value.encode() if isinstance(value,str) else value


class CloudPilot(unittest.TestCase):
    def setUp(self):
        self.store, self.results, self.jobs = Bucket(), Bucket(), {}
        self.enterContext(patch.dict(os.environ, {"K_SERVICE":"episcaf-pilot","EPISCAF_PILOT_MODE":cloud.MODE,
            "EPISCAF_RESULT_BUCKET":"results","EPISCAF_GPU_IMAGE":"worker@digest","EPISCAF_WORKER_SA":"worker"}))
        self.enterContext(patch.object(cloud, "store", return_value=self.store))
        self.enterContext(patch.object(cloud, "outputs", return_value=self.results))
        self.enterContext(patch.object(cloud, "batch_request", side_effect=self.batch))
        self.alice = self.enterContext(TestClient(app,base_url="https://testserver"))
        self.bob = self.enterContext(TestClient(app,base_url="https://testserver"))
        for user, c in [("alice", self.alice), ("bob", self.bob)]:
            token = user+"."+secrets.token_urlsafe(32)
            cloud.write(self.store.blob("testers/"+user),{"digest":cloud.pilot.digest(token)})
            self.assertEqual(c.post("/pilot/session", headers={"Authorization":"Bearer "+token}).status_code,200)

    def batch(self, method, suffix="", **kwargs):
        from requests import Response
        resp = Response()
        if method == "POST":
            job_id = kwargs["params"]["jobId"]
            resp.status_code = 409 if job_id in self.jobs else 200
            self.jobs[job_id] = kwargs["json"]
        else:
            resp.status_code = 200
            resp._content = b'{"status":{"state":"RUNNING"}}'
        return resp

    def review(self, route, structure, regions):
        return self.alice.post("/pilot/reviews",json={"route":route,"structure":structure,
            "filename":"trial.pdb","chain":"Y","regions":regions})

    def test_private_paths(self):
        for route in ["growpatch","antibody"]:
            pdb = self.alice.get("/pilot/reference.pdb?route="+route).text
            region = worker.contact_region() if route=="antibody" else "Y:20-34"
            plan = self.review(route,pdb,[region]).json()
            self.assertEqual(plan["regions"],[region])
            self.assertEqual(sum(r["backbone_trials"] for r in plan["contigs"]),16)
            self.assertEqual(plan["retained_limit"],8)
            url = "/pilot/runs/"+plan["id"]
            self.assertEqual(self.bob.post(url).status_code,404)
            self.assertEqual(self.alice.post(url).status_code,200)
            self.assertEqual(self.alice.post(url).status_code,200)
            spec = self.jobs["ep-"+plan["id"]]["taskGroups"][0]["taskSpec"]
            self.assertEqual(spec["maxRetryCount"],0)
            self.assertEqual(spec["runnables"][0]["container"]["commands"][:4],["--route",route,"--region",region])
            self.assertNotIn("ATOM",json.dumps(spec))
            self.results.blob(plan["id"]+"/retained.csv").upload_from_string("design_id\na\n")
            self.results.blob(plan["id"]+"/designs.zip").upload_from_string(b"zip fixture")
            cloud.write(self.results.blob(plan["id"]+"/done.json"),{"attempted":16,"retained":8,"completed":12345})
            with patch("time.time",return_value=9999999999):
                self.assertEqual(self.alice.get(url).json()["stage"],"completed")
            for suffix in ("","/results","/designs.zip"):
                self.assertEqual(self.bob.get(url+suffix).status_code,404)
                self.assertEqual(self.alice.get(url+suffix).status_code,200)
        self.assertEqual(len(self.jobs),2)
        self.assertEqual(self.bob.get("/pilot/runs").json(),[])
        cloud.write(self.store.blob("testers/alice"),{"digest":None})
        self.assertEqual(self.alice.get("/pilot/runs").status_code,401)

    def test_manual_gemini_handoff(self):
        from episcaf_service import handoff_job
        with patch.dict(os.environ, {"EPISCAF_PILOT_MODE":"manual-gemini"}):
            for route, region in [("growpatch","Y:20-34"), ("antibody",worker.contact_region())]:
                upload = self.alice.get("/pilot/reference.pdb?route="+route).text
                plan = self.review(route, upload, [region]).json()
                rid = plan["id"]
                url = "/pilot/runs/"+rid
                self.assertEqual(plan["total"],16)
                self.assertEqual(self.bob.post(url).status_code,404)
                for _ in range(2):
                    self.assertEqual(self.alice.post(url).json()["stage"],"awaiting Brandon’s submission")
                bundle = self.store.blob("inbox/alice/"+rid+"/reference.zip")
                data = bundle.download_as_bytes()
                self.assertEqual(bundle.generation,1)
                response = self.alice.get(url+"/request.zip")
                self.assertEqual(response.content, data)
                self.assertIn("attachment", response.headers["content-disposition"])
                self.assertEqual(self.bob.get(url+"/request.zip").status_code,404)
                cloud.write(self.results.blob(rid+"/receipt.json"),{"stage":"queued","job_id":"12345"})
                self.assertEqual(self.alice.get(url).json()["job_id"], "12345")
                with tempfile.TemporaryDirectory() as tmp:
                    base = Path(tmp)
                    inbox = base / "home" / rid
                    inbox.mkdir(parents=True)
                    (inbox/"request.zip").write_bytes(data)
                    def fake_generate(directory, selected_route, selected_region, source):
                        self.assertEqual((selected_route,selected_region,source),(route,region,"gemini-backbone"))
                        (directory/"retained.csv").write_text("source,design_id,s_p1_free\n")
                        (directory/"designs.zip").write_bytes(b"archive fixture")
                        return 0
                    with patch.object(worker,"generate",side_effect=fake_generate):
                        handoff_job.execute(inbox,base/"ws","12345")
                        with self.assertRaises(FileExistsError):
                            handoff_job.execute(inbox,base/"ws","12346")
                    self.assertEqual((base/"ws/pilot/requests"/rid/"request.zip").read_bytes(),data)
                    for name in ("terminal.json","retained.csv","designs.zip"):
                        self.results.blob(rid+"/"+name).upload_from_string((inbox/name).read_bytes())
                self.assertEqual(self.alice.get(url).json()["job_id"],"12345")
                self.assertEqual(self.alice.get(url+"/designs.zip").status_code,200)
                self.assertEqual(self.bob.get(url+"/results").status_code,404)
            plan = self.review("growpatch","unmatched upload",["Y:20-34"]).json()
            self.assertEqual(self.alice.post("/pilot/runs/"+plan["id"]).json()["stage"],"awaiting lab review")
            path = "/pilot/runs/"+plan["id"]+"/request.zip"
            with ZipFile(io.BytesIO(self.alice.get(path).content)) as files:
                self.assertEqual(json.loads(files.read("request.json"))["regions"], ["Y:20-34"])
                self.assertEqual(files.read("upload"), b"unmatched upload")
            self.assertEqual(self.bob.get(path).status_code,404)
            cloud.write(self.store.blob("testers/alice"),{"digest":None})
            self.assertEqual(self.alice.get(path).status_code,401)
            self.assertFalse(self.jobs)

    def test_lost_receipt_reconciles_and_failed_job_is_terminal(self):
        from requests import Response
        plan = self.review("growpatch",reference.reference_pdb(),["Y:20-34"]).json()
        url = "/pilot/runs/"+plan["id"]
        def lost(method, **kwargs):
            self.batch(method, **kwargs)
            raise cloud.HTTPException(503, "Lost receipt")
        with patch.object(cloud, "batch_request", side_effect=lost):
            self.assertEqual(self.alice.post(url).status_code,503)
        self.assertEqual(self.alice.post(url).json()["stage"],"queued")
        self.assertEqual(len(self.jobs),1)
        failed = Response()
        failed.status_code = 200
        failed._content = b'{"status":{"state":"FAILED","statusEvents":[{"description":"GPU task timed out"}]}}'
        with patch.object(cloud, "batch_request", return_value=failed), patch("time.time",return_value=9999999999):
            self.assertEqual(self.alice.get(url).json()["reason"], "GPU task timed out")
            self.assertEqual(self.alice.get(url+"/designs.zip").status_code,409)

    def test_budget_and_target_enforced_on_server(self):
        pdb = reference.reference_pdb()
        self.assertEqual(self.review("growpatch",pdb,["Y:15-24","Y:31-40"]).status_code,422)
        self.assertEqual(self.review("growpatch",pdb,["Y:130-140"]).status_code,422)
        changed = pdb.replace("ATOM  ","ATOM X",1)
        plan = self.review("growpatch",changed,["Y:20-34"]).json()
        self.assertEqual(self.alice.post("/pilot/runs/"+plan["id"]).json()["stage"],"awaiting lab review")
        self.assertFalse(self.jobs)
        with patch.dict(os.environ, {"K_SERVICE":"episcaf"}):
            self.assertEqual(self.alice.get("/health").json()["pilot_mode"],"staging")
            self.assertEqual(self.alice.get("/pilot/runs").status_code,404)


class ReferenceFilter(unittest.TestCase):
    def test_reference_and_contact_scores(self):
        # CPU geometry control. The RFd3 stand-in only rearranges native coordinates.
        # Every other step, including mmCIF conversion, raw occupancy and ZIP export, is real.
        atoms = {}
        for line in reference.reference_pdb("A").splitlines():
            if line.startswith("ATOM") and line[12:16].strip() in {"N","C","CA","O"}:
                atoms.setdefault(int(line[22:26])-1,[]).append(line)

        def fixture(cmd, **kwargs):
            flags = dict(arg.lstrip("+").split("=",1) for arg in cmd[2:])
            spec, = json.loads(Path(flags["inputs"]).read_text()).values()
            design, native = epi_positions_from_contig(spec["contig"].replace(",","/"))
            mapping = dict(zip(design,native))
            lines = [line[:22]+f"{i+1:4}"+line[26:] for i in range(103)
                     for line in atoms[mapping.get(i,i)]]
            pdb = io.StringIO("\n".join(lines)+"\nEND\n")
            writer = MMCIFIO()
            writer.set_structure(PDBParser(QUIET=True).get_structure("reference",pdb))
            stream = io.StringIO()
            writer.save(stream)
            out = Path(flags["out_dir"])
            out.mkdir()
            for i in range(int(flags["diffusion_batch_size"])):
                (out/f"model_{i}.cif.gz").write_bytes(gzip.compress(stream.getvalue().encode()))
        for route, region in [("growpatch","Y:20-34"),("antibody",worker.contact_region())]:
            with tempfile.TemporaryDirectory() as tmp, patch.object(worker.subprocess,"run",fixture):
                directory = Path(tmp)
                self.assertEqual(worker.generate(directory,route,region),8)
                with ZipFile(directory/"designs.zip") as z:
                    self.assertEqual(len([k for k in z.namelist() if k.endswith(".pdb")]),8)
                    scores = list(csv.DictReader(z.read("backbone_scores.csv").decode().splitlines()))
                    self.assertEqual(len(scores),16)
                    self.assertEqual(sum(r["kept"]=="True" for r in scores),8)
                    self.assertTrue(all(float(r["s_p1_free"])>=0 for r in scores))


if __name__ == "__main__":
    unittest.main()
