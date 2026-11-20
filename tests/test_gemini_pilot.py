"""Offline Gemini acceptance with fake Slurm receipts and checked-in lysozyme coordinates."""
import csv
import gzip
import io
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch
from zipfile import ZipFile

from Bio.PDB import PDBParser, MMCIFIO
from fastapi.testclient import TestClient

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from episcaf_service import gemini, lysozyme_job as job, pilot, reviewed
from episcaf_service.app import app


class GeminiPilot(unittest.TestCase):
    def setUp(self):
        temp = tempfile.TemporaryDirectory()
        self.addCleanup(temp.cleanup)
        self.root = Path(temp.name)
        private = self.root / "artifacts"
        private.mkdir(mode=0o700)
        env = patch.dict(os.environ, {"EPISCAF_PILOT_MODE":gemini.BACKEND, "K_SERVICE":"",
                                      "EPISCAF_WS":str(self.root / "legacy"),
                                      "EPISCAF_GEMINI_RUNS":str(private),
                                      "EPISCAF_PILOT_DB":str(self.root / "runs.sqlite")})
        env.start()
        self.addCleanup(env.stop)
        self.alice = self.enterContext(TestClient(app))
        self.bob = self.enterContext(TestClient(app))
        self.token = pilot.issue("alice")
        self.signin(self.alice, self.token)
        self.signin(self.bob, pilot.issue("bob"))
        self.pdb = self.alice.get("/pilot/reference.pdb").text

    def signin(self, client, token):
        self.assertEqual(client.post("/pilot/session", headers={"Authorization":"Bearer "+token}).status_code, 200)

    def review(self, pdb, regions, chain="Y"):
        return self.alice.post("/pilot/reviews", json=dict(filename="my-upload.pdb", structure=pdb,
                                                           chain=chain, regions=regions))

    def launch(self, pdb, regions, chain="Y"):
        response = self.review(pdb, regions, chain)
        self.assertEqual(response.status_code, 200, response.text)
        plan = response.json()
        url = "/pilot/runs/" + plan["id"]
        first = self.alice.post(url)
        self.assertEqual(first.status_code, 200, first.text)
        self.assertEqual(first.json(), self.alice.post(url).json())
        return plan, url

    def test_reference_budget_job_receipt_and_failure(self):
        data = self.alice.post("/enumerate?chain=Y&disco=0", content=self.pdb)
        regions = [p["epitope"] for p in data.json()["patches"][:3]]
        self.assertEqual(len(regions), 3)
        plan, url = self.launch(self.pdb, regions)
        self.assertEqual(plan["total"], 24)
        self.assertEqual(plan["sequence_predictions"], 192)
        self.assertEqual(plan["backend"], gemini.BACKEND)
        self.assertEqual(self.alice.get(url).json()["stage"], "requested")
        with patch.object(subprocess, "run", return_value=subprocess.CompletedProcess([],0,"712;gemini\n","")) as submit:
            gemini.dispatch(False)
            pilot.initialize()
            gemini.dispatch(False)
            self.assertEqual(submit.call_count, 1)
            self.assertEqual(submit.call_args.args[0][:2], ["sbatch","--parsable"])
        directory = gemini.run_dir(plan["id"])
        prepared = job.prepare(directory, json.loads((directory / "selection.json").read_text())["regions"])
        self.assertEqual([r["backbone_trials"] for r in plan["contigs"]],
                         [r["backbone_trials"] for r in prepared["contigs"]])
        self.assertEqual(sum(r["backbone_trials"] for r in prepared["contigs"]), 24)
        self.assertEqual(len(list((directory / "inputs").glob("*.json"))), len(plan["contigs"]))
        self.assertEqual(self.review(self.pdb, ["Y:130-140"]).status_code, 422)
        with patch.object(gemini, "slurm", side_effect=[[], [["712","PENDING","2026-11-21T12:00:00"]]]):
            gemini.dispatch(True)
        status = self.alice.get(url).json()
        self.assertTrue(status["estimated_start"])
        self.assertIsNone(status["completed"])
        self.alice.delete("/pilot/session")
        with TestClient(app) as fresh:
            self.signin(fresh, self.token)
            self.assertEqual(fresh.get(url).json()["job_id"], "712")
            job.stage(directory, "Predicting structures")
            with patch.object(gemini, "slurm", return_value=[["712","TIMEOUT","1:0","Unknown"]]):
                gemini.dispatch(True)
            self.assertIn("Predicting structures (TIMEOUT", fresh.get(url).json()["reason"])
        self.assertEqual(self.bob.get(url).status_code, 404)
        self.assertEqual(self.alice.get("/results/"+plan["id"]).status_code, 404)

    def test_operator_delivery_for_own_upload(self):
        # A structure with a different chain remains a private lab request.
        pdb = self.pdb.replace(" Y ", " A ")
        plan, url = self.launch(pdb, ["A:20-34"], "A")
        rid = plan["id"]
        self.assertEqual(plan["backend"], "gemini-reviewed")
        with patch.object(subprocess, "run") as submit:
            gemini.dispatch(True)
            submit.assert_not_called()
        path = reviewed.export(rid)
        self.assertEqual((path / "upload.pdb").read_text(), pdb)
        self.assertEqual(json.loads((path / "request.json").read_text())["regions"], ["A:20-34"])
        scores = self.root / "control.csv"
        scores.write_text("run_id,design_id,epitope,sequence,composite\n"+rid+",control,A:20-34,ACDE,0.5\n")
        with self.assertRaises(ValueError):
            reviewed.publish(rid, scores, self.root, "713")
        self.assertEqual(self.alice.get(url).json()["stage"], "awaiting lab review")
        (self.root / "control.cif").write_text("data_control\n# offline import fixture\n")
        reviewed.publish(rid, scores, self.root, "713")
        self.assertIn("lab-reviewed", self.alice.get(url+"/results").text)
        with ZipFile(io.BytesIO(self.alice.get(url+"/designs.zip").content)) as z:
            self.assertEqual(set(z.namelist()), {"ranked.csv","models/control.cif"})
        self.assertEqual(self.bob.get("/pilot/runs").json(), [])
        for route in (url, url+"/results", url+"/results?download=true", url+"/designs.zip"):
            self.assertEqual(self.bob.get(route).status_code, 404)

    def test_receipt_recovery_and_missing_artifacts(self):
        plan, url = self.launch(self.pdb, ["Y:20-34"])
        with patch.object(subprocess, "run", side_effect=subprocess.TimeoutExpired("sbatch",90)):
            gemini.dispatch(False)
        self.assertEqual(self.alice.get(url).json()["stage"], "submitting")
        with patch.object(gemini, "slurm", return_value=[["714","ep-"+plan["id"]]]):
            gemini.dispatch(True)
        self.assertEqual(self.alice.get(url).json()["job_id"], "714")
        with patch.object(gemini, "slurm", return_value=[["714","COMPLETED","0:0","Unknown"]]):
            gemini.dispatch(True)
        self.assertIn("before the scored designs", self.alice.get(url).json()["reason"])
        with patch.dict(os.environ, {"K_SERVICE":"public-cloud-run"}):
            self.assertEqual(pilot.mode(), "staging")
            self.assertEqual(self.alice.get(url).status_code, 404)

    def test_reference_worker_stage_contracts(self):
        # Stand-ins write only deterministic controls. Exercise the real fixed-PDB
        # conversion, sequence checks and AF3 JSON emitter without GPU inference.
        from scripts.stage05_extract_metrics import epi_positions_from_contig
        directory = self.root / "worker"
        directory.mkdir()
        (directory / "selection.json").write_text(json.dumps({"regions":["Y:20-34"]}))
        native_ca = [line for line in job.reference_pdb("A").splitlines() if line[12:16].strip() == "CA"]
        call = subprocess.run
        batches = []
        def fixture_command(args, **kwargs):
            if args[0] == "rfd3":
                flags = dict(x.lstrip("+").split("=", 1) for x in args[2:])
                spec, = json.loads(Path(flags["inputs"]).read_text()).values()
                design_idx, native_idx = epi_positions_from_contig(spec["contig"].replace(",", "/"))
                names = {d:native_ca[n][17:20] for d,n in zip(design_idx,native_idx)}
                lines = [line[:17]+names.get(i,"ALA")+line[20:] for i,line in enumerate(native_ca[:103])]
                pdb = directory / "fixture.pdb"
                pdb.write_text("\n".join(lines)+"\nEND\n")
                writer = MMCIFIO()
                writer.set_structure(PDBParser(QUIET=True).get_structure("control", pdb))
                stream = io.StringIO()
                writer.save(stream)
                output = Path(flags["out_dir"])
                output.mkdir(parents=True)
                self.assertEqual(flags["n_batches"], "1")
                batches.append(int(flags["diffusion_batch_size"]))
                for i in range(batches[-1]):
                    (output / f"control_{i}.cif.gz").write_bytes(gzip.compress(stream.getvalue().encode()))
            elif "-seqs_per_struct" in args:
                for backbone in (directory / "fixed").glob("*.pdb"):
                    for seq in range(int(args[args.index("-seqs_per_struct")+1])):
                        path = directory / "03_mpnn" / f"{backbone.stem}_dldesign_{seq}.pdb"
                        path.write_bytes(backbone.read_bytes())
            elif args[0] == "bash":
                self.assertEqual(kwargs["env"]["SLURM_ARRAY_TASK_ID"], "1")
            else:
                self.assertEqual(Path(args[1]).name, "stage04_af3_emit_jsons.py")
                return call(args, **kwargs)
        with patch.dict(os.environ, {"SLURM_JOB_ID":"710"}), patch.object(subprocess, "run", fixture_command), \
                patch.object(job, "collect") as collect:
            job.execute(directory)
        self.assertEqual(sum(batches), 8)
        self.assertEqual(len(collect.call_args.args[1]), 64)
        self.assertEqual(len(list((directory / "04_af3/inputs").glob("*.json"))), 64)

    def test_reference_cpu_scores_and_download(self):
        plan, url = self.launch(self.pdb, ["Y:20-34"])
        with patch.object(subprocess, "run", return_value=subprocess.CompletedProcess([],0,"715","")):
            gemini.dispatch(False)
        directory = gemini.run_dir(plan["id"])
        # Use the native 129-residue reference for both structures and a synthetic
        # confidence matrix strictly to exercise extraction, ranking and delivery.
        design = "Lysozyme__contig0__seed0__rep0_model_0_fixed_dldesign_1"
        output = directory / "04_af3/outputs" / design.lower()
        output.mkdir(parents=True)
        for folder in ("native","03_mpnn"):
            (directory / folder).mkdir()
        native = directory / "native/Lysozyme.pdb"
        native.write_text(job.reference_pdb("A"))
        mpnn = directory / "03_mpnn" / (design+".pdb")
        mpnn.write_bytes(native.read_bytes())
        cif = MMCIFIO()
        cif.set_structure(PDBParser(QUIET=True).get_structure("control", native))
        cif.save(str(output / "control_model.cif"))
        (output / "control_confidences.json").write_text(json.dumps({"pae":[[2.0]*129]*129}))
        (directory / "designs.csv").write_text(
            "id,contig_id,contig_string,growpatch_epitope\nLysozyme,0,19-19/A20-34/95-95,A:20-34\n")
        job.collect(directory, {design:mpnn})
        with patch.object(gemini, "slurm", return_value=[["715","COMPLETED","0:0","2026-11-21T13:00:00"]]):
            gemini.dispatch(True)
        self.assertEqual(self.alice.get(url).json()["stage"], "completed")
        result, = csv.DictReader(io.StringIO(self.alice.get(url+"/results").text))
        self.assertEqual(result["design_id"], design)
        self.assertEqual(result["epitope"], "Y:20-34")
        self.assertEqual(len(result["sequence"]), 129)
        self.assertEqual(float(result["s_p1_free"]), 0.0)
        self.assertEqual(float(result["epitope_pae"]), 2.0)
        self.assertLess(float(result["overall_rmsd"]), 1e-6)
        self.assertEqual(float(result["composite"]), 0.0)
        self.assertNotIn("aomclash_pred", result)
        self.assertEqual(self.alice.get(url+"/designs.zip").status_code, 200)
        self.assertEqual(self.bob.get(url+"/designs.zip").status_code, 404)


if __name__ == "__main__":
    unittest.main()
