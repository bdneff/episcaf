"""Exercise the operator commands without a cloud account, SSH or a Slurm installation."""
import io
import json
from pathlib import Path
import subprocess
import sys
import unittest
from unittest.mock import patch
from zipfile import ZipFile

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from episcaf_service.gcloud import pilot_transfer as transfer


class Transfer(unittest.TestCase):
    def test_request_to_publication(self):
        rid = "a"*32
        source = f"gs://episcaf-requests-123/inbox/john/{rid}/reference.zip"
        out = f"gs://episcaf-results-123/{rid}/"
        bucket = {source:b"request fixture"}
        home = {}
        def commands(command, **kw):
            command = list(command)
            code, stdout = 0, b""
            if command[:3] == ["gcloud","config","get-value"]:
                stdout = b"lab"
            elif command[:3] == ["gcloud","projects","describe"]:
                stdout = b"123"
            elif command[:4] == ["gcloud","storage","objects","describe"]:
                code = int(command[4] not in bucket)
            elif command[:3] == ["gcloud","storage","cat"]:
                stdout = bucket.get(command[3],b"")
                code = int(command[3] not in bucket)
            elif command[:3] == ["gcloud","storage","cp"]:
                src, dst = command[-2:]
                if src.startswith("gs://"):
                    Path(dst).write_bytes(bucket[src])
                elif dst in bucket:
                    code = 1
                else:
                    self.assertIn("--if-generation-match=0",command)
                    bucket[dst] = Path(src).read_bytes()
            elif command[:2] == ["ssh","gemini"]:
                self.assertNotIn("sbatch", command[-1])
                self.assertIn(rid, command[-1])
                home["request.zip"] = kw["stdin"].read()
            elif command[0] == "scp":
                Path(command[2]).write_bytes(home[command[1].split("/")[-1]])
            else:
                self.fail(command)
            if code and kw.get("check"):
                raise subprocess.CalledProcessError(code,command)
            return subprocess.CompletedProcess(command,code,stdout.decode() if kw.get("text") else stdout)
        def invoke(*args):
            with patch.object(sys,"argv",["pilot_transfer",args[0],"john",rid,*args[1:]]):
                transfer.main()
        with patch.object(transfer.subprocess,"run",commands):
            invoke("fetch")
            self.assertEqual(home["request.zip"],bucket[source])
            invoke("receipt","42")
            invoke("receipt","42")
            home["terminal.json"] = json.dumps(dict(id=rid,stage="completed",job_id="42")).encode()
            home["retained.csv"] = b"design_id\nmodel_0\n"
            buf = io.BytesIO()
            with ZipFile(buf,"w") as z:
                z.writestr("retained.csv",home["retained.csv"])
                z.writestr("backbones/model_0.pdb",b"fixture")
            home["designs.zip"] = buf.getvalue()
            with self.assertRaises(ValueError):
                invoke("publish","43")
            self.assertNotIn(out+"terminal.json",bucket)
            invoke("publish","42")
            invoke("publish","42")
            self.assertEqual(bucket[out+"terminal.json"],home["terminal.json"])
            with self.assertRaises(ValueError):
                invoke("fail","42","--reason","timeout")
            bucket.pop(out+"terminal.json")
            invoke("fail","42","--reason","Slurm reported TIMEOUT")
            invoke("fail","42","--reason","Slurm reported TIMEOUT")
            self.assertEqual(json.loads(bucket[out+"terminal.json"])["stage"],"failed")


if __name__ == "__main__":
    unittest.main()
