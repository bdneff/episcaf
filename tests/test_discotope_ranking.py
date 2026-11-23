"""Exercise the live upload endpoint with a local stand-in for BioLib inference."""
import csv
import io
import os
from pathlib import Path
import sys
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

from Bio.PDB import MMCIFIO, PDBParser
from Bio.SeqUtils import seq1
from fastapi.testclient import TestClient

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from episcaf_service import disco, lysozyme_job
from episcaf_service.app import app


class Ranking(unittest.TestCase):
    def test_uploaded_multichain_mmcif(self):
        # Change author labels and chain ID so a filename lookup or naïve CSV
        # res_id -> author number join cannot pass this test.
        structure = PDBParser(QUIET=True).get_structure("test", lysozyme_job.REFERENCE)
        chain = structure[0]["Y"]
        chain.id = "Protein_X"
        for residue in list(chain)[::-1]:
            het, num, ins = residue.id
            residue.id = (het, num + 200, "A" if num == 2 and het == " " else ins)
        writer = MMCIFIO()
        writer.set_structure(structure)
        out = io.StringIO()
        writer.save(out)
        def cli(args):
            self.assertEqual(args[-1], "alphafold")
            model = PDBParser(QUIET=True).get_structure("wire", args[1])[0]
            self.assertEqual(list(model.child_dict), ["A"])
            residues = list(model["A"])
            self.assertEqual(len(residues),129)
            def save(directory):
                path = Path(directory)/"antigen_A_discotope3.csv"
                with path.open("w",newline="") as fh:
                    table = csv.writer(fh)
                    table.writerow(["res_id","residue","DiscoTope-3.0_score"])
                    for r in residues:
                        table.writerow([r.id[1],seq1(r.resname),r.id[1]])
            return SimpleNamespace(save_files=save,get_exit_code=lambda:0)
        api = TestClient(app)
        with patch.dict(os.environ,{"EPISCAF_DISCOTOPE":"biolib","EPISCAF_DISCOTOPE_SCORED_PDB":""}), patch.dict(
                sys.modules,{"biolib":SimpleNamespace(load=lambda _:SimpleNamespace(cli=cli))}):
            self.assertTrue(api.get("/health").json()["disco_available"])
            url = "/enumerate?fmt=cif&chain=Protein_X&structure_type=alphafold"
            raw = api.post(url+"&disco=0",content=out.getvalue()).json()
            ranked = api.post(url+"&disco=1&top_k=3",content=out.getvalue()).json()
            self.assertTrue(ranked["has_disco"],ranked)
            expected = {p["i"]:sum((i-201)/128 for i in p["resis"]) for p in raw["patches"]}
            self.assertGreater(len(expected),3)
            self.assertEqual([p["i"] for p in ranked["patches"]],sorted(expected,key=lambda i:(-expected[i],i))[:3])
            for p in ranked["patches"]:
                self.assertAlmostEqual(p["disco_sum"],expected[p["i"]],places=4)
                self.assertEqual(p["disco_n"],p["n_res"])
            with patch.object(disco,"_from_biolib",side_effect=RuntimeError("test failure")):
                fallback = api.post(url+"&disco=1&top_k=3",content=out.getvalue()).json()
                self.assertFalse(fallback["has_disco"])
                self.assertEqual(fallback["patches"],raw["patches"])

    def test_baked_file_is_structure_specific(self):
        with tempfile.TemporaryDirectory() as tmp:
            pdb = Path(tmp)/"source.pdb"
            pdb.write_text(lysozyme_job.reference_pdb())
            self.assertTrue(disco._from_scored_pdb(pdb,pdb,"Y"))
            other = Path(tmp)/"other.pdb"
            other.write_text(pdb.read_text().replace("LYS", "ALA"))
            self.assertIsNone(disco._from_scored_pdb(pdb,other,"Y"))
        with patch.dict(os.environ,{"EPISCAF_DISCOTOPE":"off","EPISCAF_DISCOTOPE_SCORED_PDB":""}):
            self.assertFalse(disco.available())


if __name__ == "__main__":
    unittest.main()
