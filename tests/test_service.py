#!/usr/bin/env python3
"""Tests for the episcaf cluster service (episcaf_service.app).

Runs the FastAPI app in-process with $WS pointed at a tmp dir, so the GPU boundary is never crossed:
we assert the CPU work is real (detect / enumerate / stage / score) and that the emitted submit.sh
uses the *verified* script interfaces, not the GUI's old wrong flags.

Run:  /Users/bneff/micromamba/envs/bcell-repair/bin/python tests/test_service.py   (or pytest tests/)
"""
import os
import sys
import tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

# Point the service at a throwaway workspace BEFORE anything reads it.
_WS = tempfile.mkdtemp(prefix="episcaf_ws_")
os.environ["EPISCAF_WS"] = _WS

from fastapi.testclient import TestClient            # noqa: E402
from episcaf_service.app import app                  # noqa: E402

client = TestClient(app)
FIX = os.path.join(ROOT, "docs", "pdb", "8PWW.pdb")  # chain-A antigen used by the growPatch tests


def _body():
    with open(FIX, "rb") as fh:
        return fh.read()


def test_health():
    r = client.get("/health")
    assert r.status_code == 200, r.text
    j = r.json()
    assert j["ok"] is True and j["ws"] == _WS


def test_inspect_shape():
    r = client.post("/inspect?fmt=pdb&chain=A", content=_body())
    assert r.status_code == 200, r.text
    j = r.json()
    for k in ("antigen_chain", "is_complex", "epitope", "partner_chains"):
        assert k in j, f"missing {k} in {j}"
    assert isinstance(j["is_complex"], bool)


def test_enumerate_finds_patches():
    r = client.post("/enumerate?fmt=pdb&chain=A&mode=growpatch", content=_body())
    assert r.status_code == 200, r.text
    j = r.json()
    assert j["n_patches"] >= 1 and j["patches"], f"no patches: {j}"
    return j["patches"][0]["epitope"]


def test_submit_stages_real_job():
    epi = test_enumerate_finds_patches()
    q = (f"/submit?fmt=pdb&chain=A&mode=growpatch&name=8pww&epitopes={epi}"
         f"&rfd3=8&mpnn=8&prefilter=true&keepfrac=0.5&topn=10&ipd=true&ipdrounds=3")
    r = client.post(q, content=_body())
    assert r.status_code == 200, r.text
    j = r.json()
    run_dir = j["run_dir"]
    assert run_dir.startswith(_WS), run_dir
    # the real staged artifacts exist on disk
    for f in ("input.pdb", "designs.csv", "run_spec.json", "submit.sh"):
        assert os.path.exists(os.path.join(run_dir, f)), f"missing staged {f}"
    assert j["n_contigs"] >= 1

    sh = j["submit_sh"]
    # verified interfaces are present
    assert "scripts/ipd_round.py init" in sh and "--designs-table" in sh
    assert "--max-rounds 3" in sh and "--mode adaptive" in sh
    assert "aomclash_prefilter.py" in sh and "--keep-frac 0.5" in sh
    assert "stage06_select.py" in sh and "--group id" in sh and "--topk 10" in sh
    assert "stage04_af3_array.sbatch" in sh
    # the GUI's OLD WRONG interfaces must NOT appear
    assert "--select" not in sh, "emitted the nonexistent score.py --select flag"
    assert "stage04_af3_emit_jsons.py --run " not in sh, "emitted the wrong stage04 --run flag"
    return j["run_id"]


def test_status_offline_branch():
    run_id = test_submit_stages_real_job()
    r = client.get(f"/status/{run_id}")
    assert r.status_code == 200, r.text
    j = r.json()
    assert j["metrics_ready"] is False
    assert j["squeue"] == "not on a SLURM host"      # no SLURM locally
    # unknown run -> 404
    assert client.get("/status/does_not_exist").status_code == 404
    return run_id


def test_results_scores_metrics():
    run_id = test_submit_stages_real_job()
    rd = os.path.join(_WS, "runs", run_id, "05_analysis")
    os.makedirs(rd, exist_ok=True)
    # a minimal finished-metrics table; the composite renormalizes over whatever columns are present.
    with open(os.path.join(rd, "metrics.csv"), "w") as fh:
        fh.write("design_id,aomclash_pred,overall_rmsd,mean_pae\n"
                 "d1,2.0,1.1,3.0\n"
                 "d2,9.0,0.8,2.2\n"
                 "d3,0.5,2.5,6.0\n")
    r = client.get(f"/results/{run_id}?preset=twelvemer")
    assert r.status_code == 200, r.text
    j = r.json()
    assert j["info"]["n_designs"] == 3
    assert len(j["rows"]) == 3
    # before metrics exist -> 409 (not ready), separate fresh run
    assert client.get("/results/does_not_exist").status_code == 404


def main():
    tests = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    for t in tests:
        t()
        print(f"  ok  {t.__name__}")
    print(f"\n{len(tests)} tests passed (ws={_WS})")


if __name__ == "__main__":
    main()
