"""episcaf cluster service -- the real job broker behind the GUI tool.

A FastAPI app meant to run on a cluster-access host (gemini) and be reached over a private SSH
tunnel (`ssh -L 8808:localhost:8808 <cluster>`). It does every CPU step for real -- epitope
detection, growPatch enumeration, contig staging, and scoring -- and stops at the GPU boundary:
`/submit` stages a run under $WS and hands back the exact submit script; the user runs the GPU jobs.
Claude never submits them (the Gemini boundary).

Wire convention matches the local stdlib server (`scripts/growpatch_serve.py`): the structure is the
raw request body, everything else is a query param, so the same frontend code talks to either backend
by swapping the base URL.

Run (after `pip install fastapi uvicorn`):
    uvicorn episcaf_service.app:app --host 127.0.0.1 --port 8808
See episcaf_service/README.md for the tunnel + GUI wiring.
"""
from __future__ import annotations

import os
import shutil
import subprocess
import tempfile

from fastapi import FastAPI, Query, Request, Response
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import JSONResponse

from episcaf_analysis.growpatch import read as R
from episcaf_service import scoring, staging

# $WS is persistent cluster storage; /scratch is ephemeral, never used for run artifacts.
# EPISCAF_WS overrides it (used by the local test suite to point at a tmp dir).
try:
    from configs.paths import WORKSPACE as _WS_DEFAULT
    _WS_DEFAULT = str(_WS_DEFAULT)
except Exception:
    _WS_DEFAULT = "/tgen_labs/altin/alphafold3/workspace/episcaf_v2_bneff"


def ws_root() -> str:
    return os.environ.get("EPISCAF_WS", _WS_DEFAULT)


app = FastAPI(title="episcaf scaffolding service", version="0.1.0")
# Localhost-only MVP reached cross-origin from the GUI (Pages / local server), so the browser needs
# CORS to call it at all. Blanket-permissive is acceptable while it is tunnel-private; tightening it
# is part of the deferred public-exposure hardening (private key + rate limit).
app.add_middleware(CORSMiddleware, allow_origins=["*"], allow_methods=["*"], allow_headers=["*"])


def _tmp_structure(body: bytes, fmt: str):
    tmp = tempfile.NamedTemporaryFile(suffix="." + ("cif" if fmt == "cif" else "pdb"),
                                      delete=False, mode="wb")
    tmp.write(body)
    tmp.close()
    return tmp.name


@app.get("/health")
def health():
    ws = ws_root()
    return {"ok": True, "service": "episcaf", "ws": ws, "ws_exists": os.path.isdir(ws),
            "slurm": bool(shutil.which("squeue"))}


@app.post("/inspect")
async def inspect(request: Request, fmt: str = Query("pdb"), chain: str = Query("")):
    body = await request.body()
    path = _tmp_structure(body, fmt)
    try:
        return R.detect_complex(path, chain or None)
    except Exception as e:
        return JSONResponse(status_code=400, content={"error": f"{type(e).__name__}: {e}"})
    finally:
        os.unlink(path)


@app.post("/enumerate")
async def enumerate_patches(request: Request, fmt: str = Query("pdb"), chain: str = Query("A"),
                            sasa_cut: float = Query(0.2), footprint_radius: float = Query(14.0),
                            methods: str = Query("growpatch,eyeview"), top_k: int = Query(0),
                            mode: str = Query("growpatch")):
    body = await request.body()
    path = _tmp_structure(body, fmt)
    try:
        ms = tuple(m for m in methods.split(",") if m)
        return staging.viewer_data(path, chain or "A", sasa_cut, footprint_radius, ms, top_k, mode)
    except Exception as e:
        return JSONResponse(status_code=400, content={"error": f"{type(e).__name__}: {e}"})
    finally:
        os.unlink(path)


@app.post("/submit")
async def submit(request: Request, fmt: str = Query("pdb"), chain: str = Query("A"),
                 mode: str = Query("growpatch"), name: str = Query("antigen"),
                 epitopes: str = Query(""), preset: str = Query(""),
                 variants: int = Query(4), total_len: int = Query(103),
                 rfd3: int = Query(8), mpnn: int = Query(8), seed: int = Query(0),
                 helix_bias: float = Query(0.0),
                 budget_mode: str = Query("n_designs"), budget_value: float = Query(200.0),
                 ipd: bool = Query(False), ipdrounds: int = Query(1),
                 ipd_mode: str = Query("adaptive"), prefilter: bool = Query(False),
                 keepfrac: float = Query(0.5), topn: int = Query(10), group: str = Query("id")):
    """Stage a run under $WS and return {run_id, run_dir, submit_sh, ...}. No GPU work here."""
    body = await request.body()
    structure_text = body.decode("utf-8", errors="replace")
    epis = [e for e in epitopes.split("|") if e]
    if not preset:
        preset = "antibody_softgate" if mode == "antibody" else "twelvemer"
    cfg = dict(name=name or "antigen", variants=variants, total_len=total_len, rfd3=rfd3,
               mpnn=mpnn, seed=seed, mode=mode, helix_bias=helix_bias,
               budget={"mode": budget_mode, "value": budget_value})
    settings = dict(ipd=ipd, ipdrounds=ipdrounds, ipd_mode=ipd_mode, prefilter=prefilter,
                    keepfrac=keepfrac, topn=topn, group=group)
    try:
        return staging.persist_run(ws_root(), structure_text, fmt, chain or "A", mode, epis, cfg,
                                   preset, settings)
    except Exception as e:
        return JSONResponse(status_code=400, content={"error": f"{type(e).__name__}: {e}"})


def _run_dir(run_id: str):
    # run_id comes from a staged manifest, but guard against path traversal regardless.
    if "/" in run_id or ".." in run_id:
        return None
    rd = os.path.join(ws_root(), "runs", run_id)
    return rd if os.path.isdir(rd) else None


@app.get("/status/{run_id}")
def status(run_id: str):
    rd = _run_dir(run_id)
    if rd is None:
        return JSONResponse(status_code=404, content={"error": f"no staged run {run_id!r}"})
    out = {"run_id": run_id, "run_dir": rd}
    spec_path = os.path.join(rd, "run_spec.json")
    if os.path.exists(spec_path):
        import json
        with open(spec_path) as fh:
            out["run_spec"] = json.load(fh)
    out["metrics_ready"] = os.path.exists(os.path.join(rd, "05_analysis", "metrics.csv"))
    # Live SLURM state is cluster-only; report it when squeue is on PATH, else say so.
    if shutil.which("squeue"):
        try:
            q = subprocess.run(["squeue", "--me", "--noheader", "--format=%i %T %j"],
                               capture_output=True, text=True, timeout=15)
            out["squeue"] = [l for l in q.stdout.splitlines() if run_id[:12] in l] or "no matching jobs"
        except Exception as e:
            out["squeue_error"] = f"{type(e).__name__}: {e}"
    else:
        out["squeue"] = "not on a SLURM host"
    return out


@app.get("/results/{run_id}")
def results(run_id: str, preset: str = Query("twelvemer"), topk: int = Query(0),
            fmt: str = Query("json")):
    rd = _run_dir(run_id)
    if rd is None:
        return JSONResponse(status_code=404, content={"error": f"no staged run {run_id!r}"})
    try:
        df, info = scoring.score_run(rd, preset=preset, topk=topk)
    except FileNotFoundError as e:
        return JSONResponse(status_code=409, content={"error": str(e), "metrics_ready": False})
    except Exception as e:
        return JSONResponse(status_code=400, content={"error": f"{type(e).__name__}: {e}"})
    if fmt == "csv":
        return Response(content=df.to_csv(index=False), media_type="text/csv")
    return {"info": info, "columns": list(df.columns),
            "rows": df.head(topk or 50).to_dict(orient="records")}
