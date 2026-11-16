# episcaf cluster service

The real backend behind the GUI tool. The static page (GitHub Pages / local server) can't run Python
and the design spine (RFd3 → ProteinMPNN → AlphaFold3) is GPU/SLURM work, so this small FastAPI
service runs on a cluster-access host and does everything up to the GPU boundary **for real**:
epitope detection, growPatch enumeration, contig staging, and scoring. It stops at the GPU: `/submit`
stages a run under `$WS` and hands back the exact submit script — **you** launch the GPU jobs.

## Why it's shaped this way
- **Private SSH tunnel, not public.** v1 is reachable only through `ssh -L`, so there's no TLS, CORS
  policy, API key, or rate-limiting to get right yet. Public exposure + a private key + per-key
  rate-limiting is the deferred hardening step (see memory `webtool-server-future`).
- **Stage, don't auto-submit.** The service never runs an `sbatch`. It writes the run and prints the
  command; you stay in control of GPU spend. (Claude also never submits — the Gemini boundary.)
- **Durable artifacts under `$WS`, never `/scratch`.** `/scratch` is ephemeral.

## Endpoints (wire-compatible with the local stdlib server)
Structure = raw request body; everything else = query params.

| method | path | does |
|---|---|---|
| GET  | `/health` | liveness + resolved `$WS` + whether SLURM is on PATH |
| POST | `/inspect?fmt=&chain=` | antibody-present vs lone-antigen + the contact epitope |
| POST | `/enumerate?fmt=&chain=&methods=&…` | growPatch candidate epitopes |
| POST | `/submit?fmt=&chain=&mode=&name=&epitopes=a\|b&…` | stage a run → `{run_id, run_dir, submit_sh, …}` |
| GET  | `/status/{run_id}` | staged spec, whether metrics exist, live `squeue` (on a SLURM host) |
| GET  | `/results/{run_id}?preset=&topk=&fmt=json\|csv` | binding-model + composite score of the finished metrics |

`/submit` query params mirror the GUI settings: `rfd3`, `mpnn`, `helix_bias`, `variants`, `total_len`,
`ipd`, `ipdrounds`, `ipd_mode`, `prefilter`, `keepfrac`, `topn`, `group`, `preset`
(default `antibody_softgate` for a complex, `twelvemer` for a lone antigen).

## Deploy (you run these on the cluster)
```bash
conda activate ~/rfd3/env/rfd3_py312
pip install -r episcaf_service/requirements.txt
cd <repo root>                     # imports are repo-relative
export EPISCAF_WS=/tgen_labs/altin/alphafold3/workspace/episcaf_v2_bneff   # optional; defaults to this
uvicorn episcaf_service.app:app --host 127.0.0.1 --port 8808
```
Then from your laptop open a tunnel and point the GUI at it:
```bash
ssh -L 8808:localhost:8808 <gemini-host>
```
In the GUI tool's **Backend URL** field enter `http://127.0.0.1:8808`. Uploading a structure and
clicking run now POSTs to `/submit` and shows the real `run_dir` + `submit.sh`; run that script on the
cluster, then use **Fetch results** (`/results/{run_id}`) once AF3 finishes.

## What's staged into `$WS/runs/<run_id>/`
`input.pdb` · `designs.csv` (RFd3 contigs for the chosen epitopes) · `run_spec.json` · `submit.sh`
(the verified command sequence: `ipd_round.py init` → RFd3 array → optional AOMclash prefilter →
MPNN/AF3 → `aomclash_score_designs.py` → `stage06_select.py`).

## Not in v1 (deferred)
Public endpoint, TLS/CORS policy, private API key + rate-limiting, auto-submitting GPU jobs.
