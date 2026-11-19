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
| POST | `/enumerate?fmt=&chain=&top_k=&…` | growPatch candidate epitopes, ranked by mean DiscoTope when a scorer is configured |
| POST | `/stage?fmt=&chain=&epitopes=a\|b&…` | build RFd3 contigs for the chosen epitopes → `{run_spec, contigs_csv, n_contigs}` (download; no `$WS` write) |
| POST | `/submit?fmt=&chain=&mode=&name=&epitopes=a\|b&…` | stage a run under `$WS` → `{run_id, run_dir, submit_sh, …}` |
| GET  | `/status/{run_id}` | staged spec, whether metrics exist, live `squeue` (on a SLURM host) |
| GET  | `/results/{run_id}?preset=&topk=&fmt=json\|csv` | binding-model + composite score of the finished metrics |

`/submit` query params mirror the GUI settings: `rfd3`, `mpnn`, `helix_bias`, `variants`, `total_len`,
`ipd`, `ipdrounds`, `ipd_mode`, `prefilter`, `keepfrac`, `topn`, `group`, `preset`
(default `antibody_softgate` for a complex, `twelvemer` for a lone antigen).

## Deploy (you run these on the cluster)
The service also **serves the GUI itself** (StaticFiles over `./docs`), so the whole tool runs
same-origin on one host — no CORS, no "backend URL" field, no separate web host.
```bash
conda activate ~/rfd3/env/rfd3_py312
pip install -r episcaf_service/requirements.txt
# optional: DiscoTope ranking of uploaded antigens (runs on BioLib's cloud, no local GPU)
pip install pybiolib && export EPISCAF_DISCOTOPE=biolib
cd <repo root>                     # imports are repo-relative; StaticFiles serves ./docs
export EPISCAF_WS=/tgen_labs/altin/alphafold3/workspace/episcaf_v2_bneff   # optional; defaults to this
uvicorn episcaf_service.app:app --host 127.0.0.1 --port 8808
```
Then from your laptop open a tunnel and browse to the tool over it (it is `localhost` to the browser,
so the page's connection check passes and its calls stay same-origin):
```bash
ssh -L 8808:localhost:8808 <gemini-host>
# then open in your browser:
http://localhost:8808/growpatch_tool.html
```
Uploading a structure runs real growPatch enumeration (and DiscoTope ranking, if enabled) on the
cluster; **Prepare design inputs** calls `/stage` and downloads the contigs + run spec. To persist a
run under `$WS` and get `submit.sh`, call `/submit`; run that script on the cluster, then read results
back with `/status/{run_id}` and `/results/{run_id}` once AF3 finishes.

## DiscoTope scoring (candidate ranking)
Ranking an uploaded antigen needs per-residue DiscoTope-3.0 scores, which are **not** computed in the
base service. The provider (`episcaf_service/disco.py`) is off by default — uploads then show
"DiscoTope scores unavailable" and the full candidate set, unranked. Turn it on with either:
- `EPISCAF_DISCOTOPE=biolib` (+ `pip install pybiolib`) — runs DiscoTope-3.0 on BioLib's cloud for each
  upload (needs outbound internet from the host; no local GPU).
- `EPISCAF_DISCOTOPE_SCORED_PDB=<path>` — read scores from a pre-baked DiscoTope PDB's CA B-factor
  column (for a fixed antigen, or a test fixture; `EPISCAF_DISCOTOPE_SCALE` defaults to 100).

Scores are min-max normalized within the chain and are never read from an uploaded structure's own
B-factors. The recorded examples (`docs/data/growpatch_*.js`) carry their scores already.

## What's staged into `$WS/runs/<run_id>/`
`input.pdb` · `designs.csv` (RFd3 contigs for the chosen epitopes) · `run_spec.json` · `submit.sh`
(the verified command sequence: `ipd_round.py init` → RFd3 array → optional AOMclash prefilter →
MPNN/AF3 → `aomclash_score_designs.py` → `stage06_select.py`).

### Storage split (durable records, disposable GPU outputs)
`submit.sh` begins by sending the large GPU outputs to scratch, which self-deletes when swept, while
the run record and the scored metrics stay durable:
- **Durable** (under `$RUN_DIR` on `$WS`): `input.pdb`, `designs.csv`, `run_spec.json`, `submit.sh`,
  and `05_analysis/` (`metrics.csv`, `ranked.csv`) — the things you keep.
- **Scratch** (`$EPISCAF_SCRATCH/<run_id>/`, default `/scratch/$USER/episcaf/<run_id>/`): `02_rfd3`,
  `03_mpnn`, `04_af3`, symlinked in from `$RUN_DIR` so paths still resolve.

So keep code + env in your home dir, let the service write records to `$WS`, and the bulky
RFd3/MPNN/AF3 outputs land on scratch. **Extract what you need from `05_analysis` before scratch is
purged** — only `02/03/04` are disposable (this is why the old runs needed `AF3_REMAP`).

**If the host running the service can't see `$WS`.** The tgen_labs workspace is not mounted on the
login node, so a service running there cannot write records to it (`/submit` would 400; `/enumerate`
and `/stage` are unaffected — they never touch `$WS`). Point `EPISCAF_WS` at a path visible from both
the service host and the compute nodes — your home dir:
```bash
export EPISCAF_WS=$HOME/episcaf_runs
```
Home is NFS-mounted everywhere, so the service writes the record there and the GPU job reads
`designs.csv` / writes `05_analysis` back to the same durable spot. Only the big `02/03/04` outputs go
to scratch. Copy records/metrics to tgen_labs later, from a node that can see it, if you want.

## Public hosting, locked to your group
GitHub Pages can't run the backend, and the backend already serves the UI, so you don't need Pages.
To put the whole tool on the internet for your group:

1. **Expose it with a Cloudflare Tunnel** — outbound from the host, so no cluster inbound-firewall
   change. `cloudflared tunnel --url http://127.0.0.1:8808` gives a public HTTPS URL for a quick test;
   a named tunnel on your own domain (e.g. `episcaf.yourdomain.com`) is the durable form. Confirm your
   cluster permits running `cloudflared`.
2. **Require login with Cloudflare Access** in front of the tunnel: a Zero-Trust application that
   allows only your group's emails/domain. Login happens at Cloudflare's edge before any request
   reaches the backend — real identity auth, no app code, no API key to leak.
3. **Rate limiting, two layers.** Cloudflare's own rate-limit rules at the edge, plus the app's
   `EPISCAF_RATE_LIMIT` (slowapi — `pip install slowapi`; default `30/minute` per client IP) on the
   expensive endpoints. Use both: behind a proxy the app sees Cloudflare's IP, so the *edge* rule is
   what throttles individual users; the app limit is a backstop.

The frontend enables itself whenever the backend answers `/health`, at any origin (localhost tunnel or
a hosted URL); a static host with no backend stays disabled and points at the recorded example.

**Simpler interim (no public exposure):** keep the SSH-tunnel model — group members `ssh -L 8808:...`
(cluster SSH access is itself the auth) and open `http://localhost:8808`. One shared uvicorn serves
everyone; each person runs their own tunnel. No domain, no Cloudflare, nothing public.

## Not in v1 (deferred)
Auto-submitting GPU jobs; an in-app login/token (Cloudflare Access covers group auth at the edge).
