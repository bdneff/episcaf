# episcaf cluster service

The real backend behind the GUI tool. The static page (GitHub Pages / local server) can't run Python
and the design spine (RFd3 → ProteinMPNN → AlphaFold3) is GPU/SLURM work, so this small FastAPI
service runs on a cluster-access host and does everything up to the GPU boundary **for real**:
epitope detection, growPatch enumeration, contig staging, and scoring. It stops at the GPU: `/submit`
stages a run under `$WS` and hands back a submission template — **you** launch the GPU jobs.

## Staged Gemini pilot

The hosted pilot saves a private bundle and notifies Brandon when configured. He submits from his Gemini home directory. The compute node writes the durable run to the lab workspace. [Deploy the web service, enable email and submit a request](gcloud/README.md).

## Private simulated pilot

Use the web service Python environment. From the repository root, create a private local directory outside the checkout and issue a token for each tester. Keep the same database path across service restarts. The example writes the token to a file readable by the operator.

```bash
umask 077
mkdir -p "$HOME/.episcaf-pilot"
export EPISCAF_PILOT_DB="$HOME/.episcaf-pilot/runs.sqlite"
export EPISCAF_PILOT_MODE=simulated
python -m episcaf_service.pilot issue alice > "$HOME/.episcaf-pilot/alice.token"
python -m uvicorn episcaf_service.app:app --host 127.0.0.1 --port 8808
```

Open the candidate tool at `http://127.0.0.1:8808/growpatch_tool.html`. Upload a PDB, choose its antigen chain, generate candidates, check three regions, sign in with the tester token and choose **Prepare design inputs**. The frozen review shows eight backbone trials for each region, 24 overall, split across the contigs actually generated. **Start test run** launches one simulated job for that review and opens its status on **Results**. A new review can launch another run.

The service polls the simulated scheduler every second. Return to **Results → My runs** in a fresh browser, sign in with the same token and open the completed result. The CSV records the trial budget of each selected region. The simulator produces no structural measurements. To exercise failure, restart with `EPISCAF_SIM_FAIL=1` and launch a new test run. Already submitted jobs retain their outcome.

Issue a new token for the same owner to rotate it. Run `python -m episcaf_service.pilot revoke alice` to revoke access, including existing sessions. Owner names identify people across rotations. The service compares token digests on every private request. The browser sends the sign-in header and the path-scoped HttpOnly session cookie only to the origin serving the pilot. Use HTTPS or a loopback tunnel for private access. Keep the database and its backups on a private local filesystem. It contains the reviewed uploads and results. Cloud Run (`K_SERVICE`) and the default service configuration expose staging only. The legacy staging API uses a separate artifact namespace.

### Reproduce the browser acceptance test

From the repo root, run `python tests/test_service.py` in the service environment. For the browser test, install Chrome and Playwright Core, then run the following. The test starts its own loopback services, uses a temporary database and strips the partner chains from the checked-in 8PWW fixture. DiscoTope is disabled for this test.

```bash
npm install --prefix /tmp/episcaf-browser playwright-core
PYTHON="$(command -v python)" PLAYWRIGHT_MODULE=/tmp/episcaf-browser/node_modules/playwright-core/index.mjs node tests/pilot_browser.mjs
```

### Private Gemini pilot

The advisor setup guide covers the separate private service, tester tokens, the 3HFM lysozyme worker, and private lab review of other PDB and mmCIF uploads. [Gemini setup and first run](GEMINI_PILOT.md)

The legacy generated `submit.sh` still has a commented RFd3 submission, an AF3 array placeholder, and a manual retained-backbone step. The private lysozyme worker runs its own stages within a single Slurm allocation. The bare-backbone cutoff remains an open scientific decision.

## Staging endpoints

The CPU endpoints can run on a local, hosted or cluster-access service. The private pilot above uses its own database and endpoints under `/pilot`. Legacy `/submit` writes staging artifacts under `$WS` and has no tester ownership model. Use the private pilot for tester runs.

| method | path | does |
|---|---|---|
| GET  | `/health` | liveness + resolved `$WS` + whether SLURM is on PATH |
| POST | `/inspect?fmt=&chain=` | antibody-present vs lone-antigen + the contact epitope |
| POST | `/enumerate?fmt=&chain=&top_k=&disco=&…` | growPatch candidate epitopes. `disco=0` returns them fast in generation order (unranked); `disco=1` (default) ranks by **total** DiscoTope when a scorer is configured. The GUI calls it twice — unranked first for instant display, then ranked in the background |
| POST | `/prepare` | JSON filename, chain and selected regions → contigs and the fixed trial allocation for review |
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
configured host. In staging mode **Prepare design inputs** calls `/prepare` and offers both downloads. To persist a
run under `$WS` and get `submit.sh`, call `/submit`; complete the template for your cluster, then read results
back with `/status/{run_id}` and `/results/{run_id}` once AF3 finishes.

## DiscoTope scoring (candidate ranking)
Ranking an uploaded antigen needs per-residue DiscoTope-3.0 scores, which are **not** computed in the
base service. The provider (`episcaf_service/disco.py`) is off by default — uploads then show
"DiscoTope scores unavailable" and the full candidate set, unranked. Turn it on with either:
- `EPISCAF_DISCOTOPE=biolib` (+ `pip install pybiolib`) — runs DiscoTope-3.0 on BioLib's cloud for the selected protein chain when the researcher enables ranking
  (needs outbound internet from the host; no local GPU).
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

## Always-on public deploy (no terminal to babysit)
Candidate generation needs no GPU and no cluster, so the web tool doesn't belong on Gemini at all. Host
it on a long-running platform (Render / Fly.io / Railway) for a stable public link; keep Gemini only for
the actual GPU design runs. Note: GitHub Pages can't run the Python, and Vercel's serverless timeouts
(10-60 s) are shorter than a ~1-minute enumeration, so neither fits — you need a persistent server.

The repo ships a Render blueprint (`render.yaml`) and a self-contained deps file
(`episcaf_service/requirements-web.txt`, verified in a clean venv). To deploy:
1. Push the repo to GitHub.
2. At render.com (free, sign in with GitHub) → **New +** → **Blueprint** → pick the repo. Render reads
   `render.yaml`, installs the deps, and runs `uvicorn episcaf_service.app:app --host 0.0.0.0 --port $PORT`.
3. It gives you a URL like `https://episcaf.onrender.com`; share `…/growpatch_tool.html`.

The app serves its own UI, so there's no separate frontend to host. DiscoTope ranking works here (the host
has internet + `pybiolib`), so uploads come back ranked. Notes: the free plan sleeps after ~15 min idle
(≈30-60 s cold start on the next hit; `starter` stays warm); `/submit` writes to ephemeral `/tmp` and its
`submit.sh` references cluster paths, so run real GPU jobs on Gemini, not here. The link is public — lock
it to your group later with Cloudflare Access or a password if needed.

### Google Cloud Run (recommended if you have GCP credits)
Cloud Run scales to zero (idle = $0, so the $300 free trial lasts a long time) and allows long requests
(up to 60 min), so the ~1-minute growPatch enumeration fits. The repo ships a `Dockerfile` and
`.dockerignore`; deploy entirely from **Cloud Shell** (a browser terminal — nothing to install locally):

```bash
# in the GCP console: create/select a project and enable the $300 free-trial billing, then open Cloud Shell:
gcloud config set project <PROJECT_ID>
gcloud services enable run.googleapis.com cloudbuild.googleapis.com artifactregistry.googleapis.com
git clone https://github.com/bdneff/episcaf.git && cd episcaf
gcloud run deploy episcaf \
    --source . \
    --region us-central1 \
    --allow-unauthenticated \
    --memory 2Gi --cpu 2 --timeout 300
```
Cloud Build reads the `Dockerfile`, builds the image, and deploys. It prints a URL like
`https://episcaf-xxxxx.us-central1.run.app`. Check `/health`, then share `…/growpatch_tool.html`.
`--allow-unauthenticated` makes it public for testing; lock it to your group later with Cloud Run IAM or
Identity-Aware Proxy. Set env vars with `--set-env-vars` (the Dockerfile already defaults
`EPISCAF_DISCOTOPE=biolib`, `EPISCAF_RATE_LIMIT=20/minute`).

`EPISCAF_RATE_LIMIT` is per **real** client IP: the limiter keys on the right-most `X-Forwarded-For`
entry (the IP Cloud Run's frontend observed), not the socket peer (Google's frontend) and not the
spoofable left-most entry — so one abuser can't share-bucket everyone or rotate fake IPs to bypass it.
This is an interim abuse throttle; per-user access control is the auth step (Cloud Run IAM / IAP),
still deferred. If you ever put a sanitizing load balancer in front, revisit which XFF entry to trust
(a GCLB appends `client, lb`, making the client the second-from-right).

## Not in v1 (deferred)
Real GPU submission and Slurm accounting are pending.
