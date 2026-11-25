# episcaf platform map — pages, features, endpoints, flow

**Keep this current.** Update it in the same change whenever a page, link, endpoint, feature, or the
engine changes. It's the single source of truth for *what the platform is and where things live*.

_Last updated: 2026-11-24 (RFD4 engine live; both design avenues wired to click-to-submit)._

---

## 1. What it is, in one paragraph

episcaf designs protein scaffolds that present a target **epitope** the way it sits natively, then
scores the designs. A user uploads a structure, picks an epitope by one of **two avenues**, clicks
submit, and a GPU campaign runs on the lab's Gemini cluster and returns designs **ranked by a
composite score** in **My runs**. Two layers: a static **`docs/` UI** and a FastAPI backend
(`episcaf_service/`) deployed on Cloud Run; the heavy GPU work runs on Gemini, fired by a data-mover
worker that polls the request bucket.

## 2. The two design avenues (same engine underneath)

| Avenue | Page | How the epitope is chosen | Clash scoring |
|---|---|---|---|
| **Find antibody contacts** | `scaffold.html` | Upload a complex → auto-detect antibody (`/inspect`) **or** mark antibody chains/residues manually → residues within the contact cutoff (4 Å) become the epitope (computed in-browser) | real antibody (strict + §3.3 shadow) once `antibody_chains` is wired through; **AOMclash** planned for the no-antibody case |
| **Find candidate regions** | `growpatch_tool.html` | Upload a lone antigen → growPatch enumerates candidate surface patches (optionally DiscoTope-ranked) → pick regions | no antibody → **AOMclash** estimate (planned; currently 0) |

Both build an epitope islands string (e.g. `A:205-216;A:278-292`) and submit it to the **same RFD4
engine**: RFD4-Proteina → clash pre-filter → AF2 (ColabFold) → composite → ranked CSV + design PDBs.

## 3. Pages (`docs/`)

| Page | H1 | Purpose | Links out to |
|---|---|---|---|
| `home.html` | Epitope scaffolds… | Landing. Two CTAs: **Find antibody contacts → `scaffold.html`**, **candidate regions → `growpatch_tool.html`** | scaffold, growpatch_tool |
| `scaffold.html` | Configure scaffold design | **Antibody-contacts avenue.** Upload + 3D viewer + sequence panel; auto-detect (`/inspect`) or manual antibody chain/residue selection; in-browser 4 Å contact epitope; **Run on Gemini** (token sign-in + submit to RFD4 pilot). Also a local "simulate / prepare inputs" path when connected to the local staging server. | growpatch_tool (candidate regions), results.html#myruns |
| `growpatch_tool.html` | Find candidate regions | **growPatch avenue.** Upload lone antigen → enumerate candidate patches (DiscoTope-rankable) → pick regions → review → **launch to the RFD4 pilot**. (`?route=antibody` is a legacy simplified antibody flow — superseded by scaffold.html for the entry links.) | scaffold.html, results.html#myruns |
| `results.html` | Results | **My runs** (token sign-in → live run list with status pills) + a results viewer: loads a run's ranked CSV, sortable/filterable table, downloads designs.zip + CSV. | growpatch_tool |
| `methods.html` | Methods | Docs / methods writeup. | aomclash, growpatch_tool |
| `pipeline.html` | Design workflow | Pipeline blueprint. | aomclash |
| `aomclash.html` | AOMclash | Explains AOMclash (predicted antibody obstruction). | scaffold.html |

Shared helpers: **`pilot.js`** (`pilotRequest`, `pilotSession` — token auth + My-runs refresh; used by
scaffold.html, growpatch_tool.html, results.html). `site.css` (styles).

## 4. Navigation flow

```
home.html ─┬─► scaffold.html  (antibody contacts) ─┐
           └─► growpatch_tool.html (candidate regions) ─┤──► Launch on Gemini
                                                        │       (POST /pilot/reviews → /pilot/runs)
                                                        └──► results.html?run=<id>  (My runs: status → download)
```

## 5. Backend endpoints (`episcaf_service/`)

**Public / local staging server** (`app.py`; also `scripts/growpatch_serve.py` for pure-local dev):
`GET /health`, `GET /reference.pdb`, `POST /inspect` (antibody/complex detection → `detect_complex`),
`POST /enumerate` (growPatch candidates), `POST /prepare|/stage|/submit`, `GET /status/{id}`,
`GET /results/{id}` — the local "stage inputs then run on the cluster yourself" path (scaffold.html's
simulate/prepare flow uses `/submit`).

**Private pilot** (`pilot.py`, prefix `/pilot`, token-gated via cookie):
`POST /session` / `DELETE /session` (sign in/out), `POST /enumerate`, `POST /reviews` (stage a run →
`cloud.review`), `POST /runs/{rid}` (launch → `cloud.launch`), `GET /runs` (My runs), `GET /runs/{rid}`
(status), `GET /runs/{rid}/results` (ranked CSV), `GET /runs/{rid}/designs.zip`, `GET /runs/{rid}/request.zip`.

Pilot modes (`EPISCAF_PILOT_MODE`): `manual-gemini` (**the RFD4 pilot** — any antigen+epitope → Gemini),
`cloud-lysozyme` (Google Batch RFd3 demo), `staging`/`simulated` (local/dev).

## 6. RFD4 engine (`episcaf_service/rfd4/`) — runs on Gemini

Port of John Altin's RFD4 pipeline (see `~/Downloads/RFD4_pipeline_description.pdf`). Modules:
`config` (tunable dials; composite is Brandon's), `topology` (islands → 103-mer layout), `geom`
(Kabsch), `clash` (strict + §3.3 shadow; AOMclash fallback planned), `score` (4 metrics + composite),
`generate` (antigen+islands → RFD4 YAML spec; PDB→CIF via atomworks), `fold` (pluggable oracle;
AF2/ColabFold now, AF3/Boltz later), `read` (RFD4 annotated-CIF + PDB readers), `campaign`
(orchestrate single-pass → ranked `final_results.csv` + `designs.zip`).

**Request lifecycle:** web `POST /pilot/reviews` → `cloud.review` stages a bundle (antigen upload +
epitope regions + `workflow=rfd4-scaffold-v1`) → data-mover worker (`gcloud/pilot_worker.py`) polls,
fetches, `sbatch`es `handoff_rfd4.sbatch` on gpu-a100 → `handoff_rfd4_job.py` runs `campaign` →
writes `retained.csv`/`designs.zip`/`terminal.json` → worker publishes to the result bucket → My runs.

## 7. Cluster / infra (Gemini)

- **gpu-a100**: where RFD4-Proteina + ColabFold run. **data-mover** (`gemini-data1`): the worker (GCS +
  tgen_labs). **login**: no tgen_labs. Tools/weights on **tgen_labs** (`…/episcaf_v2_bneff/`), home is
  quota-bound. Setup: `episcaf_service/gcloud/rfd4_setup.md`.
- Worker: scrontab `worker.scrontab` → `run_worker.sh` → `pilot_worker.py tick` every 5 min;
  `EPISCAF_SBATCH` defaults to `handoff_rfd4.sbatch`.

## 8. Deployment

- **Pilot (Cloud Run)**: in Cloud Shell, `cd ~/episcaf-manual-pilot && git pull --ff-only origin
  <branch> && bash episcaf_service/gcloud/deploy.sh`. Rebuilds the web image (UI + backend) + redeploys.
  Safe to run while a Gemini job is in flight (different layer).
- **Worker (Gemini)**: `cd ~/episcaf-advisor && git pull --ff-only origin <branch>`; next scrontab tick
  uses it. Branch currently: `rfd4-engine`.

## 9. Current vs past/deprecated

**Current:** RFD4-Proteina → AF2 → composite engine (validated on Gemini); both avenues wired to
click-to-submit; scaffold.html is the antibody page; results.html shows the composite-ranked CSV.

**Deprecated / superseded:** the 3HFM-locked RFd3 reference pilot (`handoff.sbatch`/`handoff_job.py`,
`cloud_backbones.py` RFd3 path) — replaced by the general RFD4 path; `growpatch_tool.html?route=antibody`
simplified antibody flow — superseded by scaffold.html for entry links. The RFd3→MPNN→AF3 manuscript
pipeline (`scripts/`, `episcaf_pipeline/`) is separate research code, not the web platform.

## 10. Known gaps / TODO

- Wire `antibody_chains` from the submit UI → real-antibody clash on the antibody avenue.
- **AOMclash** for the no-antibody avenue (currently clash = 0 without an antibody).
- Phase 3: seeded partial-diffusion rounds (budget + early-stop) for *passing* designs at scale.
- Output-presentation pass on results.html for the RFD4 composite schema.
- Optional: antigen upload > 103-residue epitope guardrails; AF3/Boltz oracle backends.
