# episcaf platform map — pages, features, endpoints, flow

**Keep this current.** Update it in the same change whenever a page, link, endpoint, feature, or the
engine changes. It's the single source of truth for *what the platform is and where things live*.

_Last updated: 2026-11-25 (RFD4 engine live; scaffold-count budget; benchmark timing in hand)._

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
| `scaffold.html` | Configure scaffold design | **Antibody-contacts avenue.** Upload + 3D viewer + sequence panel; auto-detect (`/inspect`) or manual antibody chain/residue selection; in-browser contact epitope with adjustable **cutoff** + **gap-fill** controls; **Run on Gemini** (token sign-in + submit to RFD4 pilot). Also a local "simulate / prepare inputs" path when connected to the local staging server. | growpatch_tool (candidate regions), results.html#myruns |
| `growpatch_tool.html` | Find candidate regions | **growPatch avenue.** Upload lone antigen → enumerate candidate patches (DiscoTope-rankable) → pick regions (each checked candidate shows as a
removable cell in the "Prepare scaffold-design inputs" box) → review → **launch to the RFD4 pilot**. (`?route=antibody` is a legacy simplified antibody flow — superseded by scaffold.html for the entry links.) | scaffold.html, results.html#myruns |
| `results.html` | Results | **My runs** (token sign-in → live run list with status pills) + a results viewer: loads a run's ranked CSV, sortable/filterable table, downloads designs.zip + CSV. For an RFD4 run, **clicking a design row** lazy-fetches + unzips designs.zip (JSZip) and opens a docked **3D structure viewer** (`designview.js`: 3Dmol, scaffold/native/epitope/antibody/steric-clash/AOM toggles, pan pad). | growpatch_tool |
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
`generate` (antigen+islands → RFD4 YAML spec; PDB→CIF via atomworks), `fold` (pluggable
**structure-prediction step**; AF2/ColabFold now, AF3/Boltz later), `read` (RFD4 annotated-CIF + PDB
readers), `timing` (per-stage wall-clock → GPU-hours per task type), `viewer` (per design: align the
AF2 model onto the native epitope frame → scaffold PDB + config JSON, with one shared `native.pdb`
(antigen + antibody, Fab auto-detected by contact) → `viewer/` in designs.zip; epitope RMSD reproduces
the CSV's), `campaign` (orchestrate single-pass → ranked `final_results.csv` + `designs.zip` +
`timings.json` + `viewer/`).

**Request lifecycle:** web `POST /pilot/reviews` → `cloud.review` stages a bundle (antigen upload +
epitope regions + `workflow=rfd4-scaffold-v1`) → data-mover worker (`gcloud/pilot_worker.py`) polls,
fetches, `sbatch`es `handoff_rfd4.sbatch` on gpu-a100 → `handoff_rfd4_job.py` runs `campaign` →
writes `retained.csv`/`designs.zip`/`terminal.json` → worker publishes to the result bucket → My runs.

**Scaffold-count budget (the user knob):** the submit UI asks for **N scaffolds** (`num_replicates` in
the request manifest; default **32**, UI cap **128**). N is the number carried into the
structure-prediction step — AF2, ~95% of run cost — *not* the number of RFD4 designs made. RFD4
generation is cheap, so backbone-clashing designs don't count against N: `campaign.run_campaign`
over-generates in waves (`wave_size`) applying the backbone pre-filter, until N clash-free backbones
exist or **128 RFD4 tries** are spent, then `select_scaffolds` hands AF2 exactly N — clash-free first,
filled with the lowest-clash of the rest if fewer than N pass. So cost scales with N (AF2 runs), while
clash-filtering is free. `campaign.json` records `target_scaffolds` / `n_generated` /
`n_passed_prefilter` / `n_scored`; the two helpers are pure and unit-tested (`tests/test_rfd4_budget.py`).

**Benchmark (2026-11-25, 1× A100-40GB, 32 designs serial):** 85.8 min wall = 1.43 GPU-hr. Split:
**structure prediction 95.5%** (1.37 GPU-hr, 154 s/design mean, one 1064 s outlier), generation 3.6%
(one 188 s batched call), prefilter + scoring + startup <1%. Composite recomputes to the stored value
exactly on all 32 real rows (parity gate passed on real data). → fold is the whole cost and is
per-design independent, so the SLURM-array fan-out (§10) is the high-leverage parallelization.

## 7. Cluster / infra (Gemini)

- **gpu-a100**: where RFD4-Proteina + ColabFold run. **data-mover** (`gemini-data1`): the worker (GCS +
  tgen_labs). **login**: no tgen_labs. Tools/weights on **tgen_labs** (`…/episcaf_v2_bneff/`), home is
  quota-bound. Setup: `episcaf_service/gcloud/rfd4_setup.md`.
- Worker: scrontab `worker.scrontab` → `run_worker.sh` → `pilot_worker.py tick` every 5 min;
  `EPISCAF_SBATCH` defaults to `handoff_rfd4.sbatch`.

## 8. Deployment

- **Single Cloud Run service** (`episcaf-pilot`): in Cloud Shell, `cd ~/episcaf-manual-pilot && git
  pull --ff-only origin <branch> && bash episcaf_service/gcloud/deploy.sh`. Builds + deploys the one
  service (UI + backend). Safe to run while a Gemini job is in flight (different layer). **The pilot URL
  is the only entry point** — the former public `episcaf` landing service was retired (2026-11-24).
- **Worker (Gemini)**: `cd ~/episcaf-advisor && git pull --ff-only origin <branch>`; next scrontab tick
  uses it. Branch currently: `rfd4-engine`.

## 9. Current vs past/deprecated

**Current:** RFD4-Proteina → AF2 → composite engine (validated on Gemini); both avenues wired to
click-to-submit; scaffold.html is the antibody page; results.html shows the composite-ranked CSV.

**Deprecated / superseded:** the 3HFM-locked RFd3 reference pilot (`handoff.sbatch`/`handoff_job.py`,
`cloud_backbones.py` RFd3 path) — replaced by the general RFD4 path; `growpatch_tool.html?route=antibody`
simplified antibody flow — superseded by scaffold.html for entry links; the public `episcaf` landing
service (mode=staging, the "Open the private pilot" invite) — retired, so there is one URL. The RFd3→MPNN→AF3 manuscript
pipeline (`scripts/`, `episcaf_pipeline/`) is separate research code, not the web platform.

## 10. Known gaps / TODO

In-flight this round (engine done where noted; wiring outward):
- **Scaffold-count knob** — engine done (`target_scaffolds`/budget routine, tested). TODO: surface the
  N input (default 32, max 128) in the submit UI → carry as `num_replicates` through the request bundle.
- **AF2-fold parallelization** — in-job concurrent folding is DONE: `campaign.fold_selected` folds the
  selected N across the job's GPUs (one colabfold per GPU, pinned via `CUDA_VISIBLE_DEVICES`),
  thread-safe timing, serial by default. Enable by setting `--gres=gpu:N` + `EPISCAF_FOLD_WORKERS=N` in
  `handoff_rfd4.sbatch`; the worker still tracks the one job. The ColabFold MSA-vs-model split now
  records (timing fix). **Remaining:** cross-node **SLURM array** for scaling beyond one node's GPUs
  (gen → array fold+score over the `select_scaffolds` output → gather) — needs the gpu-a100 node GPU
  count to choose in-job-N vs array. Memory `af2-parallelization-plan`.
- **Step-aware status UX — DONE (pending a cluster run to watch it)**: the campaign emits progress
  ("generating scaffolds" → "predicting structures k/N") via a `progress` callback; `handoff_rfd4_job`
  writes it to `inbox/progress.json`; the worker surfaces it as the My-runs stage while the job runs.
  results.html already renders `run.stage`. Engine + worker compile-clean; validate on a live run.
- **Per-design viewer — DONE** (`rfd4/viewer.py` emits `viewer/` into designs.zip; `docs/designview.js`
  docked 3Dmol panel wired into results.html on row-click). Verified on the benchmark: epitope RMSD
  reproduces the CSV exactly, Fab auto-detected (L+H). Memory `results-viewer-plan`.

New, from John's updated spec (commit 9a38b52, §1.1/§1.4):
- **Antibody auto-detection — DONE (contact fallback)**: when `antibody_chains` is empty, the campaign
  detects the Fab from the complex (`viewer.detect_antibody_chains`: ≥10 heavy-atom pairs <4.5 Å of the
  epitope) and scores clashes against it; a true lone antigen detects nothing → no-antibody path. So
  scored clashes match the viewer. `tests/test_rfd4_viewer.py`. **Future:** SAbDab-annotation-first +
  the VL/VH partner rule + `require_antibody_detection` (needs the SAbDab table).
- **Compact input — DONE** (`generation.compact_input`, auto for multi-chain / high-resnum epitopes):
  `generate.write_compact_antigen` feeds RFD4 a single-chain antigen with the epitope renumbered from 1,
  gap 9 between islands; scoring still uses the original native (motif map remaps by order).
  `tests/test_rfd4_compact.py`.

Carried:
- **AOMclash no-antibody term — DONE**: when no Fab is detected, the campaign reads the vendored
  occupancy map (`episcaf_analysis/aomclash`, numpy/scipy + `aom_map.pkl`, ~3 ms/design) and computes
  `aom_s_p1_free` per design (lower = more accessible); it **replaces the clash count in the composite**
  (same 0.45 weight, John §3.5), so lone-antigen designs are ranked by accessibility (aom 0.26→comp 1.02,
  146→0.37). Surfaced as a results column. Remaining: the real AOM occupancy shells in the viewer mesh
  (currently a sphere). A future calibration (isotonic s_p1_free→clash-count) could put it on a count
  scale, but the raw score ranks identically (monotone).
- Phase 3: seeded partial-diffusion rounds. **Decision logic DONE + tested** (`rfd4/rounds.py`:
  `select_seeds` / `should_stop` / `plan_children` / `dedupe_children`, dials in GenerationConfig;
  `tests/test_rfd4_rounds.py`). **Remaining:** the round loop that generates children via RFD4 partial
  diffusion (`task=existing-structure`, noise 0.3/0.5/0.8) and folds them — needs the partial-diffusion
  generation spec + a cluster run. §6 eval panel shows 9/24 epitopes only pass *with* rounds; §6.1 is a
  24-epitope validation set (8PWW: 70 passing, best composite 2.94).
- Optional: antigen upload > 103-residue epitope guardrails; AF3/Boltz predictor backends.
