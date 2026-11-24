# Claude handoff · Google Cloud demo · 2026-11-23

**Start here in Brandon’s Mac checkout** at `~/Desktop/projects/episcaf/episcaf_v2`. The implementation is partly finished and uncommitted. Cloud Shell and GitHub do not have these changes. Read root `CLAUDE.md`, `docs/HANDOFF.md`, `docs/PIPELINE.md` and `episcaf_service/gcloud/README.md`. Inspect `git status --short --branch` and the current diff before editing. At handoff, local `main` and the remote-tracking `origin/pilot/private-simulated-runs` point to `2826684`. Brandon handles deployment. No live Google Cloud or Gemini job was submitted.

## What Brandon wants to show John

John opens the private pilot and uploads an antigen alone or an antibody–antigen complex. With antigen alone, he chooses a chain, sees the structure, ranks candidate regions with DiscoTope, and picks regions. With a complex, he selects the antigen chain and inspects the antibody contact epitope. He reviews the exact count and submits once. Google Batch starts RFdiffusion3. For every selected epitope the demo returns eight backbone coordinates. John can close his browser, sign back in, and download them from Results → My runs. The demo ends there. Keep the simple form, chain choice, residue dragging and fixed design defaults. PepSeq is the broader application.

**Current coverage matters for the demo.** The web interface accepts user-provided structures. The GPU worker uses the checked-in 3HFM lysozyme reference, with antigen alone or its bound antibody. Both routes require the matching reference on chain Y. Other uploads go to lab review. Brandon asked for arbitrary uploads through generation. That end-to-end capability is outstanding. Describe the reference scope when demonstrating the implementation.

## First two blockers

1. **Finish the Cloud Build images.** Both build YAMLs now request named Docker targets. Root `Dockerfile` still has one unnamed stage, so even the manual deployment script would fail to build from this tree. Supply the lightweight `web` target and a private `web-discotope` CPU target. Keep the public/manual image lightweight. The private adapter in `episcaf_service/disco.py` expects a Python environment at `/opt/discotope/venv/bin/python` and the upstream model directory at `/opt/discotope/models`. Verify both named images build. Include a real scoring smoke test on lysozyme that checks scores map back to the original residue numbers. A successful CLI exit alone does not establish that the page ranked any candidates.
2. **Finish the private ranking route.** The latest page posts Google-hosted scoring to `/pilot/enumerate`. The server currently serves only `/enumerate`. The tester cookie is scoped to `/pilot`. Connect the private endpoint to the existing enumeration logic with server-side tester authentication and prevent anonymous inference via the old endpoint on the cloud demo. Preserve unranked scanning and the existing explicit BioLib opt-in on other deployments. Test signed-out, signed-in, and revoked-token requests against both endpoint paths. Sign-in and the latest frontend endpoint switch have yet to be tested together.

The DiscoTope source revision inspected is `35d9f2e55f97eaba2a7acefbc394db58fb9670bc`. On this Mac, temporary copies are under `/private/tmp/episcaf-discotope35` (upstream source and models.zip), `/private/tmp/episcaf-disco-verify` (disposable Python environment) and `/private/tmp/episcaf-disco-torch/hub/checkpoints` (downloaded ESM-IF1 checkpoint). Check the environment before reuse. No inference result was verified from these assets. The Linux image must build independently of them. [Upstream DiscoTope code, installation and license](https://github.com/Magnushhoie/DiscoTope-3.0/tree/35d9f2e55f97eaba2a7acefbc394db58fb9670bc)

## What is already in the diff

- `cloud.py`, `cloud_backbones.py` and `pilot.py` implement an eight-backbone budget per reference region, up to three regions per cloud run. New cloud results contain a coordinate ZIP and CSV manifest. Manual Gemini and legacy scored results have separate paths.
- `growpatch_tool.html` and `results.html` display the reviewed count, provider, job status and downloads. The private demo defaults to local ranking for the antigen route. `scaffold.html` removes tuning controls and keeps direct selection and dragging.
- `gcloud/deploy.sh cloud-batch` under `episcaf_service/` sets up CPU Cloud Run and an on-demand L4 Batch worker. The plain script invocation still selects manual Gemini. Include the currently untracked `episcaf_service/gcloud/build-web.yaml` in the eventual commit. Scripts were exercised with a stand-in for gcloud. Project permissions, quota, image startup, and wall time remain unverified.

## Verification checkpoint

- Ten cloud-backend tests and one manual-transfer test passed using local stand-ins. These cover exact per-region counts, archive contents, incomplete output, duplicate launch, owner isolation and historical downloads.
- The simulated pilot browser test passed including three regions × eight attempts, a duplicate launch, restart, fresh sign-in, failure display and owner isolation. The scaffold browser test passed for drag selection, fixed-default staging and simulation.
- The cloud browser test passed for an eight-backbone contact request, a 24-backbone antigen request and an unmatched upload sent to review **before** the private ranking route changed. The manual browser regression also passed. Storage, Batch, ranking and coordinates were stand-ins. Real Linux image builds, real local DiscoTope ranking and the first real cloud GPU result remain unverified.

Local test commands, from the checkout root on Brandon’s Mac:

```bash
export PYTHON=/Users/bneff/micromamba/envs/bcell-repair/bin/python
export PYTHONPATH=/private/tmp/episcaf-cloud-libs:/private/tmp/episcaf-pilot-gemmi
export PLAYWRIGHT_MODULE=/private/tmp/episcaf-pilot-browser/node_modules/playwright-core/index.mjs
"$PYTHON" tests/test_discotope_ranking.py
"$PYTHON" tests/test_cloud_backbones.py
"$PYTHON" tests/test_manual_transfer.py
node tests/residue_drag.mjs
node tests/pilot_browser.mjs
node tests/cloud_browser.mjs
EPISCAF_PILOT_MODE=manual-gemini node tests/cloud_browser.mjs
```

## Before handing Brandon deployment commands

Verify the two blockers, run the tests, and inspect the full staged diff. Preserve the independent edits in `manuscript/figures/FIGURES.md`, `scripts/aomclash_score_designs.py`, `scripts/stage03_mpnn_submit.py`, `scripts/fig_composite_sigmoids.py`, `scripts/fig_metric_table.py`, `tests/test_af3_score_merge.py` and `tests/test_mpnn_receipts.py`. Stage paths explicitly. Once verified, publish the pilot changes to `pilot/private-simulated-runs` and give Brandon the Cloud Shell commands in the cloud README.

Brandon has linked a credit card. His billing **account type**, remaining credit, and GPU quota have not been checked. Google requires paid billing for GPU VMs. Unexpired welcome credit can still apply after upgrade, and costs can exceed that credit. Each demo run requests one GPU with a one-hour task limit. Separate submissions can run concurrently. [Google Cloud free trial terms](https://docs.cloud.google.com/free/docs/free-cloud-features) · [Batch quotas](https://docs.cloud.google.com/batch/quotas)
