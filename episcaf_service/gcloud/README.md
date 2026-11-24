# Cloud Run pilot

Testers sign in to the private Cloud Run site. The operator selects which backend takes new requests when deploying the service. Use Google Batch for the GPU demo or manual Gemini submission for the lab workflow. The public site links to the private pilot.

## Deploy the Google Batch demo

Use the project that owns the existing `episcaf` Cloud Run service in `us-central1`. Check its billing and quotas **before building the GPU image**. A non-billable Free Trial account cannot attach a GPU to a VM or request a quota increase. After activating paid billing, unused eligible credit lasts until the original expiry date. Additional usage can charge the billing account. [Google Cloud trial terms](https://docs.cloud.google.com/free/docs/free-cloud-features#during-free-trial)

The project needs global GPU quota and regional NVIDIA L4 quota, a Batch-compatible default network, and a zone in `us-central1` with G2 capacity. The deployment identity needs permission to enable APIs, create the service accounts and grant the listed roles, build images, and deploy both web revisions. The script checks for a reported zero GPU quota before submitting the large build. [Batch setup](https://docs.cloud.google.com/batch/docs/get-started) · [Batch GPU requirements](https://docs.cloud.google.com/batch/docs/create-run-job-gpus)

Use the existing pilot checkout in Cloud Shell.

```bash
cd ~/episcaf-manual-pilot
git pull --ff-only
gcloud config get-value project
bash episcaf_service/gcloud/deploy.sh cloud-batch
```

Open the printed **Pilot** URL. Existing tester tokens work because both modes use the same private buckets. The job starts after a tester signs in, reviews a supported run and confirms submission. Cloud Build builds a small public web image, a private Cloud Run image with CPU DiscoTope, and a GPU worker from the pinned Foundry image with bundled RFdiffusion3 weights. The pilot records the resulting worker image digest in its configuration. The worker account can read the image, report its state, write logs and create result objects. The worker has no access to tester records. [Foundry container documentation](https://github.com/RosettaCommons/foundry#docker-image)

Choose a structure for the first GPU run.

- **Find candidate regions**. Download and upload the linked 3HFM lysozyme, choose chain Y and generate candidates. To rank the candidates with DiscoTope, turn ranking on before generating. On the private Google Batch demo, the CPU model scores the selected chain on Cloud Run. Select the regions you want.
- **Find antibody contacts**. Download and upload the linked 3HFM complex, choose chain Y and find its contact epitope. The contacts come from the bound antibody.

Review the total of eight backbone attempts per selected region. Open **Results → My runs** to follow the job and download the generated backbones when it completes. For two selected regions, a successful ZIP should contain 16 backbone structures. RFdiffusion3 stops at backbone generation. Other uploads support inspection and private requests for lab review.

Each launch requests one GPU task with a one-hour execution limit and zero automatic task retries. A repeated launch for the same reviewed run reuses its Batch job ID. Separate submissions can each incur charges. Build time, stored images, buckets and compute all contribute to the project bill. To inspect a run, copy its Batch job ID from My runs and use `gcloud batch jobs describe "$JOB_ID" --location=us-central1`.

To switch new requests back to the manual handoff, let active cloud jobs finish, then run `bash episcaf_service/gcloud/deploy.sh manual-gemini` in the same checkout.

## Deploy the manual Gemini pilot

Select the project that owns the existing episcaf service in us-central1. The operator needs permission to build a web image, deploy Cloud Run and provision or use the private buckets and web service account. The deployment builds a CPU web image and prints the pilot URL. Issue one token per person. Reissuing a token rotates it while preserving their runs.

```bash
git clone --single-branch --branch pilot/private-simulated-runs https://github.com/bdneff/episcaf.git episcaf-pilot
cd episcaf-pilot
gcloud config get-value project
bash episcaf_service/gcloud/deploy.sh
bash episcaf_service/gcloud/tester.sh issue john
```

For an existing pilot checkout use `git pull --ff-only` before deploying. Give John the printed pilot URL and the contents of his private token file separately. To revoke his access run `bash episcaf_service/gcloud/tester.sh revoke john`.

## Automated data-mover worker (Gemini-native)

Instead of fetching and submitting each reference run by hand, a worker on the Gemini **data mover** polls the buckets every five minutes, submits new runs to Slurm automatically, and keeps a live status in the result bucket that drives the **My runs** state machine (*queued for cluster → submitted → queued on cluster → running → completed / failed*). The data mover reaches Google Cloud Storage directly, so no laptop bridge is involved. `gcloud` runs from the `google/cloud-sdk` singularity image with a buckets-scoped service-account key.

**Gate:** approval is per tester, not per submission. Any holder of a valid token who clicks **Launch on Gemini** gets an auto-fired job within ~5 min. Each run is bounded (one GPU, one-hour cap, 16 attempts / keep ≤8) and a token revokes instantly, so revoking a tester stops further jobs.

### One-time cloud setup (Cloud Shell / laptop, in the pilot project)

Create a service account scoped to just the two buckets and download a key:

```bash
project=$(gcloud config get-value project)
number=$(gcloud projects describe "$project" --format='value(projectNumber)')
sa="episcaf-worker@$project.iam.gserviceaccount.com"
gcloud iam service-accounts create episcaf-worker --display-name="episcaf Gemini data-mover worker" || true
gcloud storage buckets add-iam-policy-binding "gs://episcaf-requests-$number" --member="serviceAccount:$sa" --role=roles/storage.objectViewer
gcloud storage buckets add-iam-policy-binding "gs://episcaf-results-$number"  --member="serviceAccount:$sa" --role=roles/storage.objectUser
gcloud iam service-accounts keys create sa-key.json --iam-account="$sa"
```

The worker reads the request bucket (to list `inbox/` and download `reference.zip`) and reads+writes the result bucket. **Security note:** `objectViewer` on the request bucket also lets the key read the hashed tester digests under `testers/`; the raw tokens are never stored, but to withhold even the digests, replace that binding with an IAM condition limited to `resource.name.startsWith("projects/_/buckets/episcaf-requests-'"$number"'/objects/inbox/")`, or stage requests in a separate inbox bucket. The key file is a long-lived credential: keep it `chmod 600`, in `$HOME` on Gemini, never on shared `/scratch` or `/tgen_labs`.

### One-time Gemini setup (on `gemini-data1`)

```bash
# Repo checkout the worker and sbatch script live in:
git clone --single-branch --branch pilot/private-simulated-runs https://github.com/bdneff/episcaf.git "$HOME/episcaf-advisor"   # or: git -C "$HOME/episcaf-advisor" pull --ff-only
mkdir -p "$HOME/.episcaf-pilot" && chmod 700 "$HOME/.episcaf-pilot"
# Copy the key here from the machine that created it, then lock it down:
#   scp sa-key.json gemini-data1:~/.episcaf-pilot/sa-key.json
chmod 600 "$HOME/.episcaf-pilot/sa-key.json"
# Pull gcloud once (no install needed):
singularity pull "$HOME/.episcaf-pilot/cloud-sdk.sif" docker://google/cloud-sdk:slim
# Fill in the bucket names:
cp "$HOME/episcaf-advisor/episcaf_service/gcloud/worker.env.example" "$HOME/.episcaf-pilot/worker.env"
chmod 600 "$HOME/.episcaf-pilot/worker.env"
# edit worker.env: set EPISCAF_REQUEST_BUCKET=episcaf-requests-$number and EPISCAF_RESULT_BUCKET=episcaf-results-$number
```

Dry-run one tick before scheduling (reads only, submits nothing):

```bash
bash "$HOME/episcaf-advisor/episcaf_service/gcloud/run_worker.sh" --dry-run
```

Then install the schedule (fires every 5 min on the `data-mover` partition; scrontab skips a tick while the previous one is still running):

```bash
scrontab "$HOME/episcaf-advisor/episcaf_service/gcloud/worker.scrontab"
scrontab -l            # confirm it is installed
```

Logs append to `~/.episcaf-pilot/worker.log`. To pause the worker, `scrontab -r`. The manual `pilot_transfer.py` flow below still works and is the fallback if the worker is paused.

## Rank candidates on an uploaded antigen

On **Find candidate regions**, upload a PDB or mmCIF, choose its antigen chain, and check **Rank by DiscoTope-3.0**. Choose the structure type and click **Generate candidates**. In the manual pilot, BioLib receives the selected protein chain. In the cloud demo, the private Cloud Run service processes it. Regions appear in generation order while scoring runs, then reorder by the sum of their within-chain normalized residue scores. If scoring fails, all candidates remain available. Scoring accepts user uploads independently of eligibility for the Gemini reference job.

## Review requests in your browser

For your own upload, open **Results → My runs → View run** to see the selected regions or **Download submitted request (ZIP)**. Existing staged requests work too. Open `request.json` inside the ZIP to inspect its manifest.

To review someone else’s request, use the download link in the operator email while signed into the Google account with bucket access. You can also browse the request bucket’s `inbox/` folder in the Cloud Storage console. To add the link to an existing email setup, replace Code.gs with the updated `notify.gs` and save it.

## Where each request goes

1. **Cloud Run** writes an immutable private bundle to `gs://episcaf-requests-PROJECT_NUMBER/inbox/OWNER/RUN_ID/`. It contains the original structure and a JSON manifest of its hash, route, chain, regions, budget and owner. For reference runs it also contains the frozen contigs, per-contig trial counts and a copy of the batch script. A tester's retry returns the same request. A reference run appears in My runs as **queued for cluster** until the data-mover worker picks it up. Other structures appear as **awaiting lab review**.
2. **Brandon’s laptop** uses his normal Google Cloud login and SSH alias to transfer the bundle into `$HOME/.episcaf-pilot/inbox/RUN_ID/request.zip` on Gemini. The transfer needs the lab network or VPN, gcloud with access to the private buckets and a current copy of this branch. The login node only needs the home directory. Inspect the manifest there with `unzip -p request.zip request.json`.
3. **The compute node** writes the reference bundle, intermediates, scores and backbone coordinates under `$WS/pilot/requests/RUN_ID/`. The configured workspace is `/tgen_labs/altin/alphafold3/workspace/episcaf_v2_bneff`. Slurm writes its log in the home inbox. The job copies the small final CSV, ZIP and completion record back to that inbox for transfer. A second job for the same ID exits before inference.
4. **Brandon** publishes the numeric Slurm receipt from his laptop and returns the result files after completion. My runs then offers the scores and structures to the owner. The manual pilot shows the last published status. Closing the browser has no effect on the job.

## Email Brandon when a tester submits

Use a private Apps Script project under a lab-approved Google account with read access to the request bucket. Set its Code.gs to the contents of `notify.gs`. In Project Settings, show the manifest and replace `appsscript.json` with the supplied manifest. Set the Script Properties `REQUEST_BUCKET` (bucket name only) and `OPERATOR_EMAIL` (your recipient). Run `notifyOperator` once to authorize it, then add a **time-driven, every 5 minutes** trigger for that function. The script emails the tester name, run ID, budget and transfer/submission commands. The original structure and token stay out of email. Stored run IDs suppress repeat notifications. An interrupted execution just after mail delivery can repeat the email. All commands refer to the same run.

If your Workspace administrator disables Apps Script, the private inbox is also available from a signed-in terminal. This command lists submitted bundles, including those submitted during a notification outage.

```bash
number=$(gcloud projects describe "$(gcloud config get-value project)" --format='value(projectNumber)')
gcloud storage ls --recursive "gs://episcaf-requests-$number/inbox/"
```

## First Gemini submission

For a fresh checkout, run this once on Gemini:

```bash
git clone --single-branch --branch pilot/private-simulated-runs https://github.com/bdneff/episcaf.git "$HOME/episcaf-advisor"
```

Once, on Gemini, install the same pilot branch in `$HOME/episcaf-advisor` (or update that checkout with `git pull --ff-only`). The batch script uses the established Mamba module and RFd3 environment. For the first test John downloads the linked 3HFM lysozyme reference in either route and selects one region on chain Y. Each job requests one GPU for up to one hour and attempts 16 backbones. It returns up to eight structures using the lowest half of scorable raw occupancy values. Other uploads produce a review bundle only. The standalone backbone filter has no validated absolute cutoff.

The notification supplies OWNER and RUN_ID. On your laptop, in the pilot checkout, run:

```bash
python3 episcaf_service/gcloud/pilot_transfer.py fetch "$OWNER" "$RUN_ID"
```

On Gemini, use the run ID from the notification. Inspect the manifest, then submit once:

```bash
cd "$HOME/.episcaf-pilot/inbox/$RUN_ID"
unzip -p request.zip request.json
sbatch --parsable "$HOME/episcaf-advisor/episcaf_service/handoff.sbatch"
```

If the submission response is lost, check Slurm before submitting again. On your laptop, record the numeric receipt and publish after successful completion:

```bash
python3 episcaf_service/gcloud/pilot_transfer.py receipt "$OWNER" "$RUN_ID" "$JOB_ID"
python3 episcaf_service/gcloud/pilot_transfer.py publish "$OWNER" "$RUN_ID" "$JOB_ID"
```

If Slurm reports failure, inspect the job log. To return the explanation use `pilot_transfer.py fail "$OWNER" "$RUN_ID" "$JOB_ID" --reason "…"`. For other uploaded structures use `pilot_transfer.py review "$OWNER" "$RUN_ID"` to retrieve the request for inspection.

Project IAM, institutional email delivery, installed RFd3 model execution and GPU runtime require verification under the lab accounts. No always-running Gemini process or Open OnDemand app is needed.
