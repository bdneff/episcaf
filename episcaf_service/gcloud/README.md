# Manual Gemini pilot

Testers sign in to episcaf on Cloud Run. Both epitope routes accept private requests. Brandon reviews the request, transfers it and submits it from his own Gemini login. The website has no cluster credentials. The public episcaf service stays in staging mode.

## Deploy the web service from Cloud Shell

Select the project that owns the existing episcaf service in us-central1. The operator needs permission to build a web image, deploy Cloud Run and provision or use the private buckets and web service account. The deployment builds a CPU web image and prints the pilot URL. Issue one token per person. Reissuing a token rotates it while preserving their runs.

```bash
git clone --single-branch --branch pilot/private-simulated-runs https://github.com/bdneff/episcaf.git episcaf-pilot
cd episcaf-pilot
gcloud config get-value project
bash episcaf_service/gcloud/deploy.sh
bash episcaf_service/gcloud/tester.sh issue john
```

For an existing pilot checkout use `git pull --ff-only` before deploying. Give John the printed pilot URL and the contents of his private token file separately. To revoke his access run `bash episcaf_service/gcloud/tester.sh revoke john`.

## Rank candidates on an uploaded antigen

On **Find candidate regions**, upload a PDB or mmCIF, choose its antigen chain, and check **Rank by DiscoTope-3.0**. Choose the structure type and click **Generate candidates**. BioLib receives the selected protein chain. Regions appear in generation order while scoring runs, then reorder by the sum of their within-chain normalized residue scores. If scoring fails, all candidates remain available. Scoring accepts user uploads independently of eligibility for the Gemini reference job.

## Review requests in your browser

For your own upload, open **Results → My runs → View run** to see the selected regions or **Download submitted request (ZIP)**. Existing staged requests work too. Open `request.json` inside the ZIP to inspect its manifest.

To review someone else’s request, use the download link in the operator email while signed into the Google account with bucket access. You can also browse the request bucket’s `inbox/` folder in the Cloud Storage console. To add the link to an existing email setup, replace Code.gs with the updated `notify.gs` and save it.

## Where each request goes

1. **Cloud Run** writes an immutable private bundle to `gs://episcaf-requests-PROJECT_NUMBER/inbox/OWNER/RUN_ID/`. It contains the original structure and a JSON manifest of its hash, route, chain, regions, budget and owner. For reference runs it also contains the frozen contigs, per-contig trial counts and a copy of the batch script. A tester's retry returns the same request. A reference run appears in My runs as **awaiting Brandon’s submission**. Other structures appear as **awaiting lab review**.
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
