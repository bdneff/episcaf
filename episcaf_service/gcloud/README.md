# Deploy the private backbone pilot from Cloud Shell

Use the Google Cloud project that owns the existing `episcaf` service in `us-central1`. Billing, GPU quota for one L4, a Batch-compatible network, and permission to build images, deploy Cloud Run, create service accounts and assign IAM roles are required. The GPU build pins the upstream Linux/amd64 Foundry image by digest and uses its bundled weights. Allow time and disk space for the initial image build. The lab project pays for builds, storage and GPU jobs. The one-hour task timeout limits a single execution, not total project spend.

From a fresh checkout of the pilot branch:

```bash
gcloud config get-value project
git clone --branch pilot/private-simulated-runs --single-branch https://github.com/bdneff/episcaf.git episcaf-cloud-pilot
cd episcaf-cloud-pilot
bash episcaf_service/gcloud/deploy.sh
bash episcaf_service/gcloud/tester.sh issue john
```

For a later update, use `git pull --ff-only` in that checkout and repeat the deploy command. The tester command prints the path of a private token file in Cloud Shell. Deliver its contents to John separately from the link. Reissuing for the same owner rotates the token. `bash episcaf_service/gcloud/tester.sh revoke john` revokes existing sessions and leaves that owner's runs in place.

The script builds the CPU web and GPU worker images in Artifact Registry, creates two private buckets and two service accounts, and deploys a separate `episcaf-pilot` Cloud Run service. The existing `episcaf` site receives the text updates and a link to the pilot. Its own launch endpoints stay in staging mode. The pilot URL is browser-accessible. Reading a run or starting a paid job requires a tester token. The worker account can create result objects but cannot read tester records. No patient data, tokens, SSH keys, or Gemini paths are required.

## First run

Open the printed pilot URL. Choose one starting point:

- **Find candidate regions** — download the linked 3HFM lysozyme PDB, upload it, choose chain Y, generate regions and select one.
- **Find antibody contacts** — download and upload the linked 3HFM complex, choose chain Y and click **Find contact residues**. The server computes the antigen residues within 4 Å of bound chains H and L and verifies the region again at review.

Sign in, choose **Prepare design inputs**, review the 16 attempts and click **Start cloud run**. Both routes use the existing relative backbone filter: retain the lowest-scoring half of the scorable trials, up to eight per run. The antibody supplies the epitope for the second route. Both filters rank raw backbone occupancy. An absolute bare-backbone quality threshold remains open. Each ZIP contains the retained backbone PDBs, their scores and the filter decisions for all attempts. ProteinMPNN sequences and AF3 predictions belong to the downstream Gemini workflow. Other uploads can be saved for lab review and do not start a cloud GPU.

Batch provisions one GPU task per launch and performs the work after the browser closes. Retries of the same launch use the same Batch job ID. Cloud Run can scale to zero. The worker writes the archive and CSV before writing its completion record. My runs reads the saved record when the researcher returns. While jobs are active, the web service checks Batch on visits, at most about once every 30 seconds per run.

To check a submission from Cloud Shell, copy its `ep-…` Batch job ID from My runs and run:

```bash
gcloud batch jobs describe "$JOB_ID" --location=us-central1
```

A successful first cloud run should produce a completed Batch job and a downloadable archive. Compare the retained row count with the number of PDBs. Local tests use native-coordinate controls and stand-ins for GPU inference and Google APIs. They cannot verify project IAM, network, GPU availability or execution of the installed checkpoint.

[Google Batch GPU jobs](https://docs.cloud.google.com/batch/docs/create-run-job-gpus) · [Rosetta Commons Foundry worker image](https://rosettacommons.org/2026/02/25/6111/)
