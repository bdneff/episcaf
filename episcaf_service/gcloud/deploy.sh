#!/usr/bin/env bash
# Run in Cloud Shell from the checkout for the lab project.
set -euo pipefail
cd "$(dirname "$0")/../.."
if (( $# > 1 )); then echo 'Usage: deploy.sh [manual-gemini|cloud-batch]' >&2; exit 2; fi
case ${1:-manual-gemini} in
  manual-gemini) mode=manual-gemini; web_account=episcaf-inbox ;;
  cloud-batch) mode=cloud-lysozyme; web_account=episcaf-web ;;
  *) echo 'Usage: deploy.sh [manual-gemini|cloud-batch]' >&2; exit 2 ;;
esac
project=$(gcloud config get-value project)
region=us-central1
# Fail before creating resources if the expected public service is absent.
gcloud run services describe episcaf --region="$region" --project="$project" --format='value(status.url)' >/dev/null
number=$(gcloud projects describe "$project" --format='value(projectNumber)')
web="$web_account@$project.iam.gserviceaccount.com"
requests="episcaf-requests-$number"
results="episcaf-results-$number"
repo="$region-docker.pkg.dev/$project/episcaf-pilot"
tag=$(git rev-parse --short=12 HEAD)
env_vars="EPISCAF_PILOT_MODE=$mode,EPISCAF_REQUEST_BUCKET=$requests,EPISCAF_RESULT_BUCKET=$results"
web_image="$repo/web:$tag"
web_resources=(--memory=2Gi --cpu=2 --timeout=300)
accounts=("$web_account")
apis=(run.googleapis.com iam.googleapis.com cloudbuild.googleapis.com artifactregistry.googleapis.com storage.googleapis.com)
if [[ $mode == cloud-lysozyme ]]; then
  accounts+=(episcaf-backbones)
  apis+=(batch.googleapis.com compute.googleapis.com logging.googleapis.com)
fi
gcloud services enable "${apis[@]}"
if [[ $mode == cloud-lysozyme ]]; then
  # Avoid starting the large paid build if the project has zero GPU quota.
  check_quota() {
    local metric=$1
    local limit
    limit=$(jq -r --arg metric "$metric" '.quotas[] | select(.metric == $metric) | .limit')
    if [[ $limit == 0 || $limit == 0.0 ]]; then
      echo "$metric is zero. Check paid billing and GPU quota before deploying." >&2
      exit 1
    elif [[ -z $limit ]]; then
      echo "Check $metric in Cloud Quotas. Compute Engine did not report its limit." >&2
    fi
  }
  gcloud compute project-info describe --format=json | check_quota GPUS_ALL_REGIONS
  gcloud compute regions describe "$region" --format=json | check_quota NVIDIA_L4_GPUS
fi
for account in "${accounts[@]}"; do
  gcloud iam service-accounts describe "$account@$project.iam.gserviceaccount.com" >/dev/null 2>&1 ||
    gcloud iam service-accounts create "$account"
done
for bucket in "$requests" "$results"; do
  gcloud storage buckets describe "gs://$bucket" >/dev/null 2>&1 ||
    gcloud storage buckets create "gs://$bucket" --location="$region" --uniform-bucket-level-access --public-access-prevention
done
gcloud artifacts repositories describe episcaf-pilot --location="$region" >/dev/null 2>&1 ||
  gcloud artifacts repositories create episcaf-pilot --repository-format=docker --location="$region"
gcloud storage buckets add-iam-policy-binding "gs://$requests" --member="serviceAccount:$web" --role=roles/storage.objectUser >/dev/null
gcloud storage buckets add-iam-policy-binding "gs://$results" --member="serviceAccount:$web" --role=roles/storage.objectViewer >/dev/null
if [[ $mode == cloud-lysozyme ]]; then
  worker="episcaf-backbones@$project.iam.gserviceaccount.com"
  project_role() {
    gcloud projects add-iam-policy-binding "$project" --member="serviceAccount:$1" --role="$2" --condition=None --quiet >/dev/null
  }
  project_role "$web" roles/batch.jobsEditor
  project_role "$worker" roles/batch.agentReporter
  project_role "$worker" roles/logging.logWriter
  gcloud iam service-accounts add-iam-policy-binding "$worker" --member="serviceAccount:$web" --role=roles/iam.serviceAccountUser --quiet >/dev/null
  gcloud artifacts repositories add-iam-policy-binding episcaf-pilot --location="$region" --member="serviceAccount:$worker" --role=roles/artifactregistry.reader --quiet >/dev/null
  gcloud storage buckets add-iam-policy-binding "gs://$results" --member="serviceAccount:$worker" --role=roles/storage.objectCreator >/dev/null
  web_image="$repo/web-discotope:$tag"
  web_resources=(--memory=8Gi --cpu=4 --concurrency=1 --timeout=900)
  gcloud builds submit . --config=episcaf_service/gcloud/build.yaml \
    --substitutions="_WEB=$repo/web:$tag,_DEMO=$web_image,_GPU=$repo/backbones:$tag"
  gpu=$(gcloud artifacts docker images describe "$repo/backbones:$tag" --format='value(image_summary.fully_qualified_digest)')
  [[ $gpu == "$repo/backbones@sha256:"* ]] || { echo 'Could not resolve the worker image digest.' >&2; exit 1; }
  env_vars+=",EPISCAF_DISCOTOPE=local,EPISCAF_PROJECT=$project,EPISCAF_REGION=$region,EPISCAF_GPU_IMAGE=$gpu,EPISCAF_WORKER_SA=$worker"
else
  env_vars+=",EPISCAF_DISCOTOPE=biolib"
  gcloud builds submit . --config=episcaf_service/gcloud/build-web.yaml --substitutions="_WEB=$web_image"
fi
gcloud run deploy episcaf-pilot --image="$web_image" --region="$region" --service-account="$web" --allow-unauthenticated \
  --min-instances=0 "${web_resources[@]}" --set-env-vars="$env_vars"
pilot_url=$(gcloud run services describe episcaf-pilot --region="$region" --format='value(status.url)')
gcloud run services update episcaf --region="$region" --image="$repo/web:$tag" --update-env-vars="EPISCAF_PILOT_MODE=staging,EPISCAF_PILOT_URL=$pilot_url"
printf 'Pilot: %s/home.html\n' "$pilot_url"
