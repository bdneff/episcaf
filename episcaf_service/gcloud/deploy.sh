#!/usr/bin/env bash
# Run from a fresh checkout in Cloud Shell with the lab project selected.
set -euo pipefail
cd "$(dirname "$0")/../.."
project=$(gcloud config get-value project)
region=us-central1
# Fail before creating any resources if Cloud Shell points at a different project.
gcloud run services describe episcaf --region="$region" --project="$project" --format='value(status.url)' >/dev/null
number=$(gcloud projects describe "$project" --format='value(projectNumber)')
web="episcaf-inbox@$project.iam.gserviceaccount.com"
requests="episcaf-requests-$number"
results="episcaf-results-$number"
repo="$region-docker.pkg.dev/$project/episcaf-pilot"
tag=$(git rev-parse --short=12 HEAD)
gcloud services enable run.googleapis.com iam.googleapis.com cloudbuild.googleapis.com artifactregistry.googleapis.com storage.googleapis.com
gcloud iam service-accounts describe "$web" >/dev/null 2>&1 ||
  gcloud iam service-accounts create episcaf-inbox
for bucket in "$requests" "$results"; do
  gcloud storage buckets describe "gs://$bucket" >/dev/null 2>&1 ||
    gcloud storage buckets create "gs://$bucket" --location="$region" --uniform-bucket-level-access --public-access-prevention
done
gcloud artifacts repositories describe episcaf-pilot --location="$region" >/dev/null 2>&1 ||
  gcloud artifacts repositories create episcaf-pilot --repository-format=docker --location="$region"
gcloud storage buckets add-iam-policy-binding "gs://$requests" --member="serviceAccount:$web" --role=roles/storage.objectUser >/dev/null
gcloud storage buckets add-iam-policy-binding "gs://$results" --member="serviceAccount:$web" --role=roles/storage.objectViewer >/dev/null
gcloud builds submit . --tag="$repo/web:$tag"
gcloud run deploy episcaf-pilot --image="$repo/web:$tag" --region="$region" --service-account="$web" --allow-unauthenticated \
  --min-instances=0 --memory=2Gi --cpu=2 --timeout=300 \
  --set-env-vars="EPISCAF_PILOT_MODE=manual-gemini,EPISCAF_REQUEST_BUCKET=$requests,EPISCAF_RESULT_BUCKET=$results,EPISCAF_DISCOTOPE=biolib"
pilot_url=$(gcloud run services describe episcaf-pilot --region="$region" --format='value(status.url)')
# Keep the public site in staging mode and keep its existing service configuration.
gcloud run services update episcaf --region="$region" --image="$repo/web:$tag" --update-env-vars="EPISCAF_PILOT_URL=$pilot_url"
printf 'Pilot: %s/home.html\n' "$pilot_url"
