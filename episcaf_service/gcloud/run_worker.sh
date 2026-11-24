#!/usr/bin/env bash
# episcaf data-mover worker launcher. scrontab runs this every ~5 min on the Gemini `data-mover`
# node. It performs one pilot_worker tick: gcloud via the google/cloud-sdk singularity image,
# authenticated with a buckets-scoped service-account key. Config and secrets live OUTSIDE git, in
# ~/.episcaf-pilot/worker.env. See episcaf_service/gcloud/README.md ("Automated data-mover worker").
set -euo pipefail

EPISCAF_HOME="${EPISCAF_HOME:-$HOME/.episcaf-pilot}"
CONF="$EPISCAF_HOME/worker.env"
[ -r "$CONF" ] || { echo "Missing $CONF (see episcaf_service/gcloud/README.md)" >&2; exit 1; }
# shellcheck disable=SC1090
source "$CONF"   # must set EPISCAF_REQUEST_BUCKET and EPISCAF_RESULT_BUCKET; may override the rest

: "${EPISCAF_REQUEST_BUCKET:?set EPISCAF_REQUEST_BUCKET in worker.env}"
: "${EPISCAF_RESULT_BUCKET:?set EPISCAF_RESULT_BUCKET in worker.env}"
REPO="${EPISCAF_REPO:-$HOME/episcaf-advisor}"
SIF="${EPISCAF_CLOUD_SDK_SIF:-$EPISCAF_HOME/cloud-sdk.sif}"
KEY="${EPISCAF_SA_KEY:-$EPISCAF_HOME/sa-key.json}"

[ -r "$SIF" ] || { echo "Missing cloud-sdk image $SIF (singularity pull it once)" >&2; exit 1; }
[ -r "$KEY" ] || { echo "Missing service-account key $KEY" >&2; exit 1; }

export CLOUDSDK_CONFIG="$EPISCAF_HOME/gcloud-config"
mkdir -p "$CLOUDSDK_CONFIG"
chmod 700 "$EPISCAF_HOME" "$CLOUDSDK_CONFIG" 2>/dev/null || true

# gcloud runs inside the cloud-sdk image. $HOME is auto-mounted, so $KEY and $CLOUDSDK_CONFIG are
# visible at their host paths. SINGULARITYENV_ passes the config dir through --cleanenv on any version.
export SINGULARITYENV_CLOUDSDK_CONFIG="$CLOUDSDK_CONFIG"
export EPISCAF_GCLOUD="singularity exec --cleanenv $SIF gcloud"

# Activate the buckets-scoped service account (local, idempotent, no network).
$EPISCAF_GCLOUD auth activate-service-account --key-file="$KEY" --quiet

export EPISCAF_REQUEST_BUCKET EPISCAF_RESULT_BUCKET
export EPISCAF_SBATCH="${EPISCAF_SBATCH:-$REPO/episcaf_service/handoff.sbatch}"
export EPISCAF_INBOX="${EPISCAF_INBOX:-$EPISCAF_HOME/inbox}"

exec python3 "$REPO/episcaf_service/gcloud/pilot_worker.py" tick "$@"
