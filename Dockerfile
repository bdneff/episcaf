# syntax=docker/dockerfile:1
# episcaf web service container (Google Cloud Run, or any container host).
#
# Two named build targets share one app base:
#   --target web            lightweight public/manual image. Candidate generation is CPU-only; DiscoTope
#                           ranking calls BioLib (EPISCAF_DISCOTOPE=biolib). This is also the DEFAULT
#                           stage, so a bare `docker build .` or `gcloud run deploy --source .` stays light.
#   --target web-discotope  private CPU image for the cloud-batch demo: adds a self-contained DiscoTope-3.0
#                           install so ranking runs locally (EPISCAF_DISCOTOPE=local) with no BioLib call.
# Build YAMLs: build-web.yaml builds `web`; build.yaml builds `web`, `web-discotope`, and the GPU worker.

############################
# Shared application base (lightweight).
############################
FROM python:3.12-slim AS appbase

# Build tools so any sdist-only dependency (e.g. freesasa) can compile if no wheel is available.
RUN apt-get update && apt-get install -y --no-install-recommends build-essential \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app

# Install deps first (cached layer) from the self-contained web requirements.
COPY episcaf_service/requirements-web.txt ./requirements-web.txt
RUN pip install --no-cache-dir -r requirements-web.txt

# App code (the big data/results/.git dirs are excluded via .dockerignore).
COPY . .

ENV EPISCAF_DISCOTOPE=biolib \
    EPISCAF_RATE_LIMIT=20/minute \
    EPISCAF_WS=/tmp/episcaf_runs \
    PORT=8080
EXPOSE 8080

# Cloud Run provides $PORT (default 8080); bind 0.0.0.0. Shell form so $PORT expands; exec so uvicorn
# is PID 1 and receives shutdown signals.
CMD exec uvicorn episcaf_service.app:app --host 0.0.0.0 --port ${PORT:-8080}

############################
# Private DiscoTope image (cloud-batch demo). CPU-only; built on the app base.
############################
FROM appbase AS web-discotope

# DiscoTope-3.0 runs in its OWN venv so its heavy torch/ESM-IF1 deps never touch the web app runtime.
# The adapter in episcaf_service/disco.py shells `/opt/discotope/venv/bin/python -m discotope3.main
# --models_dir /opt/discotope/models` (see EPISCAF_DISCOTOPE_HOME below).
#
# [UNVERIFIED - build-time]: this stage is NOT built in local development (no Docker there). Three
# things must be confirmed in the real Cloud Build and are deliberately not hard-coded here:
#   1. torch resolves to a CPU wheel (upstream requirements.txt may otherwise pull a large CUDA build);
#      pass a CPU index via --build-arg if needed.
#   2. the ESM-IF1 checkpoint DiscoTope pulls via torch.hub is available at build or run time.
#   3. the model weights archive - supplied as DISCOTOPE_MODELS_URL rather than an invented URL.
# The smoke test at the end of this stage is the gate: a broken install fails `docker build`.
ARG DISCOTOPE_REV=35d9f2e55f97eaba2a7acefbc394db58fb9670bc
ARG DISCOTOPE_REPO=https://github.com/Magnushhoie/DiscoTope-3.0
ARG DISCOTOPE_MODELS_URL=""

RUN apt-get update && apt-get install -y --no-install-recommends git curl unzip ca-certificates \
    && rm -rf /var/lib/apt/lists/*

# Pinned upstream checkout, installed with the project's OWN declared dependencies (not invented pins).
RUN python -m venv /opt/discotope/venv \
    && /opt/discotope/venv/bin/pip install --no-cache-dir --upgrade pip \
    && git clone "$DISCOTOPE_REPO" /opt/discotope/src \
    && git -C /opt/discotope/src checkout "$DISCOTOPE_REV" \
    && /opt/discotope/venv/bin/pip install --no-cache-dir -r /opt/discotope/src/requirements.txt \
    && /opt/discotope/venv/bin/pip install --no-cache-dir /opt/discotope/src

# Model weights live outside the pip package. Supply the upstream archive at build time:
#   --build-arg DISCOTOPE_MODELS_URL=<signed URL to models.zip>
RUN mkdir -p /opt/discotope/models \
    && if [ -n "$DISCOTOPE_MODELS_URL" ]; then \
         curl -fSL "$DISCOTOPE_MODELS_URL" -o /tmp/models.zip \
         && unzip -q /tmp/models.zip -d /opt/discotope/models \
         && rm /tmp/models.zip; \
       else \
         echo "DISCOTOPE_MODELS_URL not set: provide the upstream models archive at build time." >&2; \
         exit 1; \
       fi

ENV EPISCAF_DISCOTOPE=local \
    EPISCAF_DISCOTOPE_HOME=/opt/discotope

# Gate: rank the bundled lysozyme reference locally and confirm scores map back to author residues.
RUN python episcaf_service/gcloud/discotope_smoke.py

############################
# Default (last) stage: the lightweight public image. Keeps bare `docker build .` light.
############################
FROM appbase AS web
