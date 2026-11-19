# episcaf web service container (Google Cloud Run, or any container host).
# Runs the growPatch candidate tool; the FastAPI app serves its own UI, so this one image is the
# whole thing. Candidate generation is CPU-only (no GPU/cluster); DiscoTope ranking calls BioLib.
FROM python:3.12-slim

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
