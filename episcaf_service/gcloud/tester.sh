#!/usr/bin/env bash
# Issue/rotate or revoke a cloud tester. Only a digest is uploaded.
set -euo pipefail
umask 077
action=${1:?usage: tester.sh issue|revoke owner}
owner=${2:?supply a persistent tester name}
[[ "$owner" =~ ^[a-zA-Z0-9_-]{1,48}$ ]] && [[ "$action" == issue || "$action" == revoke ]] || exit 2
number=$(gcloud projects describe "$(gcloud config get-value project)" --format='value(projectNumber)')
directory="$HOME/.episcaf-pilot/$number"
mkdir -p "$directory"
temp=$(mktemp "$directory/.issue.XXXXXX")
trap 'rm -f "$temp"' EXIT
python3 - "$action" "$owner" "$temp" <<'PY' | gcloud storage cp - "gs://episcaf-requests-$number/testers/$owner"
import hashlib, json, secrets, sys
from pathlib import Path
action, owner, path = sys.argv[1:]
token = owner + "." + secrets.token_urlsafe(32)
if action == "issue":
    Path(path).write_text(token+"\n")
print(json.dumps({"digest": hashlib.sha256(token.encode()).hexdigest() if action == "issue" else None}))
PY
if [[ "$action" == issue ]]; then
  mv "$temp" "$directory/$owner.token"
  printf 'Token saved to %s/%s.token\n' "$directory" "$owner"
fi
