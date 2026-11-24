#!/usr/bin/env python3
"""Authenticated candidate enumeration for the private pilot (/pilot/enumerate).

Verifies tester authentication (signed-out and revoked -> 401, signed-in -> candidates) and that the
public /enumerate suppresses DiscoTope ranking for anonymous callers on the hosted demo while keeping
unranked candidate scanning. No GPU/cluster and no real DiscoTope are touched.

Run:  /Users/bneff/micromamba/envs/bcell-repair/bin/python tests/test_pilot_enumerate.py
"""
import os
import sys
import tempfile
from pathlib import Path
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

# Simulated pilot mode with a throwaway token DB, set before the app imports read them.
_DB = tempfile.NamedTemporaryFile(prefix="episcaf_pilot_", suffix=".db", delete=False).name
os.environ["EPISCAF_PILOT_MODE"] = "simulated"
os.environ["EPISCAF_PILOT_DB"] = _DB

from fastapi.testclient import TestClient            # noqa: E402
from episcaf_service import pilot                    # noqa: E402
from episcaf_service.app import app                  # noqa: E402

FIX = ROOT / "docs" / "pdb" / "8PWW.pdb"             # chain-A antigen used by the growPatch tests
BODY = FIX.read_bytes()
URL = "/pilot/enumerate?fmt=pdb&chain=A&disco=0"


def test_requires_tester_and_enumerates():
    pilot.initialize()
    alice = pilot.issue("alice")
    pilot.issue("bob")                               # a second owner exists; its token is never used here

    with TestClient(app) as client:
        # Signed out -> 401, no enumeration.
        r = client.post(URL, content=BODY)
        assert r.status_code == 401, (r.status_code, r.text)

        # Sign in Alice: Bearer token in the header (never the URL); the cookie persists on the client.
        s = client.post("/pilot/session", headers={"Authorization": "Bearer " + alice})
        assert s.status_code == 200 and s.json()["owner"] == "alice", s.text

        # Signed in -> candidates; unranked because no DiscoTope provider is configured here.
        r = client.post(URL, content=BODY)
        assert r.status_code == 200, r.text
        data = r.json()
        assert data["n_patches"] >= 1 and data["patches"], data
        assert data["has_disco"] is False, data

        # Revoke Alice's token -> 401 even though the browser cookie is still set.
        with pilot.database() as db:
            db.execute("DELETE FROM testers WHERE owner=?", ("alice",))
        r = client.post(URL, content=BODY)
        assert r.status_code == 401, (r.status_code, r.text)


def test_public_enumerate_blocks_anonymous_hosted_ranking():
    # On the hosted demo (local DiscoTope), an anonymous ranked request to the public /enumerate is
    # forced to unranked scanning; the signed-in path at /pilot/enumerate is the only ranked route.
    with patch.object(pilot, "mode", lambda: "cloud-lysozyme"), \
            patch.dict(os.environ, {"EPISCAF_DISCOTOPE": "local"}):
        with TestClient(app) as client:
            r = client.post("/enumerate?fmt=pdb&chain=A&disco=1", content=BODY)
            assert r.status_code == 200, r.text
            data = r.json()
            assert data["patches"], data
            assert data["has_disco"] is False, data      # ranking suppressed for the anonymous caller


if __name__ == "__main__":
    test_requires_tester_and_enumerates()
    test_public_enumerate_blocks_anonymous_hosted_ranking()
    print("pilot enumerate auth OK")
