#!/usr/bin/env python3
"""Build-time smoke test for the private `web-discotope` image.

Runs the *real* local DiscoTope-3.0 adapter (EPISCAF_DISCOTOPE=local -> /opt/discotope/venv) on the
bundled 3HFM lysozyme reference, chain Y, and asserts the returned scores map back to the original
author residue numbers of that chain. A clean CLI exit from DiscoTope alone does not prove the page
would rank anything; this checks the adapter's residue mapping end to end. Any failure exits non-zero
so `docker build` fails loudly rather than shipping an image that silently cannot rank.

Invoked from the Dockerfile `web-discotope` stage:  python episcaf_service/gcloud/discotope_smoke.py
"""
import sys

from episcaf_service import disco, lysozyme_job


def main() -> int:
    if not disco.available():
        print("DiscoTope provider not configured (expected EPISCAF_DISCOTOPE=local).", file=sys.stderr)
        return 1
    expected = {disco._label(r) for r in disco._residues(lysozyme_job.REFERENCE, "Y")}
    if not expected:
        print("No residues found on chain Y of the lysozyme reference.", file=sys.stderr)
        return 1

    scores = disco.score_cif(lysozyme_job.REFERENCE, "Y", "solved")
    if not scores:
        print("DiscoTope returned no scores for the lysozyme reference.", file=sys.stderr)
        return 1

    stray = set(scores) - expected
    if stray:
        print(f"Scores reference residues absent from chain Y: {sorted(stray)[:5]}", file=sys.stderr)
        return 1
    if len(scores) < 0.9 * len(expected):
        print(f"Only {len(scores)}/{len(expected)} residues scored; mapping looks incomplete.", file=sys.stderr)
        return 1
    if not all(0.0 <= v <= 1.0 for v in scores.values()):
        print("Normalized DiscoTope scores fell outside [0, 1].", file=sys.stderr)
        return 1

    print(f"DiscoTope smoke OK: {len(scores)}/{len(expected)} chain-Y residues scored; "
          "author residue numbers preserved.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
