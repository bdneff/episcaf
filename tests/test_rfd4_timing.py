"""ColabFold log parsing for the benchmark MSA-vs-model split. The per-model "took Ns" line is the
stable signal; MSA time is derived by the campaign as fold-stage wall - model_seconds.

Run: python tests/test_rfd4_timing.py
"""
import tempfile
from pathlib import Path

from episcaf_service.rfd4 import timing

LOG = """2026-11-25 10:00:00,000 Running colabfold 1.5.5 (cache ...)
2026-11-25 10:00:01,000 Query 1/1: scaffold_0 (length 103)
2026-11-25 10:00:44,000 Sleeping for 8s. Reason: PENDING
2026-11-25 10:01:10,456 model_1 took 25.3s (3 recycles) with pLDDT 90.1 and ptm 0.823
2026-11-25 10:01:35,789 model_2 took 24.9s (3 recycles) with pLDDT 88.0 and ptm 0.80
2026-11-25 10:02:00,000 model_3 took 23.7s (3 recycles) with pLDDT 91.2 and ptm 0.85
2026-11-25 10:02:24,000 model_4 took 24.1s (3 recycles) with pLDDT 89.0 and ptm 0.82
2026-11-25 10:02:48,000 model_5 took 23.5s (3 recycles) with pLDDT 90.5 and ptm 0.84
2026-11-25 10:03:00,000 reranking models by 'plddt' metric
"""


def test_parse_sums_model_times():
    with tempfile.TemporaryDirectory() as d:
        p = Path(d) / "log.txt"
        p.write_text(LOG)
        out = timing.parse_colabfold_log(p)
    assert out["log_found"] is True
    assert len(out["models"]) == 5, out
    assert out["models"][0] == {"name": "model_1", "seconds": 25.3}
    assert abs(out["model_seconds"] - 121.5) < 1e-6, out["model_seconds"]


def test_missing_log_is_flagged():
    out = timing.parse_colabfold_log("/no/such/log.txt")
    assert out == {"log_found": False}


def test_log_without_model_lines():
    with tempfile.TemporaryDirectory() as d:
        p = Path(d) / "log.txt"
        p.write_text("2026-11-25 Query 1/1: x (length 103)\nReranking...\n")
        out = timing.parse_colabfold_log(p)
    assert out["log_found"] is True and "model_seconds" not in out


def test_msa_split_math():
    # the campaign computes msa = elapsed - model_seconds; emulate it here.
    out = {"model_seconds": 121.5}
    elapsed = 228.4
    msa = round(max(0.0, elapsed - out["model_seconds"]), 3)
    assert msa == 106.9


if __name__ == "__main__":
    import sys
    fns = [v for k, v in sorted(globals().items()) if k.startswith("test_") and callable(v)]
    bad = 0
    for fn in fns:
        try:
            fn(); print(f"ok   {fn.__name__}")
        except Exception as e:  # noqa: BLE001
            bad += 1; print(f"FAIL {fn.__name__}: {e}")
    print(f"\n{len(fns)-bad}/{len(fns)} passed")
    sys.exit(1 if bad else 0)
