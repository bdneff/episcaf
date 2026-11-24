#!/usr/bin/env python3
"""The RFD4 condition-spec builder: contig construction, the C_CRDGRP/C_SEQ conditions, and the
hand-emitted YAML (bare boolean/int keys). Pure logic; the ingest half is validated separately
against a real motif run. Run:  python3 tests/test_rfd4_generate.py
"""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from episcaf_service.rfd4 import generate as Gen  # noqa: E402


def test_contig_two_islands():
    spec = Gen.condition_spec("ag.cif", "A:205-216;A:278-292")
    # even-split (26,25,25) interleaved with the two islands
    assert spec.contig == "E0:26-26, A205-216, E1:25-25, A278-292, E2:25-25", spec.contig


def test_contig_single_residue_island():
    spec = Gen.condition_spec("ag.cif", "A:252")       # length 1 -> single-residue token
    assert ", A252, " in (", " + spec.contig + " , "), spec.contig
    assert "A252-252" not in spec.contig


def test_conditions_structure():
    spec = Gen.condition_spec("ag.cif", "A:205-216;A:278-292")
    assert spec.crdgrp == {1: ["sel('A/*/[205-216]/[N,CA,C,O]')", "sel('A/*/[278-292]/[N,CA,C,O]')"]}
    assert spec.seq_fixed == ["sel('A/*/[205-216]')", "sel('A/*/[278-292]')"]


def test_yaml_emission():
    spec = Gen.condition_spec("/data/ag.cif", "A:205-216;A:278-292", name="scaffold")
    y = Gen.to_yaml(spec)
    # bare int/bool keys, not quoted strings
    assert "\n      1:\n" in y, y
    assert "\n      true:\n" in y, y
    assert 'contig: "E0:26-26, A205-216, E1:25-25, A278-292, E2:25-25"' in y
    assert "    C_CRDGRP:" in y and "    C_SEQ:" in y
    assert '- select: "sel(\'A/*/[205-216]/[N,CA,C,O]\')"' in y
    assert "  input: /data/ag.cif" in y
    # structural sanity: C_CRDGRP appears before C_SEQ
    assert y.index("C_CRDGRP") < y.index("C_SEQ")


def test_generate_argv():
    argv = Gen.generate_argv("/tmp/spec.yaml", "/tmp/out", 8)
    assert argv[0] == "rfd4-generate"
    assert "inputs=/tmp/spec.yaml" in argv and "out_dir=/tmp/out" in argv
    assert "num_replicates=8" in argv and "~metrics" in argv


def main():
    tests = [v for k, v in sorted(globals().items()) if k.startswith("test_") and callable(v)]
    for t in tests:
        t()
        print(f"ok  {t.__name__}")
    print(f"\n{len(tests)} passed")


if __name__ == "__main__":
    main()
