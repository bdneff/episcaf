#!/usr/bin/env python3
"""Unit tests for the IPD-loop ROUND DRIVER (episcaf_pipeline/ipd_driver.py).

These exercise the DRIVER's decisions + emitted artifacts on SYNTHETIC scored CSVs -- NOT GPU output.
No RFd3/MPNN/AF3: a tiny fake antigen PDB + a tiny designs table stand in, and we feed hand-built
scored CSVs back through `advance` to check that the right RFd3 inputs (de-novo + partial-diffusion
children) and the right round-state fall out.

Run:  python tests/test_ipd_driver.py     (or: pytest tests/)
"""
import csv
import glob
import json
import os
import re
import sys
import tempfile

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from episcaf_pipeline import ipd_driver          # noqa: E402
from episcaf_pipeline.ipd import IPDConfig        # noqa: E402

_DID = re.compile(r"^(?P<id>.+?)__contig(?P<cid>\d+)__seed(?P<seed>\d+)__rep(?P<rep>\d+)")


def _fake_antigen(path):
    # a single CA atom is enough -- the driver only references the path, never parses it
    with open(path, "w") as fh:
        fh.write("ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00           C\n")
        fh.write("END\n")


def _designs_table(path):
    rows = [
        # eA: two flank-split variants of epitope A:10-14
        dict(id="eA", contig_id=0, contig_string="5-5/A10-14/5-5", contig_length="15-15",
             growpatch_epitope="A:10-14", epitope_chunks=1),
        dict(id="eA", contig_id=1, contig_string="4-4/A10-14/6-6", contig_length="15-15",
             growpatch_epitope="A:10-14", epitope_chunks=1),
        # eB: two variants of epitope B:20-23
        dict(id="eB", contig_id=0, contig_string="6-6/B20-23/5-5", contig_length="15-15",
             growpatch_epitope="B:20-23", epitope_chunks=1),
        dict(id="eB", contig_id=1, contig_string="5-5/B20-23/6-6", contig_length="15-15",
             growpatch_epitope="B:20-23", epitope_chunks=1),
    ]
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)


def _manifest_ids(base, rnd):
    mf = os.path.join(base, f"round{rnd:02d}", "02_rfd3", "inputs_manifest.csv")
    with open(mf) as fh:
        return [r["design_id"] for r in csv.DictReader(fh)]


def _group(did):
    m = _DID.match(did)
    return m.group("id")


def _write_scored(path, rows):
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["design_id", "id", "composite", "is_pass"])
        w.writeheader()
        w.writerows(rows)


def _child_jsons(base, rnd):
    """RFd3 input JSONs in a round that are partial-diffusion children (carry partial_T)."""
    out = []
    for p in glob.glob(os.path.join(base, f"round{rnd:02d}", "02_rfd3", "inputs", "*.json")):
        payload = json.loads(open(p).read())
        inner = next(iter(payload.values()))
        if inner.get("partial_diffusion"):
            out.append((p, inner))
    return out


def _cfg():
    # small budget so the whole multi-round run writes only a handful of JSONs, but big enough that
    # children clear the per-(winner,noise) >=1 threshold once a streak is earned.
    return IPDConfig(mode="adaptive", topn_winners=2, budget_per_epitope=16,
                     children_step=0.25, children_cap=0.6, noise_levels=(5, 10), max_rounds=6)


def test_round0_init_writes_state_and_inputs():
    with tempfile.TemporaryDirectory() as base:
        antigen = os.path.join(base, "antigen.pdb"); _fake_antigen(antigen)
        dt = os.path.join(base, "designs.csv"); _designs_table(dt)
        state = ipd_driver.init_round0(base, antigen, dt, _cfg())

        assert os.path.exists(os.path.join(base, ipd_driver.STATE_FILE))
        assert set(state["epitopes"]) == {"eA", "eB"}
        assert state["round_idx"] == 0
        assert all(state["epitopes"][g]["best"] is None for g in ("eA", "eB"))
        # round 0 emits the full de-novo pool (4 contigs) and NO children
        ids = _manifest_ids(base, 0)
        assert len(ids) == 4 and not _child_jsons(base, 0)
        # a de-novo JSON is a standard stage02 motif-scaffolding input vs the antigen
        j = json.loads(open(os.path.join(base, "round00", "02_rfd3", "inputs", ids[0] + ".json")).read())
        inner = next(iter(j.values()))
        assert inner["input"].endswith("antigen.pdb") and "partial_diffusion" not in inner
        assert inner["select_fixed_atoms"]            # BKBN fixed atoms from the contig
        print("  ok: test_round0_init_writes_state_and_inputs")


def test_advance_denovo_only_when_no_streak():
    """First advance: no prior best -> streak 0 -> de-novo only (children must earn their budget)."""
    with tempfile.TemporaryDirectory() as base:
        antigen = os.path.join(base, "antigen.pdb"); _fake_antigen(antigen)
        dt = os.path.join(base, "designs.csv"); _designs_table(dt)
        ipd_driver.init_round0(base, antigen, dt, _cfg())

        # score round 0: both epitopes present, no success
        r0 = [dict(design_id=d, id=_group(d), composite=0.50, is_pass="") for d in _manifest_ids(base, 0)]
        scored0 = os.path.join(base, "scored0.csv"); _write_scored(scored0, r0)
        res = ipd_driver.advance(base, scored0)

        assert res["round_idx"] == 1 and not res["complete"]
        for g in ("eA", "eB"):
            assert res["summary"][g]["status"] == "continue"
            assert res["summary"][g]["children_fraction"] == 0.0
            assert res["summary"][g]["n_denovo"] == 16 and res["summary"][g]["n_child"] == 0
        assert not _child_jsons(base, 1)                      # no children yet
        assert len(_manifest_ids(base, 1)) == 32              # 16 de-novo x 2 epitopes
        st = ipd_driver.load_state(base)
        assert st["epitopes"]["eA"]["best"] == 0.5 and st["epitopes"]["eA"]["streak"] == 0
        print("  ok: test_advance_denovo_only_when_no_streak")


def test_advance_improving_emits_children_and_increments_streak():
    with tempfile.TemporaryDirectory() as base:
        antigen = os.path.join(base, "antigen.pdb"); _fake_antigen(antigen)
        dt = os.path.join(base, "designs.csv"); _designs_table(dt)
        ipd_driver.init_round0(base, antigen, dt, _cfg())

        # round 0 scored -> round 1 (de-novo only, streak 0, best=0.5)
        r0 = [dict(design_id=d, id=_group(d), composite=0.50, is_pass="") for d in _manifest_ids(base, 0)]
        s0 = os.path.join(base, "s0.csv"); _write_scored(s0, r0)
        ipd_driver.advance(base, s0)

        # round 1 scored: eA IMPROVES (0.80 > 0.50) -> streak 1 -> children;
        #                 eB PLATEAUS (0.50, no gain)  -> stop, emits nothing.
        r1 = []
        for d in _manifest_ids(base, 1):
            g = _group(d)
            r1.append(dict(design_id=d, id=g, composite=(0.80 if g == "eA" else 0.50), is_pass=""))
        s1 = os.path.join(base, "s1.csv"); _write_scored(s1, r1)
        res = ipd_driver.advance(base, s1)

        assert res["round_idx"] == 2
        assert res["summary"]["eB"]["status"] == "stop" and res["summary"]["eB"]["reason"] == "plateau"
        a = res["summary"]["eA"]
        assert a["status"] == "continue" and a["children_fraction"] == 0.25
        # cf 0.25 x budget 16 = 4 child slots = 2 winners x 2 noise, 1 each; 12 de-novo
        assert a["n_child"] == 4 and a["n_denovo"] == 12

        children = _child_jsons(base, 2)
        assert len(children) == 4
        for _, inner in children:
            assert inner["partial_T"] in (5, 10)                      # FLAGGED field present
            assert inner["input"].endswith(".pdb") and "round01" in inner["input"]  # parent backbone
            assert inner["select_fixed_atoms"]                        # epitope held fixed
        assert {inner["partial_T"] for _, inner in children} == {5, 10}

        # all round-2 emitted designs belong to eA only (eB stopped)
        assert {_group(d) for d in _manifest_ids(base, 2)} == {"eA"}

        st = ipd_driver.load_state(base)
        assert st["epitopes"]["eA"]["streak"] == 1 and st["epitopes"]["eA"]["best"] == 0.80
        assert st["epitopes"]["eB"]["status"] == "stopped"
        print("  ok: test_advance_improving_emits_children_and_increments_streak")


def test_success_stops_and_emits_nothing():
    with tempfile.TemporaryDirectory() as base:
        antigen = os.path.join(base, "antigen.pdb"); _fake_antigen(antigen)
        dt = os.path.join(base, "designs.csv"); _designs_table(dt)
        ipd_driver.init_round0(base, antigen, dt, _cfg())
        r0 = [dict(design_id=d, id=_group(d), composite=0.50, is_pass="") for d in _manifest_ids(base, 0)]
        s0 = os.path.join(base, "s0.csv"); _write_scored(s0, r0)
        ipd_driver.advance(base, s0)

        # round 1 scored: eA has a 4-filter PASS, eB plateaus -> BOTH stop -> round 2 emits nothing
        r1 = []
        for d in _manifest_ids(base, 1):
            g = _group(d)
            r1.append(dict(design_id=d, id=g, composite=0.9 if g == "eA" else 0.50,
                           is_pass="true" if g == "eA" else ""))
        s1 = os.path.join(base, "s1.csv"); _write_scored(s1, r1)
        res = ipd_driver.advance(base, s1)

        assert res["summary"]["eA"]["reason"] == "success"
        assert res["summary"]["eA"]["n_denovo"] == 0 and res["summary"]["eA"]["n_child"] == 0
        assert res["complete"] is True
        assert _manifest_ids(base, 2) == []            # empty manifest, nothing to generate
        assert not _child_jsons(base, 2)
        print("  ok: test_success_stops_and_emits_nothing")


if __name__ == "__main__":
    for name, fn in sorted(globals().items()):
        if name.startswith("test_") and callable(fn):
            fn()
    print("all ipd_driver tests passed")
