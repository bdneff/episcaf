#!/usr/bin/env python3
"""Unit tests for the IPD loop controller (episcaf_pipeline/ipd.py) -- pure logic, no cluster."""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from episcaf_pipeline import ipd   # noqa: E402


def rows(*specs):
    # specs: (design_id, group, composite, is_pass)
    return [dict(design_id=d, id=g, composite=c, is_pass=p) for d, g, c, p in specs]


def test_best_in_group_and_winners():
    cfg = ipd.IPDConfig(topn_winners=2)
    rs = rows(("d1", "eA", 0.9, ""), ("d2", "eA", 0.7, ""), ("d3", "eA", 0.95, ""),
              ("d4", "eB", 0.5, "true"))
    assert ipd.best_in_group(rs, cfg) == {"eA": 0.95, "eB": 0.5}
    w = ipd.select_winners(rs, cfg)
    assert [r["design_id"] for r in w["eA"]] == ["d3", "d1"]      # top-2 by composite (high)
    assert [r["design_id"] for r in w["eB"]] == ["d4"]


def test_winners_put_passers_first():
    cfg = ipd.IPDConfig(topn_winners=2)
    rs = rows(("hi", "e", 0.99, ""), ("passlo", "e", 0.60, "true"), ("mid", "e", 0.80, ""))
    w = ipd.select_winners(rs, cfg)
    assert w["e"][0]["design_id"] == "passlo"                    # a success outranks a higher score
    assert [r["design_id"] for r in w["e"]] == ["passlo", "hi"]


def test_allocate_splits_budget():
    # 64 budget, 50% children, 2 winners x 2 noise = 4 slots -> 8 each -> 32 children, 32 de-novo
    n_denovo, per_child = ipd.allocate(64, 0.5, n_winners=2, noise_levels=(5, 10))
    assert (n_denovo, per_child) == (32, 8)
    assert ipd.allocate(64, 0.0, 2, (5, 10)) == (64, 0)          # all de-novo
    assert ipd.allocate(10, 1.0, 0, (5, 10)) == (10, 0)          # no winners -> all de-novo


def test_plan_round_success_stops():
    cfg = ipd.IPDConfig()
    plan, _ = ipd.plan_round(rows(("d1", "e", 0.8, "true")), {}, 0, cfg)
    assert plan["e"]["status"] == "stop" and plan["e"]["reason"] == "success"


def test_plan_round_fixed_mode_children_and_denovo():
    cfg = ipd.IPDConfig(mode="fixed", topn_winners=2, budget_per_epitope=64,
                        children_fraction=0.5, noise_levels=(5, 10))
    rs = rows(("d1", "e", 0.80, ""), ("d2", "e", 0.70, ""))
    plan, state = ipd.plan_round(rs, {"e": {"best": 0.5, "streak": 0}}, 0, cfg)   # 0.80 > 0.50
    p = plan["e"]
    assert p["status"] == "continue" and p["n_denovo"] == 32
    assert len(p["children"]) == 4 and all(c["count"] == 8 for c in p["children"])  # 2 winners x 2 noise
    assert {c["noise"] for c in p["children"]} == {5, 10}
    assert {c["parent"] for c in p["children"]} == {"d1", "d2"}
    assert state["e"]["best"] == 0.80 and state["e"]["streak"] == 1


def test_adaptive_children_budget_ramps_with_improvement_streak():
    cfg = ipd.IPDConfig(mode="adaptive", topn_winners=2, budget_per_epitope=100,
                        children_step=0.25, children_cap=0.6, noise_levels=(5, 10))
    rs = rows(("d1", "e", 0.80, ""), ("d2", "e", 0.70, ""))
    # round 0: no prior -> streak 0 -> de-novo only (children earn their budget later)
    p0, s0 = ipd.plan_round(rs, {}, 0, cfg)
    assert p0["e"]["children_fraction"] == 0.0 and p0["e"]["children"] == [] and p0["e"]["n_denovo"] == 100
    assert s0["e"]["streak"] == 0
    # round 1: improved over round 0 -> streak 1 -> cf 0.25
    p1, s1 = ipd.plan_round(rows(("d3", "e", 0.90, "")), s0, 1, cfg)
    assert p1["e"]["children_fraction"] == 0.25 and s1["e"]["streak"] == 1 and p1["e"]["children"]
    # round 2: still improving -> streak 2 -> cf 0.50
    p2, s2 = ipd.plan_round(rows(("d4", "e", 0.95, "")), s1, 2, cfg)
    assert p2["e"]["children_fraction"] == 0.50 and s2["e"]["streak"] == 2
    # streak keeps climbing but the cap holds
    assert ipd.adaptive_fraction(5, cfg) == 0.6


def test_plan_round_plateau_and_maxrounds_stop():
    cfg = ipd.IPDConfig(stop_eps=0.01, max_rounds=4)
    rs = rows(("d1", "e", 0.80, ""))
    assert ipd.plan_round(rs, {"e": {"best": 0.80}}, 1, cfg)[0]["e"]["reason"] == "plateau"
    assert ipd.plan_round(rs, {"e": {"best": 0.5}}, 3, cfg)[0]["e"]["reason"] == "max_rounds"


if __name__ == "__main__":
    for name, fn in sorted(globals().items()):
        if name.startswith("test_") and callable(fn):
            fn()
            print(f"  ok: {name}")
    print("all ipd controller tests passed")
