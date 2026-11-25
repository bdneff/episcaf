"""Scaffold-budget routine: over-generate cheap RFD4 backbones to find N clash-free ones, cap the
tries, and always hand the structure-prediction step exactly N (best-by-clash when short).

Tests the two pure helpers (select_scaffolds, wave_size) that drive run_campaign's generation loop.
No data, no subprocess — run: python tests/test_rfd4_budget.py
"""
from episcaf_service.rfd4 import campaign
from episcaf_service.rfd4.config import DEFAULT


def _cand(i, clashes, passes):
    return campaign.Candidate(f"scaffold_{i}", design=None, backbone_clashes=clashes, passes=passes)


def test_select_prefers_passers_then_lowest_clash():
    # 3 passers (0 clashes) + 3 failers (varying clashes); ask for 5 -> all passers + 2 best failers.
    pool = [_cand(0, 0, True), _cand(1, 0, True), _cand(2, 0, True),
            _cand(3, 9, False), _cand(4, 4, False), _cand(5, 6, False)]
    sel = campaign.select_scaffolds(pool, 5)
    ids = [c.design_id for c in sel]
    assert len(sel) == 5, ids
    assert set(ids[:3]) == {"scaffold_0", "scaffold_1", "scaffold_2"}, ids   # passers first
    assert ids[3:] == ["scaffold_4", "scaffold_5"], ids                      # then clash 4 then 6, not 9
    assert "scaffold_3" not in ids                                           # worst (9) left out


def test_select_enough_passers_ignores_failers():
    pool = [_cand(i, 0, True) for i in range(6)] + [_cand(99, 1, False)]
    sel = campaign.select_scaffolds(pool, 4)
    assert len(sel) == 4 and all(c.passes for c in sel)


def test_select_returns_fewer_only_when_underfilled():
    pool = [_cand(0, 0, True), _cand(1, 2, False)]
    sel = campaign.select_scaffolds(pool, 10)         # asked 10, only 2 generated
    assert len(sel) == 2


def test_select_fallback_all_failed():
    # none pass the filter -> still return N, the lowest-clash ones (Brandon: "take the best n").
    pool = [_cand(i, c, False) for i, c in enumerate([8, 2, 5, 3, 9])]
    sel = campaign.select_scaffolds(pool, 3)
    assert [c.backbone_clashes for c in sel] == [2, 3, 5]


def test_wave_first_wave_is_target():
    assert campaign.wave_size(32, 0, 0, 8, 128) == 32


def test_wave_stops_when_target_met():
    assert campaign.wave_size(32, 32, 40, 8, 128) == 0
    assert campaign.wave_size(32, 33, 40, 8, 128) == 0


def test_wave_stops_at_try_cap():
    assert campaign.wave_size(32, 10, 128, 8, 128) == 0      # budget spent, even though short
    assert campaign.wave_size(32, 10, 130, 8, 128) == 0


def test_wave_sizes_from_survival_rate():
    # 16 passers from 32 generated -> survival 0.5; need 16 more -> ceil(16/0.5*1.2)=39, within budget.
    assert campaign.wave_size(32, 16, 32, 8, 128) == 39


def test_wave_clamps_to_remaining_budget():
    # need many, survival low, but only 10 tries left -> clamp to 10.
    assert campaign.wave_size(128, 1, 118, 8, 128) == 10


def test_wave_at_least_one_job():
    # tiny shortfall with high survival would round below a job; floor at designs_per_job.
    assert campaign.wave_size(32, 31, 31, 8, 128) == 8


def test_config_defaults():
    g = DEFAULT.generation
    assert g.target_scaffolds == 32 and g.max_scaffolds == 128 and g.max_generation_tries == 128


if __name__ == "__main__":
    import sys
    fns = [v for k, v in sorted(globals().items()) if k.startswith("test_") and callable(v)]
    failed = 0
    for fn in fns:
        try:
            fn()
            print(f"ok   {fn.__name__}")
        except Exception as e:  # noqa: BLE001
            failed += 1
            print(f"FAIL {fn.__name__}: {e}")
    print(f"\n{len(fns) - failed}/{len(fns)} passed")
    sys.exit(1 if failed else 0)
