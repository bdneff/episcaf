"""Phase-3 round control: seed selection, stopping, child planning, de-dup (John §1.3/§2.3).
Pure decision logic; the cluster runs the round loop on top. Run: python tests/test_rfd4_rounds.py
"""
import numpy as np

from episcaf_service.rfd4 import rounds
from episcaf_service.rfd4.config import DEFAULT


def _d(composite, offset):
    # a 103-CA design displaced by `offset` Å along x, so CA RMSD between two is |offset diff|.
    ca = np.zeros((103, 3)); ca[:, 0] = offset
    return {"composite": composite, "ca": ca, "id": f"d{offset}"}


def test_ca_rmsd_is_positionwise():
    a = np.zeros((10, 3)); b = np.zeros((10, 3)); b[:, 0] = 2.0
    assert abs(rounds.ca_rmsd(a, b) - 2.0) < 1e-9


def test_select_seeds_diverse_and_above_threshold():
    # two high scorers 0.2 A apart (too close) + two more spread out; expect the best + the far ones.
    scored = [_d(2.0, 0.0), _d(1.9, 0.2), _d(1.5, 3.0), _d(1.2, 6.0), _d(0.4, 9.0)]
    seeds = rounds.select_seeds(scored, DEFAULT)
    ids = [s["id"] for s in seeds]
    assert "d0.0" in ids                       # best kept
    assert "d0.2" not in ids                   # within 1.0 A of d0.0 -> dropped
    assert "d9.0" not in ids                   # composite 0.4 <= 0.5 -> never seeds
    assert len(seeds) <= DEFAULT.generation.n_seeds_max


def test_select_seeds_caps_at_n_max():
    scored = [_d(2.0 - 0.1 * i, 5.0 * i) for i in range(10)]   # all far apart, all passing
    seeds = rounds.select_seeds(scored, DEFAULT)
    assert len(seeds) == DEFAULT.generation.n_seeds_max        # 5


def test_should_stop_on_success_designs():
    assert rounds.should_stop([1.9, 2.0, 1.85], round_idx=1, cfg=DEFAULT)      # 3 >= 1.8
    assert not rounds.should_stop([1.9, 2.0, 1.2], round_idx=1, cfg=DEFAULT)   # only 2 pass


def test_should_stop_on_max_rounds():
    assert rounds.should_stop([0.5, 0.6], round_idx=5, cfg=DEFAULT)
    assert not rounds.should_stop([0.5, 0.6], round_idx=4, cfg=DEFAULT)


def test_plan_children_one_per_noise_capped():
    seeds = [_d(1.0, i) for i in range(3)]
    pairs = rounds.plan_children(seeds, DEFAULT)
    assert len(pairs) == 3 * len(DEFAULT.generation.noise_levels)              # 3 seeds x 3 noise
    assert {p[1] for p in pairs} == set(DEFAULT.generation.noise_levels)


def test_dedupe_drops_near_copies():
    kept = [_d(2.0, 0.0)]
    children = [_d(1.0, 0.1), _d(1.0, 5.0), _d(1.0, 5.05)]   # 0.1 near parent; 5.05 near 5.0
    out = rounds.dedupe_children(children, kept, DEFAULT)
    offs = sorted(c["ca"][0, 0] for c in out)
    assert offs == [5.0]                                     # 0.1 (near parent) and 5.05 (near 5.0) dropped


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
