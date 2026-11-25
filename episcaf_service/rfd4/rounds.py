"""Seeded partial-diffusion round control (John §1.3/§2.3) — the decision logic for Phase 3.

Three pure, tested pieces: pick the next round's seeds, decide when to stop, and plan + de-duplicate
the children. The round LOOP that generates children (RFD4 partial diffusion, task=existing-structure)
and folds them runs on the cluster and layers on top of these; it is the remaining Phase-3 work.

CA RMSD here is position-wise (same-length designs in the shared native epitope frame, as
campaign.evaluate_design already produces), not a re-alignment.
"""
from __future__ import annotations

import numpy as np

from episcaf_service.rfd4.config import DEFAULT, EngineConfig


def ca_rmsd(a: np.ndarray, b: np.ndarray) -> float:
    """Position-wise CA RMSD between two equal-length CA coordinate arrays (Nx3)."""
    a, b = np.asarray(a, float), np.asarray(b, float)
    n = min(len(a), len(b))
    if n == 0:
        return float("inf")
    return float(np.sqrt(((a[:n] - b[:n]) ** 2).sum(1).mean()))


def select_seeds(scored, cfg: EngineConfig = DEFAULT):
    """Seeds for the next round: the best-composite designs that clear `seed_min_composite` and sit at
    least `seed_min_separation` CA RMSD from every seed already chosen, up to `n_seeds_max`. `scored` is
    an iterable of objects/dicts with `.composite` and `.ca` (Nx3 CA coords). Pooled over all rounds —
    pass every design scored so far. Returns the chosen items (a diverse, high-scoring set)."""
    g = cfg.generation
    def comp(d):
        return d["composite"] if isinstance(d, dict) else d.composite
    def ca(d):
        return d["ca"] if isinstance(d, dict) else d.ca
    ranked = sorted((d for d in scored if comp(d) is not None and comp(d) > g.seed_min_composite),
                    key=comp, reverse=True)
    seeds = []
    for d in ranked:
        if len(seeds) >= g.n_seeds_max:
            break
        if all(ca_rmsd(ca(d), ca(s)) >= g.seed_min_separation for s in seeds):
            seeds.append(d)
    return seeds


def should_stop(composites, round_idx: int, cfg: EngineConfig = DEFAULT) -> bool:
    """Stop once `success_designs` designs reach `success_composite`, or after `max_rounds`. `round_idx`
    is 1-based (the round just completed). `composites` is every composite scored so far."""
    c, g = cfg.composite, cfg.generation
    n_pass = sum(1 for x in composites if x is not None and x >= c.success_composite)
    return n_pass >= c.success_designs or round_idx >= g.max_rounds


def plan_children(seeds, cfg: EngineConfig = DEFAULT):
    """One child per (seed, noise level) — a spread of partial-diffusion noise around each seed — capped
    at `af2_cap_per_round`. Returns (seed, noise) pairs in seed-major order."""
    g = cfg.generation
    pairs = [(s, noise) for s in seeds for noise in g.noise_levels]
    return pairs[:g.af2_cap_per_round]


def dedupe_children(children, kept, cfg: EngineConfig = DEFAULT):
    """Drop a child within `dedupe_rmsd` CA RMSD of any already-kept design (or an earlier kept child),
    so a round does not re-score near-copies of its parents. `children`/`kept` carry `.ca` (or ['ca'])."""
    g = cfg.generation
    def ca(d):
        return d["ca"] if isinstance(d, dict) else d.ca
    keep = list(kept)
    out = []
    for ch in children:
        if all(ca_rmsd(ca(ch), ca(k)) >= g.dedupe_rmsd for k in keep):
            out.append(ch)
            keep.append(ch)
    return out
