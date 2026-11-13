"""IPD loop controller: pick winners, spend the next round's budget, decide when to stop.

The IPD ("iterative protein design") loop runs rounds of generate -> AOMclash backbone pre-filter ->
MPNN -> AF3 -> score, seeding each round from the best of the last. The generate/MPNN/AF3 steps are
GPU/cluster work; THIS module is the pure controller that decides, from a round's scored designs,
what the next round should be. It runs locally, is deterministic, and is unit-tested; a driver
(cluster-side) emits each round's inputs from `plan_round` and feeds the next round's scores back.

Each round's budget per epitope is split between:
  * de-novo  -- fresh RFd3 contigs (exploration of new scaffolds), and
  * children -- PARTIAL-DIFFUSION neighbours of the round's winners: re-noise a good backbone part way
                and denoise with the epitope fixed, exploring the neighbourhood of something that
                worked (exploitation).
How much budget goes to children is the exploration/exploitation dial -- John's open false-discovery
question (are partial-diffusion children genuinely better, or just best-of-N in a lucky basin?). The
default is **adaptive**: start de-novo-only and let each epitope EARN its children budget by a streak
of round-over-round improvement (children_fraction = children_step x streak, capped). So an epitope
whose winners' neighbourhood keeps paying off gets exploited harder, one that stops improving stays
exploratory -- the data decides per epitope instead of us guessing. `mode="fixed"` uses a flat
`children_fraction` instead.

Stopping: an epitope stops when it has a design meeting `success` (e.g. a four-filter pass), or its
best score has not improved by more than `stop_eps` over the last round (a plateau), or `max_rounds`
is reached.
"""
from __future__ import annotations

from dataclasses import dataclass, field


@dataclass
class IPDConfig:
    score_key: str = "composite"
    better: str = "high"                 # 'high' or 'low' -- orientation of score_key
    group_key: str = "id"                # the epitope/target key
    pass_key: str = "is_pass"            # truthy column marking a success (optional; "" to ignore)
    topn_winners: int = 3                # winners carried forward per epitope
    budget_per_epitope: int = 64         # designs to spawn next round per still-open epitope
    mode: str = "adaptive"               # 'adaptive' (children earn budget via improvement) or 'fixed'
    children_fraction: float = 0.5       # fixed mode: flat share spent on partial-diffusion children
    children_step: float = 0.25          # adaptive: children_fraction = step x improvement streak ...
    children_cap: float = 0.6            # ... capped here
    noise_levels: tuple = (5, 10)        # partial-diffusion re-noise steps per winner (RFd3 partial_T)
    stop_eps: float = 0.0                # min score improvement to count as progress
    max_rounds: int = 4


def adaptive_fraction(streak, cfg):
    """Children budget earned by a streak of round-over-round improvement (0 at streak 0 -> de-novo)."""
    return min(cfg.children_cap, cfg.children_step * max(0, streak))


def _num(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


def _truthy(v):
    return str(v).strip().lower() in ("1", "true", "yes", "pass", "t")


def best_in_group(rows, cfg):
    """{group: best score} for the configured score_key/orientation (ignores unscored rows)."""
    out = {}
    for r in rows:
        g, s = r.get(cfg.group_key), _num(r.get(cfg.score_key))
        if s is None:
            continue
        if g not in out or (s > out[g] if cfg.better == "high" else s < out[g]):
            out[g] = s
    return out


def select_winners(rows, cfg):
    """Top `topn_winners` rows per group: successes first (if pass_key), then best score."""
    from collections import defaultdict
    groups = defaultdict(list)
    for r in rows:
        if _num(r.get(cfg.score_key)) is not None:
            groups[r.get(cfg.group_key)].append(r)
    sign = -1 if cfg.better == "high" else 1
    winners = {}
    for g, rs in groups.items():
        rs.sort(key=lambda r: (0 if (cfg.pass_key and _truthy(r.get(cfg.pass_key))) else 1,
                               sign * _num(r.get(cfg.score_key))))
        winners[g] = rs[:cfg.topn_winners]
    return winners


def has_success(rows, cfg):
    """{group: True} where any design is a success (pass_key truthy)."""
    out = {}
    if not cfg.pass_key:
        return out
    for r in rows:
        if _truthy(r.get(cfg.pass_key)):
            out[r.get(cfg.group_key)] = True
    return out


def allocate(budget, children_fraction, n_winners, noise_levels):
    """Split a round budget into (n_denovo, child_specs-per-(winner,noise)).

    children get `round(budget*children_fraction)` slots, spread over winners x noise levels (at least
    one each if any budget and winners); the rest is de-novo. Returns (n_denovo, per_child_count)."""
    n_child_total = round(budget * children_fraction)
    slots = max(1, n_winners * len(noise_levels)) if n_winners else 0
    per_child = (n_child_total // slots) if slots else 0
    n_child = per_child * slots
    return budget - n_child, per_child


def plan_round(cur_rows, prev_state, round_idx, cfg):
    """Decide the next round per epitope from this round's scored designs.

    `prev_state` is {group: {"best": score, "streak": int}} from the PREVIOUS round (or {} on round 0).
    Returns (plan, next_state). plan is {group: {status, reason, winners, children_fraction, n_denovo,
    children}} with children a list of {parent, noise, count} partial-diffusion specs; a stopped
    epitope gets status='stop'. next_state carries each group's updated best + improvement streak, to
    pass back in next round (this is what makes the adaptive children budget work).
    """
    prev_state = prev_state or {}
    cur_best = best_in_group(cur_rows, cfg)
    succ = has_success(cur_rows, cfg)
    winners = select_winners(cur_rows, cfg)
    plan, next_state = {}, {}
    for g, cb in cur_best.items():
        ps = prev_state.get(g, {})
        pb = ps.get("best")
        improved_over_pb = pb is not None and (cb > pb + cfg.stop_eps if cfg.better == "high"
                                               else cb < pb - cfg.stop_eps)
        streak = (ps.get("streak", 0) + 1) if improved_over_pb else 0   # 0 on round 0 (no prior)
        next_state[g] = {"best": cb, "streak": streak}

        def stop(reason):
            plan[g] = dict(status="stop", reason=reason, winners=winners[g],
                           children_fraction=0.0, n_denovo=0, children=[])
        if succ.get(g):
            stop("success"); continue
        if round_idx + 1 >= cfg.max_rounds:
            stop("max_rounds"); continue
        if pb is not None and not improved_over_pb:          # a real plateau (needs a prior round)
            stop("plateau"); continue

        cf = cfg.children_fraction if cfg.mode == "fixed" else adaptive_fraction(streak, cfg)
        n_denovo, per_child = allocate(cfg.budget_per_epitope, cf, len(winners[g]), cfg.noise_levels)
        children = ([dict(parent=w.get("design_id", w.get(cfg.group_key)), noise=nz, count=per_child)
                     for w in winners[g] for nz in cfg.noise_levels] if per_child else [])
        plan[g] = dict(status="continue", reason="improving", winners=winners[g],
                       children_fraction=round(cf, 3), n_denovo=n_denovo, children=children)
    return plan, next_state
