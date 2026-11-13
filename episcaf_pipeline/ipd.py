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
`children_fraction` is that exploration/exploitation dial and is left to the caller -- it is the open
question from John's reports (are partial-diffusion children finding genuinely better designs, or just
best-of-N in a lucky basin? his own false-discovery caveat), so the controller does not hard-code it.

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
    children_fraction: float = 0.5       # share of that budget spent on partial-diffusion children
    noise_levels: tuple = (5, 10)        # partial-diffusion re-noise steps per winner (RFd3 partial_T)
    stop_eps: float = 0.0                # min score improvement to count as progress
    max_rounds: int = 4


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


def plan_round(cur_rows, prev_best, round_idx, cfg):
    """Decide the next round per epitope from this round's scored designs.

    `prev_best` is {group: best score} from the PREVIOUS round (or {} on round 0). Returns
    {group: {status, reason, winners, n_denovo, children}} where children is a list of
    {parent, noise, count} partial-diffusion specs. A stopped epitope gets status='stop'.
    """
    cur_best = best_in_group(cur_rows, cfg)
    succ = has_success(cur_rows, cfg)
    winners = select_winners(cur_rows, cfg)
    plan = {}
    for g, cb in cur_best.items():
        if succ.get(g):
            plan[g] = dict(status="stop", reason="success", winners=winners[g],
                           n_denovo=0, children=[])
            continue
        pb = prev_best.get(g)
        improved = pb is None or (cb > pb + cfg.stop_eps if cfg.better == "high"
                                  else cb < pb - cfg.stop_eps)
        if round_idx + 1 >= cfg.max_rounds:
            plan[g] = dict(status="stop", reason="max_rounds", winners=winners[g],
                           n_denovo=0, children=[])
        elif not improved:
            plan[g] = dict(status="stop", reason="plateau", winners=winners[g],
                           n_denovo=0, children=[])
        else:
            n_denovo, per_child = allocate(cfg.budget_per_epitope, cfg.children_fraction,
                                           len(winners[g]), cfg.noise_levels)
            children = ([dict(parent=w.get("design_id", w.get(cfg.group_key)), noise=nz, count=per_child)
                         for w in winners[g] for nz in cfg.noise_levels] if per_child else [])
            plan[g] = dict(status="continue", reason="improving", winners=winners[g],
                           n_denovo=n_denovo, children=children)
    return plan
