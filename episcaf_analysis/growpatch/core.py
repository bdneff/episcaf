"""growPatch core: the pure seed-and-grow patch enumerator (numpy only, no I/O).

Faithful port of the algorithm in John Altin's `segmentation/growPatch.py`
(AltinLab/epitope-universe @ upstream/main, 668de1d). Only the structure reading and the
freesasa surface step are lifted out (into `read.py`), so the same algorithm runs on our
inputs; the geometry/graph math here is byte-for-byte John's and is cross-checked against his
functions on a real antigen chain in `tests/test_growpatch_core.py`.

What growPatch does, with no antibody: it tiles an antigen's accessible surface with candidate
epitope "patches", each shaped like a plausible antibody footprint, so that whatever the real
epitope is, one of the candidates covers it. The steps here (given a surface set):

  1. adjacency: residues touch if any two heavy atoms are within `contact` (4.5 A).
  2. seed-and-grow: from every surface residue, greedily accrete the surface residue nearest the
     seed that stays spatially connected to the patch (heavy-atom distance <= `link`, 5.0 A),
     until the patch holds `size` raw residues (PATCH_SIZE=16).
  3. gap-10 fill: merge sequence gaps <= `gap` (10) and span-fill, but only over residues that
     exist in the structure, so unresolved loops break islands. This is the candidate epitope.
  4. keep candidates that are footprint-plausible: >= `min_res` (12) residues, <= `max_islands`
     (3) sequence islands, after dropping islands shorter than `min_island_len` (3) that hold no
     protected residue.
  5. Jaccard NMS: drop near-duplicate candidates (fewer islands, then larger, wins).

Residues are identified as `(chain, slot)` keys throughout, so a patch can span several antigen
chains. A "slot" is the residue's position in its chain's numbering with every insertion-coded
residue given a slot of its own (the chain<->slot<->label bookkeeping lives in `read.py`; the core
only ever sees integer slots). `present` ({chain: slots that exist with a complete backbone})
is what makes gap-fill skip unresolved loops.

Everything here is deterministic and order-agnostic: the NMS sort breaks ties on the sorted
residue list, so a run does not depend on set-iteration order once chain ids are in the keys.
"""
from __future__ import annotations

import numpy as np

# --- fitted defaults (growPatch.py @ upstream/main; see FEATURES/PARAMS note in the package) ---
SASA_CUT = 0.2          # relSASA surface cutoff (applied in read.py; here only as a documented default)
PATCH_SIZE = 16         # raw residues a seed grows to before gap-fill
GAP = 10                # sequence gaps <= GAP are bridged in gap-fill
MIN_RES = 12            # min residues in the kept, gap-filled candidate
MAX_ISLANDS = 3         # max sequence islands (summed over chains) in a candidate
MIN_ISLAND_LEN = 3      # gap-filled islands shorter than this are dropped unless protected
CONTACT = 4.5           # heavy-atom distance (A) for the residue adjacency graph
GROW_LINK = 5.0         # heavy-atom distance (A) a new residue must be within to join a patch
NMS_JACCARD = 0.8       # two candidates with Jaccard >= this are near-duplicates


# ------------------------------------------------------------------ geometry / graph

def min_atom_dist(a, b):
    """Minimum heavy-atom distance between two residues' atom arrays (port of min_atom_dist)."""
    return np.sqrt(((a[:, None, :] - b[None, :, :]) ** 2).sum(-1)).min()


def build_adjacency(residues, contact=CONTACT):
    """Residue contact graph: i~j if any two heavy atoms are within `contact` (port).

    `residues` is a list of dicts with "atoms" (N,3 heavy-atom coords) and "rep" (a representative
    point, the side-chain centroid). A 22 A representative-point prefilter skips the obvious
    non-contacts. Returns a list of neighbour-index sets, one per residue.
    """
    n = len(residues)
    adj = [set() for _ in range(n)]
    reps = np.array([r["rep"] for r in residues])
    rr = np.sqrt(((reps[:, None, :] - reps[None, :, :]) ** 2).sum(-1))
    for i in range(n):
        for j in range(i + 1, n):
            if rr[i, j] > 22.0:
                continue
            if min_atom_dist(residues[i]["atoms"], residues[j]["atoms"]) <= contact:
                adj[i].add(j)
                adj[j].add(i)
    return adj


def grow_patch(seed_idx, surf_idx, residues, adj, size=PATCH_SIZE, link=GROW_LINK):
    """Grow one patch from `seed_idx` (port of grow_patch).

    Greedily accrete the surface residue nearest the seed (smallest squared representative-point
    distance) that is adjacent to the current patch and within `link` heavy-atom distance of some
    patch member, until the patch holds `size` residues or nothing more can join. `surf_idx` is a
    set of surface residue indices. Returns a frozenset of residue indices.
    """
    seed_rep = residues[seed_idx]["rep"]
    patch = {seed_idx}
    while len(patch) < size:
        frontier = set()
        for p in patch:
            frontier |= (adj[p] & surf_idx)
        frontier -= patch
        ok = []
        for f in frontier:
            fa = residues[f]["atoms"]
            if any(min_atom_dist(fa, residues[p]["atoms"]) <= link for p in patch):
                ok.append(f)
        if not ok:
            break
        nxt = min(ok, key=lambda f: np.sum((residues[f]["rep"] - seed_rep) ** 2))
        patch.add(nxt)
    return frozenset(patch)


# ------------------------------------------------------------------ sets of (chain, slot) keys

def _by_chain(keys):
    out = {}
    for c, n in keys:
        out.setdefault(c, set()).add(n)
    return out


def islands(keys):
    """Sequence islands of a set of (chain, slot) keys -> [(chain, first, last), ...] (port).

    An island never crosses a chain; chains are listed in sorted order.
    """
    segs = []
    for c, ns in sorted(_by_chain(keys).items()):
        xs = sorted(ns)
        start = prev = xs[0]
        for x in xs[1:]:
            if x == prev + 1:
                prev = x
            else:
                segs.append((c, start, prev))
                start = prev = x
        segs.append((c, start, prev))
    return segs


def gap_fill(keys, gap=GAP, present=None):
    """Merge sequence gaps <= `gap` and span-fill each segment, chain by chain (port of gap_fill).

    A gap is never bridged between chains. With `present` ({chain: slots that exist in the
    structure with a complete backbone}) the fill skips the rest, so unresolved loops and dropped
    residues are never reported as epitope and a patch's islands break there.
    """
    out = set()
    for c, ns in _by_chain(keys).items():
        xs = sorted(ns)
        keep = (lambda n: True) if present is None else (lambda n, p=present.get(c, ()): n in p)
        start = prev = xs[0]
        for x in xs[1:]:
            if x - prev <= gap + 1:
                prev = x
            else:
                out |= {(c, n) for n in range(start, prev + 1) if keep(n)}
                start = prev = x
        out |= {(c, n) for n in range(start, prev + 1) if keep(n)}
    return out


def trim_short_islands(keys, min_len=MIN_ISLAND_LEN, protect=frozenset()):
    """Drop gap-filled islands shorter than `min_len` unless they hold a protected residue (port).

    Removes stray singletons that a spatial patch drags in far from the footprint core. In
    production enumeration `protect` is empty; the benchmark passes the reference epitope so a real
    epitope residue is never trimmed.
    """
    keep = set()
    for c, a, b in islands(keys):
        seg = {(c, n) for n in range(a, b + 1)}
        if (b - a + 1) >= min_len or (seg & protect):
            keep |= seg
    return keep


def jaccard(a, b):
    """|a ∩ b| / |a ∪ b| (port); 0 for two empty sets."""
    a, b = set(a), set(b)
    return len(a & b) / len(a | b) if (a or b) else 0.0


# ------------------------------------------------------------------ enumeration

def enumerate_raw_patches(residues, adj, surf_idx, size=PATCH_SIZE, link=GROW_LINK, min_raw=3):
    """Grow a patch from every surface residue; keep the distinct raw patches with >= `min_raw`.

    Returns a set of frozensets of residue indices (the raw, pre-gap-fill patches).
    """
    surf_idx = frozenset(surf_idx)
    raw = set()
    for s in surf_idx:
        p = grow_patch(s, surf_idx, residues, adj, size, link)
        if len(p) >= min_raw:
            raw.add(p)
    return raw


def nms_dedupe(cands, nms_jaccard=NMS_JACCARD):
    """Jaccard non-max suppression on gap-filled candidate sets (port of the main() NMS block).

    `cands` is a list of (raw_keys, gf_keys, n_islands). Fewer islands wins, then larger, then the
    sorted residue list (a deterministic tiebreak so the result is independent of set-iteration /
    hash-seed order). Keeps a candidate if its Jaccard to every already-kept candidate is
    < `nms_jaccard`. Returns the kept list in the same tuple shape.
    """
    ordered = sorted(cands, key=lambda t: (t[2], -len(t[1]), sorted(t[1])))
    kept = []
    for raw_r, gf_r, nisl in ordered:
        if all(jaccard(gf_r, k[1]) < nms_jaccard for k in kept):
            kept.append((raw_r, gf_r, nisl))
    return kept


def candidate_patches(residues, adj, surf_idx, idx_key, present,
                      size=PATCH_SIZE, link=GROW_LINK, gap=GAP,
                      min_res=MIN_RES, max_islands=MAX_ISLANDS, min_island_len=MIN_ISLAND_LEN,
                      protect=frozenset(), nms_jaccard=NMS_JACCARD):
    """Full candidate enumeration for one (surface set, patch size): the heart of growPatch.

    Grows raw patches from every surface seed, gap-10 fills each into a candidate epitope, keeps
    the footprint-plausible ones (>= `min_res` residues, <= `max_islands` islands after trimming
    short unprotected islands), and Jaccard-dedupes. Returns the kept list of
    (raw_keys, gf_keys, n_islands), deduped and NMS-ordered.

    Inputs:
      residues, adj  -- from build_adjacency over the residue list
      surf_idx       -- residue indices whose relSASA >= the cutoff (from read.py; in a test, any set)
      idx_key        -- {residue index: (chain, slot)} to map raw patches back to keys
      present        -- {chain: slots with a complete backbone} for gap-fill
      protect        -- reference residues never trimmed (empty in production)
    """
    raw_patches = enumerate_raw_patches(residues, adj, surf_idx, size, link)
    cands = []
    for p in raw_patches:
        raw_r = {idx_key[i] for i in p}
        gf_r = gap_fill(raw_r, gap, present)
        gf_r = trim_short_islands(gf_r, min_island_len, protect)
        if not gf_r:
            continue
        nisl = len(islands(gf_r))
        if len(gf_r) >= min_res and nisl <= max_islands:
            cands.append((raw_r, gf_r, nisl))
    return nms_dedupe(cands, nms_jaccard)


def greedy_set_cover(target, patch_sets):
    """Greedy minimum set-cover of `target` keys by named patch sets (port; benchmark/validation).

    `patch_sets` is a list of (name, key_set). Returns (chosen, uncovered): the patches picked, each
    covering the most still-uncovered target residues, and whatever target residues no patch covers.
    """
    remaining, chosen, pool = set(target), [], list(patch_sets)
    while remaining:
        best = max(pool, key=lambda ps: len(ps[1] & remaining), default=None)
        if best is None or not (best[1] & remaining):
            break
        chosen.append(best)
        remaining -= best[1]
        pool.remove(best)
    return chosen, remaining
