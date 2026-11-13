"""Score a design's antibody-clash propensity: the deployable AOMclash scorer.

This is the glue that makes AOMclash *do* something: given an antigen, an epitope, and a design
structure, it places the design in the epitope frame, reads the occupancy map, and returns
`s_p1_free` -- the single carved-occupancy feature we ship as the clash score (see
results/aomclash/feature_set_decision.md). It is the learned, finer-grained successor to the
native-aware cylinder (a continuous occupancy-weighted sum instead of a geometric head-count), and
replaces the cylinder as the no-antibody accessibility metric.

Orchestration ported from John Altin's `aom_score.py` / `aom_features.target_context`
(AltinLab/epitope-universe @ upstream/feat/clash-method), using our `features.py` and `aom_map.py`:

    build_target(antigen, epitope, map_path)   # frame + native-antigen KD-tree + the smoothed map
    score_design(target, design_atoms, epi_pos)  # -> dict incl. s_p1_free (None if atoms missing)

`epi_pos` is the design chain's residue ids that correspond to the epitope, in native-epitope order;
the pipeline stage that runs this knows it from the contig (the same epitope mapping stage05 uses).
`native_epi_pos(epitope)` is the default for a structure that keeps native epitope numbering.

The map read excludes a target's own SAbDab complexes when `exclude_key`/`exclude_pdb` are given
(leave-one-antigen-out, for scoring a design whose antigen is in the map); omit them for a new antigen.
"""
from __future__ import annotations

import os
import re
from dataclasses import dataclass

import numpy as np
from scipy.spatial import cKDTree

from . import aom_map
from . import features as F

_RES_RE = re.compile(r"^([A-Za-z0-9]+):(-?\d+)(?:-(-?\d+))?$")
DEFAULT_MAP = os.path.join(os.path.dirname(__file__), "aom_map.pkl")


def parse_epitope(spec):
    """'A:205-216;A:278-292' or chainless '205-216,278' (one chain) -> [(chain, a, b), ...]."""
    isl, cur = [], None
    for tok in str(spec).replace(";", ",").split(","):
        tok = tok.strip()
        if not tok:
            continue
        if ":" in tok:
            m = _RES_RE.match(tok)
            if not m:
                raise ValueError(f"bad epitope token {tok!r}")
            cur = m.group(1)
            a, b = int(m.group(2)), int(m.group(3) or m.group(2))
        else:
            m = re.match(r"^(-?\d+)(?:-(-?\d+))?$", tok)
            if not m or cur is None:
                raise ValueError(f"epitope token {tok!r} needs a chain (e.g. 'A:205-216')")
            a, b = int(m.group(1)), int(m.group(2) or m.group(1))
        isl.append((cur, a, b))
    return isl


def native_epi_pos(epitope):
    """Default epi_pos: the native epitope residue ids in order (for structures keeping native numbering)."""
    isl = parse_epitope(epitope) if isinstance(epitope, str) else epitope
    return [s for _, a, b in isl for s in range(a, b + 1)]


@dataclass
class Target:
    """Everything about one (antigen, epitope) that a design is scored against."""
    isl: list
    ref: np.ndarray            # native epitope backbone (N,CA,C,O per residue, island order)
    o: np.ndarray
    Rf: np.ndarray
    ag_tree: cKDTree
    epi_L: np.ndarray
    P1: np.ndarray
    P2: np.ndarray
    target_feats: dict
    n_complexes: int


def build_target(antigen, epitope, map_path=DEFAULT_MAP, exclude_key=None, exclude_pdb=None):
    """Frame the epitope and load the occupancy map (port of aom_features.target_context).

    `antigen` is a PDB path or already-read atoms [(name, chain, resid, x, y, z)]; `epitope` is a spec
    string or [(chain, a, b)]. Antibody chains in the antigen file are ignored (only epitope chains are
    used for the frame and the clash KD-tree). Returns a Target, or None if the epitope cannot be framed.
    """
    atoms = F.read_pdb_atoms(antigen) if isinstance(antigen, str) else list(antigen)
    isl = parse_epitope(epitope) if isinstance(epitope, str) else list(epitope)
    ep_chains = {c for c, _, _ in isl}

    ca = [(c, s, np.array(xyz)) for n, c, s, *xyz in atoms if n == "CA"]
    epi_ca = np.array([p for c, s, p in ca for ch, a, b in isl if c == ch and a <= s <= b])
    if len(epi_ca) < 5:
        return None
    ag_ca = np.array([p for c, s, p in ca if c in ep_chains])
    ag_atoms = np.array([xyz for _, c, _, *xyz in atoms if c in ep_chains], dtype=float)
    n_isl = len(isl)

    tf = F.target_features(epi_ca, ag_ca, ag_atoms, n_isl)
    if tf is None:
        return None
    o, Rf, ag_tree, epi_L, target_feats = tf

    ob = {(c, s, n): np.array(xyz) for n, c, s, *xyz in atoms if n in ("N", "CA", "C", "O")}
    try:
        ref = np.array([ob[(ch, s, n)] for ch, a, b in isl for s in range(a, b + 1)
                        for n in ("N", "CA", "C", "O")])
    except KeyError:
        return None

    field, n = aom_map.mean_field(aom_map.load_map(map_path), exclude_key, exclude_pdb)
    P1, P2 = aom_map.smoothed_maps(field)
    return Target(isl=isl, ref=ref, o=o, Rf=Rf, ag_tree=ag_tree, epi_L=epi_L,
                  P1=P1, P2=P2, target_feats=target_feats, n_complexes=n)


def score_design(target, design, epi_pos, segs=()):
    """Features of one design against `target`, including `s_p1_free` (port of the design_feats call).

    `design` is a PDB path or read atoms; its epitope must be chain 'A' at residue ids `epi_pos`
    (native-epitope order). Returns the feature dict, or None if the design lacks the epitope backbone
    or any scaffold atoms. `s_p1_free` is the shipped clash score.
    """
    atoms = F.read_pdb_atoms(design) if isinstance(design, str) else list(design)
    segs = list(segs) or [1]      # segs drive only the seg_* shape features; s_p1_free is independent,
                                  # so a placeholder is fine when the caller only wants the clash score
    feat = F.design_feats(atoms, list(epi_pos), segs, target.ref, target.o, target.Rf,
                          target.ag_tree, target.epi_L, target.P1, target.P2)
    if feat is None:
        return None
    return {**target.target_feats, **feat}


def aomclash_score(target, design, epi_pos, segs=()):
    """Just the shipped clash score `s_p1_free` for a design (None if it cannot be scored)."""
    feat = score_design(target, design, epi_pos, segs)
    return None if feat is None else float(feat["s_p1_free"])
