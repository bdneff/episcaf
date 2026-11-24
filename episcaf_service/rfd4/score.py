"""The four structure metrics and the composite score (section 4). Composite is Brandon's design.

composite = base + pass_bonus * gate, with s(x;m,k) = 1/(1+exp(k*(x-m))) (lower-is-better on each
term). The base is a weighted blend; the gate is a product of steep sigmoids that only approaches 1
when all four thresholds are comfortably met, so a full pass scores ~1.8-3.0 and anything that misses
a threshold stays below 1. Ported from the §4.3 Python reference; verified against the §4.2 worked
examples in tests/test_rfd4_scoring.py (the parity gate for scoring).
"""
from __future__ import annotations

import math

import numpy as np

from episcaf_service.rfd4.config import DEFAULT, CompositeConfig
from episcaf_service.rfd4.geom import kabsch_rmsd


def sigmoid(x: float, midpoint: float, steepness: float) -> float:
    """1/(1+exp(k*(x-m))), clipped like John's R (exp arg to [-60,60]) to avoid overflow."""
    z = max(-60.0, min(60.0, steepness * (x - midpoint)))
    return 1.0 / (1.0 + math.exp(z))


def _is_missing(v) -> bool:
    return v is None or (isinstance(v, float) and math.isnan(v))


def composite_breakdown(af_clashes, epitope_rmsd, overall_rmsd, pae,
                        cfg: CompositeConfig = DEFAULT.composite):
    """Return {base, gate, composite} or None if any input is missing (NA-safe like the R)."""
    if any(_is_missing(v) for v in (af_clashes, epitope_rmsd, overall_rmsd, pae)):
        return None
    t = cfg.epitope_rmsd_target
    s = sigmoid

    base = (cfg.w_clash * s(af_clashes, *cfg.base_clash)
            + cfg.w_epitope_rmsd * s(epitope_rmsd, t, 4.0 / t)
            + cfg.w_overall_rmsd * s(overall_rmsd, *cfg.base_overall_rmsd)
            + cfg.w_pae * s(pae, *cfg.base_pae))

    k = cfg.gate_steepness
    gate = (s(epitope_rmsd, t, k / t)
            * s(overall_rmsd, cfg.gate_overall_rmsd_mid, k)
            * s(pae, cfg.gate_pae_mid, k)
            * s(af_clashes, cfg.gate_clash_mid, k))

    return {"base": base, "gate": gate, "composite": base + cfg.pass_bonus * gate}


def composite(af_clashes, epitope_rmsd, overall_rmsd, pae,
              cfg: CompositeConfig = DEFAULT.composite):
    """Composite score, or None if any metric is missing."""
    parts = composite_breakdown(af_clashes, epitope_rmsd, overall_rmsd, pae, cfg)
    return None if parts is None else parts["composite"]


def passes(composite_value, cfg: CompositeConfig = DEFAULT.composite) -> bool:
    return composite_value is not None and composite_value >= cfg.success_composite


# --- the two RMSD measures (both Kabsch-superposed) ---

def epitope_backbone_rmsd(design_epitope_bb: np.ndarray, native_epitope_bb: np.ndarray) -> float:
    """epitopeRMSD: AF2 model's epitope backbone (N,CA,C,O) vs the native PDB's, Kabsch-aligned."""
    return kabsch_rmsd(np.asarray(native_epitope_bb, float), np.asarray(design_epitope_bb, float))


def overall_ca_rmsd(af2_ca: np.ndarray, design_ca: np.ndarray) -> float:
    """overallRMSD: AF2 model CA vs the RFD4 design CA over all residues, Kabsch-aligned."""
    return kabsch_rmsd(np.asarray(design_ca, float), np.asarray(af2_ca, float))
