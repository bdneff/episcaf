"""Tunable dials for the RFD4 engine — the single place thresholds and weights live.

Mirrors the project convention that scorer weights are dials, not code (cf.
episcaf_analysis/presets.py). Values are ported verbatim from John's campaign_args.R /
RFD4_pipeline_description.pdf so our numbers match his reference (the parity gate). The
COMPOSITE block is Brandon's scoring design; keep the provenance note if you retune it.
"""
from __future__ import annotations

from dataclasses import dataclass, field

# --- fixed campaign geometry (tamarind.R top block) ---
PROTEIN_LENGTH = 103          # total designed length, epitope residues included


# --- clash filter (section 3; clash_method.R) ---
@dataclass(frozen=True)
class ClashConfig:
    af_clash_cutoff: float = 4.0        # a scaffold heavy atom clashes within this of the antibody
    shadow: float = 1.5                 # clash_antigen_shadow_angstrom: within this of a native
                                        #   antigen atom -> the position is antibody-compatible
    guard: float = 2.5                  # clash_overlap_guard_angstrom: closer than this to the
                                        #   antibody is a genuine overlap, shadow rule does not excuse it
    backbone_threshold_antibody: int = 3   # drop a backbone with more clashing residues (pre-AF2)
    backbone_threshold_aom: int = 7        # cutoff for the --clash-method aom fallback


# --- composite score (section 4). Brandon's design. s(x;m,k)=1/(1+exp(k*(x-m))), lower-is-better. ---
@dataclass(frozen=True)
class CompositeConfig:
    epitope_rmsd_target: float = 0.5    # t: centres the epitope terms
    success_composite: float = 1.8      # a design "passes" at/above this
    success_designs: int = 3            # stop a campaign once this many pass
    pass_bonus: float = 2.0             # composite = base + pass_bonus * gate
    gate_steepness: float = 12.0        # steepness of every gate sigmoid

    # base term weights (sum to 1.0): accessibility(AFclashes), epitope_rmsd, overall_rmsd, pae
    w_clash: float = 0.45
    w_epitope_rmsd: float = 0.25
    w_overall_rmsd: float = 0.20
    w_pae: float = 0.10

    # base-term sigmoid (midpoint, steepness); epitope terms derive from the target t
    base_clash: tuple[float, float] = (6.0, 0.5)
    base_overall_rmsd: tuple[float, float] = (2.0, 4.0)
    base_pae: tuple[float, float] = (5.0, 0.8)
    # epitope base: (t, 4.0/t); gate epitope: (t, gate_steepness/t). Computed in score.py from target.

    # gate midpoints for the non-epitope terms (steepness = gate_steepness)
    gate_overall_rmsd_mid: float = 2.0
    gate_pae_mid: float = 5.0
    gate_clash_mid: float = 0.5


# --- generation defaults. PILOT SCALE: deliberately small to bound per-click GPU cost; John's
# research scales (125 jobs/wave, up to 5000 designs) are far larger. Tune per the pilot budget. ---
@dataclass(frozen=True)
class GenerationConfig:
    designs_per_job: int = 8            # RFD4-Proteina returns 8 designs/job (numDesigns)
    sampling_effort: str = "thinking_high"   # thinking_low|medium|high
    unindex_motif: bool = True          # model places the islands (C_CRDGRP, not absolute coords)
    fix_motif_sequence: bool = True     # keep the native epitope sequence (C_SEQ)
    refold_designs: bool = False        # RFD4's own refolder off; AF2 is the structure-prediction step
    pilot_denovo_jobs: int = 4          # round-1 de novo jobs for the pilot (4*8 = 32 designs)
    noise_levels: tuple[float, ...] = (0.3, 0.5, 0.8)   # partial-diffusion children (phase 3)
    max_rounds: int = 5                 # phase 3


# --- structure-prediction step (a pluggable foundation-model seam; §1.4) ---
# The backend is swappable: AF2/ColabFold (validate the port against John here — the composite is
# AF2-CALIBRATED), then AF3 (better; already on Gemini; non-commercial only) or Boltz-2 (Apache,
# commercial-safe). Swapping the backend needs the composite dials above re-validated for that
# backend's PAE/confidence scale — a tuning pass, not a code change.
@dataclass(frozen=True)
class PredictorConfig:
    backend: str = "af2"                # "af2" | "af3" | "boltz"
    # AF2/ColabFold settings (John's §1.4): 5 models, 3 recycles, no relax, MSA on, no templates
    num_models: int = 5
    num_recycles: int = 3
    num_relax: int = 0
    use_msa: bool = True
    template_mode: str = "none"
    choose_best: bool = True            # keep rank_001


@dataclass(frozen=True)
class EngineConfig:
    protein_length: int = PROTEIN_LENGTH
    clash: ClashConfig = field(default_factory=ClashConfig)
    composite: CompositeConfig = field(default_factory=CompositeConfig)
    generation: GenerationConfig = field(default_factory=GenerationConfig)
    predictor: PredictorConfig = field(default_factory=PredictorConfig)


DEFAULT = EngineConfig()
