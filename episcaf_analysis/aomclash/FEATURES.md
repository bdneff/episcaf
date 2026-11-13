# AOMclash features

**The shipped scorer is one feature: `s_p1_free`** — the antibody-occupancy map read at every
scaffold atom, summed over the atoms outside native-antigen space (`free`). It is isotonic-calibrated
onto the AF3 clash count (`model.py:fit_prior` / `score_prior`); for ranking designs within a target
the calibration is cosmetic, so the raw sum ranks identically. See its row under **Occupancy** below.

This file is the reference for **what every computed feature is** — the extractor emits 111 columns,
and we document all of them so the set can still be changed deliberately without re-reading the math
in `features.py`. The richer **39-feature gradient-boosted model** (`model.py:FEATURES_39`) was tested
and **rejected** — it beats the single feature by only ~0.02 at real overfitting risk — and is kept
only as a head-to-head baseline. What we tried and the numbers live in
`results/aomclash/explorations/` (not front and center); this page just defines the inputs.

- Changing the model's inputs is editing the group lists in `model.py` (`OCCUPANCY`,
  `OUTWARD_COUNTS`, …). Adding a *new* kind of feature is editing `design_feats`/`target_features`
  here, then adding its name to a group in `model.py`. Nothing else keys off the list.

## The frame everything is measured in

Every per-design feature is computed after putting the design into the **epitope frame**
(`features.py:frame_from`, `design_feats`):

- **origin** = centroid of the native epitope Cα atoms.
- **z (outward normal)** = origin minus the centroid of antigen Cα within 25 Å → points *away* from
  the antigen body, into the space an antibody would occupy. `z > 0` is "the antibody side."
- **x** = the epitope's main in-plane axis; **y** = z × x.
- The design is Kabsch-aligned onto the native epitope first, so a design's scaffold atoms land in
  the same frame as the native antigen and the occupancy map.

Quantities each atom carries (`design_feats`, lines ~177–191):

| symbol | meaning |
|---|---|
| `p1` | occupancy map read at the atom — **P(antibody here)**, map smoothed ~2 Å (`P1`) |
| `p2` | same, map smoothed ~5 Å (`P2`) — coarser |
| `dag` | distance from the atom to the nearest **native antigen** heavy atom |
| `free` | `dag > 4.0` — atom is **not** inside native-antigen space (not explained by mimicking the real antigen) |
| `dist` | `‖L‖` — distance from the epitope-frame origin |
| `rr` | radial distance in the x–y plane |
| `z` | height along the outward axis |

Per-residue roll-ups: `pmax` = max `p1` over a residue's atoms; `r_sum` terms use the
log-complement trick to get "probability at least one atom of the residue clashes."

---

## The 39 model features

### Occupancy (15) — readouts of the antibody-occupancy map
The core signal: how much of the scaffold sits where antibodies usually are.

| name | definition | reads as |
|---|---|---|
| `s_p1` | Σ `p1` over all scaffold atoms | total antibody-occupancy the scaffold sits in |
| `s_p2` | Σ `p2` | same at coarser (5 Å) smoothing |
| `s_p1_free` | Σ `p1·free` | **the carved score** — occupancy over only free atoms; the ~0.87 map-only baseline |
| `s_p1_near` | Σ `p1·(~free)` | occupancy over atoms that *are* in native-antigen space (the carved-out part) |
| `n_p05` | #atoms `p1 > 0.05` | atoms in any non-trivial occupancy |
| `n_p15` | #atoms `p1 > 0.15` | " |
| `n_p30` | #atoms `p1 > 0.30` | " |
| `n_p50` | #atoms `p1 > 0.50` | atoms in high-occupancy (likely-clash) voxels |
| `n_p15_free` | #atoms `p1 > 0.15 & free` | same, free only |
| `n_p30_free` | #atoms `p1 > 0.30 & free` | " |
| `r_sum` | Σ_res `1 − Π(1 − p1)` | per-residue P(≥1 atom clashes), summed — residue-level expected clash count |
| `r_sum_free` | same, using `p1·free` | carved version of `r_sum` |
| `r_n30` | #residues `pmax > 0.3` | residues with a moderately occupied atom |
| `r_n50` | #residues `pmax > 0.5` | " |
| `r_n70` | #residues `pmax > 0.7` | residues with a strongly occupied atom |

### Outward-side counts (4) — raw geometric head-counts (unweighted cousin of occupancy)
No map; just count scaffold atoms in nested balls/caps on the antibody side.

| name | definition | reads as |
|---|---|---|
| `b_near20` | #atoms `dist < 20` | scaffold mass near the epitope, **any** direction (a sphere, not z-gated) |
| `b_out15` | #atoms `dist < 15 & z > 0` | mass in a tight outward cap |
| `b_out25` | #atoms `dist < 25 & z > 0` | mass in a wider outward cap |
| `b_out35_hi` | #atoms `dist < 35 & z > 4` | mass up high on the outward axis (where the antibody body sits) |

### Antigen contact (4) — proximity to the native antigen
| name | definition | reads as |
|---|---|---|
| `ag3` | #atoms `dag < 3` | scaffold mass overlapping native antigen (hard clash with the real structure) |
| `ag4` | #atoms `dag < 4` | " |
| `ag6` | #atoms `dag < 6` | " |
| `ag_free_out25` | #atoms `dist < 25 & z > 0 & free` | outward-cap atoms **not** in antigen space — mass that genuinely sticks into antibody territory, not explained by mimicking the antigen |

### Scaffold shape (10) — geometry of the design alone, no map
| name | definition | reads as |
|---|---|---|
| `z_mean` | mean atom height along z | how far the scaffold bulk sits on the antibody side |
| `z_p90` | 90th-percentile height | " (tail) |
| `z_max` | max height | tallest scaffold point above the epitope |
| `r_mean` | mean x–y radial distance | lateral spread |
| `ca_rg` | radius of gyration of scaffold Cα | overall compactness |
| `n_scaf` | # scaffold residues | size of the scaffold |
| `n_seg` | # scaffold segments (contiguous runs) | how fragmented the scaffold is |
| `seg_max` | longest segment length | longest continuous run |
| `seg_sum` | total residues across segments | " (total) |
| `apr` | atoms per residue = `len(X)/n_scaf` | mean heavy atoms/residue — proxy for side-chain bulk / composition |

### Epitope descriptors (6) — per-target, identical for every design of one epitope
From `target_features`. **Highest value per feature**: dropping this whole group was the single
biggest ablation hit (within-campaign ρ 0.919 → 0.891).

| name | definition | reads as |
|---|---|---|
| `frame_ratio` | eigenvalue ratio of the epitope's in-plane spread | elongated vs round epitope footprint |
| `n_epi` | # epitope residues | epitope size |
| `n_isl` | # epitope islands | discontiguity |
| `ag_out25` | #native-antigen atoms `z > 0 & dist < 25` | how much native antigen sits on the outward side — context the antibody must also clear |
| `ag_out25_hi` | of those, #with `z > 5` | native antigen mass high on the outward side |
| `ag_near15` | #native-antigen atoms `dist < 15` | local native density around the epitope |

---

## Dropped from the model (72 computed, not trained on)

Kept in the extractor so the full bank is one `feats=None` away and the cylinder stays available for
a head-to-head, but excluded from `FEATURES_39`. Ablation: dropping the 70 histograms costs ~0.007
within-campaign ρ (0.919 → 0.912); the 2 cylinder counts add ~0.003. A plateau, not a floor.

| names | count | what they are |
|---|---|---|
| `h0`…`h34` | 35 | 2-D histogram of **all** scaffold atoms over 5 radial rings × 7 height bands (`RB`×`ZB`). Spatial distribution of scaffold mass — mostly redundant with the occupancy sums. |
| `hf0`…`hf34` | 35 | same histogram over **free** atoms only |
| `cyl_aware`, `cyl_plain` | 2 | the native-aware cylinder count and the plain count — Brandon's old geometric surrogate, kept for the cylinder-vs-AOMclash comparison, not an input |

## Changing the set — the honest caveats

- **The shipped model is the single `s_p1_free` feature.** The notes below apply to the rejected
  39-feature alternative, kept only as a head-to-head baseline (`results/aomclash/explorations/`).
- **Features are highly correlated** (the 70 histograms ≈ the occupancy sum; the `n_p*` ladders are
  near-monotone in `s_p1`). Prune by **group**, judged under the grouped CV in `model.py`, never by
  per-feature importance — a per-feature plot is misleading here.
- **39 is a safe plateau, not the minimum.** A forward-selection pass for the smallest set that
  still holds ρ ≈ 0.9 is a later trim; it needs the assembled AF3 feature matrix (cluster-side, not
  built yet). See `results/aomclash/feature_set_decision.md`.
- The ablation numbers above are from **John's AF2 campaign** (`aom_ablation.py`). Re-confirm the
  plateau on our AF3 labels when the feature matrix exists before treating 39 as final.
