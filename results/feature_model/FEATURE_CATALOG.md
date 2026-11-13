# Binding-model features: catalog, candidates, and the keep rule

The interpretable ML model that ranks scaffolded designs by predicted antibody binding. This is the
feature roadmap: what we already use, what new **interpretable** per-epitope/antigen features are worth
adding, and the rule for keeping them. Grounded in the DP4 binding ML (`results/dp4_ml/dp4_ml_findings.tex`)
and a survey of the sibling projects (file leads below are reuse candidates — verify before trusting).

## The principle (from the DP4 result)

Binding is **additive and mildly non-linear** — a GAM / scorecard matched the best model; forbidding
interactions cost nothing. So the model stays an **interpretable additive logistic regression / EBM
on a small set of features**, not a black box. Signal is **secondary-structure- and native-context-
dominated** (AUC ~0.73 design-only → ~0.83 with native context); **dynamics (B-factor/ENM) was the
weakest block and confounded by the bound state.** New features must clear a real bar, not 0.0001 AUC.

**Keep rule.** A candidate is kept only if, computed on the DP4 join and evaluated under the existing
grouped-by-epitope CV (`scripts/dp4_binder_ml.py`), it (a) is not near-redundant with a feature we
already have (check correlation), and (b) adds **meaningful marginal AUC** as an interpretable term.
Otherwise it is cut. Prefer features that are causal *before* an antibody exists (apo/static), not
bound-state quantities that are downstream of binding (the dynamics confound).

## Baseline — features already computed/used

Design: `epitope_rmsd`, `overall_rmsd`, `mean/epitope/scaffold_pae`, `ptm`, `plddt`,
`fraction_disordered`, `af3_n_clash_res`, `aomclash s_p1_free`, (legacy) native-aware cylinder.
Epitope/antigen context: helix/strand/loop fractions, island count/size, residue count, Rg,
composition fractions (hydrophobic/charged/aromatic/Gly/Pro), relSASA, ΔSASA, B-factor(z), GNM/ENM
flexibility. The native-context group gave the DP4 lift; the dynamics subset was weakest.

## Candidate NEW interpretable features (ranked; add-first at the top)

**Tier 1 — cheap, interpretable, orthogonal, reuse exists:**
1. **Residue depth** — distance of a residue from the molecular surface (burial beyond what SASA
   sees; SASA saturates, depth keeps ranking). New. `Bio.PDB.ResidueDepth` (MSMS) or KD-tree to the
   surface points from `epitope-universe/src/epitope_universe/sasa.py`.
2. **Protrusion / contact-number** — few Cα neighbours in a ~10 Å ball ⇒ sticks out ⇒ antibody-
   accessible (the classic Thornton protrusion index). **Already coded**:
   `scripts/webtool_epitope_propensity.py`; the `growpatch/eyeview.py` apex-depth is a richer cousin.
   Promote to a model feature.
3. **Epitope crowding by the rest of the antigen** — native-antigen mass on the antibody-approach
   side the paratope must clear. **Already computed** as AOMclash target features `ag_out25`,
   `ag_out25_hi`, `ag_near15` — and dropping that group was AOMclash's single biggest ablation hit.
   Promote these into the binding model too. (Native-context, the DP4 lift direction.)
4. **Glycosylation proximity** — distance to the nearest N-linked sequon (N-X-S/T, X≠P); glycan-
   shielded surface is poorly accessible. New, trivial (sequence regex + Cα distance), fully interpretable.
5. **Electrostatic potential / net charge over the patch** — charge complementarity determinant.
   New as a static feature; cheap per-residue Coulomb/formal-charge sum, or APBS+pdb2pqr for PB.

**Tier 2 — higher value, more work:**
6. **Sequence conservation / Shannon entropy** (per residue, from an MSA) — conserved surface patches
   are common functional epitopes. New; the one candidate needing external data (jackhmmer/ConSurf).
   Highest-value truly-new signal.
7. **Surface curvature / planarity** of the patch (antibodies favour flat/convex caps). Partial today
   (`frame_ratio`, `Rg`); fit a plane/quadric to surface points. New-ish, interpretable.
8. **Hydrophobic moment / amphipathicity** — spatial organization of hydrophobicity, not just the
   fraction (which we have). Cheap add-on to the composition code.
9. **Contact order / sequence dispersion of the epitope** — mean sequence separation of spatially
   contiguous patch residues (local vs assembled-from-far); complements island count. Trivial.

## Defer / treat carefully
- **Bound-state / MD-heavy, confounded** (MLCE coupling, gRINN interaction energy, MM-GBSA
  decomposition, holo interaction energy, apo-MD RMSF/DCCM/FRESEAN): interpretable but downstream of
  binding or expensive; the DP4 finding argues against leaning on dynamics. Revisit only with apo MD
  (separately gated — see [[dynamics-features-research]]).
- **Native antibody-interaction type counts** (arpeggio H-bond/ionic/aromatic tallies; cached in
  `epitope-universe/scripts/arpeggio_precompute.py`): a native-context label like the ones that
  helped, but it describes the *bound native complex*, so it's only available for known-antibody
  targets, not a lone antigen — not usable as a predictor on a growPatch candidate.
- **DiscoTope-3.0 / BepiPred-3.0 propensity** (baking path exists: `scripts/webtool_discotope_bake.py`):
  a strong single feature but a black box — keep as an optional ranker/term, not a core interpretable one.

## Plan
Compute the Tier-1 set on the DP4 join (`data/dp4_binding/dp4_im0276_joined.csv`, local), check
orthogonality + marginal AUC under grouped-epitope CV, keep the ones that clear the bar, and fold the
survivors into both the composite and the logistic binding model as additive, inspectable terms.
Tier-2 (esp. conservation) next. Sibling projects seen: `episcaf_v3` (energetics), `minibinder`,
`md`, `IMMREP26`, `collab_learned_geometry`, `tcr*` — feature-method reuse noted in the survey.
