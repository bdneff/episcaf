# episcaf_v3 — designing from an *energetic* definition of the epitope

This is a new design avenue for the episcaf project. It is a subdirectory of the `episcaf`
repository (whose root holds the shipped **v2** method — the pipeline, scoring, and manuscript that
produced the DP4 library), sitting alongside the other components such as `dp4_8vdl/`. It does
**not** replace v2 yet: v2 is the method that produced the DP4 library and is the current source of
truth. v3 is where we test whether a different, more principled way
of defining the epitope — and of constraining the designs — produces better scaffolds. If it does,
it replaces the old strategy; if it does not, we have learned why, on the record.

## The thesis in a paragraph
v2 defines the epitope **geometrically**: the residues an antibody contacts (a heavy-atom distance
cutoff), and it hands RFdiffusion a contig that fixes the identity and backbone of *all* of them.
That treats every contact as equally load-bearing, which they are not. v3 defines the epitope
**energetically**: it asks which residues actually carry the antigen–antibody binding free energy,
and lets the design constraints follow the energetics — fix identity only where it is load-bearing,
fix shape where a residue's job is structural, and free everything else. The bet is that this gives
the generator more room to build a foldable scaffold without giving up the interaction that matters.

## How this is built
Reproducible, informed, one decision at a time. Every methodological choice — the energy method,
the MD protocol, the classification thresholds, the contig strategy, the filters — is a dated entry
in `docs/DECISIONS.md` with the options considered, the choice, the rationale, and the check that
backs it. Nothing becomes a "standard" until it is an entry with provenance. The manuscript
(`manuscript/main.pdf`, built with `tectonic main.tex`) is the living narrative record, read top to
bottom, exactly as in v2.

## Layout (grown as we build, not pre-populated)
- `manuscript/` — the living record (`main.tex` → `main.pdf`); opens by explaining v2 and why v3.
- `docs/DECISIONS.md` — the decision log; the reproducibility backbone.
- `energetics/` — (to build) MD + per-residue interaction-energy decomposition + Cat 1/2/3 labels.
- `contigs/` — (to build) the three-tier contig builder.
- `filters/` — (to build) the RFdiffusion backbone pre-filter (orientation + antibody clash, pre-MPNN).

## Relationship to the other projects
- **The repository root (the v2 method)** — the shipped method and data, one level up (`..` from
  here); v3 reuses its structures, scoring, and AF3 stages by relative path rather than forking them
  where sensible.
- **The sibling `bcell_epitope` project** (a separate repository under `projects/`) — prior art for
  the MD + per-residue interaction-energy machinery (its `campaign2_holo_epitopeness` computes almost
  exactly the interface decomposition v3 needs). It is a **template and reference**, not a protocol we
  inherit; v3's energy and design choices are decided here, in `docs/DECISIONS.md`.
