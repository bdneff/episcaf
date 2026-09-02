# episcaf_v3 — designing from an energetic definition of the epitope

A new design avenue for the episcaf project. It is a subdirectory of the `episcaf` repository (whose
root holds the shipped **v2** method — the pipeline, scoring, and manuscript that produced the DP4
library), next to the other components like `dp4_8vdl/`. It does not replace v2 yet: v2 is the
current method and the source of truth. v3 tests whether a different way of defining the epitope, and
of constraining the designs, gives better scaffolds. If it does, it takes over; if it does not, we
will have written down why.

## The idea, in a paragraph
v2 defines the epitope by geometry: the residues an antibody contacts (a heavy-atom distance cutoff),
and it hands RFdiffusion a contig that fixes the identity and backbone of all of them. But those
contacts do not contribute equally — a few carry most of the binding energy and the rest mostly fill
space. v3 defines the epitope by energy instead. It asks which residues actually carry the
antigen–antibody binding, and constrains the design accordingly: fix identity only for those
residues, fix shape (but not identity) for the ones whose job is structural, and leave the rest free.
Fixing fewer residues should give RFdiffusion more room to fold a good scaffold while still holding
the residues that do the binding.

## How this is built
One decision at a time, and we write each one down. Every methodological choice — the energy method,
the MD protocol, the classification thresholds, the contig strategy, the filters — is a dated entry
in `docs/DECISIONS.md`: the options, the choice, the reason, and how to check it. A choice is not
settled until it is an entry there. The manuscript (`manuscript/main.pdf`, built with `tectonic
main.tex`) is the running write-up, read top to bottom, as in v2.

## Layout (built out as we go)
- `manuscript/` — the running write-up (`main.tex` → `main.pdf`); opens by explaining v2 and why v3.
- `docs/DECISIONS.md` — the decision log, where every choice is recorded.
- `energetics/` — MD + per-residue interaction-energy decomposition + Cat 1/2/3 labels (to build).
- `contigs/` — the three-tier contig builder (to build).
- `filters/` — the RFdiffusion backbone pre-filter (orientation + antibody clash, before ProteinMPNN) (to build).

## Relationship to the other projects
- **The repository root (v2)** — the shipped method and data, one level up (`..`); v3 reuses its
  structures, scoring, and AlphaFold3 stages by relative path instead of forking them.
- **The `bcell_epitope` project** (a separate repository under `projects/`) — prior art for the MD
  and per-residue interaction-energy machinery; its `campaign2_holo_epitopeness` computes almost the
  same interface decomposition v3 needs. We use it as a template for how to run the simulations, not a
  protocol to copy wholesale — v3's energy and design choices are decided here, in `docs/DECISIONS.md`.
