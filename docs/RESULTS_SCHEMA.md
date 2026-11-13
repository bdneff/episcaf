# episcaf results CSV — the output schema

The useful output file: one row per design, with the AlphaFold3 metrics and the composite score, so a
run can be ranked and filtered (the LatentX-style deliverable). It is what `stage05` produces and what
`docs/results.html` reads. Metric columns are optional — the viewer shows whatever is present — but
these are the canonical names.

## Identity
| column | meaning |
|---|---|
| `design_id` | unique design: `<target>__contig<c>__seed<s>__rep<r>` |
| `id` / `target` | the epitope/target this design presents |
| `epitope` | epitope spec (`A:298-307;B:37-47`) |
| `source` | how the epitope was chosen: `growpatch`, `antibody` (known Ab contact), `known` |
| `contig_id`, `seed`, `rep` | the generation coordinates |

## Fidelity — does it present the epitope in its native geometry?
| column | meaning | good |
|---|---|---|
| `epitope_chunk_rmsd` | epitope backbone RMSD to native (Å) | ≤ 1 |
| `overall_rmsd` | whole-window backbone RMSD to the MPNN model (Å) | ≤ 2 |

## Confidence — AlphaFold3
| column | meaning | good |
|---|---|---|
| `mean_pae` | mean predicted aligned error (Å) | < 5 |
| `epitope_pae` | PAE within the epitope | low |
| `scaffold_pae` | PAE within the scaffold | low |
| `ptm` | predicted TM-score | high |
| `plddt` | mean pLDDT, if present | high |

## Accessibility / clash — can an antibody reach the epitope?
| column | meaning | good |
|---|---|---|
| `af3_n_clash_res` | real AF3 clashing residues vs the native antibody (antibody arm — ground truth) | 0 |
| `aomclash_pred` | **AOMclash** predicted clash on the no-antibody arm (the learned, finer-grained successor to the cylinder) | low |
| `cylinder_native_aware`, `cylinder_ca_clashes`, `native_in_cylinder` | legacy cylinder surrogate — kept for before/after, superseded by `aomclash_pred` |

## Score / selection
| column | meaning |
|---|---|
| `composite` | the weighted score (dials in `episcaf_analysis/presets.py`) |
| `rank_in_group` | rank within the target/epitope |
| `is_pass` | DP3 four-filter: `overall_rmsd ≤ 2 ∧ epitope_chunk_rmsd ≤ 1 ∧ mean_pae < 5 ∧ (af3_n_clash_res == 0` on the antibody arm, else low `aomclash_pred)` |

The four-filter (`overall_rmsd ≤ 2`, `epitope_chunk_rmsd ≤ 1`, `mean_pae < 5`, clash == 0) is the
in-silico pass used throughout; the composite ranks within it.
