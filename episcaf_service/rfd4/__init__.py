"""RFD4-Proteina epitope-scaffolding engine.

A Python port of John Altin's RFD4 pipeline (AltinLab/epitope-universe, branch
feat/rfd4-generator, efficientScaffolding/) for running on Gemini/Slurm instead of the
Tamarind cloud broker. See ~/.claude/plans/jolly-questing-perlis.md and
RFD4_pipeline_description.pdf for the spec. The composite score is Brandon's design.

Single-pass (round 1) first; the seeded partial-diffusion round loop is added on top.
"""
