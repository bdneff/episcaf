"""Build the RFD4-Proteina input and drive `rfd4-generate`.

The input is a YAML condition spec (RFD4 compiles it to an annotated CIF). For our unindexed
multi-island motif scaffolding we mirror the pattern that the quickstart binder spec uses and that
we confirmed runs on Gemini:

    <name>:
      input: <single-chain antigen structure>
      contig: "E0:g0-g0, <island0>, E1:g1-g1, <island1>, ..., E_n:gn-gn"
      conditions:
        C_CRDGRP: {1: [<island backbone atoms>...]}   # hold the epitope's internal geometry as ONE
                                                       #   rigid group -> native inter-island
                                                       #   arrangement preserved, placement free (unindexed)
        C_SEQ:    {true: [<island residues>...]}       # keep the native epitope sequence

`contig` crops the input to just the epitope islands (preserved segments) and inserts the scaffold
as fixed-length expandable segments (E0..En), total length = PROTEIN_LENGTH. RFD4 honours the total
and, because the motif is unindexed, chooses the islands' placement; the output CIF's per-residue
provenance records where each landed (parsed by the ingest step).

Emitted by hand (no PyYAML dep) because the spec needs bare YAML booleans/ints as mapping keys
(`true:`, `1:`), which JSON cannot express. INGEST (parsing the generated CIF back into
coords+sequence+epitope positions) is a TODO pinned to a real motif run's output layout.
"""
from __future__ import annotations

from dataclasses import dataclass

from episcaf_service.rfd4.config import DEFAULT, PROTEIN_LENGTH, GenerationConfig
from episcaf_service.rfd4.topology import Island, Topology, plan_topology


def _residue_range(island: Island) -> str:
    return f"[{island.start}]" if island.length == 1 else f"[{island.start}-{island.end}]"


def residue_select(island: Island) -> str:
    """AtomWorks selection for a whole island: sel('<chain>/*/<range>')."""
    return f"sel('{island.chain}/*/{_residue_range(island)}')"


def backbone_select(island: Island) -> str:
    """AtomWorks selection for an island's backbone atoms: sel('<chain>/*/<range>/[N,CA,C,O]')."""
    return f"sel('{island.chain}/*/{_residue_range(island)}/[N,CA,C,O]')"


def contig_string(topology: Topology) -> str:
    """Interleave fixed-length expandable scaffold segments with the preserved epitope islands."""
    segs = topology.scaffold_segments
    islands = topology.islands
    tokens: list[str] = []
    for i, island in enumerate(islands):
        g = segs[i]
        tokens.append(f"E{i}:{g}-{g}")
        rng = f"{island.start}-{island.end}" if island.length > 1 else f"{island.start}"
        tokens.append(f"{island.chain}{rng}")
    g = segs[-1]
    tokens.append(f"E{len(islands)}:{g}-{g}")
    return ", ".join(tokens)


@dataclass
class ConditionSpec:
    name: str
    input_path: str
    contig: str
    crdgrp: dict            # {group_id(int): [select_str, ...]}
    seq_fixed: list         # [select_str, ...] held under C_SEQ true  (empty -> omit C_SEQ)


def condition_spec(antigen_input: str, epitope, cfg: GenerationConfig = DEFAULT.generation,
                   protein_length: int | None = None, name: str = "scaffold") -> ConditionSpec:
    """Build the motif-scaffolding condition spec for an antigen file + epitope islands."""
    top = epitope if isinstance(epitope, Topology) else plan_topology(epitope, protein_length or PROTEIN_LENGTH)
    islands = top.islands
    crdgrp = {1: [backbone_select(i) for i in islands]} if cfg.unindex_motif else {}
    # indexed fallback would fix absolute backbone coords instead; unindex is the default/validated path.
    seq_fixed = [residue_select(i) for i in islands] if cfg.fix_motif_sequence else []
    return ConditionSpec(name=name, input_path=antigen_input, contig=contig_string(top),
                         crdgrp=crdgrp, seq_fixed=seq_fixed)


def _yaml_key(k) -> str:
    if isinstance(k, bool):          # bool before int: True/False are ints in Python
        return "true" if k else "false"
    if isinstance(k, int):
        return str(k)
    return str(k)


def to_yaml(spec: ConditionSpec) -> str:
    """Emit the condition spec as YAML with bare boolean/int mapping keys."""
    lines = [f"{spec.name}:",
             f"  input: {spec.input_path}",
             f'  contig: "{spec.contig}"',
             "  conditions:"]
    if spec.crdgrp:
        lines.append("    C_CRDGRP:")
        for gid, selects in spec.crdgrp.items():
            lines.append(f"      {_yaml_key(gid)}:")
            for sel in selects:
                lines.append(f'        - select: "{sel}"')
    if spec.seq_fixed:
        lines.append("    C_SEQ:")
        lines.append(f"      {_yaml_key(True)}:")
        for sel in spec.seq_fixed:
            lines.append(f'        - select: "{sel}"')
    return "\n".join(lines) + "\n"


def write_spec(spec: ConditionSpec, path) -> str:
    text = to_yaml(spec)
    with open(path, "w") as fh:
        fh.write(text)
    return text


def generate_argv(spec_path: str, out_dir: str, num_replicates: int,
                  hf_ckpt: str = "nvidia/RFD4-Proteina", drop_metrics: bool = True) -> list[str]:
    """The `rfd4-generate` command (metrics off; we score with AF2 ourselves).

    TODO: map GenerationConfig.sampling_effort (thinking_*) onto RFD4's recycle_strategy preset
    (cookbook uses e.g. recycle_strategy=high/front_loaded) once confirmed against a real run.
    """
    argv = ["rfd4-generate", f"inputs={spec_path}", f"out_dir={out_dir}",
            f"num_replicates={num_replicates}", f"hf_ckpt_path={hf_ckpt}"]
    if drop_metrics:
        argv.append("~metrics")
    return argv
