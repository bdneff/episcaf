"""Epitope-island parsing and scaffold topology.

Two jobs, both ported from tamarind.R (sample_rfd3_topology ~2357/2386, position_spans,
rfd4_topology_from_map ~2692):

  plan_topology(islands)      the NOMINAL layout used to BUILD the RFD4 input: the epitope
                              held in native order, the (103 - epitope) scaffold split as
                              evenly as possible over (n_islands + 1) gaps, extras to the
                              leading gaps. RFD4 honours only the total length (unindexed), so
                              per-gap sizes are a starting point, not a constraint.

  topology_from_map(map_str)  the ACTUAL per-design layout, parsed back from RFD4's
                              motif_residue_map output, since an unindexed run lets the model
                              reorder and place the islands. A map that does not add up to the
                              epitope length, or that runs past the protein, is rejected
                              (BadResidueMap) and the design is dropped.
"""
from __future__ import annotations

import re
from dataclasses import dataclass

from episcaf_service.rfd4.config import PROTEIN_LENGTH


class BadResidueMap(ValueError):
    """RFD4's motif_residue_map is inconsistent with the epitope / protein length."""


@dataclass(frozen=True)
class Island:
    chain: str
    start: int
    end: int

    @property
    def length(self) -> int:
        return self.end - self.start + 1


@dataclass(frozen=True)
class Topology:
    islands: tuple[Island, ...]           # in native order
    scaffold_segments: tuple[int, ...]    # n_islands + 1 gap lengths, in layout order
    epitope_positions: tuple[int, ...]    # 1-based positions in the 1..L design the epitope occupies
    island_order: tuple[int, ...]         # index into `islands` for each island as laid out, in order
    protein_length: int

    @property
    def scaffold_positions(self) -> tuple[int, ...]:
        epi = set(self.epitope_positions)
        return tuple(p for p in range(1, self.protein_length + 1) if p not in epi)

    @property
    def total_epitope_length(self) -> int:
        return len(self.epitope_positions)


_ISLAND_RE = re.compile(r"^\s*([A-Za-z0-9]+)\s*:\s*(-?\d+)(?:\s*-\s*(-?\d+))?\s*$")


def parse_islands(epitope: str) -> tuple[Island, ...]:
    """Parse 'A:205-216;A:278-292' (or a single 'A:252') into ordered Islands.

    Semicolon-separated; each island is '<chain>:<start>[-<end>]' in the PDB's own numbering.
    """
    if not epitope or not epitope.strip():
        raise ValueError("empty epitope specification")
    islands: list[Island] = []
    for token in epitope.split(";"):
        if not token.strip():
            continue
        m = _ISLAND_RE.match(token)
        if not m:
            raise ValueError(f"cannot parse epitope island {token!r} (want '<chain>:<start>[-<end>]')")
        chain, start_s, end_s = m.group(1), m.group(2), m.group(3)
        start = int(start_s)
        end = int(end_s) if end_s is not None else start
        if end < start:
            raise ValueError(f"island {token!r} has end < start")
        islands.append(Island(chain=chain, start=start, end=end))
    if not islands:
        raise ValueError(f"no islands parsed from {epitope!r}")
    return tuple(islands)


def even_split(total: int, n: int) -> tuple[int, ...]:
    """Split `total` into `n` parts as evenly as possible, extra residues to the leading parts.

    Ported from tamarind.R: rep(total %/% n, n); extra = total %% n; parts[1:extra] += 1.
    e.g. even_split(76, 3) -> (26, 25, 25).
    """
    if n <= 0:
        raise ValueError("n must be positive")
    base, extra = divmod(total, n)
    return tuple(base + (1 if i < extra else 0) for i in range(n))


def plan_topology(islands, protein_length: int = PROTEIN_LENGTH) -> Topology:
    """Nominal de-novo topology: even-split scaffold around the epitope, native island order."""
    if isinstance(islands, str):
        islands = parse_islands(islands)
    islands = tuple(islands)
    total_epitope = sum(i.length for i in islands)
    scaffold_length = protein_length - total_epitope
    if scaffold_length < len(islands) + 1:
        raise ValueError(
            f"epitope length {total_epitope} leaves no room for a scaffold in {protein_length} residues"
        )
    segments = even_split(scaffold_length, len(islands) + 1)

    epitope_positions: list[int] = []
    pos = 0
    for k, island in enumerate(islands):
        pos += segments[k]                      # leading gap
        epitope_positions.extend(range(pos + 1, pos + island.length + 1))
        pos += island.length
    pos += segments[-1]                          # trailing gap
    assert pos == protein_length, (pos, protein_length)

    return Topology(
        islands=islands,
        scaffold_segments=segments,
        epitope_positions=tuple(epitope_positions),
        island_order=tuple(range(len(islands))),
        protein_length=protein_length,
    )


_MAP_SIDE_RE = re.compile(r"^\s*([A-Za-z]*)\s*(\d+)(?:-(\d+))?\s*$")


def _parse_range(side: str) -> tuple[int, int]:
    m = _MAP_SIDE_RE.match(side)
    if not m:
        raise BadResidueMap(f"cannot parse residue range {side!r}")
    start = int(m.group(2))
    end = int(m.group(3)) if m.group(3) is not None else start
    if end < start:
        raise BadResidueMap(f"range {side!r} has end < start")
    return start, end


def topology_from_map(map_str: str, total_epitope_length: int,
                      protein_length: int = PROTEIN_LENGTH) -> dict:
    """Parse an RFD4 motif_residue_map like 'R312-328>R28-44;R362-364>R83-85'.

    Left of '>' is the native (source) range, right is where it landed in the design (dest).
    Returns dict(epitope_positions, scaffold_segments, island_order_by_dest). Raises BadResidueMap
    if the mapped positions do not sum to the epitope length or run past the protein.
    """
    entries = []
    for part in map_str.split(";"):
        part = part.strip()
        if not part:
            continue
        if ">" not in part:
            raise BadResidueMap(f"map entry {part!r} missing '>'")
        src_s, dst_s = part.split(">", 1)
        src = _parse_range(src_s)
        dst = _parse_range(dst_s)
        if (dst[1] - dst[0]) != (src[1] - src[0]):
            raise BadResidueMap(f"map entry {part!r}: source and dest lengths differ")
        entries.append((src, dst))
    if not entries:
        raise BadResidueMap("empty motif_residue_map")

    entries.sort(key=lambda e: e[1][0])          # by dest start = layout order
    epitope_positions: list[int] = []
    for _src, dst in entries:
        epitope_positions.extend(range(dst[0], dst[1] + 1))

    if len(epitope_positions) != total_epitope_length:
        raise BadResidueMap(
            f"mapped {len(epitope_positions)} residues, expected {total_epitope_length}"
        )
    if epitope_positions[-1] > protein_length or epitope_positions[0] < 1:
        raise BadResidueMap(f"mapped positions run outside 1..{protein_length}")

    # scaffold_segments = leading gap, gaps between islands, trailing gap
    starts = [dst[0] for _s, dst in entries]
    ends = [dst[1] for _s, dst in entries]
    segments = [starts[0] - 1]
    for i in range(1, len(entries)):
        segments.append(starts[i] - ends[i - 1] - 1)
    segments.append(protein_length - ends[-1])
    if any(g < 0 for g in segments):
        raise BadResidueMap(f"overlapping islands in map {map_str!r}")

    return {
        "epitope_positions": tuple(epitope_positions),
        "scaffold_segments": tuple(segments),
        "island_order_by_dest": tuple(range(len(entries))),
    }
