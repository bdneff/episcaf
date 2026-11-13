"""growPatch I/O: read an antigen structure and its surface, and run the enumerator on it.

This is the adapter between a structure file and the pure algorithm in `core.py`. It ports the
reading and relSASA parts of John Altin's `segmentation/growPatch.py` (AltinLab/epitope-universe @
upstream/main), faithfully enough that `enumerate_patches` reproduces his candidate set
(cross-checked in `tests/test_growpatch_core.py`).

Two jobs:
  1. read the antigen chain(s) into the residue list `core` consumes -- standard residues with a
     complete N/CA/C/O backbone (RFdiffusion3 motifs and the EpiScaf clash pre-filter need all four
     for every fixed residue), each carrying its heavy-atom coords and a side-chain-centroid `rep`;
  2. compute per-residue relative SASA with freesasa, normalized by the Tien et al. 2013 theoretical
     max-ASA table, so the surface set is `{residue : relSASA >= SASA_CUT}`.

Residues are keyed `(chain, slot)` (see `core`). A "slot" gives every insertion-coded residue its
own integer position (170, 170A, 170B -> slots 170, 171, 172); the slot<->label bookkeeping lives on
the `Antigen` object (not in module globals, so several structures can be read in one process). Only
formatting and the viewer selections need labels; the algorithm sees only integer slots.

relSASA engine note: John uses freesasa + Tien 2013, and we match that exactly for parity. The b-cell
project computes the identical Tien-2013 relSASA with Biopython `Bio.PDB.SASA.ShrakeRupley` instead,
which would drop the freesasa dependency (biopython is already everywhere freesasa is not). That is a
measured swap for later, not a silent one -- see `PARAMS.md`; the two engines' areas are not
byte-identical, so the surface set at a fixed cut could shift and would need re-checking.
"""
from __future__ import annotations

import re
import string
from dataclasses import dataclass, field

import numpy as np

from . import core
from . import eyeview as _eyeview

# Tien et al. 2013 theoretical max ASA (A^2) -- the relSASA denominator (verbatim from growPatch.py;
# independently reproduced in bcell_epitope/.../recompute_apo_baselines.py).
MAXASA = {
    "ALA": 129.0, "ARG": 274.0, "ASN": 195.0, "ASP": 193.0, "CYS": 167.0,
    "GLN": 225.0, "GLU": 223.0, "GLY": 104.0, "HIS": 224.0, "ILE": 197.0,
    "LEU": 201.0, "LYS": 236.0, "MET": 224.0, "PHE": 240.0, "PRO": 159.0,
    "SER": 155.0, "THR": 172.0, "TRP": 285.0, "TYR": 263.0, "VAL": 174.0,
}

BACKBONE = {"N", "CA", "C", "O"}

_RANGE_RE = re.compile(r"^(-?\d+[A-Za-z]*)(?:-(-?\d+[A-Za-z]*))?$")
_RES_RE = re.compile(r"^(-?\d+)([A-Za-z]*)$")


# ------------------------------------------------------------------ slot <-> label bookkeeping

def _ins_ord(icode):
    """insertion code -> rank ('' 0, 'A' 1, ..., 'AA' 27), as in reference/sabdb.R (port)."""
    v = 0
    for ch in icode.strip().upper():
        v = v * 26 + (ord(ch) - 64)
    return v


def _slot_table(ch):
    """Slots of every standard residue of a Bio.PDB chain (port of growPatch._slot_table).

    Hetero residues are not numbered and take no slot. Returns
    ({(number, icode): slot}, {slot: label}, {label: slot}, [(number, ins_ord)] of inserted residues).
    """
    rs = sorted(((r.id[1], r.id[2].strip()) for r in ch if r.id[0] == " "),
                key=lambda t: (t[0], _ins_ord(t[1])))
    slot_of, shift = {}, 0
    labels, slots, ins = {}, {}, []
    for num, ic in rs:
        if ic:
            shift += 1
            ins.append((num, _ins_ord(ic)))
        slot = num + shift
        slot_of[(num, ic)] = slot
        labels[slot] = f"{num}{ic}"
        slots[f"{num}{ic}"] = slot
    return slot_of, labels, slots, ins


# ------------------------------------------------------------------ the antigen

@dataclass
class Antigen:
    """An antigen's tiled chains, read into what `core` needs, plus the label bookkeeping.

    residues  -- list of dicts {chain, resno(=slot), key, resname, atoms (N,3 heavy), rep}; the exact
                 shape `core.build_adjacency` / `core.candidate_patches` expect.
    by_key    -- {(chain, slot): residue index}
    idx_key   -- {residue index: (chain, slot)}
    present   -- {chain: frozenset(slots with a complete backbone)} -- gap-fill skips the rest
    labels    -- {chain: {slot: label}} for formatting ('170B'); the number itself if unknown
    _slots    -- {chain: {label: slot}}  _ins -- {chain: [(number, ins_ord)]}  (for parse_epitope)
    relsasa   -- {(chain, slot): relSASA}; filled by read_antigen
    """
    chains: list
    residues: list
    by_key: dict
    idx_key: dict
    present: dict
    labels: dict
    _slots: dict = field(default_factory=dict)
    _ins: dict = field(default_factory=dict)
    relsasa: dict = field(default_factory=dict)

    @property
    def multi(self):
        return len(self.chains) > 1

    def label(self, chain, slot):
        return self.labels.get(chain, {}).get(slot, str(slot))

    def rs(self, key):
        return self.relsasa.get(key, 0.0)

    def surface_indices(self, sasa_cut):
        """Residue indices whose relSASA >= sasa_cut (the surface set the enumerator seeds from)."""
        return frozenset(i for i, r in enumerate(self.residues) if self.rs(r["key"]) >= sasa_cut)

    # -- formatting (labels back from slots) --
    def fmt_islands(self, keys):
        """'205-216; 274-292' for one chain, 'A:298-307; B:37-47' for several."""
        pre = (lambda c: f"{c}:") if self.multi else (lambda c: "")
        return "; ".join(
            f"{pre(c)}{self.label(c, a)}-{self.label(c, b)}" if a != b else f"{pre(c)}{self.label(c, a)}"
            for c, a, b in core.islands(keys))

    def fmt_epitope(self, keys):
        """EpiScaf contig string, always chain-prefixed: 'A:298-307;B:37-47'."""
        return ";".join(
            f"{c}:{self.label(c, a)}-{self.label(c, b)}" if a != b else f"{c}:{self.label(c, a)}"
            for c, a, b in core.islands(keys))

    def vmd_selection(self, keys):
        """A VMD atomselect string highlighting a patch: 'chain A and resid 298 to 307 ...'.

        Insertion-coded residues use the bare number (VMD `resid` is numeric); a patch that spans an
        insertion is still selected by its number range. For the headless renderer / interactive view.
        """
        parts = []
        for c, a, b in core.islands(keys):
            na = re.match(r"-?\d+", self.label(c, a)).group()
            nb = re.match(r"-?\d+", self.label(c, b)).group()
            rng = f"resid {na} to {nb}" if a != b else f"resid {na}"
            parts.append(f"(chain {c} and {rng})")
        return " or ".join(parts)

    # -- optional: parse a reference/known epitope spec to (chain, slot) keys (benchmark / protect) --
    def _slot_of_label(self, chain, num, icode=""):
        label = f"{num}{icode}"
        known = self._slots.get(chain, {}).get(label)
        if known is not None:
            return known
        key = (num, _ins_ord(icode))
        return num + sum(1 for t in self._ins.get(chain, ()) if t < key) + (1 if icode else 0)

    def parse_epitope(self, spec):
        """Residue spec -> set of (chain, slot) (port of growPatch.parse_resno_spec).

        Bare list for a single chain ('278-280,288-289'); chain-prefixed for several
        ('A:298-299,303-307;B:37,39-40'). Insertion codes allowed ('60B-61'). Used to protect a known
        epitope from short-island trimming, and by the parity cross-check.
        """
        out, cur = set(), (self.chains[0] if len(self.chains) == 1 else None)
        for tok in spec.replace(";", ",").split(","):
            tok = tok.strip()
            if not tok:
                continue
            if ":" in tok:
                cur, tok = (x.strip() for x in tok.split(":", 1))
                if cur not in self.chains:
                    raise ValueError(f"epitope names chain {cur!r}, not in chains {self.chains}")
            elif cur is None:
                raise ValueError(f"residue token {tok!r} has no chain; use 'A:298-307;B:37-47'")
            m = _RANGE_RE.match(tok)
            if not m:
                raise ValueError(f"bad residue token: {tok!r}")
            (na, ia), (nb, ib) = (_RES_RE.match(m.group(1)).groups(),
                                  _RES_RE.match(m.group(2) or m.group(1)).groups())
            na, nb = int(na), int(nb)
            if not ia and not ib:
                out |= {(cur, self._slot_of_label(cur, n)) for n in range(na, nb + 1)}
            else:
                out |= {(cur, n) for n in range(self._slot_of_label(cur, na, ia),
                                                self._slot_of_label(cur, nb, ib) + 1)}
        return out


# ------------------------------------------------------------------ reading

def parse_chains(spec):
    """'A' / 'A,B' / 'A+B' / 'A;B' -> ['A', 'B'] (order kept, duplicates dropped)."""
    out = []
    for c in re.split(r"[,+;\s]+", str(spec).strip()):
        if c and c not in out:
            out.append(c)
    if not out:
        raise ValueError("no chain given")
    return out


def _read_chain(model, chain_id):
    """Standard, complete-backbone residues of one chain -> (residue dicts, label maps, dropped)."""
    if chain_id not in model:
        raise KeyError(f"chain {chain_id} not in structure")
    ch = model[chain_id]
    slot_of, labels, slots, ins = _slot_table(ch)
    residues, incomplete = [], []
    for r in ch:
        if r.id[0] != " " or r.resname not in MAXASA:
            continue
        slot = slot_of[(r.id[1], r.id[2].strip())]
        if not BACKBONE <= {a.get_name() for a in r}:
            incomplete.append(f"{r.id[1]}{r.id[2].strip()}")
            continue
        heavy = np.array([a.coord for a in r if a.element != "H"], dtype=float)
        if heavy.shape[0] == 0:
            continue
        sc = np.array([a.coord for a in r
                       if a.element != "H" and a.get_name() not in ("N", "C", "O")], dtype=float)
        rep = sc.mean(axis=0) if sc.shape[0] else heavy.mean(axis=0)
        residues.append(dict(chain=chain_id, resno=slot, key=(chain_id, slot),
                             resname=r.resname, atoms=heavy, rep=rep))
    residues.sort(key=lambda d: d["resno"])
    return residues, labels, slots, ins, incomplete


def _relsasa_group(cif_path, group):
    """freesasa on the chains in `group` together -> {(chain, slot): relSASA} (port)."""
    import os
    import tempfile
    import freesasa
    from Bio.PDB import MMCIFParser, PDBIO, Select
    from Bio.PDB.Chain import Chain
    from Bio.PDB.Model import Model
    from Bio.PDB.Structure import Structure

    model = next(iter(MMCIFParser(QUIET=True).get_structure("x", cif_path)))

    # the temp file is PDB (one-character chain ids): write the group under single-char names and
    # read back under the real names, so a chain called 'BBB' does not break.
    codes = string.ascii_uppercase + string.ascii_lowercase + string.digits
    if len(group) > len(codes):
        raise ValueError(f"cannot measure more than {len(codes)} chains together")
    rename = {c: codes[i] for i, c in enumerate(group)}
    slots_c, new, nm = {}, Structure("x"), Model(0)
    new.add(nm)
    for c in group:
        if c not in model:
            continue
        slots_c[c] = _slot_table(model[c])[2]
        nc = Chain(rename[c])
        for r in model[c]:
            if r.id[0] == " " and r.resname in MAXASA:
                nc.add(r.copy())
        nm.add(nc)

    class NoH(Select):
        def accept_atom(self, a):
            return a.element != "H"

    tmp = tempfile.NamedTemporaryFile(suffix=".pdb", delete=False, mode="w")
    tmp.close()
    io = PDBIO()
    io.set_structure(new)
    io.save(tmp.name, NoH())
    areas = freesasa.calc(freesasa.Structure(tmp.name)).residueAreas()
    os.unlink(tmp.name)

    out = {}
    for c in group:
        if rename[c] not in areas:
            continue
        for label, ra in areas[rename[c]].items():
            slot = slots_c[c].get(label.strip())
            mx = MAXASA.get(ra.residueType)
            if slot is not None and mx:
                out[(c, slot)] = ra.total / mx
    return out


def compute_relsasa(cif_path, chains, scope="antigen", context=()):
    """Relative SASA -> {(chain, slot): relSASA} (port of growPatch.compute_relsasa).

    scope 'antigen' (default): all `chains` together, so inter-chain interface residues read as
        buried. `context` adds further antigen chains to the calculation that are not tiled (e.g. a
        trimer copy carrying no epitope) -- they bury interface residues but are not reported.
    scope 'chain': each chain in isolation (the convention of the SASA-filtered epitope table), so
        interface residues read as exposed.
    """
    if scope == "chain":
        out = {}
        for c in chains:
            out.update(_relsasa_group(cif_path, [c]))
        return out
    group = list(chains) + [c for c in context if c not in chains]
    return {k: v for k, v in _relsasa_group(cif_path, group).items() if k[0] in chains}


def read_antigen(cif_path, chains, sasa_scope="antigen", context=(), relsasa=True, verbose=False):
    """Read antigen chain(s) into an `Antigen`, with relSASA computed (unless relsasa=False).

    `chains` is a chain id, a list, or a spec string ('A,B'). Raises if a chain has no standard,
    complete-backbone residue. The returned object is everything the enumerator needs.
    """
    from Bio.PDB import MMCIFParser

    chains = parse_chains(chains) if isinstance(chains, str) else list(chains)
    model = next(iter(MMCIFParser(QUIET=True).get_structure("x", cif_path)))

    residues, labels, slots_m, ins_m = [], {}, {}, {}
    for c in chains:
        rc, lab_c, slot_c, ins_c, incomplete = _read_chain(model, c)
        if not rc:
            raise ValueError(f"chain {c}: no standard residues with a complete N/CA/C/O backbone")
        labels[c], slots_m[c], ins_m[c] = lab_c, slot_c, ins_c
        if verbose:
            print(f"chain {c}: {len(rc)} standard residues "
                  f"({lab_c.get(rc[0]['resno'])}-{lab_c.get(rc[-1]['resno'])})"
                  + (f"; dropped {len(incomplete)} incomplete-backbone" if incomplete else ""))
        residues += rc

    by_key = {r["key"]: i for i, r in enumerate(residues)}
    idx_key = {i: r["key"] for i, r in enumerate(residues)}
    present = {c: frozenset(n for (cc, n) in by_key if cc == c) for c in chains}
    ag = Antigen(chains=chains, residues=residues, by_key=by_key, idx_key=idx_key,
                 present=present, labels=labels, _slots=slots_m, _ins=ins_m)
    if relsasa:
        ctx = parse_chains(context) if isinstance(context, str) and context.strip() else list(context)
        ag.relsasa = compute_relsasa(cif_path, chains, sasa_scope, ctx)
    return ag


# ------------------------------------------------------------------ the production entry

def enumerate_patches(cif_path, chains, sasa_cut=core.SASA_CUT, patch_size=core.PATCH_SIZE,
                      sasa_scope="antigen", context=(), gap=core.GAP, min_res=core.MIN_RES,
                      max_islands=core.MAX_ISLANDS, min_island_len=core.MIN_ISLAND_LEN,
                      contact=core.CONTACT, grow_link=core.GROW_LINK, nms_jaccard=core.NMS_JACCARD,
                      protect_epitope=None, antigen=None, verbose=False):
    """Enumerate candidate epitope patches for an antigen -- the no-antibody production entry.

    Reads the structure (or takes a prebuilt `antigen`), builds the surface set at `sasa_cut`, runs
    the seed-and-grow enumerator (`core.candidate_patches`), and returns the kept candidates as
    dicts, NMS-ordered (fewer islands, then larger):

        {patch, n_raw, n_res, n_islands, epitope, islands, keys, vmd_selection}

    where `epitope` is the EpiScaf contig string ('A:298-307;B:37-47'). `protect_epitope` (a spec
    string or a set of (chain, slot) keys) shields a known epitope's residues from short-island
    trimming; production leaves it None (nothing protected). Pass `antigen` to reuse a read + relSASA
    across several (sasa_cut, patch_size) without re-reading the file.
    """
    ag = antigen if antigen is not None else read_antigen(
        cif_path, chains, sasa_scope, context, relsasa=True, verbose=verbose)

    protect = frozenset()
    if protect_epitope:
        protect = (ag.parse_epitope(protect_epitope) if isinstance(protect_epitope, str)
                   else frozenset(protect_epitope))

    adj = core.build_adjacency(ag.residues, contact=contact)
    surf_idx = ag.surface_indices(sasa_cut)
    kept = core.candidate_patches(
        ag.residues, adj, surf_idx, ag.idx_key, ag.present,
        size=patch_size, link=grow_link, gap=gap, min_res=min_res, max_islands=max_islands,
        min_island_len=min_island_len, protect=protect, nms_jaccard=nms_jaccard)

    out = []
    for k, (raw_r, gf_r, nisl) in enumerate(kept):
        out.append(dict(
            patch=k, n_raw=len(raw_r), n_res=len(gf_r), n_islands=nisl,
            epitope=ag.fmt_epitope(gf_r), islands=ag.fmt_islands(gf_r),
            keys=sorted(gf_r), vmd_selection=ag.vmd_selection(gf_r)))
    return out, ag


# ------------------------------------------------------------------ combined entry (both strategies)

def enumerate_candidates(cif_path, chains, methods=("growpatch", "eyeview"),
                         sasa_cut=core.SASA_CUT, sasa_scope="antigen", context=(),
                         gap=core.GAP, min_res=core.MIN_RES, max_islands=core.MAX_ISLANDS,
                         min_island_len=core.MIN_ISLAND_LEN,
                         patch_size=core.PATCH_SIZE, grow_link=core.GROW_LINK, contact=core.CONTACT,
                         footprint_radius=14.0, depth=None, n_dirs=150,
                         nms_jaccard=core.NMS_JACCARD, top_k=None,
                         protect_epitope=None, antigen=None, verbose=False):
    """Combined no-antibody epitope enumerator -- the production entry.

    Runs the chosen `methods` and unions their candidates, then a final Jaccard NMS collapses patches
    both strategies found (tagged as found by both -- a confidence signal). growPatch (seed-and-grow)
    and eyeview (antibody's-eye-view direction caps) are complementary: on the conformational
    benchmark the union reaches recall@0.8 ~0.92 vs ~0.80 for either alone (see
    `results/growpatch/ENUMERATOR_DECISION.md`).

    All the dials are user knobs:
      surface        : sasa_cut, sasa_scope, context
      patch shape    : patch_size + grow_link (growpatch), footprint_radius + depth (eyeview);
                       depth defaults to footprint_radius
      footprint      : gap (gap-fill size), min_res, max_islands, min_island_len
      exhaustiveness : nms_jaccard (lower -> less dedup -> more, more-overlapping candidates),
                       top_k (hard cap on the number returned, in NMS rank order)

    Returns (candidates, antigen). Each candidate adds `methods` (sorted list of the strategies that
    produced it) to the usual fields {patch, n_res, n_islands, epitope, islands, keys, vmd_selection}.
    """
    ag = antigen if antigen is not None else read_antigen(
        cif_path, chains, sasa_scope, context, relsasa=True, verbose=verbose)
    if depth is None:
        depth = footprint_radius
    protect = frozenset()
    if protect_epitope:
        protect = (ag.parse_epitope(protect_epitope) if isinstance(protect_epitope, str)
                   else frozenset(protect_epitope))

    pool = []
    if "growpatch" in methods:
        gp, _ = enumerate_patches(cif_path, chains, antigen=ag, sasa_cut=sasa_cut,
                                  patch_size=patch_size, gap=gap, min_res=min_res,
                                  max_islands=max_islands, min_island_len=min_island_len,
                                  contact=contact, grow_link=grow_link, nms_jaccard=nms_jaccard,
                                  protect_epitope=protect)
        pool += [(p, "growpatch") for p in gp]
    if "eyeview" in methods:
        ev = _eyeview.enumerate_eyeview(ag, sasa_cut=sasa_cut, n_dirs=n_dirs,
                                        footprint_radius=footprint_radius, depth=depth, gap=gap,
                                        min_res=min_res, max_islands=max_islands,
                                        min_island_len=min_island_len, nms_jaccard=nms_jaccard,
                                        protect=protect)
        pool += [(p, "eyeview") for p in ev]

    # cross-method Jaccard NMS: fewer islands, then larger, then sorted keys (deterministic)
    ordered = sorted(pool, key=lambda pm: (pm[0]["n_islands"], -pm[0]["n_res"], pm[0]["keys"]))
    kept = []
    for p, method in ordered:
        ks = set(p["keys"])
        dup = next((q for q in kept if core.jaccard(ks, set(q["keys"])) >= nms_jaccard), None)
        if dup is None:
            p = dict(p)
            p["_methods"] = {method}
            kept.append(p)
        else:
            dup["_methods"].add(method)
    if top_k:
        kept = kept[:top_k]

    out = [dict(patch=k, methods=sorted(p["_methods"]), n_res=p["n_res"], n_islands=p["n_islands"],
                epitope=p["epitope"], islands=p["islands"], keys=p["keys"],
                vmd_selection=p["vmd_selection"]) for k, p in enumerate(kept)]
    return out, ag


# ------------------------------------------------------------------ CLI

def _tag_from(cif_path, chains):
    import os
    m = re.search(r"pdb_0*([0-9a-z]{4})", os.path.basename(cif_path))
    return (m.group(1) if m else "epi") + "_" + "".join(chains)


def _main(argv=None):
    """CLI: enumerate candidate epitope patches for an antigen, optionally sweeping SASA x size."""
    import argparse
    import csv
    import os

    ap = argparse.ArgumentParser(description=_main.__doc__)
    ap.add_argument("--cif", required=True, help="antigen structure (mmCIF)")
    ap.add_argument("--chain", default="A", help="antigen chain(s) to tile: 'A' or 'A,B'")
    ap.add_argument("--sasa-scope", choices=("chain", "antigen"), default="antigen")
    ap.add_argument("--sasa-context-chains", default="",
                    help="further antigen chains in the relSASA but not tiled (e.g. 'C,G')")
    ap.add_argument("--sasa-cuts", default=str(core.SASA_CUT),
                    help="comma-separated relSASA surface cutoffs to sweep")
    ap.add_argument("--patch-sizes", default=str(core.PATCH_SIZE),
                    help="comma-separated raw patch sizes to sweep")
    ap.add_argument("--protect-epitope", default=None,
                    help="known-epitope spec protected from short-island trimming (benchmark only)")
    ap.add_argument("--gap", type=int, default=core.GAP)
    ap.add_argument("--min-res", type=int, default=core.MIN_RES)
    ap.add_argument("--max-islands", type=int, default=core.MAX_ISLANDS)
    ap.add_argument("--min-island-len", type=int, default=core.MIN_ISLAND_LEN)
    ap.add_argument("--contact", type=float, default=core.CONTACT)
    ap.add_argument("--grow-link", type=float, default=core.GROW_LINK)
    ap.add_argument("--nms-jaccard", type=float, default=core.NMS_JACCARD)
    ap.add_argument("--outdir", default=".")
    ap.add_argument("--tag", default="")
    a = ap.parse_args(argv)

    chains = parse_chains(a.chain)
    ag = read_antigen(a.cif, chains, a.sasa_scope, a.sasa_context_chains, relsasa=True, verbose=True)

    sasa_cuts = [float(x) for x in a.sasa_cuts.split(",")]
    sizes = [int(x) for x in a.patch_sizes.split(",")]
    grid, rows = [], []
    for cut in sasa_cuts:
        n_surf = len(ag.surface_indices(cut))
        for size in sizes:
            patches, _ = enumerate_patches(
                a.cif, chains, sasa_cut=cut, patch_size=size, gap=a.gap, min_res=a.min_res,
                max_islands=a.max_islands, min_island_len=a.min_island_len, contact=a.contact,
                grow_link=a.grow_link, nms_jaccard=a.nms_jaccard, protect_epitope=a.protect_epitope,
                antigen=ag)
            grid.append(dict(sasa_cut=cut, patch_size=size, n_surface=n_surf,
                             n_candidates=len(patches)))
            for p in patches:
                rows.append(dict(sasa_cut=cut, patch_size=size, **{k: p[k] for k in (
                    "patch", "n_raw", "n_res", "n_islands", "epitope", "islands", "vmd_selection")}))
            print(f"SASA_cut={cut}  patch_size={size}  surface={n_surf}  candidates={len(patches)}")

    os.makedirs(a.outdir, exist_ok=True)
    tag = a.tag or _tag_from(a.cif, chains)
    p1 = os.path.join(a.outdir, f"growpatch_patches_{tag}.csv")
    if rows:
        with open(p1, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
            w.writeheader()
            w.writerows(rows)
    p2 = os.path.join(a.outdir, f"growpatch_grid_{tag}.csv")
    with open(p2, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(grid[0].keys()))
        w.writeheader()
        w.writerows(grid)
    print(f"wrote {p1}\nwrote {p2}")


if __name__ == "__main__":
    _main()
