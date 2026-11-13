#!/usr/bin/env python3
"""Build an RFd3 designs table from growPatch candidate epitopes -- the no-antibody arm's input.

growPatch (episcaf_analysis/growpatch) enumerates candidate conformational epitopes on an antigen
surface when there is no antibody to tell us where the epitope is. This script turns each candidate
into a motif-scaffolding contig that fixes ALL of its islands together -- presenting the whole
discontinuous epitope as one unit, the way an antibody sees it -- with free scaffold flanks and
inter-island linkers, to a fixed total length. The output is the designs table
`stage01_compile_contigs` consumes (columns id, contig_id, contig_string, contig_length), so growPatch
candidates flow into the same RFd3 -> MPNN -> AF3 pipeline as the known-antibody arm.

This REPLACES the old tiled-12-mer no-antibody approach: tiling was sequence-linear and could only
propose contiguous stretches, so it could never present a conformational epitope; a growPatch contig
fixes residues that are close in 3D but scattered in sequence (and across chains).

Contig format (RFd3 after utils.contig_to_rfd3 turns '/' into ','):
    free / {chain}{a}-{b} / linker / {chain}{c}-{d} / free       e.g. "12-12/A298-307/6-6/A330-340/9-9"
The free spans sum to (total_len - sum of island spans); V variants randomize the split (seeded,
reproducible), moving the epitope around in the construct. Motif residues keep their native
coordinates from the input antigen PDB, so the islands' relative 3D arrangement is preserved.

Usage (one antigen):
    python -m episcaf_pipeline.build_growpatch_designs --cif antigen.cif --chain A \
        --pdb-id 8pww --total-len 103 --variants 4 --out results/growpatch/designs_8pww.csv
Batch: --manifest a CSV with columns pdb_id,chain,cif (one antigen per row).
"""
from __future__ import annotations

import argparse
import csv
import os
import random
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from episcaf_analysis.growpatch import read as R   # noqa: E402

_ISL = re.compile(r"^([A-Za-z0-9]+):(-?\d+)(?:-(-?\d+))?$")
TOTAL_LEN = 103          # PepSeq max construct length


def parse_epitope(epitope):
    """'A:298-307;B:37-47' -> [(chain, a, b), ...] sorted by (chain, a)."""
    out = []
    for tok in epitope.split(";"):
        tok = tok.strip()
        if not tok:
            continue
        m = _ISL.match(tok)
        if not m:
            raise ValueError(f"bad island token {tok!r} in {epitope!r}")
        ch, a, b = m.group(1), int(m.group(2)), int(m.group(3) or m.group(2))
        out.append((ch, a, b))
    return sorted(out)


def free_splits(budget, n_free, variants, rng):
    """`variants` ways to split `budget` residues into `n_free` free spans, each >= 1.

    V=1 is the even split (deterministic). V>1 draws distinct random compositions (seeded), which
    slide the epitope around in the construct. Falls back to fewer if the budget can't make V
    distinct compositions.
    """
    if budget < n_free:
        return []
    if variants <= 1:
        base, rem = divmod(budget, n_free)
        return [[base + (1 if i < rem else 0) for i in range(n_free)]]
    seen, out, tries = set(), [], 0
    max_tries = variants * 50
    while len(out) < variants and tries < max_tries:
        tries += 1
        if n_free == 1:
            comp = (budget,)
        else:
            cuts = sorted(rng.sample(range(1, budget), n_free - 1))
            comp = tuple(b - a for a, b in zip([0] + cuts, cuts + [budget]))
        if comp not in seen:
            seen.add(comp)
            out.append(list(comp))
    return out


def contigs_for(islands, total_len, variants, rng):
    """Contig strings presenting all `islands` (list of (chain,a,b)) with free flanks/linkers."""
    spans = [b - a + 1 for _, a, b in islands]
    budget = total_len - sum(spans)
    n_free = len(islands) + 1
    out = []
    for v, free in enumerate(free_splits(budget, n_free, variants, rng)):
        segs = []
        for i, (ch, a, b) in enumerate(islands):
            segs.append(f"{free[i]}-{free[i]}")
            segs.append(f"{ch}{a}-{b}")
        segs.append(f"{free[-1]}-{free[-1]}")
        out.append((v, "/".join(segs)))
    return out


def build_rows(pdb_id, chains, candidates, total_len, variants, seed):
    """designs rows for one antigen's growPatch candidates."""
    rng = random.Random(seed)
    rows, skipped = [], []
    for p in candidates:
        islands = parse_epitope(p["epitope"])
        cons = contigs_for(islands, total_len, variants, rng)
        if not cons:
            skipped.append((p["epitope"], sum(b - a + 1 for _, a, b in islands)))
            continue
        base_id = f"{pdb_id.lower()}{''.join(chains)}_gp{p['patch']}"
        for cid, contig_string in cons:
            rows.append({
                "id": base_id,
                "contig_id": cid,
                "contig_string": contig_string,
                "contig_length": f"{total_len}-{total_len}",
                "epitope_chunks": p["n_islands"],
                "source": "growpatch",
                "pdb_id": pdb_id.lower(),
                "chains": "".join(chains),
                "growpatch_epitope": p["epitope"],
                "growpatch_methods": "+".join(p["methods"]),
                "growpatch_n_res": p["n_res"],
                "variant": cid,
            })
    return rows, skipped


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--cif", help="one antigen structure (mmCIF)")
    ap.add_argument("--chain", default="A", help="antigen chain(s): 'A' or 'A,B'")
    ap.add_argument("--pdb-id", default="", help="id for the design rows (default: from --cif name)")
    ap.add_argument("--manifest", help="CSV with columns pdb_id,chain,cif (batch mode)")
    ap.add_argument("--total-len", type=int, default=TOTAL_LEN)
    ap.add_argument("--variants", type=int, default=4,
                    help="flank-split contigs per candidate (V); V=1 is the centered split")
    ap.add_argument("--sasa-cut", type=float, default=0.2)
    ap.add_argument("--footprint-radius", type=float, default=14.0)
    ap.add_argument("--methods", default="growpatch,eyeview",
                    help="enumerators to union (comma-separated)")
    ap.add_argument("--top-k", type=int, default=0, help="cap candidates per antigen (0 = all)")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--out", required=True)
    a = ap.parse_args(argv)

    if a.manifest:
        jobs = [(r["pdb_id"], r["chain"], r["cif"]) for r in csv.DictReader(open(a.manifest))]
    elif a.cif:
        jobs = [(a.pdb_id or os.path.basename(a.cif).split(".")[0], a.chain, a.cif)]
    else:
        ap.error("give --cif or --manifest")

    methods = tuple(m.strip() for m in a.methods.split(",") if m.strip())
    all_rows, total_skipped, n_antigens = [], 0, 0
    for pdb_id, chain, cif in jobs:
        from episcaf_analysis.growpatch.read import parse_chains
        chains = parse_chains(chain)
        cands, _ = R.enumerate_candidates(cif, chain, methods=methods, sasa_cut=a.sasa_cut,
                                          footprint_radius=a.footprint_radius,
                                          top_k=(a.top_k or None))
        rows, skipped = build_rows(pdb_id, chains, cands, a.total_len, a.variants, a.seed)
        all_rows += rows
        total_skipped += len(skipped)
        n_antigens += 1
        print(f"{pdb_id} [{chain}]: {len(cands)} candidates -> {len(rows)} contig rows"
              + (f"  ({len(skipped)} too big for {a.total_len}-mer)" if skipped else ""))

    if not all_rows:
        sys.exit("no contigs built")
    cols = ["id", "contig_id", "contig_string", "contig_length", "epitope_chunks", "source",
            "pdb_id", "chains", "growpatch_epitope", "growpatch_methods", "growpatch_n_res", "variant"]
    os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
    if a.out.endswith(".parquet"):
        import pandas as pd                                   # stage01_compile_contigs reads parquet
        pd.DataFrame(all_rows, columns=cols).to_parquet(a.out, index=False)
    else:
        with open(a.out, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=cols)
            w.writeheader()
            w.writerows(all_rows)
    n_design = len(all_rows)
    print(f"\nwrote {a.out}: {n_antigens} antigen(s), {n_design} contig rows "
          f"(x8 RFD3 x8 MPNN = {n_design*64:,} designs at full fanout)"
          + (f"; {total_skipped} candidates too big for the construct" if total_skipped else ""))


if __name__ == "__main__":
    main()
