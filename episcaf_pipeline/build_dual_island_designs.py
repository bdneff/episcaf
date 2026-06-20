#!/usr/bin/env python3
"""
build_dual_island_designs.py  --  emit the per-island RFD3 design ledger for John's run.

John's plan: take the 46 dual-island (`epitope_chunks==2`) mAb epitopes, scaffold each of
their two islands *individually* as a 103-mer (the PepSeq max length), composite-score, and
keep the top 5 per island. This script builds the `designs.parquet` that the pipeline's
stage02/stage03 consume, one row per scaffolded island.

For each epitope we read its native dual-island contig from `dp2.parquet`. An *island* is one
contig `A`-segment (e.g. `A35-85`). For each island with span >= 2 we emit a single-island
construct:

    <N>/A<a>-<b>/<C>        with  N + (b-a+1) + C = 103

i.e. the full island A-segment is the fixed motif, flanked by a scaffold budget (103 - span)
split into an N- and C-flank, to a fixed 103-residue total. We emit `--variants` V such
contigs per island by randomly sampling V distinct N-flank lengths (seeded, reproducible),
which sweeps where the island sits along the 103-mer -- the single-island analog of Lawson's
inter-island-gap sweep (see manuscript sec:perisland). All V contigs of an island, and both
islands of an epitope, get distinct `contig_id`s under the same ledger `id`, so the
per-epitope antigen PDB (`<id>.pdb`, configs.paths ABDB_CLEANED_PDB_DIR) still resolves and
stage02 keeps them all.

Run size: 87 islands x V x 8 RFD3 x 8 MPNN = 5,568*V designs. Default V=20 -> 1,740 contigs,
111,360 designs (order of Lawson's 151k DP3 set). The RFD3 step is 1,740 x 8 seeds = 13,920
array tasks (submit in chunks under SLURM MaxArraySize; the driver does this).

Decisions baked in (flag if they change):
  - Total length 103 (PepSeq max). Lawson's native dual-island contigs were 104; the change
    to 103 is unexplained but expected immaterial (see manuscript sec:perisland).
  - Fixed atoms = the antibody-contact residues *within* each island (the island's slice of
    `epitope_resindices`), matching Lawson's RFD3 setup -- the whole A-segment is the motif,
    only contacts are held. Set `--fix-whole-island` to instead fix the entire A-segment span.
  - Size-1 islands are skipped (cannot be presented alone); the partner island is still built.
  - Flanks are fixed-length (e.g. `43-43`), not ranges; V=20 random splits per island (--seed
    for reproducibility). V=1 keeps a single centered split (the minimal 5,568-design run).

Seeds/reps (the array size) are NOT applied here -- stage02 expands them:
    python -m episcaf_pipeline init   --dataset results/dual_island_designs.parquet --run_dir runs/<run>
    python -m episcaf_pipeline stage02 --run_dir runs/<run> --seeds 0,1,2,3 --reps 1
    python -m episcaf_pipeline stage03 --run_dir runs/<run> \
        --pdb_dir /tgen_labs/altin/alphafold3/workspace/episcaf-experiments/data/abdb/complex_pdbfiles/cleaned
    sbatch --array=1-<rows> episcaf_pipeline/hpc/sbatch/rfd3_array.sbatch runs/<run>

Usage:
  python3 -m episcaf_pipeline.build_dual_island_designs \
      --dp2 ../known_antigen/analysis/full_run/dp2.parquet \
      --out results/dual_island_designs.parquet
"""
from __future__ import annotations
import argparse
import random
import re
from pathlib import Path

import pandas as pd

A_SEGMENT = re.compile(r"A(\d+)-(\d+)")
TOTAL_LEN = 103          # PepSeq max construct length
MIN_ISLAND_SIZE = 2      # a single residue cannot be presented on its own


def flank_splits(budget: int, variants: int, rng: random.Random) -> list[int]:
    """Return `variants` N-flank lengths for an island over a scaffold `budget` (C = budget-N).

    V=1 keeps the single centered split (so the minimal run is unchanged). For V>1 we sample
    N at random (distinct) over [1, budget-1] -- a flank as small as 1, matching Lawson's
    minimum -- which randomizes where the island sits along the construct. The rng is seeded
    by the caller, so the ledger is reproducible. Capped at budget-1 distinct placements.
    """
    if variants <= 1:
        return [budget // 2]
    pool = list(range(1, budget))                  # possible N values, C = budget - N
    return sorted(rng.sample(pool, min(variants, len(pool))))


def build_rows(df: pd.DataFrame, total_len: int, fix_whole_island: bool,
               variants: int, seed: int) -> list[dict]:
    rng = random.Random(seed)
    dual = df[df.epitope_chunks == 2].drop_duplicates("id")
    rows: list[dict] = []
    for _, r in dual.iterrows():
        segs = [(int(a), int(b)) for a, b in A_SEGMENT.findall(r["contig_string"])]
        if len(segs) != 2:
            raise ValueError(f"{r['id']}: expected 2 A-segments, got {segs} "
                             f"from {r['contig_string']!r}")
        contacts_all = sorted(int(x) for x in r["epitope_resindices"])
        cid = 0  # running contig id within this epitope (island x flank-variant)
        for isl_idx, (a, b) in enumerate(segs):
            span = b - a + 1
            if span < MIN_ISLAND_SIZE:
                continue  # skip size-1 island; partner still emitted
            budget = total_len - span
            if budget < 2:
                raise ValueError(f"{r['id']} island A{a}-{b}: span {span} leaves "
                                 f"no room in a {total_len}-mer")
            # 0-based original-antigen indices fixed as backbone for this island
            if fix_whole_island:
                fixed = list(range(a - 1, b))
            else:
                fixed = [x for x in contacts_all if a - 1 <= x <= b - 1]
            for v_idx, n_flank in enumerate(flank_splits(budget, variants, rng)):
                c_flank = budget - n_flank
                rows.append({
                    "id": r["id"],
                    "contig_id": cid,
                    "contig_string": f"{n_flank}-{n_flank}/A{a}-{b}/{c_flank}-{c_flank}",
                    "contig_length": f"{total_len}-{total_len}",
                    "epitope_resindices": fixed,
                    "epitope_chunks": 1,            # now a single-island target
                    "island_segment": f"A{a}-{b}",  # audit only (stripped by stage02 SAFE_COLS)
                    "island_index": isl_idx,
                    "island_size": span,
                    "variant": v_idx,
                    "n_flank": n_flank,
                    "c_flank": c_flank,
                    "n_fixed": len(fixed),
                })
                cid += 1
    return rows


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dp2", type=Path,
                    default=Path("../known_antigen/analysis/full_run/dp2.parquet"))
    ap.add_argument("--out", type=Path, default=Path("results/dual_island_designs.parquet"))
    ap.add_argument("--total-len", type=int, default=TOTAL_LEN)
    ap.add_argument("--variants", type=int, default=20,
                    help="random flank-split contigs per island (V). "
                         "V=1 -> 5,568 designs; V=20 -> ~111k; V=27 -> ~150k (Lawson scale)")
    ap.add_argument("--seed", type=int, default=0,
                    help="seed for the random flank splits (keeps the ledger reproducible)")
    ap.add_argument("--fix-whole-island", action="store_true",
                    help="fix the entire A-segment span instead of only its contact residues")
    args = ap.parse_args()

    df = pd.read_parquet(args.dp2)
    rows = build_rows(df, args.total_len, args.fix_whole_island, args.variants, args.seed)
    out = pd.DataFrame(rows)

    args.out.parent.mkdir(parents=True, exist_ok=True)
    out.to_parquet(args.out, index=False)
    # Full ledger as a tracked CSV (parquet is gitignored): epitope_resindices is serialized
    # as a python-literal list so the driver can rebuild the parquet on the cluster with no
    # dp2 dependency (pull-and-go). This CSV is the shippable source of the run.
    csv_out = args.out.with_suffix(".csv")
    out.to_csv(csv_out, index=False)

    n_epi = out["id"].nunique()
    n_islands = out.drop_duplicates(["id", "island_index"]).shape[0]
    n_contigs = len(out)
    print(f"wrote {args.out}  ({n_contigs} contig rows, {n_islands} islands, {n_epi} epitopes)")
    print(f"wrote {csv_out}  (tracked, full ledger)")
    print(f"  total length   : {args.total_len}   variants/island: {args.variants} (seed {args.seed})")
    print(f"  fixed atoms    : {'whole A-segment' if args.fix_whole_island else 'contacts within island'}")
    print(f"  contigs        : {n_islands} islands x {args.variants} variants = {n_contigs}")
    print(f"  designs total  : {n_contigs} contigs x 8 RFD3 x 8 MPNN = {n_contigs*64:,}")
    print(f"  RFD3 array     : {n_contigs} contigs x (seeds x reps) tasks, set at stage02")


if __name__ == "__main__":
    main()
