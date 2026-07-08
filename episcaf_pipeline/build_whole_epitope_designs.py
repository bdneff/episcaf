#!/usr/bin/env python3
"""
build_whole_epitope_designs.py  --  emit the C1 (whole-epitope) RFD3 design ledger at 103.

C1 presents each mAb epitope *whole* -- all its islands held together in native geometry, the
way the antibody sees it -- as the DP3 comparator. Our first C1 pool *reproduced* Lawson's DP3
run and so inherited his **104**-residue contigs; the assay length is **103**, and the DP3
experience (truncating those 104mers to 103 gave generally weaker signal) means we do not want
to trim after the fact. This script rebuilds C1 natively at 103 by taking Lawson's exact
whole-epitope contigs and dropping a *single scaffold residue*, so the epitope, its islands, and
Lawson's inter-island-spacing sweep are all preserved -- only one flank (or, for an epitope whose
islands are flush at both termini, one inter-island spacer residue) is one shorter.

Why edit Lawson's contigs rather than regenerate them: C1's whole point is to be the DP3
comparator, so we keep his island placements and spacing sweep intact and change only the length.
This is the whole-epitope analog of `build_dual_island_designs.py` (C2), which instead *generated*
fresh 103-mer single-island contigs.

For each epitope we read every unique whole-epitope contig from `dp2.parquet`
(e.g. `15-15/A1-16/57-57/A81-81/15-15`, always summing to 104), shorten it to 103, and fix the
epitope **contact** residues (`epitope_resindices`, antigen-frame) as the motif backbone -- exactly
Lawson's setup, and exactly what `build_dual_island_designs.py` fixes. The A-segments are antigen
residue numbers and are untouched by the length edit, so the fixed atoms carry over unchanged.

The 104->103 length edit (`shorten_to_103`): drop one residue from the **larger terminal scaffold
flank** (keeps both termini >=1 when possible); if both terminal flanks are already 0 -- the
`3ux9_1P` island-flush-at-both-ends case -- drop it from the largest interior scaffold spacer
instead, so no epitope residue is ever lost. Native 103 generation therefore dissolves the
"cannot trim 3ux9" edge case that the post-hoc 104->103 assembly trim hit.

Emits the same minimal ledger schema stage01/stage02 consume (`id, contig_id, contig_string,
contig_length, epitope_resindices, epitope_chunks`), so the run is driven exactly like the
dual-island one:

    python -m episcaf_pipeline init   --dataset results/whole_epitope_designs.parquet --run_dir runs/<run>
    python -m episcaf_pipeline stage01 --run_dir runs/<run> --seeds 0 --reps 1
    python -m episcaf_pipeline stage02 --run_dir runs/<run> --pdb_dir <cleaned_antigen_pdbs>
    sbatch --array=1-<rows> episcaf_pipeline/hpc/sbatch/rfd3_array.sbatch runs/<run>

(or just `bash scripts/run_whole_epitope_rfd3.sh`). RFD3 emits 8 backbones/contig (n_batches=1),
so 1 task per contig -> ~2.2k RFD3 array tasks (56 mAbs x ~41 contigs), 64 designs/contig.

Usage:
  python3 -m episcaf_pipeline.build_whole_epitope_designs \
      --dp2 ../known_antigen/analysis/full_run/dp2.parquet \
      --drop-targets 2h32,4xwo,7a3t \
      --out results/whole_epitope_designs.parquet
"""
from __future__ import annotations
import argparse
import re
from pathlib import Path

import pandas as pd

TOTAL_LEN = 103          # PepSeq max construct length (was 104 in Lawson's DP3 contigs)
_TOK = re.compile(r"^(A?)(\d+)-(\d+)$")


def _parse(contig: str) -> list[tuple[bool, int, int]]:
    """Contig string -> [(is_island, a, b)] tokens. Scaffold flanks/spacers are 'N-N' (a==b);
    islands are 'Aa-b'."""
    toks = []
    for t in contig.split("/"):
        m = _TOK.match(t)
        if not m:
            raise ValueError(f"unparseable contig token {t!r} in {contig!r}")
        toks.append((m.group(1) == "A", int(m.group(2)), int(m.group(3))))
    return toks


def _contig_len(contig: str) -> int:
    return sum((b - a + 1) if isl else a for isl, a, b in _parse(contig))


def _emit(toks: list[tuple[bool, int, int]]) -> str:
    out = []
    for isl, a, b in toks:
        out.append(f"A{a}-{b}" if isl else f"{a}-{a}")
    return "/".join(out)


def shorten_to_103(contig: str, total_len: int = TOTAL_LEN) -> tuple[str, str]:
    """Drop one *scaffold* residue from a (total_len+1)-mer contig, preserving every epitope
    (A-segment) residue and the inter-island spacing. Returns (new_contig, where).

    Prefer the larger terminal scaffold flank (ties -> C-terminal); only if both terminal flanks
    are already 0 (island flush at both ends, e.g. 3ux9) shorten the largest interior scaffold
    spacer. Idempotent-safe: raises if the contig is already <= total_len or has no scaffold to cut.
    """
    toks = _parse(contig)
    cur = sum((b - a + 1) if isl else a for isl, a, b in toks)
    if cur <= total_len:
        raise ValueError(f"contig {contig!r} already length {cur} (<= {total_len})")
    if cur != total_len + 1:
        raise ValueError(f"contig {contig!r} length {cur}, expected {total_len + 1}")

    n_idx, c_idx = 0, len(toks) - 1
    n_is_flank = not toks[n_idx][0]
    c_is_flank = not toks[c_idx][0]
    n_len = toks[n_idx][1] if n_is_flank else -1
    c_len = toks[c_idx][1] if c_is_flank else -1

    # terminal flanks first, larger one (ties -> C-terminal), but only if it has a residue to give
    cand = []
    if c_is_flank and c_len >= 1:
        cand.append((c_len, c_idx, "C"))
    if n_is_flank and n_len >= 1:
        cand.append((n_len, n_idx, "N"))
    if cand:
        cand.sort(key=lambda x: (-x[0], x[2] == "N"))  # larger first; on tie prefer C (False<True)
        _, idx, where = cand[0]
    else:
        # both terminal flanks are 0 (flush) -> shorten the largest interior scaffold spacer
        interior = [(toks[i][1], i) for i in range(1, len(toks) - 1) if not toks[i][0] and toks[i][1] >= 1]
        if not interior:
            raise ValueError(f"no scaffold residue to drop in {contig!r}")
        interior.sort(key=lambda x: -x[0])
        idx, where = interior[0][1], f"spacer@{interior[0][1]}"

    isl, a, b = toks[idx]
    toks[idx] = (False, a - 1, a - 1)
    new = _emit(toks)
    assert _contig_len(new) == total_len, (contig, new, _contig_len(new))
    return new, where


def _drop(id_: str, drop: set[str]) -> bool:
    stem = str(id_).lower().split("_")[0]
    return stem in drop or any(str(id_).lower().startswith(d) for d in drop)


def build_rows(df: pd.DataFrame, total_len: int, drop: set[str]) -> list[dict]:
    rows: list[dict] = []
    u = df.drop_duplicates(["id", "contig_string"]).sort_values(["id", "contig_string"])
    for id_, grp in u.groupby("id", sort=True):
        if _drop(id_, drop):
            continue
        # epitope contact residues (antigen-frame) are per-epitope, identical across its contigs
        contacts = sorted(int(x) for x in grp.iloc[0]["epitope_resindices"])
        chunks = int(grp.iloc[0]["epitope_chunks"])
        for cid, (_, r) in enumerate(grp.iterrows()):
            new_contig, where = shorten_to_103(r["contig_string"], total_len)
            rows.append({
                "id": id_,
                "contig_id": cid,
                "contig_string": new_contig,
                "contig_length": f"{total_len}-{total_len}",
                "epitope_resindices": contacts,
                "epitope_chunks": chunks,
                "orig_contig_string": r["contig_string"],  # audit (stripped by stage02 SAFE_COLS)
                "length_edit": where,                       # audit
                "n_fixed": len(contacts),                   # audit
            })
    return rows


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dp2", type=Path,
                    default=Path("../known_antigen/analysis/full_run/dp2.parquet"))
    ap.add_argument("--out", type=Path, default=Path("results/whole_epitope_designs.parquet"))
    ap.add_argument("--total-len", type=int, default=TOTAL_LEN)
    ap.add_argument("--drop-targets", default="2h32,4xwo,7a3t",
                    help="comma-separated epitope stems to exclude (the 56-mAb set). "
                         "'' keeps all 59.")
    args = ap.parse_args()

    drop = {d.strip().lower() for d in args.drop_targets.split(",") if d.strip()}
    df = pd.read_parquet(args.dp2)
    rows = build_rows(df, args.total_len, drop)
    out = pd.DataFrame(rows)

    args.out.parent.mkdir(parents=True, exist_ok=True)
    out.to_parquet(args.out, index=False)
    # Tracked CSV ledger (parquet is gitignored): epitope_resindices serialized as a python-literal
    # list so the cluster driver rebuilds the parquet with NO dp2 dependency (pull-and-go), exactly
    # like the dual-island ledger.
    csv_out = args.out.with_suffix(".csv")
    out.to_csv(csv_out, index=False)

    n_epi = out["id"].nunique()
    n_contigs = len(out)
    lens = out["contig_string"].map(_contig_len)
    where = out["length_edit"].str.replace(r"@\d+", "", regex=True).value_counts().to_dict()
    print(f"wrote {args.out}  ({n_contigs} contigs, {n_epi} epitopes)")
    print(f"wrote {csv_out}  (tracked, full ledger)")
    print(f"  dropped stems  : {sorted(drop)}")
    print(f"  contig lengths : {sorted(lens.unique())}  (all == {args.total_len}: {bool((lens == args.total_len).all())})")
    print(f"  length edit    : {where}")
    print(f"  designs total  : {n_contigs} contigs x 8 RFD3 x 8 MPNN = {n_contigs * 64:,}")
    print(f"  RFD3 array     : {n_contigs} tasks (1/contig, 8 backbones each)")


if __name__ == "__main__":
    main()
