#!/usr/bin/env python3
"""Join DP4 PepSeq assay binding data to the design metrics already in dp4_library.csv.

STATUS: pre-staged and waiting for data. The DP4 library was ordered 2026-08-05 (~4-week
synthesis + assay turnaround), so the assay CSV does not exist yet. This script is the turnkey
bridge from "what binds" to "what the design looks like", written and self-checked now so that
the moment the data lands we run one command instead of writing first-draft code under pressure.

It is the DP4 analogue of `scripts/build_dp3_binding_join.py`, but MUCH simpler: DP3 had to
round-trip through dp2.parquet to fetch AF3 metrics and run a separate cylinder job. In DP4 every
metric AND the composite score are already columns in `data/libraries/dp4_library.csv`
(epitope_rmsd, overall_rmsd, mean_pae, af3_clashes, cylinder_clashes, composite, rank_in_group,
is_global_pass, island_index, category, design_ID, target). So the join is a single merge on the
`library_member` key -- no parquet, no cluster round-trip.

WHAT IS KNOWN NOW (verified against the shipped files):
  - dp4_library.csv: 36,000 rows, `library_member` unique + contiguous DP4_1..DP4_36000.
  - `target` carries the cognate id for every antibody arm as `<pdb>_<N>P` (e.g. 2qqn_0P),
    exactly like DP3's `Target`; 8VDL uses `8VDL_epitope`/`8VDL_hotspots`; the EPCR minibinders
    use `fold_pfemp1_epcr_*` (note the known cosmetic truncation `..._mod`, John 2026-08-05).

WHAT IS UNKNOWN UNTIL THE DATA ARRIVES (the ONE thing to confirm, then this runs):
  - the DP4 assay CSV's exact column layout: (a) does it keep a `library_member` column (then the
    join is the trivial route below); (b) the NoAb baseline column name; (c) the per-antibody
    intensity column names. DP3 used `NoAb...` for baseline and `X<pdb>_...` for each antibody;
    those regexes are the defaults here and are overridable with --noab-re / --ab-re, or you can
    supply an explicit --ab-map file (one `pdb,column_name` line per antibody) if DP4 names differ.

Run NOW (no assay data needed -- proves the library half of the join is sound):
    python scripts/build_dp4_binding_join.py --self-check

Run WHEN THE DATA LANDS:
    python scripts/build_dp4_binding_join.py --assay data/dp4_binding/<run1>.csv \
                                             --assay data/dp4_binding/<run2>.csv
    # -> results/dp4_binding_metrics.csv

See docs/DP4_RESULTS_ANALYSIS.md for the full analysis plan this join feeds.
"""
from __future__ import annotations
import argparse
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
LIBRARY = ROOT / "data" / "libraries" / "dp4_library.csv"
OUT = ROOT / "results" / "dp4_binding_metrics.csv"

# Design/scoring columns pulled from dp4_library.csv onto every assayed row. These are the
# columns the binding analysis regresses on -- the whole point of the join.
LIB_METRIC_COLS = [
    "category", "target", "design_ID", "island_index",
    "epitope_rmsd", "overall_rmsd", "epitope_pae", "scaffold_pae", "mean_pae",
    "ptm", "af3_clashes", "cylinder_clashes",
    "composite", "rank_in_group", "is_global_pass",
]

# Assay column conventions (DP3 defaults; override if DP4 names differ).
DEFAULT_AB_RE = r"^X([0-9a-z]{4})_"   # per-antibody intensity col; capture group = 4-char PDB id
DEFAULT_NOAB_RE = r"^NoAb"            # the NoAb baseline column

log1p10 = lambda x: np.log10(1.0 + x)  # John's transform (log10(1+intensity))


# --------------------------------------------------------------------------- library side
def load_library() -> pd.DataFrame:
    lib = pd.read_csv(LIBRARY, low_memory=False)
    missing = [c for c in ["library_member", "sequence", "designedSequence"] + LIB_METRIC_COLS
               if c not in lib.columns]
    if missing:
        sys.exit(f"ERROR: dp4_library.csv is missing expected columns: {missing}")
    return lib


def cognate_pdb_factory(known_abs: set[str]):
    """target '2qqn_0P' -> '2qqn' if that antibody was assayed, else None.

    Antibody arms store `<pdb>_<N>P`; 8VDL stores `8VDL_epitope/_hotspots`; minibinders store
    `fold_pfemp1_epcr_*`. Only ids present in the assay's antibody columns resolve to a cognate."""
    def cognate(target: str):
        pdb = str(target).split("_")[0].lower()
        return pdb if pdb in known_abs else None
    return cognate


# --------------------------------------------------------------------------- assay side
def detect_ab_columns(df: pd.DataFrame, ab_re: str) -> dict[str, str]:
    """{antibody pdb id -> intensity column} from column names matching ab_re."""
    pat = re.compile(ab_re)
    out: dict[str, str] = {}
    for c in df.columns:
        m = pat.match(c)
        if m:
            out[m.group(1).lower()] = c
    return out


def detect_noab_column(df: pd.DataFrame, noab_re: str) -> str:
    pat = re.compile(noab_re)
    cols = [c for c in df.columns if pat.match(c)]
    if len(cols) != 1:
        sys.exit(f"ERROR: expected exactly one NoAb column matching {noab_re!r}, got {cols}. "
                 f"Pass --noab-re to match the DP4 baseline column.")
    return cols[0]


def load_ab_map_file(path: Path) -> dict[str, str]:
    """Explicit `pdb,column_name` mapping, one per line -- an escape hatch if the DP4 antibody
    column names don't match a single regex."""
    out: dict[str, str] = {}
    for ln in path.read_text().splitlines():
        ln = ln.strip()
        if not ln or ln.startswith("#"):
            continue
        pdb, col = (x.strip() for x in ln.split(",", 1))
        out[pdb.lower()] = col
    return out


# --------------------------------------------------------------------------- self-check
def self_check() -> None:
    """Validate everything the join depends on on the library side, with no assay data.
    This is the 'verify against data, don't assume' gate we can run today."""
    lib = load_library()
    n = len(lib)
    print(f"[self-check] dp4_library.csv: {n:,} rows, {len(lib.columns)} columns")

    # library_member unique + contiguous DP4_1..DP4_N (the master join key)
    lm = lib["library_member"]
    assert lm.is_unique, "library_member is not unique"
    nums = sorted(int(x.split("_")[1]) for x in lm)
    assert nums == list(range(1, n + 1)), "library_member is not contiguous 1..N"
    print(f"[self-check] library_member OK: unique + contiguous DP4_1..DP4_{n:,}")

    # every metric/scoring column present and typed
    print(f"[self-check] all {len(LIB_METRIC_COLS)} metric/scoring columns present: "
          f"{', '.join(LIB_METRIC_COLS)}")

    # cognate-antibody resolvability: pretend the 8 DP3 antibodies were re-assayed and show how
    # many rows in each arm would resolve to a cognate. (Illustrative -- real ids come from data.)
    demo_abs = {"6o9i", "8cz8", "8jnk", "6xxv", "5fhx", "7ox3", "8db4", "8pww"}
    cog = cognate_pdb_factory(demo_abs)
    lib["_demo_cognate"] = lib["target"].map(cog)
    print("\n[self-check] per-category coverage (n rows | with a resolvable cognate, "
          "using the 8 DP3 antibodies as a stand-in target set):")
    g = (lib.groupby("category")
            .agg(n=("library_member", "size"),
                 demo_cognate=("_demo_cognate", lambda s: s.notna().sum()))
            .sort_values("n", ascending=False))
    print(g.to_string())
    print("\n[self-check] PASS -- the library side of the join is sound. When the DP4 assay CSV "
          "arrives, run without --self-check (see the module docstring).")


# --------------------------------------------------------------------------- main join
def build_join(assay_paths: list[Path], ab_re: str, noab_re: str,
               ab_map_file: Path | None, lib_key: str, seq_key: str) -> None:
    lib = load_library()

    frames = [pd.read_csv(p, low_memory=False) for p in assay_paths]
    for p, df in zip(assay_paths, frames):
        print(f"[join] {p.name}: {len(df):,} rows, {len(df.columns)} cols")

    # antibody + NoAb columns per run
    ab_map: dict[str, tuple[str, str]] = {}   # pdb -> (intensity col, its run's NoAb col)
    for df in frames:
        noab = detect_noab_column(df, noab_re)
        abcols = (load_ab_map_file(ab_map_file) if ab_map_file
                  else detect_ab_columns(df, ab_re))
        if not abcols:
            sys.exit(f"ERROR: no antibody columns detected (ab_re={ab_re!r}). Inspect the assay "
                     f"header and pass --ab-re or --ab-map.")
        for pdb, col in abcols.items():
            if col in df.columns:
                ab_map[pdb] = (col, noab)
    print(f"[join] antibodies detected: {sorted(ab_map)}")

    # stack all runs' assay rows, keyed by library_member (or sequence fallback)
    key = lib_key
    if all(lib_key in df.columns for df in frames):
        print(f"[join] joining on library key '{lib_key}'")
    elif all(seq_key in df.columns for df in frames):
        key = seq_key
        print(f"[join] '{lib_key}' absent from assay; falling back to sequence key '{seq_key}'")
    else:
        sys.exit(f"ERROR: assay CSV has neither '{lib_key}' nor '{seq_key}'. "
                 f"Columns: {list(frames[0].columns)[:20]}...")

    # keep each run's key + its assay columns; outer-merge runs on the key
    merged = None
    for df in frames:
        noab = detect_noab_column(df, noab_re)
        abcols = list((load_ab_map_file(ab_map_file) if ab_map_file
                       else detect_ab_columns(df, ab_re)).values())
        sub = df[[key, noab] + abcols].copy()
        merged = sub if merged is None else merged.merge(sub, on=key, how="outer")

    # attach the library metrics/composite
    lib_key_col = lib_key if key == lib_key else seq_key_to_lib_col(seq_key)
    out = merged.merge(
        lib[[lib_key_col] + [c for c in LIB_METRIC_COLS]].rename(columns={lib_key_col: key}),
        on=key, how="left")
    n_lib_hit = out["composite"].notna().sum()
    print(f"[join] {n_lib_hit:,}/{len(out):,} assay rows matched a library design")

    # cognate antibody per row -> log-enrichment
    cog = cognate_pdb_factory(set(ab_map))
    out["cognate_ab"] = out["target"].map(cog)
    sig, base = [], []
    for _, r in out.iterrows():
        ab = r["cognate_ab"]
        if ab is None or ab not in ab_map:
            sig.append(np.nan); base.append(np.nan); continue
        ab_col, noab_col = ab_map[ab]
        sig.append(r.get(ab_col, np.nan)); base.append(r.get(noab_col, np.nan))
    out["cognate_ab_signal"] = sig
    out["cognate_noab"] = base
    out["cognate_log_ab"] = log1p10(out["cognate_ab_signal"])
    out["cognate_log_noab"] = log1p10(out["cognate_noab"])
    out["cognate_log_enrichment"] = out["cognate_log_ab"] - out["cognate_log_noab"]

    OUT.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(OUT, index=False)
    print(f"\n[join] wrote {OUT}  ({len(out):,} rows, {len(out.columns)} cols)")

    # quick orientation
    s = out[out["cognate_ab"].notna()]
    if len(s):
        g = (s.groupby(["category", "cognate_ab"])
               .agg(n=("cognate_log_enrichment", "size"),
                    med_enrich=("cognate_log_enrichment", "median"),
                    med_composite=("composite", "median"))
               .sort_values("n", ascending=False))
        print(f"\n[join] cognate rows by arm x antibody (top 15):\n{g.head(15).to_string()}")


def seq_key_to_lib_col(seq_key: str) -> str:
    # assay sequence column may be the bare scaffold (designedSequence) or the 103-mer (sequence)
    return "designedSequence" if seq_key == "designedSequence" else "sequence"


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--self-check", action="store_true",
                    help="validate the library side with no assay data (run this today)")
    ap.add_argument("--assay", action="append", default=[], type=Path,
                    help="an assay CSV (repeat for multiple runs, as DP3 had IM226 + IM229)")
    ap.add_argument("--ab-re", default=DEFAULT_AB_RE,
                    help=f"regex for per-antibody intensity columns (default {DEFAULT_AB_RE!r})")
    ap.add_argument("--noab-re", default=DEFAULT_NOAB_RE,
                    help=f"regex for the NoAb baseline column (default {DEFAULT_NOAB_RE!r})")
    ap.add_argument("--ab-map", type=Path, default=None,
                    help="explicit `pdb,column` map file (escape hatch if names don't match a regex)")
    ap.add_argument("--lib-key", default="library_member",
                    help="assay column that carries the library key (default library_member)")
    ap.add_argument("--seq-key", default="designedSequence",
                    help="fallback assay column with the peptide sequence (default designedSequence)")
    args = ap.parse_args()

    if args.self_check:
        self_check()
        return
    if not args.assay:
        sys.exit("No --assay given. Run `--self-check` today, or pass the DP4 assay CSV(s) once "
                 "the data lands. See the module docstring.")
    for p in args.assay:
        if not p.exists():
            sys.exit(f"ERROR: assay file not found: {p}")
    build_join(args.assay, args.ab_re, args.noab_re, args.ab_map, args.lib_key, args.seq_key)


if __name__ == "__main__":
    main()
