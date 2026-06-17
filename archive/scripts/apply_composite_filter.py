#!/usr/bin/env python3
"""
apply_composite_filter.py

Turn per-design AF3 metrics into a shortlist:
  1. GATE   - hard sanity cuts on the reliable metrics (drop designs that aren't
              even presenting the epitope). The cylinder is NOT a gate.
  2. SCORE  - percentile-normalize each term over the gate-passing pool
              (robust to the heavy tails these metrics have), orient so higher
              is better, weighted sum -> composite in [0, 1].
  3. SELECT - within each epitope, keep the best prediction per design, then the
              top-K designs by composite. Keeps min(K, n_passing); thin epitopes
              are flagged low-yield rather than backfilled with failures.

Scoring is GLOBAL (percentiles over the whole gate-passing pool) but selection
is PER-EPITOPE, so every target gets representation while the absolute score
still says which epitopes produced genuinely good candidates.

--validate mode (known-antibody set only): defines a ground-truth "pass" using
the REAL antibody clash filter (af3_n_clash_res == 0) plus the gate metrics, then
measures how well the cylinder-based composite recovers those passers — both
WITH and WITHOUT the cylinder term, so you can see if the proxy is pulling its
weight.

Example
-------
    # production shortlist
    python scripts/apply_composite_filter.py \
        --metrics_csv runs/run_rfd3_mpnn/04_filter/metrics_cylinder_full.csv \
        --out_csv     runs/run_rfd3_mpnn/04_filter/shortlist.csv

    # validate on the known-antibody set (ground truth present)
    python scripts/apply_composite_filter.py \
        --metrics_csv runs/run_rfd3_mpnn/04_filter/metrics_cylinder_full.csv \
        --out_csv     runs/run_rfd3_mpnn/04_filter/shortlist.csv \
        --validate
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

# ---------------------------------------------------------------------------
# Tunables — edit here.
# ---------------------------------------------------------------------------
# weight + direction ("low" = smaller is better, "high" = larger is better).
# Column aliases are resolved at load time (first present wins).
WEIGHTS = {
    "epitope_rmsd": dict(weight=0.35, direction="low",
                         aliases=["epitope_chunk_rmsd"]),
    "pae":          dict(weight=0.25, direction="low",
                         aliases=["epitope_pae", "mean_pae"]),
    "overall_rmsd": dict(weight=0.15, direction="low",
                         aliases=["overall_rmsd"]),
    "ptm":          dict(weight=0.10, direction="high",
                         aliases=["ptm"]),
    "cylinder":     dict(weight=0.15, direction="low",
                         aliases=["cylinder_ca_clashes"]),
}

# hard gates: metric -> (op, threshold). Tunable production cut. Its only job is
# to drop designs that aren't presenting the epitope at all; scoring does the
# rest. Kept looser than PASS_CRITERIA on every axis so it can't silently remove
# a genuine passer before it can be ranked.
GATES = {
    "epitope_rmsd": ("<", 2.5),
}

# FIXED definition of a genuinely good design, used ONLY by --validate as ground
# truth. This is the established four-filter success criterion (overall RMSD <= 2,
# epitope RMSD <= 1, mean PAE < 5) plus the REAL antibody clash (af3_n_clash_res
# == 0, added in validate()). Independent of GATES so gate/weight tuning stays
# comparable. No pTM term — it isn't part of the official criterion.
PASS_CRITERIA = {
    "epitope_rmsd": ("<=", 1.0),
    "overall_rmsd": ("<=", 2.0),
    "pae":          ("<", 5.0),
}

GROUP_COL_ALIASES = ["id"]                                   # the "epitope"
DEDUPE_COL_ALIASES = ["token", "assay_scaffolded_epitope_id"]  # one design
TRUE_CLASH_COUNT_ALIASES = ["af3_n_clash_res"]
# pass definition must match the OFFICIAL filter (mean PAE), even though the
# composite score uses epitope_pae as its (better) ranking signal.
PASS_PAE_ALIASES = ["mean_pae", "epitope_pae"]
STATUS_OK = {"status": "ok", "af3_clash_status": "ok"}

OPS = {"<": np.less, "<=": np.less_equal, ">": np.greater, ">=": np.greater_equal}


# ---------------------------------------------------------------------------
def first_present(cols, candidates):
    for c in candidates:
        if c in cols:
            return c
    return None


def resolve_columns(df):
    """Map logical metric names -> actual column names present in df."""
    resolved = {}
    for name, spec in WEIGHTS.items():
        col = first_present(df.columns, spec["aliases"])
        if col is not None:
            resolved[name] = col
    return resolved


def apply_status_filter(df):
    for col, ok in STATUS_OK.items():
        if col in df.columns:
            before = len(df)
            df = df[df[col] == ok].copy()
            print(f"  status filter {col}=={ok!r}: {len(df):,}/{before:,}")
    return df


def apply_gates(df, colmap):
    mask = pd.Series(True, index=df.index)
    for name, (op, thr) in GATES.items():
        col = colmap.get(name)
        if col is None:
            print(f"  [gate] skip {name}: no column present")
            continue
        vals = pd.to_numeric(df[col], errors="coerce")
        mask &= OPS[op](vals, thr) & vals.notna()
    return df[mask].copy()


def percentile_score(df, colmap, weights):
    """Per-term percentile in [0,1], oriented so higher=better, weighted sum."""
    total_w = sum(weights[n] for n in weights)
    composite = pd.Series(0.0, index=df.index)
    term_cols = {}
    for name, w in weights.items():
        col = colmap[name]
        vals = pd.to_numeric(df[col], errors="coerce")
        pct = vals.rank(pct=True)                       # 1.0 = largest raw value
        if WEIGHTS[name]["direction"] == "low":
            pct = 1.0 - pct
        term = pct * (w / total_w)
        composite = composite + term
        term_cols[f"score_{name}"] = term
    out = df.copy()
    for k, v in term_cols.items():
        out[k] = v
    out["composite"] = composite
    return out


def select_topk(df, group_col, dedupe_col, k):
    """Best prediction per design, then top-k designs per epitope."""
    d = df
    if dedupe_col and dedupe_col in d.columns:
        d = (d.sort_values("composite", ascending=False)
               .drop_duplicates(dedupe_col, keep="first"))
    d = d.sort_values(["__grp__", "composite"], ascending=[True, False])
    d["rank_in_epitope"] = d.groupby("__grp__").cumcount() + 1
    return d[d["rank_in_epitope"] <= k].copy()


def active_weights(colmap):
    """Weights for terms whose column is actually present; renormalized to 1."""
    w = {n: WEIGHTS[n]["weight"] for n in WEIGHTS if n in colmap}
    s = sum(w.values())
    return {n: v / s for n, v in w.items()}, s


# ---------------------------------------------------------------------------
def run_selection(df, colmap, weights, group_col, dedupe_col, k):
    df = df.copy()
    df["__grp__"] = df[group_col].astype(str)
    scored = percentile_score(df, colmap, weights)
    selected = select_topk(scored, group_col, dedupe_col, k)
    return scored, selected


def compute_is_pass(frame, colmap):
    """Ground-truth pass = official four-filter criterion + real clash == 0.
    Uses mean_pae for the PAE criterion to match the official filter, even
    though the composite score ranks on epitope_pae."""
    gt = pd.Series(True, index=frame.index)
    for name, (op, thr) in PASS_CRITERIA.items():
        col = (first_present(frame.columns, PASS_PAE_ALIASES) if name == "pae"
               else colmap.get(name))
        if col is None:
            continue
        vals = pd.to_numeric(frame[col], errors="coerce")
        gt &= OPS[op](vals, thr) & vals.notna()
    clash_col = first_present(frame.columns, TRUE_CLASH_COUNT_ALIASES)
    if clash_col is not None:
        gt &= (pd.to_numeric(frame[clash_col], errors="coerce") == 0)
    else:
        gt &= False
    return gt


def validate(full_df, selected, colmap, group_col, top_k):
    """Recall/precision of the composite shortlist vs ground-truth passers.
    Ground truth is computed over the FULL status-ok pool (not the gated set),
    so a true passer dropped by the gate correctly counts as a recall miss."""
    if first_present(full_df.columns, TRUE_CLASH_COUNT_ALIASES) is None:
        print("  [validate] no true-clash column; cannot validate.")
        return

    truepass = full_df[compute_is_pass(full_df, colmap)]
    sel_idx = set(selected.index)
    tp = sum(i in sel_idx for i in truepass.index)
    n_truepass = len(truepass)
    n_selected = len(selected)
    recall = tp / n_truepass if n_truepass else float("nan")
    precision = tp / n_selected if n_selected else float("nan")

    # epitope-level coverage: of epitopes with any true passer, how many keep one?
    tp_groups = set(truepass[group_col].astype(str))
    sel_tp_groups = set(selected.loc[selected.index.isin(truepass.index),
                                     group_col].astype(str))
    cov = len(sel_tp_groups) / len(tp_groups) if tp_groups else float("nan")

    # row-level recall is CAPPED by top-k when passers concentrate in few
    # epitopes: you can keep at most top_k passers per epitope. Report the ceiling
    # so a low recall isn't mistaken for a bad ranker.
    cap = sum(min(len(grp), top_k) for _, grp in truepass.groupby(group_col))
    ceiling = cap / n_truepass if n_truepass else float("nan")
    # per-epitope recall: mean over passer-epitopes of (its passers kept / min(n,k))
    per_ep = []
    sel_by_grp = selected.groupby(selected[group_col].astype(str))
    for g, grp in truepass.groupby(truepass[group_col].astype(str)):
        kept = sum(i in sel_idx for i in grp.index)
        per_ep.append(kept / min(len(grp), top_k))
    per_ep_recall = float(np.mean(per_ep)) if per_ep else float("nan")

    print(f"  ground-truth passers (official filter + real clash==0): {n_truepass:,}")
    print(f"  shortlist size:                              {n_selected:,}")
    print(f"  row-level recall   (passers kept):  {recall:.3f}  "
          f"(ceiling given top_k: {ceiling:.3f})")
    print(f"  per-epitope recall (capacity-adjusted): {per_ep_recall:.3f}")
    print(f"  row-level precision(shortlist good): {precision:.3f}")
    print(f"  epitope coverage   (targets with a passer kept): {cov:.3f} "
          f"({len(sel_tp_groups)}/{len(tp_groups)})")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--metrics_csv", required=True)
    ap.add_argument("--out_csv", required=True)
    ap.add_argument("--scored_csv", default=None,
                    help="optional: write ALL scored designs (composite + term "
                         "scores + is_pass) for plotting/diagnostics")
    ap.add_argument("--top_k", type=int, default=15)
    ap.add_argument("--validate", action="store_true")
    args = ap.parse_args()

    print("Loading metrics ...")
    df = pd.read_csv(args.metrics_csv, low_memory=False)
    print(f"  {len(df):,} rows")
    df = apply_status_filter(df)

    group_col = first_present(df.columns, GROUP_COL_ALIASES)
    dedupe_col = first_present(df.columns, DEDUPE_COL_ALIASES)
    if group_col is None:
        raise SystemExit(f"no epitope/group column among {GROUP_COL_ALIASES}")
    print(f"  group (epitope) = {group_col}   dedupe (design) = {dedupe_col}")

    colmap = resolve_columns(df)
    weights, raw_sum = active_weights(colmap)
    print(f"  active score terms: {colmap}")
    if "cylinder" not in colmap:
        print("  [warn] no cylinder_ca_clashes column — scoring WITHOUT the "
              "cylinder term (weights renormalized). Add it for the real run.")

    gated = apply_gates(df, colmap)
    print(f"  gate-passing designs: {len(gated):,}/{len(df):,}")

    # drop rows missing any active score term
    need = [colmap[n] for n in weights]
    gated = gated.dropna(subset=need)

    scored, selected = run_selection(gated, colmap, weights,
                                     group_col, dedupe_col, args.top_k)

    out = Path(args.out_csv)
    out.parent.mkdir(parents=True, exist_ok=True)
    # The shortlist carries EVERY column from the metrics table (all raw metrics
    # incl. af3_n_clash_res and cylinder_ca_clashes) plus the per-term percentile
    # scores, the composite, and rank_in_epitope -- so you can select by composite
    # and still inspect any individual metric per design.
    sel_out = selected.drop(columns=["__grp__"]).copy()
    # Attach the ground-truth pass flag when the real-clash column is present, so
    # the shortlist itself shows which picks meet the official four-filter
    # criterion (incl. real antibody clash). On no-antibody production sets that
    # column is absent; inspect cylinder_ca_clashes (also in the CSV) instead.
    clash_col = first_present(selected.columns, TRUE_CLASH_COUNT_ALIASES)
    if clash_col is not None:
        sel_out["is_pass"] = compute_is_pass(selected, colmap).astype(int)
    sel_out.to_csv(out, index=False)

    n_groups = selected["__grp__"].nunique()
    print(f"\nShortlist: {len(selected):,} designs across {n_groups:,} epitopes "
          f"(top {args.top_k} each) -> {out}")
    thin = (selected.groupby("__grp__").size() < args.top_k).sum()
    print(f"  low-yield epitopes (< {args.top_k} survivors): {thin:,}")

    # quick clash readout on what we actually selected
    cyl_col = colmap.get("cylinder")
    if cyl_col is not None:
        cz = pd.to_numeric(selected[cyl_col], errors="coerce")
        print(f"  cylinder clash on shortlist: {int((cz == 0).sum()):,}/{len(selected):,} "
              f"clash-free  (median {cz.median():.0f}, max {cz.max():.0f})")
    if clash_col is not None:
        rc = pd.to_numeric(selected[clash_col], errors="coerce")
        print(f"  real clash on shortlist:     {int((rc == 0).sum()):,}/{len(selected):,} "
              f"clash-free  (is_pass column written)")

    if args.scored_csv:
        full = scored.copy()
        full["is_pass"] = compute_is_pass(full, colmap).astype(int)
        sp = Path(args.scored_csv)
        sp.parent.mkdir(parents=True, exist_ok=True)
        full.drop(columns=["__grp__"]).to_csv(sp, index=False)
        print(f"  scored table ({len(full):,} designs, {int(full['is_pass'].sum())} "
              f"passers) -> {sp}")

    if args.validate:
        print("\n=== VALIDATION (known-antibody ground truth) ===")
        print("[with cylinder term]" if "cylinder" in colmap
              else "[cylinder term ABSENT — this is the no-cylinder baseline]")
        validate(df, selected, colmap, group_col, args.top_k)

        if "cylinder" in colmap:
            print("\n[ablation: WITHOUT cylinder term]")
            w_no, _ = active_weights({k: v for k, v in colmap.items()
                                      if k != "cylinder"})
            cmap_no = {k: v for k, v in colmap.items() if k != "cylinder"}
            scored_no, selected_no = run_selection(gated, cmap_no, w_no,
                                                   group_col, dedupe_col, args.top_k)
            validate(df, selected_no, cmap_no, group_col, args.top_k)


if __name__ == "__main__":
    main()
