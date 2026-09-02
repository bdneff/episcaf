#!/usr/bin/env python3
"""Cross-reference our antibody:antigen structures against SKEMPI 2.0 (Decision 1 evidence).

Decision 1 of v3 is picking a pilot complex with experimentally measured binding energetics, so the
per-residue energy method can be validated against known hot spots before it drives any design.
SKEMPI 2.0 (Jankauskaite et al. 2019) is the standard database of measured binding-energy changes on
mutation. This script asks two things:
  (a) do any of OUR scaffolded structures appear in SKEMPI (so we could validate on our own set)?
  (b) which antibody:antigen (AB/AG) complexes in SKEMPI have full alanine scans, as pilot candidates?

Result on file (2026-09-02): NONE of our 59 appear in SKEMPI 2.0 -- our set is almost all recent
(2020s) PDBs and SKEMPI 2.0 is a 2019 database. So the pilot comes from the AB/AG fallback list; the
recommended pick is 3HFM (HyHEL-10 / hen egg lysozyme): a textbook hot-spot system with 18 alanine
measurements, and lysozyme is already the system bcell_epitope runs MD on. See docs/DECISIONS.md.

SKEMPI 2.0 is not vendored (1.5 MB; *.csv is gitignored). Download it first:
    curl -sSL -o skempi_v2.csv https://life.bsc.es/pid/skempi2/database/download/skempi_v2.csv

Run:
    python energetics/skempi_overlap.py --skempi /path/to/skempi_v2.csv
"""
from __future__ import annotations
import argparse
import csv
import math
from collections import defaultdict
from pathlib import Path

# repo root = episcaf_v2/ ; this file is episcaf_v2/episcaf_v3/energetics/skempi_overlap.py
REPO = Path(__file__).resolve().parents[2]
DEFAULT_FASTA = REPO / "data" / "sequences" / "dp3_mab_antigens.fasta"
R_KCAL = 1.987204e-3  # gas constant, kcal/mol/K


def our_pdb_codes(fasta: Path) -> set[str]:
    """The 4-char PDB codes of our scaffolded structures (record ids like '7ox3_0P')."""
    codes = set()
    for line in fasta.read_text().splitlines():
        if line.startswith(">"):
            codes.add(line[1:].split("_")[0].strip().lower())
    return codes


def ddg_bind(aff_mut: str, aff_wt: str, temp: str) -> float | None:
    """ΔΔG = RT ln(Kd_mut / Kd_wt); +ve = binding weakened (a hot spot, for an Ala mutation)."""
    try:
        m, w = float(aff_mut), float(aff_wt)
        T = float(temp) if temp else 298.0
        return R_KCAL * T * math.log(m / w)
    except (ValueError, ZeroDivisionError):
        return None


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--skempi", type=Path, required=True, help="path to skempi_v2.csv (see docstring)")
    ap.add_argument("--fasta", type=Path, default=DEFAULT_FASTA,
                    help="our antigen record fasta (default: repo data/sequences/dp3_mab_antigens.fasta)")
    ap.add_argument("--top", type=int, default=12, help="how many AB/AG alanine-scan candidates to list")
    args = ap.parse_args()

    ours = our_pdb_codes(args.fasta)
    print(f"our structures: {len(ours)} PDB codes from {args.fasta.name}")

    rows = list(csv.reader(args.skempi.open(), delimiter=";"))
    idx = {h: i for i, h in enumerate(rows[0])}
    def col(r, name): return r[idx[name]] if idx[name] < len(r) else ""

    our_hits = defaultdict(lambda: {"n": 0, "ala": 0, "type": set(), "p1": "", "p2": "", "ddg": []})
    abag_ala = defaultdict(lambda: {"ala": 0, "p1": "", "p2": "", "ddg": []})

    for r in rows[1:]:
        if not r:
            continue
        pdb = col(r, "#Pdb").split("_")[0].lower()
        mut = col(r, "Mutation(s)_cleaned")
        htype = col(r, "Hold_out_type")
        single_ala = mut and ("," not in mut) and mut[-1:] == "A"
        d = ddg_bind(col(r, "Affinity_mut_parsed"), col(r, "Affinity_wt_parsed"), col(r, "Temperature"))
        if pdb in ours:
            h = our_hits[pdb]
            h["n"] += 1; h["type"].add(htype); h["p1"] = col(r, "Protein 1"); h["p2"] = col(r, "Protein 2")
            if single_ala:
                h["ala"] += 1
                if d is not None:
                    h["ddg"].append(d)
        if htype == "AB/AG" and single_ala and d is not None:
            a = abag_ala[pdb]
            a["ala"] += 1; a["p1"] = col(r, "Protein 1"); a["p2"] = col(r, "Protein 2"); a["ddg"].append(d)

    print("\n=== our structures found in SKEMPI ===")
    if not our_hits:
        print("  NONE of our PDB codes appear in SKEMPI 2.0 (expected: our set is mostly 2020s PDBs;")
        print("  SKEMPI 2.0 is a 2019 database). Pilot comes from the AB/AG fallback list below.")
    else:
        for pdb, h in sorted(our_hits.items()):
            dd = h["ddg"]; rng = f"{min(dd):+.2f}..{max(dd):+.2f}" if dd else "n/a"
            print(f"  {pdb}: {h['n']} mutations ({h['ala']} Ala), type={sorted(h['type'])}, "
                  f"ddG(Ala) {rng} kcal/mol | {h['p1']} / {h['p2']}")

    print(f"\n=== AB/AG complexes with alanine scans (top {args.top} by #Ala) -- pilot candidates ===")
    for pdb, a in sorted(abag_ala.items(), key=lambda kv: -kv[1]["ala"])[:args.top]:
        dd = a["ddg"]; rng = f"{min(dd):+.2f}..{max(dd):+.2f}" if dd else "n/a"
        nhot = sum(1 for x in dd if x > 1.0)
        print(f"  {pdb}: {a['ala']} Ala, ddG {rng} kcal/mol, {nhot} hot spots(>1) | {a['p1']} / {a['p2']}")


if __name__ == "__main__":
    main()
