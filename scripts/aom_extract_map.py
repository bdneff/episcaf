#!/usr/bin/env python3
"""Extract the antibody-occupancy map from John Altin's committed AOMclash model into a slim artifact.

John's `aomclash/model/aomclash_model.pkl` (AltinLab/epitope-universe @ upstream/feat/clash-method)
bundles a trained 109-feature gradient-boosted model AND the occupancy grids of every SAbDab complex
(bit-packed in `fields_packed`). We ship only the SINGLE carved-occupancy feature `s_p1_free`
(see results/aomclash/feature_set_decision.md), so we keep the MAP and drop his model: this writes
`episcaf_analysis/aomclash/aom_map.pkl` with just the packed grids + the per-complex keys needed for
leave-one-antigen-out. The map is antibody-geometry only (AF-independent), so it is the reproducible
scorer artifact our port needs; we can also rebuild it from SAbDab and cross-check.

  python scripts/aom_extract_map.py --model <aomclash_model.pkl> --out episcaf_analysis/aomclash/aom_map.pkl
"""
from __future__ import annotations

import argparse
import pickle
import warnings


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--model", required=True, help="John's aomclash_model.pkl")
    ap.add_argument("--out", default="episcaf_analysis/aomclash/aom_map.pkl")
    a = ap.parse_args(argv)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")            # sklearn version mismatch on his GBM (we drop it)
        art = pickle.load(open(a.model, "rb"))

    slim = {
        "fields_packed": art["fields_packed"],     # (n_complexes, ceil(field_length/8)) uint8
        "field_length": int(art["field_length"]),  # 35937 = 33^3
        "field_side": int(art["field_side"]),       # 33
        "keys": art["keys"],                         # biological_antigen_key per complex (LOO exclusion)
        "pdbs": art["pdbs"],                         # pdb id per complex
        "suggested_cutoffs": art.get("suggested_cutoffs"),   # {'aom':7.0,...} reference thresholds
        "source": "AltinLab/epitope-universe @ upstream/feat/clash-method aomclash/model/aomclash_model.pkl",
        "training": {k: art.get("training", {}).get(k) for k in ("date", "complexes")}
                    if isinstance(art.get("training"), dict) else None,
    }
    with open(a.out, "wb") as fh:
        pickle.dump(slim, fh, protocol=4)
    n = slim["fields_packed"].shape[0]
    print(f"wrote {a.out}: {n} complexes, {slim['field_side']}^3 grid, cutoffs={slim['suggested_cutoffs']}")


if __name__ == "__main__":
    main()
