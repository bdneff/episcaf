#!/usr/bin/env python3
"""Unit tests for the interpretable binding model: Tier-1 feature extractors + the fit/predict path.

Pure-synthetic tests for the extractor math (no data), plus an integration check on one local native
CIF (epi_size reproduction) and the fit/predict/save/load round-trip on the real DP4 table.

Run:  <env>/bin/python tests/test_binding_model.py
"""
import sys
import tempfile
import warnings
from pathlib import Path

import numpy as np

warnings.filterwarnings("ignore")
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from episcaf_analysis.binding_model import native_features as NF   # noqa: E402
from episcaf_analysis.binding_model import data as D               # noqa: E402
from episcaf_analysis.binding_model.model import BindingModel      # noqa: E402

PASS = []


def check(name, cond):
    PASS.append(cond)
    print(f"  [{'OK' if cond else 'FAIL'}] {name}")


# ----------------------------- extractor math (synthetic) -----------------------------
def test_contact_number():
    # 3 colinear CA at 5 A spacing; radius 10 catches immediate neighbours only
    ca = np.array([[0, 0, 0], [5, 0, 0], [10, 0, 0], [15, 0, 0]], float)
    mask = np.array([True, True, True, True])
    # cn: endpoints have 1 nbr (<10, excl self: 5 and 10? for atom0: 5->yes,10->no at 10.0 strict),
    # atom0 nbrs within <10: atom1(5). =1 ; atom1 nbrs: atom0(5),atom2(5)=2 ; mean over all
    cn = NF.contact_number(ca, mask, radius=10.0)
    # expected per-atom: [1,2,2,1] -> mean 1.5
    check("contact_number colinear mean", abs(cn - 1.5) < 1e-9)
    # a protruding (isolated) residue lowers CN: add a far atom
    ca2 = np.vstack([ca, [100, 0, 0]])
    mask2 = np.array([False, False, False, False, True])
    check("contact_number isolated==0", NF.contact_number(ca2, mask2, 10.0) == 0.0)


def test_charge_features():
    # K(+1) R(+1) D(-1) E(-1) A(0) -> net 0, abs 4
    names = ["LYS", "ARG", "ASP", "GLU", "ALA"]
    ca = np.array([[0, 0, 0], [10, 0, 0], [20, 0, 0], [30, 0, 0], [40, 0, 0]], float)
    net, absq, coulomb = NF.charge_features(names, ca)
    check("charge net==0", abs(net) < 1e-9)
    check("charge abs==4", abs(absq - 4.0) < 1e-9)
    # coulomb: all pairs among 4 charged; sign mix -> finite, and two like/two unlike present
    check("charge coulomb finite", np.isfinite(coulomb))
    # all-positive patch -> net positive, coulomb strictly positive (repulsive like charges)
    net2, abs2, c2 = NF.charge_features(["LYS", "ARG"], np.array([[0, 0, 0], [5, 0, 0]], float))
    check("charge all-pos net==2", abs(net2 - 2.0) < 1e-9)
    check("charge all-pos coulomb>0", c2 > 0)


def test_glyco_distance():
    # antigen chain with one sequon NAT at residues 10-12 (N at 10), epitope residue 10 sits on it
    def mk(resid, resname, ca_xyz):
        return (resid, resname, {"CA": np.array(ca_xyz, float)})
    chain = [mk(9, "ALA", [0, 0, 0]), mk(10, "ASN", [5, 0, 0]),
             mk(11, "ALA", [10, 0, 0]), mk(12, "THR", [15, 0, 0]),
             mk(13, "GLY", [20, 0, 0])]
    chains = {"A": chain}
    # epitope = residue 13 (dist 15 A from the sequon N at resid10 [5,0,0]; under the 30 A cap)
    d_far = NF.glyco_distance(chains, "A", [13])
    check("glyco distance ~15", abs(d_far - 15.0) < 1e-6)
    # epitope = residue 10 (the sequon Asn itself) -> distance 0
    d0 = NF.glyco_distance(chains, "A", [10])
    check("glyco on-sequon distance 0", abs(d0) < 1e-6)
    # no sequon -> cap
    chain2 = [mk(1, "ALA", [0, 0, 0]), mk(2, "GLY", [5, 0, 0])]
    check("glyco no-sequon cap", NF.glyco_distance({"A": chain2}, "A", [1], cap=30.0) == 30.0)


def test_residue_depth_runs():
    # a dense SOLID ball (2 A grid within r<=10): the centre atom is genuinely buried (SASA 0),
    # a boundary atom is solvent-exposed. Residue depth of the buried atom must exceed the surface one.
    g = np.arange(-10, 10.01, 2.0)
    grid = np.array([[x, y, z] for x in g for y in g for z in g])
    xyz = grid[np.linalg.norm(grid, axis=1) <= 10.0]
    elem = ["C"] * len(xyz)
    core = xyz[np.argmin(np.linalg.norm(xyz, axis=1))][None]        # nearest to centre
    surf = xyz[np.argmax(np.linalg.norm(xyz, axis=1))][None]        # on the boundary
    d_core = NF.residue_depth(xyz, elem, core)
    d_surf = NF.residue_depth(xyz, elem, surf)
    check("residue_depth core>surface", d_core > d_surf)
    check("residue_depth finite", np.isfinite(d_core))


# ----------------------------- integration: one real CIF -----------------------------
def test_extract_one_epi_size():
    cif = ROOT / "data/dp4_binding/structures/6OAN.cif"
    if not cif.exists():
        check("extract_one 6oan (SKIP: no CIF)", True); return
    feat, size = NF.extract_one(cif, "A", n_isl=2)
    check("extract_one 6oan epi_size==32", size == 32)
    check("extract_one all features present", all(k in feat for k in NF.NEW_FEATURES if k != "epitope"))
    check("extract_one features finite", all(np.isfinite(feat[k]) for k in
          ["epi_depth", "epi_cn", "ag_out25", "epi_net_charge", "epi_coulomb"]))


# ----------------------------- fit / predict / save / load -----------------------------
def test_fit_predict_path():
    df = D.load_table()
    check("table has 2720 rows", len(df) == 2720)
    check("label present", D.LABEL_COL in df.columns)
    check("all baseline features present", all(f in df.columns for f in D.BASELINE_FEATURES))
    m = BindingModel().fit(df)
    pred = m.predict(df)
    check("predict has binding_pred", "binding_pred" in pred.columns)
    check("predict keyed by design_id", "design_id" in pred.columns)
    check("predict len matches", len(pred) == len(df))
    p = pred["binding_pred"].to_numpy()
    check("binding_pred in [0,1]", (p >= 0).all() and (p <= 1).all())
    # grouped-CV AUC reproduces the published ~0.81 (tolerance for shuffle noise)
    cv = m.cross_val_auc(df, repeats=3)
    check(f"grouped-CV ROC-AUC ~0.81 (got {cv['roc_mean']:.3f})", 0.79 <= cv["roc_mean"] <= 0.84)
    # save / load round-trip gives identical predictions
    with tempfile.TemporaryDirectory() as td:
        path = m.save(Path(td) / "m.pkl")
        m2 = BindingModel.load(path)
        p2 = m2.predict(df)["binding_pred"].to_numpy()
        check("save/load identical preds", np.allclose(p, p2))


def main():
    print("extractor math:");        test_contact_number(); test_charge_features()
    test_glyco_distance(); test_residue_depth_runs()
    print("integration (real CIF):"); test_extract_one_epi_size()
    print("fit/predict path:");       test_fit_predict_path()
    n = len(PASS); k = sum(PASS)
    print(f"\n{k}/{n} checks passed")
    sys.exit(0 if k == n else 1)


if __name__ == "__main__":
    main()
