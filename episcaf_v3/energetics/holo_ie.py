#!/usr/bin/env python
"""Per-residue antigen-antibody interaction energy from a holo MD trajectory (Decision D2).

For each ANTIGEN residue, the time-averaged nonbonded interaction energy with the ANTIBODY
(reaction-field Coulomb + Lorentz-Berthelot LJ), and separately with SOLVENT. This is the quantity
`bcell_epitope` calls "epitopeness"; it is an interaction energy, NOT a binding free energy (it omits
desolvation and entropy), so treat it as a ranking signal. The v3 test (on 3HFM, where SKEMPI gives
experimental alanine ddG) is whether this ranking recovers the known hot spots. If it does, the cheap
signal suffices; if it misses hot spots driven by desolvation/water, that is the evidence we need a
fuller (MM-GBSA / alanine-scan) treatment.

Adapted directly from bcell_epitope campaign2 `analysis/scripts/holo_ie.py` (same reaction-field +
Lorentz-Berthelot physics; nonbonded params read from the GROMACS topology via OpenMM). The one
change: the antigen can be the LAST N residues, not the first -- our 3HFM build order is L,H,Y, so
lysozyme (antigen) is last. `ag_res_idx` in the output is 1..N in antigen order, i.e. the antigen's
own residue numbering when its residues are complete and in order (3HFM lysozyme = 1..129).

Run on Gemini in an env with OpenMM + MDAnalysis (bcell's `epitope-energy` conda env):
    conda activate epitope-energy
    python holo_ie.py --workdir /scratch/bneff/episcaf_run/episcaf_v3/energetics/md/3hfm/out \
                      --antigen last:129 \
                      --includedir /home/bneff/.conda/envs/grinn/share/gromacs/top \
                      --out holo_ie_mean.csv
"""
import argparse, os, warnings
import numpy as np
warnings.filterwarnings("ignore")


def parse_antigen(spec, n_res):
    """'last:129' | 'first:129' | 'range:430-558' -> set of 0-based residue indices (antigen)."""
    kind, _, val = spec.partition(":")
    if kind == "last":
        n = int(val); return list(range(n_res - n, n_res))
    if kind == "first":
        n = int(val); return list(range(0, n))
    if kind == "range":
        a, b = val.split("-"); return list(range(int(a) - 1, int(b)))  # 1-based inclusive
    raise SystemExit(f"bad --antigen spec {spec!r}")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--workdir", required=True, help="dir with md.tpr, md.xtc, topol.top, md.gro")
    ap.add_argument("--antigen", default="last:129", help="antigen residues: last:N | first:N | range:A-B")
    ap.add_argument("--includedir", default="/home/bneff/.conda/envs/grinn/share/gromacs/top",
                    help="HOST-accessible GROMACS top dir with amber99sb-ildn.ff, so OpenMM can resolve "
                         "the force-field include. NOT the containerized module path "
                         "(/usr/local/gromacs/...): GROMACS runs in a Singularity container on Gemini, "
                         "so that path is invisible to the host conda env. Use a conda env's top dir "
                         "(bcell uses the grinn env's); the force field is standard, params are identical.")
    ap.add_argument("--stride-ps", type=float, default=100.0, help="sample every ~this many ps")
    ap.add_argument("--out", default="holo_ie_mean.csv", help="per-residue mean output (in workdir)")
    ap.add_argument("--perframe", default="holo_ie_perframe.csv", help="per-frame output (in workdir)")
    args = ap.parse_args()
    os.chdir(args.workdir)

    from openmm.app import GromacsTopFile, GromacsGroFile
    from openmm import NonbondedForce, unit, app
    gro = GromacsGroFile("md.gro")
    top = GromacsTopFile("topol.top", periodicBoxVectors=gro.getPeriodicBoxVectors(),
                         includeDir=args.includedir)
    system = top.createSystem(nonbondedMethod=app.NoCutoff, constraints=None)
    nb = [f for f in system.getForces() if isinstance(f, NonbondedForce)][0]
    N = nb.getNumParticles()
    Q = np.empty(N); SIG = np.empty(N); EPS = np.empty(N)
    for i in range(N):
        q, s, e = nb.getParticleParameters(i)
        Q[i] = q.value_in_unit(unit.elementary_charge)
        SIG[i] = s.value_in_unit(unit.nanometer)
        EPS[i] = e.value_in_unit(unit.kilojoule_per_mole)

    SOLVENT = {"HOH", "WAT", "SOL", "TIP3", "NA", "CL", "NA+", "CL-", "K", "MG"}
    resname = []; resid = []; is_solvent = []
    for chain in top.topology.chains():
        for res in chain.residues():
            for a in res.atoms():
                resname.append(res.name); resid.append(res.id); is_solvent.append(res.name in SOLVENT)
    is_solvent = np.array(is_solvent)
    prot_idx = np.where(~is_solvent)[0]
    sol_idx = np.where(is_solvent)[0]

    from collections import OrderedDict
    res_atoms = OrderedDict(); prev = None; ridx = -1
    for i in prot_idx:
        key = (resid[i], resname[i])
        if key != prev:
            ridx += 1; prev = key; res_atoms[ridx] = [key, []]
        res_atoms[ridx][1].append(i)
    n_prot_res = len(res_atoms)

    ag_res = parse_antigen(args.antigen, n_prot_res)
    ab_res = [ri for ri in range(n_prot_res) if ri not in set(ag_res)]
    ab_idx = np.concatenate([np.array(res_atoms[ri][1]) for ri in ab_res])
    print(f"{n_prot_res} protein residues | antigen {len(ag_res)} (first {res_atoms[ag_res[0]][0]}, "
          f"last {res_atoms[ag_res[-1]][0]}) | antibody {len(ab_res)}", flush=True)

    KE = 138.935458; RC = 1.2; eps_in = 1.0; eps_rf = 78.0
    krf = (1.0 / RC**3) * (eps_rf - eps_in) / (2 * eps_rf + eps_in)
    crf = 1.0 / RC + krf * RC**2

    import MDAnalysis as mda
    u = mda.Universe("md.tpr", "md.xtc")
    dt_ps = u.trajectory.dt
    stride = max(1, int(round(args.stride_ps / dt_ps)))
    frames = list(range(0, len(u.trajectory), stride))
    print(f"traj {len(u.trajectory)} frames, dt {dt_ps} ps, stride {stride} -> {len(frames)} sampled",
          flush=True)

    Qsol = Q[sol_idx]; Ssol = SIG[sol_idx]; Esol = EPS[sol_idx]
    Qab = Q[ab_idx]; Sab = SIG[ab_idx]; Eab = EPS[ab_idx]

    def pair_energy(pos_a, q_a, s_a, e_a, pos_b, q_b, s_b, e_b, box):
        d = pos_a[:, None, :] - pos_b[None, :, :]
        d -= box * np.round(d / box)
        r = np.sqrt((d * d).sum(2)); r = np.clip(r, 0.06, None)
        within = r <= RC
        qq = q_a[:, None] * q_b[None, :]
        coul = np.where(within, KE * qq * (1.0 / r + krf * r * r - crf), 0.0).sum()
        sij = 0.5 * (s_a[:, None] + s_b[None, :]); eij = np.sqrt(e_a[:, None] * e_b[None, :])
        sr6 = (sij / r)**6
        lj = np.where(within, 4 * eij * (sr6 * sr6 - sr6), 0.0).sum()
        return float(coul), float(lj)

    per = []
    for k, fno in enumerate(frames):
        u.trajectory[fno]
        pos = u.atoms.positions / 10.0
        box = u.dimensions[:3] / 10.0
        t_ns = u.trajectory.time / 1000.0
        sol_pos = pos[sol_idx]; ab_pos = pos[ab_idx]
        for j, ri in enumerate(ag_res):
            (rid, rn), idxs = res_atoms[ri]; idxs = np.array(idxs)
            rpos = pos[idxs]; rq = Q[idxs]; rs = SIG[idxs]; re = EPS[idxs]
            com = rpos.mean(0)
            dab = ab_pos - com; dab -= box * np.round(dab / box)
            nab = np.where(np.linalg.norm(dab, axis=1) < (RC + 1.2))[0]
            c_ab, l_ab = (pair_energy(rpos, rq, rs, re, ab_pos[nab], Qab[nab], Sab[nab], Eab[nab], box)
                          if len(nab) else (0.0, 0.0))
            ds = sol_pos - com; ds -= box * np.round(ds / box)
            ns = np.where(np.linalg.norm(ds, axis=1) < (RC + 1.2))[0]
            c_s, l_s = (pair_energy(rpos, rq, rs, re, sol_pos[ns], Qsol[ns], Ssol[ns], Esol[ns], box)
                        if len(ns) else (0.0, 0.0))
            per.append((round(t_ns, 3), j + 1, rid, rn, c_ab, l_ab, c_ab + l_ab, c_s, l_s, c_s + l_s))
        if (k + 1) % 25 == 0:
            print(f"  {k+1}/{len(frames)} frames", flush=True)

    import csv
    with open(args.perframe, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["t_ns", "ag_res_idx", "resid", "resname",
                    "ab_coul", "ab_lj", "ab_total", "sol_coul", "sol_lj", "sol_total"])
        w.writerows(per)

    arr = {}
    for row in per:
        arr.setdefault((row[1], row[2], row[3]), []).append(row[4:])
    with open(args.out, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["ag_res_idx", "resid", "resname", "ab_coul", "ab_lj", "ab_total",
                    "sol_coul", "sol_lj", "sol_total", "n_frames"])
        for (ai, rid, rn), vals in sorted(arr.items()):
            m = np.array(vals).mean(0)
            w.writerow([ai, rid, rn] + [round(x, 4) for x in m] + [len(vals)])
    print(f"DONE: {len(per)} residue-frame rows, {len(arr)} antigen residues -> {args.out}", flush=True)


if __name__ == "__main__":
    main()
