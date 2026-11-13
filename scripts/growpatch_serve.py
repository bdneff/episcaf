#!/usr/bin/env python3
"""Local growPatch tool: upload a PDB, run the real enumerator, view the candidate epitopes.

GitHub Pages is static and cannot run Python, and growPatch needs Python (numpy, Bio.PDB, freesasa),
so the interactive "run it on your own antigen" tool is served locally by this script rather than on
Pages. It serves docs/growpatch_tool.html (the upload + 3D/flat viewer) and answers POST /enumerate by
running episcaf's actual, validated enumerator on the uploaded structure -- no reimplementation, the
same code the pipeline uses.

Run:
    /Users/bneff/micromamba/envs/bcell-repair/bin/python scripts/growpatch_serve.py --port 8765
then open http://localhost:8765/ and upload a PDB.

Note: input validation (is this antigen suitable for episcaf?) is deferred -- a later gate.
"""
from __future__ import annotations

import argparse
import functools
import http.server
import json
import os
import re
import sys
import tempfile
import urllib.parse

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from episcaf_analysis.growpatch import read as R   # noqa: E402
from episcaf_pipeline import build_growpatch_designs as BGD   # noqa: E402

DOCS = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "docs")


def _resi(label):
    m = re.match(r"-?\d+", str(label))
    return int(m.group()) if m else None


def viewer_data(cif_path, chain, sasa_cut, footprint_radius, methods, top_k, mode="growpatch"):
    """Residues + candidate patches for the viewer (like docs/data/growpatch_*.js).

    mode 'growpatch' -> enumerate candidate epitopes (lone antigen).
    mode 'antibody'  -> the epitope is already defined by the bound partner; return it as one patch,
                        no enumeration (growPatch unnecessary).
    """
    ag = R.read_antigen(cif_path, chain, relsasa=True)
    residues = [dict(resi=_resi(ag.label(c, s)), resn=r["resname"],
                     rel=round(float(ag.rs(r["key"])), 3),
                     surf=bool(ag.rs(r["key"]) >= sasa_cut),
                     p=[round(float(x), 2) for x in r["rep"]])
                for r in ag.residues for (c, s) in [r["key"]]]
    extra = {}
    if mode == "antibody":
        det = R.detect_complex(cif_path, chain)
        epi = det["epitope"]
        resis = sorted({_resi(ag.label(c, s)) for (c, s) in ag.parse_epitope(epi)}) if epi else []
        patches = ([dict(i=0, methods=["antibody-contact"], epitope=epi, n_res=len(resis),
                         n_isl=len([t for t in epi.split(";") if t]), resis=resis)] if resis else [])
        extra = dict(mode="antibody", partner_chains=det["partner_chains"])
    else:
        cands, _ = R.enumerate_candidates(cif_path, chain, antigen=ag, methods=methods,
                                          sasa_cut=sasa_cut, footprint_radius=footprint_radius,
                                          top_k=(top_k or None))
        patches = [dict(i=p["patch"], methods=p["methods"], epitope=p["epitope"],
                        n_res=p["n_res"], n_isl=p["n_islands"],
                        resis=sorted({_resi(ag.label(c, s)) for (c, s) in p["keys"]}))
                   for p in cands]
        extra = dict(mode="growpatch")
    return dict(chain=chain, sasa_cut=sasa_cut, n_residues=len(residues),
                n_surface=sum(r["surf"] for r in residues), n_patches=len(patches),
                n_growpatch=sum("growpatch" in p["methods"] for p in patches),
                n_eyeview=sum("eyeview" in p["methods"] for p in patches),
                n_both=sum(len(p["methods"]) == 2 for p in patches),
                residues=residues, patches=patches, **extra)


def stage_run(cif_path, chain, epitopes, cfg):
    """Build the RFd3 contigs for the chosen epitopes + a run spec the cluster consumes.

    Returns (run_spec dict, contigs_csv str, n_contigs, skipped). The heavy run (RFd3->MPNN->AF3) is
    cluster/GPU work; this stages the job, it does not run it.
    """
    import csv as _csv
    import datetime
    import io
    chains = R.parse_chains(chain)
    cands = [dict(patch=i, epitope=e, n_islands=len([t for t in e.split(";") if t]),
                  n_res=0, methods=[]) for i, e in enumerate(epitopes)]
    rows, skipped = BGD.build_rows(cfg["name"], chains, cands, cfg["total_len"],
                                   cfg["variants"], cfg["seed"])
    cols = ["id", "contig_id", "contig_string", "contig_length", "epitope_chunks", "source",
            "pdb_id", "chains", "growpatch_epitope", "growpatch_methods", "growpatch_n_res", "variant"]
    buf = io.StringIO()
    w = _csv.DictWriter(buf, fieldnames=cols)
    w.writeheader()
    w.writerows({k: r.get(k) for k in cols} for r in rows)
    n_epi = len(epitopes)
    run_spec = {
        "created": datetime.datetime.now().isoformat(timespec="seconds"),
        "antigen": cfg["name"], "chain": chain, "mode": cfg.get("mode", "growpatch"),
        "n_epitopes": n_epi, "epitopes": epitopes,
        "generation": {"flank_variants": cfg["variants"], "total_len": cfg["total_len"],
                       "rfd3_designs_per_contig": cfg["rfd3"], "mpnn_seqs_per_backbone": cfg["mpnn"]},
        "budget": cfg["budget"],
        "n_contigs": len(rows),
        "designs_at_full_fanout": len(rows) * cfg["rfd3"] * cfg["mpnn"],
        "note": "Stage only. Run RFd3->MPNN->AF3 on the cluster; score + view the metrics CSV in "
                "the Results tab. accessibility metric = AOMclash (aomclash_pred).",
    }
    return run_spec, buf.getvalue(), len(rows), skipped


class Handler(http.server.SimpleHTTPRequestHandler):
    def _send(self, code, body, ctype="application/json"):
        data = body if isinstance(body, bytes) else body.encode()
        self.send_response(code)
        self.send_header("Content-Type", ctype)
        self.send_header("Content-Length", str(len(data)))
        self.end_headers()
        self.wfile.write(data)

    def do_GET(self):
        if self.path == "/" or self.path.startswith("/?"):
            try:
                with open(os.path.join(DOCS, "growpatch_tool.html"), "rb") as fh:
                    self._send(200, fh.read(), "text/html; charset=utf-8")
            except FileNotFoundError:
                self._send(404, b"growpatch_tool.html not found")
            return
        return super().do_GET()                        # static files from docs/ (pdb, data, ...)

    def do_POST(self):
        parsed = urllib.parse.urlparse(self.path)
        if parsed.path not in ("/enumerate", "/inspect", "/stage"):
            self._send(404, json.dumps({"error": "unknown endpoint"}))
            return
        q = urllib.parse.parse_qs(parsed.query)
        body = self.rfile.read(int(self.headers.get("Content-Length", 0)))
        fmt = q.get("fmt", ["pdb"])[0]
        chain = q.get("chain", [""])[0] or None
        tmp = tempfile.NamedTemporaryFile(suffix="." + ("cif" if fmt == "cif" else "pdb"),
                                          delete=False, mode="wb")
        try:
            tmp.write(body)
            tmp.close()
            if parsed.path == "/inspect":
                self._send(200, json.dumps(R.detect_complex(tmp.name, chain)))
                return
            if parsed.path == "/stage":
                epitopes = [e for e in q.get("epitopes", [""])[0].split("|") if e]
                if not epitopes:
                    self._send(400, json.dumps({"error": "no epitopes selected"}))
                    return
                try:
                    cfg = dict(name=(q.get("name", ["antigen"])[0] or "antigen"),
                               variants=int(q.get("variants", ["4"])[0]),
                               total_len=int(q.get("total_len", ["103"])[0]),
                               rfd3=int(q.get("rfd3", ["8"])[0]), mpnn=int(q.get("mpnn", ["8"])[0]),
                               seed=int(q.get("seed", ["0"])[0]), mode=q.get("mode", ["growpatch"])[0],
                               budget={"mode": q.get("budget_mode", ["n_designs"])[0],
                                       "value": float(q.get("budget_value", ["200"])[0])})
                except ValueError as e:
                    self._send(400, json.dumps({"error": f"bad parameter: {e}"}))
                    return
                spec, csv_text, n, skipped = stage_run(tmp.name, chain or "A", epitopes, cfg)
                self._send(200, json.dumps(dict(run_spec=spec, contigs_csv=csv_text,
                                                n_contigs=n, skipped=len(skipped))))
                return
            try:
                sasa_cut = float(q.get("sasa_cut", ["0.2"])[0])
                footprint = float(q.get("footprint_radius", ["14"])[0])
                top_k = int(q.get("top_k", ["0"])[0])
                methods = tuple(m for m in q.get("methods", ["growpatch,eyeview"])[0].split(",") if m)
            except ValueError as e:
                self._send(400, json.dumps({"error": f"bad parameter: {e}"}))
                return
            mode = q.get("mode", ["growpatch"])[0]
            data = viewer_data(tmp.name, chain or "A", sasa_cut, footprint, methods, top_k, mode)
            self._send(200, json.dumps(data))
        except Exception as e:                          # a bad upload must not crash the server
            self._send(400, json.dumps({"error": f"{type(e).__name__}: {e}"}))
        finally:
            os.unlink(tmp.name)

    def log_message(self, fmt, *args):                  # quieter console
        sys.stderr.write("  " + (fmt % args) + "\n")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--port", type=int, default=8765)
    ap.add_argument("--host", default="127.0.0.1")
    a = ap.parse_args(argv)
    handler = functools.partial(Handler, directory=DOCS)
    srv = http.server.ThreadingHTTPServer((a.host, a.port), handler)
    print(f"growPatch tool on http://{a.host}:{a.port}/  (serving {DOCS}; Ctrl-C to stop)")
    try:
        srv.serve_forever()
    except KeyboardInterrupt:
        print("\nstopped")


if __name__ == "__main__":
    main()
