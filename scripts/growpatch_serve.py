#!/usr/bin/env python3
"""Local growPatch tool: upload a PDB, run the real enumerator, view the candidate epitopes.

GitHub Pages is static and cannot run Python, and growPatch needs Python (numpy, Bio.PDB, freesasa),
so the interactive "run it on your own antigen" tool is served locally by this script rather than on
Pages. It serves docs/growpatch_tool.html (the upload + 3D/flat viewer) and answers POST /enumerate by
running episcaf's actual, validated enumerator on the uploaded structure -- no reimplementation, the
same code the pipeline uses.

Run:
    python scripts/growpatch_serve.py
then open http://localhost:8765/ and upload a PDB.

Note: input validation (is this antigen suitable for episcaf?) is deferred -- a later gate.
"""
from __future__ import annotations

import argparse
import functools
import http.server
import json
import os
import sys
import tempfile
import urllib.parse

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from episcaf_analysis.growpatch import read as R   # noqa: E402
# viewer_data / stage_run live in episcaf_service.staging (one source of truth, shared with the
# cluster FastAPI service); _resi is re-exported so this module's callers keep working unchanged.
from episcaf_service.staging import _resi, stage_run, viewer_data   # noqa: E402,F401

DOCS = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "docs")


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
                with open(os.path.join(DOCS, "index.html"), "rb") as fh:  # redirects to home.html
                    self._send(200, fh.read(), "text/html; charset=utf-8")
            except FileNotFoundError:
                self._send(404, b"index.html not found")
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
                               helix_bias=float(q.get("helix_bias", ["0"])[0]),
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
                methods = tuple(m for m in q.get("methods", ["growpatch"])[0].split(",") if m)
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
