"""Restore exact archived reference files needed to reproduce the non-SKEMPI audit."""
import gzip
import hashlib
import json
from pathlib import Path

directory = Path(__file__).resolve().parent / 'reference_inputs'
for name, source in json.loads((directory / 'sources.json').read_text()).items():
    data = gzip.decompress((directory / (name + '.gz')).read_bytes())
    assert hashlib.sha256(data).hexdigest() == source['sha256'], name
    target = directory / name
    if target.exists() and target.read_bytes() != data:
        raise SystemExit(f'Refusing to replace different reference bytes: {target}')
    target.write_bytes(data)
    print('Verified/restored', name)
