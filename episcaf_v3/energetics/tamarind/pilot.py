"""Frozen 3HFM StaB-ddG pilot: prepare, validate, submit once, status, download.

Run from any directory. Credentials remain outside the repository.
"""
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import urllib.request
import urllib.parse

HERE = Path(__file__).resolve().parent
RUN = HERE / '3hfm_pilot'
PDB = HERE.parent / 'md/3hfm/structures/3hfm.pdb'
BASE = 'https://app.tamarind.bio/api/'


def save(name, value):
    (RUN / name).write_text(json.dumps(value, indent=2) + '\n')


def api(path, payload=None):
    key = os.environ.get('TAMARIND_API_KEY') or (
        Path.home() / '.config/tamarind/api_key').read_text().strip()
    request = urllib.request.Request(
        BASE + path,
        data=None if payload is None else json.dumps(payload).encode(),
        headers={'x-api-key': key, 'Content-Type': 'application/json'})
    with urllib.request.urlopen(request, timeout=60) as response:
        raw = response.read().decode()
        try:
            return response.status, json.loads(raw)
        except json.JSONDecodeError:
            return response.status, raw


def main():
    global RUN
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=['prepare', 'verify', 'validate', 'submit', 'status', 'download'])
    parser.add_argument('--run-dir', type=Path, default=RUN)
    parser.add_argument('--job-name', default='episcaf-v3-3hfm-stabddg-K96A-R73A-20260908')
    args = parser.parse_args()
    RUN = args.run_dir.resolve()
    RUN.mkdir(exist_ok=True, parents=True)
    if args.action == 'prepare':
        if (RUN / 'request.json').exists():
            raise SystemExit('Frozen request already exists; will not overwrite.')
        pdb = PDB.read_text()
        residues = {(line[21], int(line[22:26]), line[26]): line[17:20]
                    for line in pdb.splitlines() if line.startswith('ATOM  ')}
        assert {k[0] for k in residues} == {'H', 'L', 'Y'}
        assert residues['Y', 96, ' '] == 'LYS'
        assert residues['Y', 73, ' '] == 'ARG'
        save('request.json', {
            'jobName': args.job_name,
            'type': 'stabddg',
            'settings': {'pdbFile': pdb, 'binder1Chains': ['H', 'L'],
                         'binder2Chains': ['Y'], 'mutations': ['KY96A', 'RY73A'],
                         'mcSamples': 20, 'seed': 42}})
        save('input_manifest.json', {'source': str(PDB.relative_to(HERE.parent)),
             'sha256': hashlib.sha256(PDB.read_bytes()).hexdigest(),
             'mutations': ['KY96A', 'RY73A'], 'separate_single_mutants': True})
        print('Prepared two verified single mutations; seed 42, MC samples 20.')
        return
    request = json.loads((RUN / 'request.json').read_text())
    if args.action == 'verify':
        manifest = json.loads((RUN / 'input_manifest.json').read_text())
        assert hashlib.sha256(request['settings']['pdbFile'].encode()).hexdigest() == manifest['sha256']
        assert hashlib.sha256(PDB.read_bytes()).hexdigest() == manifest['sha256']
        assert request['settings']['mutations'] == manifest['mutations']
        with (HERE.parent / 'skempi_3hfm_ddg.csv').open() as handle:
            rows = {int(row['resid']): row for row in csv.DictReader(
                line for line in handle if not line.startswith('#'))}
        print('Frozen inline PDB and source SHA256 match:', manifest['sha256'])
        for resid in (96, 73):
            row = rows[resid]
            print(f"Y:{row['resname']}{resid}A experimental ddG = {row['ddg_kcal_mol']} kcal/mol; n={row['n_meas']}")
        print('Settings:', {k: v for k, v in request['settings'].items() if k != 'pdbFile'})
    elif args.action == 'validate':
        _, schema = api('tools/stabddg/schema')
        save('schema.json', schema)
        _, result = api('validate-job', request)
        save('validation.json', result)
        print(json.dumps({k: v for k, v in result.items() if k != 'normalized'}))
        if not result.get('valid') or result.get('unrecognized_settings'):
            raise SystemExit('Validation failed; no submission.')
    elif args.action == 'submit':
        validation = json.loads((RUN / 'validation.json').read_text())
        assert validation['valid'] and not validation.get('unrecognized_settings')
        if (RUN / 'submission_attempt.json').exists():
            raise SystemExit('Submission already attempted. Check status; do not resubmit blindly.')
        request['settings'] = validation['normalized']
        save('submitted_request.json', request)
        save('submission_attempt.json', {'jobName': request['jobName']})
        _, result = api('submit-job', request)
        save('submission.json', result)
        print(result)
    elif args.action == 'status':
        _, result = api('jobs?' + urllib.parse.urlencode({'jobName': request['jobName']}))
        save('status.json', result)
        print(json.dumps({k: result.get(k) for k in
              ['JobName', 'JobStatus', 'Created', 'Started', 'Completed', 'WeightedHours']}))
    else:
        status, url = api('result', {'jobName': request['jobName']})
        if status == 202:
            print('Archive preparing; retry download later.')
            return
        assert isinstance(url, str) and url.startswith('https://')
        # Presigned download needs no API key; do not save or print its URL.
        with urllib.request.urlopen(url, timeout=60) as response:
            (RUN / 'results.zip').write_bytes(response.read())
        print('Saved results.zip')


if __name__ == '__main__':
    main()
