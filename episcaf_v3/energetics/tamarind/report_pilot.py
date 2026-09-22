"""Rebuild the pilot comparison and LaTeX table from saved Tamarind results."""
import csv
import argparse
from datetime import datetime
import hashlib
import io
import json
from pathlib import Path
from zipfile import ZipFile

HERE = Path(__file__).resolve().parent
RUN = HERE / '3hfm_pilot'


def main():
    global RUN
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run-dir', type=Path, default=RUN)
    parser.add_argument('--table-name', default='tamarind_pilot_table.tex')
    args = parser.parse_args()
    RUN = args.run_dir.resolve()
    request = json.loads((RUN / 'request.json').read_text())
    status = json.loads((RUN / 'status.json').read_text())
    manifest = json.loads((RUN / 'input_manifest.json').read_text())
    assert status['JobStatus'] == 'Complete'
    assert status['JobName'] == request['jobName']
    with (HERE.parent / 'skempi_3hfm_ddg.csv').open() as f:
        experiment = {int(r['resid']): r for r in csv.DictReader(
            line for line in f if not line.startswith('#'))}
    with ZipFile(RUN / 'results.zip') as archive:
        output = archive.read('output.csv')
        rows = list(csv.DictReader(io.StringIO(output.decode())))
        assert {r['Mutation'] for r in rows} == set(request['settings']['mutations'])
        assert len(rows) == len(request['settings']['mutations'])
        settings = json.loads(archive.read('settings.json'))
        for key in ('binder1Chains', 'binder2Chains', 'mutations', 'mcSamples', 'seed'):
            assert settings[key] == request['settings'][key], key
        # Ensure the returned input matches the exact submitted bytes.
        assert hashlib.sha256(archive.read('input.pdb')).hexdigest() == manifest['sha256']
        processed = archive.read('pdbs/input.pdb').decode().splitlines()
        y_residues = {(int(l[22:26]), l[17:20]) for l in processed
                      if l.startswith('ATOM  ') and l[21] == 'Y'}
        assert (96, 'LYS') in y_residues and (73, 'ARG') in y_residues
        # Save small original evidence without rewriting content.
        for name in ('output.csv', 'input.csv', 'output.log'):
            (RUN / name).write_bytes(archive.read(name))
    comparison = []
    for row in rows:
        mutation = row['Mutation']
        assert row['RenumberedMutation'] == mutation
        assert mutation[1] == 'Y' and mutation[-1] == 'A'
        reference = experiment[int(mutation[2:-1])]
        prediction = float(row['Prediction'])
        measured = float(reference['ddg_kcal_mol'])
        comparison.append({'mutation': mutation, 'predicted_kcal_mol': prediction,
                           'experimental_kcal_mol': measured,
                           'signed_error_kcal_mol': prediction - measured})
    with (RUN / 'comparison.csv').open('w') as f:
        writer = csv.DictWriter(f, fieldnames=list(comparison[0]))
        writer.writeheader()
        writer.writerows(comparison)
    fmt = '%Y-%m-%d %H:%M:%S'
    elapsed = lambda a, b: (datetime.strptime(status[b], fmt) -
                            datetime.strptime(status[a], fmt)).total_seconds()
    summary = {'job_name': status['JobName'], 'execution_seconds': elapsed('Started', 'Completed'),
               'queue_seconds': elapsed('Created', 'Started'),
               'weighted_hours': status.get('WeightedHours'),
               'results_sha256': hashlib.sha256((RUN / 'results.zip').read_bytes()).hexdigest(),
               'processed_antigen_atom_residues': len(y_residues)}
    (RUN / 'summary.json').write_text(json.dumps(summary, indent=2) + '\n')
    table = [r'\begin{tabular}{lrrr}', r'\toprule',
             r'Mutation & Experiment & StaB-ddG & Prediction $-$ experiment \\', r'\midrule']
    for row in comparison:
        label = row['mutation'][0] + row['mutation'][2:]
        table.append(f"{label} & {row['experimental_kcal_mol']:+.2f} & "
                     f"{row['predicted_kcal_mol']:+.2f} & {row['signed_error_kcal_mol']:+.2f} " + r'\\')
    table += [r'\bottomrule', r'\end{tabular}']
    target = HERE.parent.parent / 'manuscript/figures' / args.table_name
    target.write_text('\n'.join(table) + '\n')
    print(json.dumps({'comparison': comparison, 'run': summary}, indent=2))


if __name__ == '__main__':
    main()
