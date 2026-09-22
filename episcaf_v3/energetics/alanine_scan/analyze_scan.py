"""Pair WT/mutant GB binding energies by frame; report block-size sensitivity.

CSV grammar follows upstream GMXMMPBSA/output_file.py and amber_outputs.py.
No real alanine-scan CSV has yet been returned for this project.
"""
import argparse
import csv
import json
import math
from pathlib import Path
import statistics


def parse(path):
    state = None
    delta = False
    header = None
    result = {'wt': {}, 'mutant': {}}
    with path.open(newline='') as f:
        for row in csv.reader(f):
            if not row:
                header = None
                continue
            first = row[0].strip()
            if 'GENERALIZED BORN:' in first:
                state = 'mutant' if 'Mutant' in first else 'wt'
                delta = False
                header = None
            elif len(row) == 1:
                delta = first == 'Delta Energy Terms'
                header = None
            elif first == 'Frame #':
                header = row if delta and state else None
            elif header:
                frame = int(first)
                if frame in result[state]:
                    raise ValueError(f'Duplicate frame {frame} in {state}')
                values = {k.strip(): float(v) for k, v in zip(header[1:], row[1:])}
                if 'TOTAL' not in values or not all(math.isfinite(v) for v in values.values()):
                    raise ValueError('Missing TOTAL or nonfinite energy')
                result[state][frame] = values
    if not result['wt'] or result['wt'].keys() != result['mutant'].keys():
        raise ValueError('Missing or mismatched WT/mutant binding-energy frames')
    return result


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--run', type=Path, required=True, help='One completed mutation directory')
    args = ap.parse_args()
    manifest = json.loads((args.run / 'run_manifest.json').read_text())
    data = parse(args.run / 'FRAME_ENERGIES.csv')
    frames = sorted(data['wt'])
    if frames != manifest['intended_frames']:
        raise ValueError('Observed frames differ from the frozen plan; do not silently truncate.')
    records = []
    for frame in frames:
        wt, mutant = data['wt'][frame], data['mutant'][frame]
        if wt.keys() != mutant.keys():
            raise ValueError('WT/mutant energy terms differ')
        record = {'frame': frame, 'wt_binding_total': wt['TOTAL'],
                  'mutant_binding_total': mutant['TOTAL'],
                  'ddg_kcal_mol': mutant['TOTAL'] - wt['TOTAL']}
        record.update({f'dd_{k}': mutant[k] - wt[k] for k in wt if k != 'TOTAL'})
        records.append(record)
    values = [r['ddg_kcal_mol'] for r in records]
    blocks = []
    for size in (5, 10, 20, 25):
        if len(values) % size:
            continue
        means = [statistics.mean(values[i:i+size]) for i in range(0, len(values), size)]
        blocks.append({'frames_per_block': size, 'n_blocks': len(means),
                       'block_means': means,
                       'sem_of_block_means': statistics.stdev(means) / math.sqrt(len(means))})
    summary = {'mutation': manifest['task']['name'], 'n_frames': len(values),
               'mean_ddg_kcal_mol': statistics.mean(values),
               'experimental_ddg_kcal_mol': manifest['task']['experimental_ddg_kcal_mol'],
               'block_sensitivity': blocks,
               'limitations': 'Block SEM assumes sufficiently independent blocks; no independent-replica or systematic-error uncertainty. CSV energies rounded to 0.01 kcal/mol upstream. Mutant topology/selection must be checked before interpreting.'}
    with (args.run / 'paired_ddg.csv').open('w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=list(records[0]))
        writer.writeheader()
        writer.writerows(records)
    (args.run / 'paired_summary.json').write_text(json.dumps(summary, indent=2) + '\n')
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()
