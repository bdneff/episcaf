"""Run ONE staged task on the user's compute machine; never submit a scheduler job."""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess

HERE = Path(__file__).resolve().parent


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--source', type=Path, required=True)
    ap.add_argument('--out', type=Path, required=True)
    ap.add_argument('--task', type=int, required=True)
    args = ap.parse_args()
    source, output = args.source.resolve(), args.out.resolve()
    if source == output or source in output.parents:
        raise SystemExit('Output must be outside the existing MD directory.')
    manifest = json.loads((HERE / 'manifest.json').read_text())
    if not 0 <= args.task < len(manifest['tasks']):
        raise SystemExit('Task index out of range.')
    task = manifest['tasks'][args.task]
    for filename in ('md.tpr', 'md_whole.xtc', 'topol.top', 'gbsa_index.ndx'):
        if not (source / filename).is_file():
            raise SystemExit(f'Missing input: {source / filename}')
    groups = [line.strip()[1:-1].strip() for line in
              (source / 'gbsa_index.ndx').read_text().splitlines() if line.strip().startswith('[')]
    if groups.count('antibody') != 1 or groups.count('antigen') != 1:
        raise SystemExit('Index must have unique groups named antibody and antigen.')
    indices = [str(groups.index(name)) for name in ('antibody', 'antigen')]
    config = HERE / 'configs' / f'{task["name"]}.in'
    reference = HERE.parent / 'md/3hfm/structures/3hfm.pdb'
    assert digest(config) == task['config_sha256']
    assert digest(reference) == manifest['pdb_sha256']
    work = output / task['name']
    work.mkdir(parents=True, exist_ok=False)  # refuse to overwrite any earlier run
    # Preserve relative topology includes. Installed force-field includes remain tied
    # to the recorded environment; generated Amber topologies are retained as evidence.
    for path in source.rglob('*'):
        if path.is_file() and path.suffix in ('.top', '.itp'):
            target = work / path.relative_to(source)
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(path, target)
    shutil.copy2(config, work / 'mmpbsa.in')
    shutil.copy2(HERE / 'run_scan.py', work / 'runner_source.py')
    shutil.copy2(HERE / 'analyze_scan.py', work / 'analysis_source.py')
    shutil.copy2(HERE / 'manifest.json', work / 'plan_manifest.json')
    shutil.copy2(source / 'gbsa_index.ndx', work / 'gbsa_index.ndx')
    (work / 'reference.pdb').write_text('\n'.join(
        l for l in reference.read_text().splitlines() if l.startswith(('ATOM  ', 'TER   '))) + '\nEND\n')
    hashes = {str(p): digest(p) for p in [source / 'md.tpr', source / 'md_whole.xtc']}
    hashes.update({str(p.relative_to(work)): digest(p) for p in work.rglob('*') if p.is_file()})
    command = ['gmx_MMPBSA', '-O', '-i', 'mmpbsa.in', '-cs', str(source / 'md.tpr'),
               '-ct', str(source / 'md_whole.xtc'), '-ci', 'gbsa_index.ndx', '-cg', *indices,
               '-cp', 'topol.top', '-cr', 'reference.pdb', '-o', 'FINAL_RESULTS.dat',
               '-eo', 'FRAME_ENERGIES.csv', '-nogui']
    (work / 'run_manifest.json').write_text(json.dumps(
        {'task': task, 'source': str(source), 'hashes': hashes, 'command': command,
         'intended_frames': manifest['frames']}, indent=2) + '\n')
    with (work / 'versions.txt').open('w') as f:
        for cmd in (['gmx', '--version'], ['gmx_MMPBSA', '--version'], ['conda', 'list', '--explicit']):
            subprocess.run(cmd, stdout=f, stderr=subprocess.STDOUT, check=True)
    with (work / 'run.log').open('w') as f:
        subprocess.run(command, cwd=work, stdout=f, stderr=subprocess.STDOUT, check=True)
    (work / 'CALCULATION_COMPLETE').write_text('Output parsing and mapping checks still required.\n')


if __name__ == '__main__':
    main()
