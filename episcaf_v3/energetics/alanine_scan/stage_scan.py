"""Generate frozen 3HFM paired WT/alanine GB inputs from the measured antigen scan."""
import csv
import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
ENERGY = HERE.parent


def main():
    pdb = ENERGY / 'md/3hfm/structures/3hfm.pdb'
    residues = {(l[21], int(l[22:26]), l[26]): l[17:20]
                for l in pdb.read_text().splitlines() if l.startswith('ATOM  ')}
    with (ENERGY / 'skempi_3hfm_ddg.csv').open() as f:
        rows = list(csv.DictReader(l for l in f if not l.startswith('#')))
    tasks = []
    (HERE / 'configs').mkdir(exist_ok=True)
    for task, row in enumerate(rows):
        resid = int(row['resid'])
        assert residues['Y', resid, ' '] == row['resname']
        name = f"Y{resid}_ALA"
        config = f'''3HFM paired WT/alanine endpoint estimate: Y:{resid}
&general
 sys_name = "3HFM_{name}",
 startframe = 1, endframe = 2000, interval = 20,
 temperature = 300.0, PBRadii = 3,
 keep_files = 2,
/
&gb
 igb = 5, saltcon = 0.150, intdiel = 1.0, extdiel = 78.5,
 surften = 0.0072, surfoff = 0.0,
/
&alanine_scanning
 mutant_res = "Y:{resid}", mutant = "ALA",
 mutant_only = 0, cas_intdiel = 0,
/
'''
        (HERE / 'configs' / f'{name}.in').write_text(config)
        tasks.append({'task': task, 'name': name, 'resid': resid, 'resname': row['resname'],
                      'experimental_ddg_kcal_mol': float(row['ddg_kcal_mol']),
                      'config_sha256': hashlib.sha256(config.encode()).hexdigest()})
    manifest = {'pdb_sha256': hashlib.sha256(pdb.read_bytes()).hexdigest(),
                'frames': list(range(1, 2001, 20)), 'tasks': tasks,
                'observable': '(Gcomplex-Gantibody-Gantigen)_mut - (Gcomplex-Gantibody-Gantigen)_WT',
                'status': 'Staged only; no cluster execution or mutant topology verification yet.'}
    (HERE / 'manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
    print(f'Generated {len(tasks)} verified residue selections; {len(manifest["frames"])} matched frames each.')


if __name__ == '__main__':
    main()
