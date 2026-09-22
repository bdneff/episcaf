"""Freeze the eight 6A77 measurements from Sheng et al. Table S1 (no API calls)."""
from pathlib import Path
import hashlib
import json
import pandas as pd

HERE = Path(__file__).resolve().parent
SOURCE = HERE / 'reference_inputs'
RUN = HERE / '6a77_non_skempi'


def main():
    RUN.mkdir(exist_ok=True)
    if (RUN / 'request.json').exists():
        raise SystemExit('Frozen run already exists.')
    table = pd.read_excel(SOURCE / 'sheng2023_Table_1.xlsx', sheet_name='Table S1')
    rows = table[table.PDBID == '6A77'].copy()
    assert len(rows) == 8 and rows.mutation.is_unique
    text = (SOURCE / '6a77.pdb').read_text()
    residues = {(l[21], l[22:26].strip() + l[26].strip()): l[17:20]
                for l in text.splitlines() if l.startswith('ATOM  ')}
    code = dict(zip(['GLU','LEU','PRO','SER','TYR'], 'ELPSY'))
    for m in rows.mutation:
        assert code[residues[m[1], m[2:-1]]] == m[0], m
    audit = {}
    for filename, separator in [('skempi_v2.csv', ';'), ('stabddg_filtered_skempi.csv', ',')]:
        data = pd.read_csv(SOURCE / filename, sep=separator)
        hits = data['#Pdb'].str.split('_').str[0].str.upper().eq('6A77').sum()
        assert hits == 0
        audit[filename] = {'pdb_hits': int(hits), 'sha256': hashlib.sha256((SOURCE / filename).read_bytes()).hexdigest()}
    rows[['PDBID', 'mutation', 'Experimental DDG', 'Location']].to_csv(RUN / 'experimental.csv', index=False)
    request = {'jobName': 'episcaf-v3-stabddg-6a77-nonskempi-20260908', 'type': 'stabddg',
               'settings': {'pdbFile': text, 'binder1Chains': ['H', 'L'],
                            'binder2Chains': ['A'], 'mutations': rows.mutation.tolist(),
                            'mcSamples': 20, 'seed': 42}}
    (RUN / 'request.json').write_text(json.dumps(request, indent=2) + '\n')
    manifest = {'pdb': '6A77', 'pdb_sha256': hashlib.sha256(text.encode()).hexdigest(),
                'experimental_source': 'https://doi.org/10.3389/fimmu.2023.1190416',
                'supplement_url': 'https://www.ebi.ac.uk/europepmc/webservices/rest/PMC10235760/supplementaryFiles',
                'table_member': 'Table_1.xlsx',
                'table_sha256': hashlib.sha256((SOURCE / 'sheng2023_Table_1.xlsx').read_bytes()).hexdigest(),
                'audit': audit,
                'scope': 'Exact PDB absent from published SKEMPI tables. Homolog-level and hosted checkpoint provenance unresolved.',
                'selection': 'All eight 6A77 entries, before predictions; other literature systems require additional mapping/background checks.',
                'mutated_partner': 'antibody', 'no_threshold_fit': True}
    (RUN / 'input_manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
    print(rows[['mutation', 'Experimental DDG']].to_string(index=False))
    print('Native identities verified; no 6A77 entries in either SKEMPI table.')


if __name__ == '__main__':
    main()
