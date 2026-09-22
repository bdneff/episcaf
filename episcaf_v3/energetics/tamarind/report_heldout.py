"""Verify 6A77 service remapping and compare its eight predictions to Table S1."""
from pathlib import Path
from zipfile import ZipFile
import hashlib
import json
import io
import pandas as pd
from scipy.stats import spearmanr, pearsonr
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
RUN = HERE / '6a77_non_skempi'


def residues(text):
    result = {}
    for line in text.splitlines():
        if line.startswith('ATOM  '):
            key = (line[21], line[22:26].strip() + line[26].strip())
            result[key] = line[17:20]
    return result


def main():
    request = json.loads((RUN / 'request.json').read_text())
    manifest = json.loads((RUN / 'input_manifest.json').read_text())
    with ZipFile(RUN / 'results.zip') as z:
        assert hashlib.sha256(z.read('input.pdb')).hexdigest() == manifest['pdb_sha256']
        original = residues(z.read('input.pdb').decode())
        processed = residues(z.read('pdbs/input.pdb').decode())
        mapping = {}
        for chain in ['H', 'L', 'A']:
            a = [(key, value) for key, value in original.items() if key[0] == chain]
            b = [(key, value) for key, value in processed.items() if key[0] == chain]
            assert [v for k, v in a] == [v for k, v in b], chain
            mapping.update({ka: kb for (ka, va), (kb, vb) in zip(a, b)})
        settings = json.loads(z.read('settings.json'))
        for key in ['binder1Chains', 'binder2Chains', 'mutations', 'seed', 'mcSamples']:
            assert settings[key] == request['settings'][key]
        predictions = pd.read_csv(io.BytesIO(z.read('output.csv')))
        for row in predictions.itertuples():
            old, new = row.Mutation, row.RenumberedMutation
            assert old[0] == new[0] and old[-1] == new[-1]
            assert mapping[old[1], old[2:-1]] == (new[1], new[2:-1])
        assert set(predictions.Mutation) == set(request['settings']['mutations'])
        assert predictions.Mutation.is_unique
        for name in ['output.csv', 'input.csv', 'output.log']:
            (RUN / name).write_bytes(z.read(name))
    exp = pd.read_csv(RUN / 'experimental.csv')
    data = predictions.merge(exp, left_on='Mutation', right_on='mutation', validate='one_to_one')
    data['error_kcal_mol'] = data.Prediction - data['Experimental DDG']
    data.to_csv(RUN / 'comparison.csv', index=False)
    x, y = data['Experimental DDG'], data.Prediction
    stats = {'n': len(data), 'spearman': float(spearmanr(x,y).correlation),
             'pearson': float(pearsonr(x,y)[0]), 'mae_kcal_mol': float((y-x).abs().mean()),
             'rmse_kcal_mol': float(((y-x)**2).mean()**0.5),
             'results_sha256': hashlib.sha256((RUN / 'results.zip').read_bytes()).hexdigest(),
             'mapping_verified': True}
    (RUN / 'statistics.json').write_text(json.dumps(stats, indent=2) + '\n')
    pilot = pd.read_csv(HERE / '3hfm_all13/comparison.csv')
    fig, axes = plt.subplots(1, 2, figsize=(11, 5))
    for ax, xx, yy, title in [(axes[0], pilot.experimental_kcal_mol, pilot.predicted_kcal_mol,
                              '3HFM · SKEMPI antigen scan (n=13)'),
                             (axes[1], x, y, '6A77 · non-SKEMPI antibody scan (n=8)')]:
        ax.scatter(xx, yy, color='#0072B2', s=55, zorder=3)
        ax.plot([-2.5,7],[-2.5,7], '--', color='gray', lw=1)
        ax.axhline(0, color='gray', lw=0.7)
        ax.axvline(0, color='gray', lw=0.7)
        ax.set_xlim(-2.5,7); ax.set_ylim(-2.5,7)
        ax.set_title(title + f'\nSpearman ρ={spearmanr(xx,yy).correlation:+.2f}', fontsize=11)
        ax.set_xlabel('Experimental ΔΔG (kcal/mol)')
        ax.set_ylabel('StaB-ddG prediction (kcal/mol)')
        ax.spines[['top','right']].set_visible(False)
        ax.grid(alpha=0.15)
    for _, row in data.iterrows():
        axes[1].annotate(row.Mutation, (row['Experimental DDG'], row.Prediction),
                         xytext=(4,4), textcoords='offset points', fontsize=8)
    fig.text(0.5,0.015,'Positive = weaker binding. Same seed/settings; no calibration or threshold fitting.\n'
             '6A77 is absent by PDB ID from published SKEMPI tables; hosted checkpoint and homology audit remain open.',
             ha='center', fontsize=9)
    fig.tight_layout(rect=(0,0.1,1,1))
    for extension in ('png','pdf'):
        fig.savefig(HERE.parent.parent / 'manuscript/figures' / f'tamarind_6a77_check.{extension}', dpi=180)
    print(json.dumps(stats, indent=2))
    print(data[['Mutation','Prediction','Experimental DDG']].to_string(index=False))


if __name__ == '__main__':
    main()
