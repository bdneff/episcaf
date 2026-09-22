"""Plot saved StaB-ddG and simulation signals against the same experimental scan."""
from pathlib import Path
import argparse
import json
import pandas as pd
from scipy.stats import spearmanr, pearsonr
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
ENERGY = HERE.parent
OUT = ENERGY.parent / 'manuscript/figures'


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run-dir', type=Path, default=HERE / '3hfm_pilot')
    args = parser.parse_args()
    run = args.run_dir.resolve()
    exp = pd.read_csv(ENERGY / 'skempi_3hfm_ddg.csv', comment='#')
    ie = pd.read_csv(ENERGY / 'md/3hfm/holo_ie_mean.csv')
    assert ie.n_frames.eq(201).all()
    gb = pd.read_csv(ENERGY / 'md/3hfm/mmgbsa_perres.csv')
    ml = pd.read_csv(run / 'comparison.csv')
    ml['ag_res_idx'] = ml.mutation.str[2:-1].astype(int)
    fig, axes = plt.subplots(2, 2, figsize=(11, 8.5), sharex=True)
    methods = [
        (f'StaB-ddG · {len(ml)} single mutants', ml, 'predicted_kcal_mol', 1,
         'Predicted binding ΔΔG (kcal/mol)'),
        ('Raw interaction energy · 201 frames', ie, 'ab_total', -1,
         'Favorable interaction energy, −E (kJ/mol)'),
        ('Lennard-Jones only · 201 frames', ie, 'ab_lj', -1,
         'Favorable LJ energy, −E (kJ/mol)'),
        ('MM-GBSA · 20 frames', gb, 'mmgbsa_dg', -1,
         'Favorable decomposition, −ΔG (kcal/mol)')]
    report = []
    points = []
    colors = {96: '#0072B2', 73: '#D55E00'}
    for ax, (title, data, column, sign, ylabel) in zip(axes.flat, methods):
        assert data.ag_res_idx.is_unique and exp.resid.is_unique
        joined = data.merge(exp, left_on='ag_res_idx', right_on='resid', validate='one_to_one')
        x, y = joined.ddg_kcal_mol, sign * joined[column]
        ax.scatter(x, y, color='#9aa3ad', s=44, edgecolors='white', zorder=3)
        for residue, color in colors.items():
            selected = joined[joined.ag_res_idx == residue].iloc[0]
            yy = sign * selected[column]
            ax.scatter(selected.ddg_kcal_mol, yy, color=color, s=90, edgecolors='white', zorder=4)
            ax.annotate('K96A' if residue == 96 else 'R73A',
                        (selected.ddg_kcal_mol, yy), xytext=(-6, 10),
                        textcoords='offset points', ha='right' if residue == 96 else 'left',
                        color=color, fontsize=10, weight='bold')
        if column == 'predicted_kcal_mol':
            ax.plot([-0.9, 7.2], [-0.9, 7.2], '--', color='#68727c', lw=1, label='Perfect prediction')
            ax.set_ylim(min(-1, y.min() - 0.5), max(7.4, y.max() + 0.5))
            ax.legend(loc='upper left', frameon=False, fontsize=9)
            rho = float(spearmanr(x, y).correlation) if len(joined) > 2 else None
            note = (f'n = {len(joined)}; Spearman ρ = {rho:+.2f}' if rho is not None
                    else 'n = 2; selected diagnostic pair')
        else:
            rho = float(spearmanr(x, y).correlation)
            note = f'n = {len(joined)}; Spearman ρ = {rho:+.2f}'
            ax.margins(y=0.2)
        ax.set_title(title + '\n' + note, loc='left', fontsize=11)
        ax.set_ylabel(ylabel, fontsize=10)
        ax.set_xlim(-0.9, 7.2)
        ax.axhline(0, color='#d6dbe0', lw=0.8)
        ax.axvline(0, color='#d6dbe0', lw=0.8)
        ax.grid(alpha=0.15)
        ax.spines[['top', 'right']].set_visible(False)
        metric = {'method': column, 'n': len(joined), 'spearman': rho}
        if column == 'predicted_kcal_mol' and len(joined) > 2:
            metric.update({'pearson': float(pearsonr(x, y)[0]),
                           'mae_kcal_mol': float((y - x).abs().mean()),
                           'rmse_kcal_mol': float(((y - x)**2).mean()**0.5)})
        report.append(metric)
        for _, row in joined.iterrows():
            points.append({'method': column, 'resid': int(row.ag_res_idx),
                           'experimental_ddg_kcal_mol': row.ddg_kcal_mol,
                           'plotted_signal': sign * row[column], 'axis': ylabel})
    for ax in axes[1]:
        ax.set_xlabel('Experimental alanine ΔΔG (kcal/mol)\nPositive = mutation weakens binding')
    fig.suptitle('3HFM / HyHEL-10–lysozyme: learned ΔΔG versus simulation signals', fontsize=14)
    fig.text(0.5, 0.015, 'Blue: K96A   Orange: R73A   Gray: other measured antigen residues\n'
             'Simulation panels show sign-flipped wild-type energy contributions, not predicted mutation ΔΔG.\n'
             f'{len(ml)} mutations have ML predictions; training overlap has not been established.',
             ha='center', fontsize=9, color='#48515a')
    fig.tight_layout(rect=(0, 0.095, 1, 0.95))
    for suffix in ('png', 'pdf'):
        fig.savefig(OUT / f'tamarind_vs_simulation_3hfm.{suffix}', dpi=180)
    pd.DataFrame(points).to_csv(run / 'plot_points.csv', index=False)
    (run / 'plot_statistics.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
