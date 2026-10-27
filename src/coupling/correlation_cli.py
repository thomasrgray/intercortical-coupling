"""Batch drivers, audit trail and figures for spike/transition correlations."""

from coupling.array_store import save_arrays
import argparse
import json
from pathlib import Path
import subprocess
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from .phase_coherence import read_session_list, resolve_path, digest_file, write_json
from .correlation import write_csv
from . import mahmood_repository as repository

REFERENCE_SHA = '9c54acb57e227dc9d13597750165b3ebc7be68e9'


def plot(rows, out, kind):
    labels = sorted({(r['stimulus'], r['stratum']) for r in rows})
    fig, axes = plt.subplots(len(labels), 1, figsize=(10, max(3, 2.8*len(labels))), squeeze=False)
    for ax, (stimulus, stratum) in zip(axes[:, 0], labels):
        subset = [r for r in rows if (r['stimulus'], r['stratum']) == (stimulus, stratum) and r['status'] == 'ok']
        if kind == 'spike':
            groups = sorted({(r['window'], r['comparison']) for r in subset})
            for i, group in enumerate(groups):
                values = [r['rho'] for r in subset if (r['window'], r['comparison']) == group]
                ax.scatter(np.full(len(values), i), values, s=8, alpha=.3)
                ax.scatter(i, np.median(values), color='black', marker='_', s=100)
            ax.set_xticks(range(len(groups)), [' / '.join(g) for g in groups], rotation=45, ha='right')
        else:
            for i, r in enumerate(subset):
                ax.plot([i, i], [r['null_low'], r['null_high']], color='.75', linewidth=5)
                ax.scatter(i, r['rho'], c='tab:blue')
            ax.set_xticks(range(len(subset)), [f"{r['gc_transition']}→{r['pc_transition']}" for r in subset])
            ax.set_xlabel('GC → PC transition index (gray: 95% shuffle interval)')
        ax.axhline(0, color='.5', lw=.5)
        ax.set_ylim(-1.05, 1.05)
        ax.set_ylabel('Spearman rho')
        ax.set_title(f'{stimulus} | {stratum}')
    fig.suptitle(f'{kind.capitalize()} correlations — descriptive session summary')
    fig.tight_layout()
    for ext in ('png', 'pdf'):
        fig.savefig(out/f'summary.{ext}', dpi=180)
    plt.close(fig)


def plot_fractions(rows, path):
    groups = sorted({(r.get('session', ''), r['stratum']) for r in rows})
    fig, axes = plt.subplots(len(groups), 1, figsize=(10, max(3.5, 3.5*len(groups))), squeeze=False)
    for ax, (session, stratum) in zip(axes[:, 0], groups):
        selected = [r for r in rows if (r.get('session', ''), r['stratum']) == (session, stratum)]
        for i, row in enumerate(selected):
            ax.plot([i, i], [row['shuffle_low_5'], row['shuffle_high_95']], c='.35', lw=3)
            ax.scatter(i, row['fraction_significant'], c='tab:red' if row['comparison']=='GC-PC' else 'tab:blue', s=20)
        ax.set_xticks(range(len(selected)), [r['stimulus']+' / '+r['comparison'] for r in selected], rotation=45, ha='right')
        ax.axhline(.05, c='.7', ls='--', lw=1)
        ax.set_ylabel('Fraction with Spearman p < .05')
        ax.set_ylim(bottom=0)
        ax.set_title(f'{session} {stratum}: dots=data; gray=5–95% shuffle interval')
    fig.tight_layout()
    for ext in ('png', 'pdf', 'svg'):
        fig.savefig(path.with_suffix('.'+ext), dpi=180)
    plt.close(fig)


def plot_transition_population(rows, path):
    groups = sorted({r['stratum'] for r in rows})
    fig, axes = plt.subplots(1, len(groups), figsize=(5*len(groups), 4), squeeze=False)
    for ax, stratum in zip(axes[0], groups):
        matrix = np.full((3, 3), np.nan)
        for row in rows:
            if row['stratum'] == stratum:
                a, b = row['gc_transition'], row['pc_transition']
                matrix[a, b] = row['fraction']
                if row['p_binomial_bonferroni'] < .05:
                    ax.text(b, a, '*', ha='center', va='center', color='red')
        im = ax.imshow(matrix, vmin=0, vmax=1, cmap='viridis')
        ax.set(xlabel='PC transition', ylabel='GC transition', title=stratum, xticks=range(3), yticks=range(3))
        fig.colorbar(im, ax=ax, label='Fraction of sessions above shuffle 90th percentile')
    fig.tight_layout()
    for ext in ('png', 'pdf', 'svg'):
        fig.savefig(path.with_suffix('.'+ext), dpi=180)
    plt.close(fig)


def route_contrasts(rows, mapping, kind):
    """Descriptive paired differences, never a test treating pairs as animals."""
    result = []
    identity = ('stratum', 'window', 'unit_a', 'unit_b', 'comparison') if kind == 'spike' else ('stratum', 'gc_transition', 'pc_transition', 'comparison')
    for odor, labels in mapping.items():
        retro, ortho = labels
        lookup = {tuple(r[k] for k in identity): r for r in rows if r['stimulus'] == retro and r['status'] == 'ok'}
        for row in rows:
            key = tuple(row[k] for k in identity)
            if row['stimulus'] == ortho and row['status'] == 'ok' and key in lookup:
                base = lookup[key]
                result.append(dict(odor=odor, retro=retro, ortho=ortho,
                                   **{k: row[k] for k in identity},
                                   rho_retro=base['rho'], rho_ortho=row['rho'],
                                   delta_ortho_minus_retro=row['rho']-base['rho'],
                                   n_retro=base['n'], n_ortho=row['n']))
    return result


def run(kind, list_path, config_path, out, repeats=None, seed=0, check=False, profile='repository'):
    if profile not in ('repository', 'adapted'):
        raise ValueError('Unknown analysis profile')
    repeats = repeats if repeats is not None else (5000 if kind == 'transition' and profile == 'repository' else 10000)
    from . import spike_coupling, transition_coupling
    module = spike_coupling if kind == 'spike' else transition_coupling
    if repeats < 1 or seed < 0:
        raise ValueError('shuffles must be positive and seed nonnegative')
    cfg = json.loads(Path(config_path).read_text())
    sessions = read_session_list(list_path)
    overrides = {resolve_path(k, Path(config_path).resolve().parent): v for k, v in cfg.get('sessions', {}).items()}
    if set(overrides)-set(sessions):
        raise ValueError('Config contains sessions absent from session list')
    if len({s.name for s in sessions}) != len(sessions):
        raise ValueError('Session basenames must be unique for output paths')
    out = Path(out)
    out.mkdir(parents=True, exist_ok=False)
    write_json(out/'config.json', dict(**cfg, resolved_profile=profile, resolved_shuffles=repeats))
    reports, resolved = [], []
    for session in sessions:
        spec = {**cfg.get('defaults', {}), **overrides.get(session, {})}
        if profile == 'repository':
            spec = dict(spec, analysis_profile=profile, min_trials=3)
            if kind == 'spike':
                spec.update(windows_s={'evoked': [0, 2]}, detrend='first_difference', single_units_only=False)
                spec.pop('rolling_s', None)
        else:
            spec = dict(spec, analysis_profile=profile)
        report = dict(session=str(session), ready=False, profile=profile, resolved_config=spec)
        try:
            if not spec.get('preprocessing') or 'REPLACE_ME' in str(spec['preprocessing']):
                raise ValueError('Declare preprocessing, sorting/QC, alignment and trial order in config')
            if int(spec.get('min_trials', 8)) < 3:
                raise ValueError('min_trials must be at least three')
            data = module.load(session, spec)
            for odor, labels in spec.get('route_pairs', {}).items():
                if len(labels) != 2 or labels[0] == labels[1] or not set(labels).issubset(set(data['stimulus'])):
                    raise ValueError(f'Invalid route_pairs for {odor}; need two present labels [retro, ortho]')
            if kind == 'spike':
                spike_coupling.windows(spec, data['time_s'])
                if profile == 'repository' and not np.allclose(np.diff(data['time_s']), .001):
                    raise ValueError('Repository spike profile requires 1 ms bins')
                if spec.get('detrend', 'first_difference') not in ('first_difference', 'none'):
                    raise ValueError('Unknown detrend method')
            elif profile == 'repository' and (any(data[r]['mode_ms'].shape[1] != 3 for r in ('GC', 'PC')) or data['GC']['metadata']['bin_width_ms'] != data['PC']['metadata']['bin_width_ms']):
                raise ValueError('Repository transition profile requires matching-bin four-state regional fits')
            elif int(spec.get('posterior_draws', 200)) < 1:
                raise ValueError('posterior_draws must be positive')
            report.update(ready=True, trials=len(data['trial_id']), accepted_trials=int(data['good_trial'].sum()),
                          sources=[str(p) for p in data['sources']])
            resolved.append((session, spec))
            del data
        except (ValueError, KeyError, OSError, ImportError, TypeError) as error:
            report['error'] = str(error)
        reports.append(report)
        print(f"{session.name}: {'ready' if report['ready'] else report['error']}", flush=True)
    write_json(out/'preflight.json', reports)
    if not all(r['ready'] for r in reports):
        raise ValueError(f'Preflight failed; see {out}/preflight.json')
    if check:
        return
    all_rows, all_fractions, all_intra = [], [], []
    status = {'state': 'running', 'completed_sessions': []}
    write_json(out/'status.json', status)
    try:
        for index, (session, spec) in enumerate(resolved):
            folder = out/session.name/('fig08_spike_coupling' if kind == 'spike' else 'fig07_transition_coupling')
            folder.mkdir(parents=True)
            print(f'{index+1}/{len(resolved)}: {session.name} | profile={profile} | shuffles={repeats}', flush=True)
            data = module.load(session, spec)
            retention = dict(trial_id=data['trial_id'].astype(str).tolist(),
                             good_trial=data['good_trial'].tolist(),
                             stimulus=data['stimulus'].astype(str).tolist(),
                             stratum=data['stratum'].astype(str).tolist())
            if kind == 'spike':
                retention.update(unit_id=data['unit_id'].astype(str).tolist(),
                                 region=data['region'].astype(str).tolist(),
                                 good_unit=data['good_unit'].tolist())
            else:
                retention.update(models={r: data[r]['metadata'] for r in ('GC', 'PC')})
            write_json(folder/'retention.json', retention)
            session_seed = seed if profile == 'repository' else seed+index
            if profile == 'repository':
                engine = repository.spike_analyze if kind == 'spike' else repository.transition_analyze
                rows, nulls = engine(data, spec, repeats, session_seed)
                if kind == 'spike':
                    fractions, fraction_nulls = repository.spike_fractions(rows, nulls, session_seed)
                    write_csv(folder/'significant_fractions.csv', fractions)
                    save_arrays(folder/'fraction_nulls.h5', **fraction_nulls)
                    plot_fractions(fractions, folder/'significant_fractions')
                    intra = repository.intra_region_pooled(rows)
                    if intra:
                        write_csv(folder/'intra_region_pooled.csv', intra)
                        all_intra.extend(dict(session=session.name, animal_id=spec.get('animal_id', 'unknown'), **r) for r in intra)
                    all_fractions.extend(dict(session=session.name, animal_id=spec.get('animal_id', 'unknown'), **r) for r in fractions)
                else:
                    setup, moving, setup_nulls = repository.transition_setup(data, session_seed)
                    write_csv(folder/'matched_transitions.csv', setup)
                    if moving:
                        write_csv(folder/'uncertainty_ranked_windows.csv', moving)
                    save_arrays(folder/'matched_nulls.h5', **setup_nulls)
            else:
                rows, nulls = module.analyze(data, spec, repeats, session_seed)
            write_csv(folder/'correlations.csv', rows)
            contrasts = route_contrasts(rows, spec.get('route_pairs', {}), kind)
            if contrasts:
                write_csv(folder/'route_contrasts.csv', contrasts)
            save_arrays(folder/'nulls.h5', **nulls)
            if kind == 'transition':
                save_arrays(folder/'posterior_summary.h5', trial_id=data['trial_id'],
                                    good_trial=data['good_trial'], stimulus=data['stimulus'], stratum=data['stratum'],
                                    gc_mode_ms=data['GC']['mode_ms'], pc_mode_ms=data['PC']['mode_ms'],
                                    gc_variance_ms2=data['GC']['variance_ms2'], pc_variance_ms2=data['PC']['variance_ms2'])
            revision = subprocess.run(['git', 'rev-parse', 'HEAD'], cwd=Path(__file__).resolve().parent,
                                      capture_output=True, text=True).stdout.strip() or 'unavailable'
            write_json(folder/'provenance.json', dict(config=spec, source_revision=revision,
                reference_revision=REFERENCE_SHA, seed=session_seed, shuffles=repeats,
                analysis_profile=profile, matched_transition_shuffles=10000 if kind == 'transition' and profile == 'repository' else None,
                inputs=[dict(path=str(p), sha256=digest_file(p)) for p in data['sources']]))
            plot(rows, folder, kind)
            all_rows.extend(dict(session=session.name, animal_id=spec.get('animal_id', 'unknown'), **r) for r in rows)
            status['completed_sessions'].append(session.name)
            write_json(out/'status.json', status)
        write_csv(out/'all_sessions.csv', all_rows)
        if all_fractions:
            write_csv(out/'all_session_significant_fractions.csv', all_fractions)
            plot_fractions(all_fractions, out/'all_session_significant_fractions')
        if all_intra:
            write_csv(out/'all_session_intra_region_pooled.csv', all_intra)
        if kind == 'transition' and profile == 'repository':
            population = repository.transition_population(all_rows)
            write_csv(out/'transition_population.csv', population)
            plot_transition_population(population, out/'transition_population')
        status['state'] = 'complete'
    except Exception as error:
        status.update(state='failed', error=str(error))
        raise
    finally:
        write_json(out/'status.json', status)


def main(kind):
    parser = argparse.ArgumentParser(description=f'Mahmood-style {kind} correlations; read-only source recordings')
    parser.add_argument('session_list')
    parser.add_argument('config')
    parser.add_argument('output')
    parser.add_argument('--check', action='store_true')
    parser.add_argument('--profile', choices=['repository', 'adapted'], default='repository')
    parser.add_argument('--shuffles', type=int, default=None, help='default: spike 10000; repository transition all-to-all 5000')
    parser.add_argument('--seed', type=int, default=0)
    args = parser.parse_args()
    run(kind, args.session_list, args.config, args.output, args.shuffles, args.seed, args.check, args.profile)
