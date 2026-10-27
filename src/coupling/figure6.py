"""Export compact Figure 6 panels from completed analyses, without recomputation."""

from coupling.array_store import load_arrays, save_arrays
import argparse
import hashlib
import json
from pathlib import Path
import tempfile
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

DEFAULT_OUTPUT = Path('figures/fig06_lfp_coupling/exports')
STYLE = {'font.size': 7, 'axes.titlesize': 8, 'axes.labelsize': 7,
         'axes.spines.top': False, 'axes.spines.right': False,
         'axes.linewidth': .6, 'lines.linewidth': 1,
         'svg.fonttype': 'none', 'pdf.fonttype': 42}
COLORS = {'ortho': '#0072B2', 'retro': '#D55E00'}
PAIRS = (('BA', 'oBA', 'BA'), ('IA', 'oIA', 'IA'), ('route', 'ortho', 'retro'))


def read_arrays(path):
    with load_arrays(path, allow_pickle=False) as data:
        return {key: data[key] for key in data.files}


def grid(data):
    t, f = data['time_s'], data['frequency_hz']
    for name, v in [('time', t), ('frequency', f)]:
        if v.ndim != 1 or len(v) < 2 or not np.isfinite(v).all() or not (np.diff(v) > 0).all():
            raise ValueError(f'{name} grid must be finite and strictly increasing')
    return t, f


def check_pair(a, b):
    for key in ('time_s', 'frequency_hz', 'band', 'band_limits_hz', 'selected_electrode_id', 'alignment_label'):
        if not np.array_equal(a[key], b[key]):
            raise ValueError(f'Cannot compare conditions with different {key}')
    if a['bootstrap_band'].shape != b['bootstrap_band'].shape:
        raise ValueError('Cannot compare different bootstrap shapes')


def save(fig, folder, name):
    try:
        for ext in ('pdf', 'svg', 'png'):
            fig.savefig(folder/f'{name}.{ext}', dpi=300, bbox_inches='tight')
    finally:
        plt.close(fig)


def heatmap(z, t, f, folder, name, title, label, xlabel, limits, cmap):
    fig, ax = plt.subplots(figsize=(2.5, 2), layout='constrained')
    im = ax.pcolormesh(t, f, z, shading='auto', cmap=cmap,
                       vmin=limits[0], vmax=limits[1], rasterized=True)
    ax.axvline(0, color='black', lw=.6, ls='--')
    ax.set(title=title, xlabel=xlabel, ylabel='Frequency (Hz)')
    fig.colorbar(im, ax=ax, label=label, pad=.03)
    save(fig, folder, name)


def phase_panels(source, folder, sources):
    metadata = json.loads((source/'conditions.json').read_text())
    sources.append(source/'conditions.json')
    if (source/'provenance.json').exists():
        sources.append(source/'provenance.json')
    results = {}
    for label in ('oBA', 'BA', 'oIA', 'IA', 'ortho', 'retro'):
        row = metadata.get(label, {})
        if row.get('status') != 'completed':
            continue
        path = (source/row['directory']/'phase_coherence.h5').resolve()
        if source not in path.parents:
            raise ValueError('Condition directory escapes the session output')
        r = read_arrays(path)
        t, f = grid(r)
        draws = r['bootstrap_band']
        if draws.ndim != 3 or draws.shape[0] < 2 or draws.shape[1:] != (len(r['band']), len(t)):
            raise ValueError(f'Invalid bootstrap dimensions: {label}')
        if not np.isfinite(draws).all() or (draws < 0).any() or (draws > 1+1e-12).any():
            raise ValueError(f'Invalid phase coherence: {label}')
        if r['coherence'].shape != (len(f), len(t)) or not np.isfinite(r['coherence']).all():
            raise ValueError(f'Invalid coherence map: {label}')
        results[label] = r
        sources.append(path)
    available = [(name, a, b) for name, a, b in PAIRS if a in results and b in results]
    if not available:
        raise ValueError('No completed matched ortho/retro conditions; run coupling-phase-figure first')
    for name, a, b in available:
        ra, rb = results[a], results[b]
        check_pair(ra, rb)
        t, f = grid(ra)
        xlabel = str(ra['alignment_label'].item())
        for label, r, route in ((a, ra, 'ortho'), (b, rb, 'retro')):
            heatmap(r['coherence'], t, f, folder, f'phase_{name}_{route}',
                    f'{name} · {route}', 'Phase coherence', xlabel, (0, 1), 'viridis')
        heatmap(ra['coherence']-rb['coherence'], t, f, folder, f'phase_{name}_difference',
                f'{name} · ortho − retro', 'Δ phase coherence', xlabel, (-1, 1), 'RdBu_r')
        values = {'time_s': t, 'frequency_hz': f, 'band': ra['band'],
                  'ortho_coherence': ra['coherence'], 'retro_coherence': rb['coherence']}
        for i, band in enumerate(ra['band']):
            fig, ax = plt.subplots(figsize=(2.6, 1.9), layout='constrained')
            for route, r in (('ortho', ra), ('retro', rb)):
                draws = r['bootstrap_band'][:, i]
                lo, hi = np.percentile(draws, [2.5, 97.5], axis=0)
                mean = draws.mean(axis=0)
                ax.plot(t, mean, color=COLORS[route], label=route.capitalize())
                ax.fill_between(t, lo, hi, color=COLORS[route], alpha=.18, linewidth=0)
                values[f'{route}_band_{i}_mean_lo_hi'] = np.stack([mean, lo, hi])
            ax.axvline(0, color='black', lw=.6, ls='--')
            ax.set(title=f'{name} · {band}', xlabel=xlabel, ylabel='Phase coherence', ylim=(0, 1))
            ax.legend(frameon=False, fontsize=6)
            save(fig, folder, f'phase_{name}_band_{i}')
        save_arrays(folder/f'phase_{name}_values.h5', **values)
    return {'matched_comparisons': [x[0] for x in available],
            'omitted_comparisons': [x[0] for x in PAIRS if x not in available],
            'interval': 'Pointwise 95% within-session bootstrap percentile interval; descriptive, not a route significance test',
            'route_weighting': 'Equal within-odor coherence means; no pooling of odor phase vectors'}


def power_panels(source, folder, sources):
    r = read_arrays(source)
    t, f = grid(r)
    regions, stimuli = r['region'].tolist(), r['stimulus'].tolist()
    z = r['z_power']
    if z.shape != (len(regions), len(stimuli), len(f), len(t)):
        raise ValueError('Power axes must be region × stimulus × frequency × time')
    if set(regions) != {'GC', 'PC'} or len(regions) != 2 or len(set(stimuli)) != len(stimuli):
        raise ValueError('Require unique stimuli and GC/PC power regions')
    if not {'oBA', 'BA', 'oIA', 'IA'}.issubset(stimuli):
        raise ValueError('Missing matched power stimuli')
    if np.isinf(z).any():
        raise ValueError('Infinite power values')
    sources.append(source)
    provenance = source.parent/'provenance.json'
    if not provenance.exists():
        raise ValueError('Power provenance.json required to label event alignment')
    sources.append(provenance)
    p = json.loads(provenance.read_text())
    event = p['preprocessing']['event_alignment']
    xlabel = f'Time from {event} (s)'
    for i, region in enumerate(regions):
        for name, a, b in PAIRS[:2]:
            for route, label in [('ortho', a), ('retro', b)]:
                panel = z[i, stimuli.index(label)]
                if not np.isfinite(panel).any():
                    raise ValueError(f'No finite power for {region}/{label}')
                heatmap(panel, t, f, folder, f'power_{region}_{name}_{route}',
                        f'{region} · {name} · {route}', 'Power (baseline z)', xlabel, (-3, 3), 'RdBu_r')
    save_arrays(folder/'power_values.h5', **r)
    return {'normalization': 'Per-trial prestimulus power z-score; equal trials within channel, then equal channels',
            'color_limits': [-3, 3], 'missing_values': 'NaN remains blank; no zero imputation', 'alignment': event}


def export(output, phase_session=None, power=None):
    """Stage all panels, then publish a complete folder; never overwrite a run."""
    output = Path(output).resolve()
    if output.exists():
        raise FileExistsError(output)
    if phase_session is None and power is None:
        raise ValueError('Supply a phase session and/or session_power.h5')
    output.parent.mkdir(parents=True, exist_ok=True)
    sources, metadata = [], {'schema_version': 1, 'scope': 'Single-session descriptive panels',
                            'style': 'Compact editable panels inspired by Mahmood 2023; not journal certification'}
    with tempfile.TemporaryDirectory(prefix='.figure6-', dir=output.parent) as tmp:
        stage = Path(tmp)/'panels'; stage.mkdir()
        with plt.rc_context(STYLE):
            if phase_session is not None:
                metadata['phase'] = phase_panels(Path(phase_session).resolve(), stage, sources)
            if power is not None:
                metadata['power'] = power_panels(Path(power).resolve(), stage, sources)
        metadata['sources'] = []
        for path in sources:
            digest = hashlib.sha256()
            with path.open('rb') as handle:
                for chunk in iter(lambda: handle.read(1024*1024), b''):
                    digest.update(chunk)
            metadata['sources'].append({'path': str(path), 'sha256': digest.hexdigest()})
        metadata['files'] = sorted(p.name for p in stage.iterdir())
        (stage/'manifest.json').write_text(json.dumps(metadata, indent=2))
        if output.exists():
            raise FileExistsError(output)
        stage.rename(output)
    return output


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--phase-session', type=Path, help='One completed session folder containing conditions.json')
    parser.add_argument('--power', type=Path, help='Completed session_power.h5 (with sibling provenance.json)')
    parser.add_argument('--output', type=Path, default=DEFAULT_OUTPUT, help='New export directory; never overwritten')
    args = parser.parse_args()
    print(export(args.output, args.phase_session, args.power))


if __name__ == '__main__':
    main()
