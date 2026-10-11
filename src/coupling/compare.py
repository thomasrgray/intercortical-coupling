"""Publication exports of two saved session power panels; descriptive comparisons."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def select_panel(data, selector):
    """Resolve REGION:STIMULUS against labels rather than hard-coded array order."""
    parts = selector.split(':')
    if len(parts) != 2:
        raise ValueError('Panel must be REGION:STIMULUS, e.g. GC:IA')
    indices = []
    for key, label in zip(('region', 'stimulus'), parts):
        matches = np.flatnonzero(data[key] == label)
        if len(matches) != 1:
            raise ValueError(f'Unknown or duplicate {key}: {label}')
        indices.append(int(matches[0]))
    idx = tuple(indices)
    return data['z_power'][idx], data['n_channels'][idx]


def compare(source, output, left='GC:IA', right='GC:oIA', clim=3.,
            difference=False, difference_clim=1., title=None):
    source, output = Path(source), Path(output)
    if output.exists():
        raise FileExistsError(output)
    for limit in (clim, difference_clim):
        if not np.isfinite(limit) or limit <= 0:
            raise ValueError('Color limits must be positive and finite')
    with np.load(source, allow_pickle=False) as archive:
        required = {'z_power', 'n_channels', 'frequency_hz', 'time_s', 'region',
                    'stimulus', 'session_id', 'animal_id', 'settings_json'}
        if not required.issubset(archive.files):
            raise ValueError(f'Missing fields: {sorted(required-set(archive.files))}')
        data = {key: archive[key] for key in required}
    f, t = data['frequency_hz'], data['time_s']
    expected = (len(data['region']), len(data['stimulus']), len(f), len(t))
    if data['z_power'].shape != expected or data['n_channels'].shape != expected:
        raise ValueError('Power/count axes must be region, stimulus, frequency, time')
    if (len(f) < 2 or len(t) < 2 or not np.isfinite(f).all()
            or not np.isfinite(t).all() or np.any(np.diff(f) <= 0)
            or np.any(np.diff(t) <= 0)):
        raise ValueError('Require finite increasing frequency and time grids')
    a, ca = select_panel(data, left)
    b, cb = select_panel(data, right)
    if not np.isfinite(a).any() or not np.isfinite(b).any():
        raise ValueError('Selected panel has no finite power')
    delta = a-b
    settings = json.loads(str(data['settings_json'].item()))
    output.mkdir(parents=True)
    panels = [(left, a, ca), (right, b, cb)]
    if difference:
        panels.append((f'{left} − {right}', delta, None))
    style = {'font.size': 9, 'axes.titlesize': 10, 'axes.labelsize': 9,
             'pdf.fonttype': 42, 'ps.fonttype': 42, 'svg.fonttype': 'none'}
    with plt.rc_context(style):
        fig, axes = plt.subplots(1, len(panels), figsize=(7.2 if not difference else 10.5, 3.5),
                                 sharex=True, sharey=True, layout='constrained')
        meshes = []
        for i, (ax, (label, z, counts)) in enumerate(zip(axes, panels)):
            is_delta = i == 2
            limit = difference_clim if is_delta else clim
            mesh = ax.pcolormesh(t*1000, f, np.ma.masked_invalid(z), shading='auto',
                                 cmap='RdBu_r' if is_delta else 'turbo',
                                 vmin=-limit, vmax=limit, rasterized=True)
            meshes.append(mesh)
            suffix = ''
            if counts is not None:
                valid = counts[np.isfinite(z)]
                low, high = int(valid.min()), int(valid.max())
                count_label = str(low) if low == high else f'{low}–{high}'
                suffix = f'\n{count_label} contributing channels'
            ax.set_title(label.replace('PC:', 'pPC:').replace(':', ' · ') + suffix)
            ax.axvline(0, color='black', lw=.8, ls='--')
            ax.set(xlim=(-1000, 1750), ylim=(0, 70), xlabel='Time from TTL onset (ms)')
            ax.set_xticks([-1000, 0, 1000, 1750])
            ax.text(0, 1.22, chr(65+i), transform=ax.transAxes, weight='bold', fontsize=12)
        axes[0].set_ylabel('Frequency (Hz)')
        fig.colorbar(meshes[0], ax=list(axes[:2]), label='Baseline-normalized power (z)', shrink=.85)
        if difference:
            fig.colorbar(meshes[2], ax=axes[2], label='Difference in normalized power (Δz)', shrink=.85)
        fig.suptitle(title or f'{data["animal_id"].item()} · single-session power', fontsize=11)
        for extension in ('png', 'pdf', 'svg'):
            fig.savefig(output/f'comparison.{extension}', dpi=600)
        plt.close(fig)
    np.savez_compressed(output/'comparison_values.npz', left_z=a, right_z=b,
                        difference_z=delta, frequency_hz=f, time_s=t,
                        left_n_channels=ca, right_n_channels=cb)
    metadata = dict(source=str(source.resolve()), input_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
                    left=left, right=right, difference_direction='left minus right',
                    color_limit=clim, difference_color_limit=difference_clim,
                    session_id=str(data['session_id'].item()), animal_id=str(data['animal_id'].item()),
                    settings=settings, statistical_inference=False)
    (output/'comparison.json').write_text(json.dumps(metadata, indent=2))
    caption = (f'Single-session LFP power from {metadata["animal_id"]} ({metadata["session_id"]}): '
               f'{left} and {right}. Power was normalized per trial/channel/frequency against '
               f'the {settings.get("baseline_s")} s baseline before averaging trials and channels. '
               'Both panels use the same color scale. Time zero marks TTL onset, not measured sensory arrival. '
               'Channel counts report contributors, not independent animals. ')
    if difference:
        caption += f'The difference panel shows {left} minus {right} in normalized power units. '
    caption += 'Descriptive sample from one animal; no significance test or coupling estimate.\n'
    (output/'caption.txt').write_text(caption)
    return output


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('session_power_npz')
    parser.add_argument('output_directory')
    parser.add_argument('--left', default='GC:IA')
    parser.add_argument('--right', default='GC:oIA')
    parser.add_argument('--clim', type=float, default=3.)
    parser.add_argument('--difference', action='store_true')
    parser.add_argument('--difference-clim', type=float, default=1.)
    parser.add_argument('--title')
    args = parser.parse_args()
    print(compare(args.session_power_npz, args.output_directory, args.left, args.right,
                  args.clim, args.difference, args.difference_clim, args.title))


if __name__ == '__main__':
    main()
