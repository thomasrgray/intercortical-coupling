"""Mahmood 2023 phase coherence with read-only inputs and portable output paths.

Independent implementation of the estimator in updated/lfp_coherence_aggregate.py.
See docs/mahmood_phase_coherence.md for paper/code differences and input contracts.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess

import numpy as np
from scipy.signal import detrend, stft
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

REFERENCE_SHA = '9c54acb57e227dc9d13597750165b3ebc7be68e9'
BANDS = {'repository': [('Theta', 3, 8), ('Alpha', 8, 13), ('Beta', 13, 30),
                        ('Low gamma', 30, 70), ('High gamma', 70, 100)],
         'paper': [('Delta', 1, 4), ('Theta', 4, 8), ('Mu', 8, 13),
                   ('Beta', 12, 30), ('Gamma', 30, 100)]}


def write_json(path, value):
    path.write_text(json.dumps(value, indent=2) + '\n')


def resolve_path(value, base):
    path = Path(os.path.expandvars(str(value))).expanduser()
    return (base / path).resolve()


def read_session_list(path):
    path = Path(path).resolve()
    rows = [line.strip() for line in path.read_text().splitlines()
            if line.strip() and not line.lstrip().startswith('#')]
    paths = [resolve_path(row, path.parent) for row in rows]
    if not paths or len(set(paths)) != len(paths):
        raise ValueError('Session list is empty or contains duplicate resolved paths')
    for p in paths:
        if not p.is_dir():
            raise ValueError(f'Session directory does not exist: {p}')
    return paths


def only_file(folder, patterns):
    files = sorted({p for pattern in patterns for p in folder.glob(pattern)})
    if len(files) != 1:
        raise ValueError(f'{folder}: expected one {patterns}, found {len(files)}; set an explicit path')
    return files[0]


def flatten_ids(value):
    if isinstance(value, list):
        return [i for group in value for i in flatten_ids(group)]
    return [int(value)]


def load_h5(source, session, spec):
    import h5py
    info_path = resolve_path(spec['info'], session) if 'info' in spec else only_file(session, ['*.info'])
    info = json.loads(info_path.read_text())
    layout = info['electrode_layout']
    aliases = spec.get('region_keys', {'GC': 'gc', 'PC': 'pc'})
    mapping = {region: flatten_ids(layout[key]) for region, key in aliases.items()}
    if set(mapping) != {'GC', 'PC'}:
        raise ValueError('region_keys must map GC and PC to exact .info layout keys')
    if set(mapping['GC']) & set(mapping['PC']):
        raise ValueError('GC and PC physical electrode mappings overlap')
    fs, onset = float(spec['fs_hz']), float(spec['onset_s'])
    dins = spec['stimulus_by_din']
    if not dins or len(set(dins.values())) != len(dins):
        raise ValueError('Provide unique stimulus labels for the DINs to include')
    arrays, labels, event_ids = [], [], []
    with h5py.File(source, 'r') as h:
        if 'Parsed_LFP' not in h or 'Parsed_LFP_channels' not in h:
            raise ValueError('Missing /Parsed_LFP or /Parsed_LFP_channels; extract LFPs with blech_clust first (raw DAT is not a parsed session)')
        electrodes = np.asarray(h['Parsed_LFP_channels']).reshape(-1)
        nodes = {}
        for name in h['Parsed_LFP']:
            match = re.fullmatch(r'dig_in_(\d+)_LFPs', name)
            if match:
                number = int(match[1])
                if number in nodes:
                    raise ValueError(f'Duplicate parsed DIN {number}')
                nodes[number] = name
        for din, label in sorted(dins.items(), key=lambda item: int(item[0])):
            if int(din) not in nodes:
                raise ValueError(f'Missing /Parsed_LFP/dig_in_{din}_LFPs')
            a = np.asarray(h['Parsed_LFP'][nodes[int(din)]], dtype=float)
            if a.ndim != 3 or a.shape[0] != len(electrodes):
                raise ValueError('Parsed LFP axes must be channel, trial, sample')
            arrays.append(a)
            labels.extend([label] * a.shape[1])
            event_ids.extend([f'DIN{din}:trial{i}' for i in range(a.shape[1])])
        x = np.concatenate(arrays, axis=1)
        # Follow the paper's channel eligibility when requested; never infer units from row order.
        eligible = spec.get('eligible_electrode_ids')
        if eligible is None:
            if 'unit_descriptor' not in h:
                raise ValueError('Need /unit_descriptor or explicit eligible_electrode_ids from reviewed channel QC')
            units = np.asarray(h['unit_descriptor'])
            if not units.dtype.names or 'electrode_number' not in units.dtype.names:
                raise ValueError('unit_descriptor lacks electrode_number')
            if 'single_unit' in units.dtype.names:
                units = units[units['single_unit'].astype(bool)]
            eligible = np.unique(units['electrode_number']).tolist()
    region = []
    for electrode in electrodes:
        hits = [name for name, ids in mapping.items() if electrode in ids]
        if len(hits) != 1:
            raise ValueError(f'Parsed physical electrode {electrode} has no unique GC/PC mapping')
        region.append(hits[0])
    provenance = spec.get('preprocessing')
    if not provenance or 'REPLACE_ME' in str(provenance):
        raise ValueError('Declare preprocessing (filtering, reference, extraction and QC) in config')
    return dict(x=x, time_s=np.arange(x.shape[-1])/fs-onset,
                electrode_id=electrodes, region=np.array(region), stimulus=np.array(labels),
                event_id=np.array(event_ids), good_channel=np.isin(electrodes, eligible),
                good_trial=np.ones(x.shape[:2], dtype=bool),
                provenance={'declared_preprocessing': provenance, 'info_path': str(info_path),
                            'info_sha256': digest_file(info_path), 'eligible_electrode_ids': eligible,
                            'fs_hz': fs, 'onset_s': onset, 'stimulus_by_din': dins})


def load_npz(source):
    from .power import validate
    with np.load(source, allow_pickle=False) as d:
        provenance = validate(d)
        result = {key: d[key] for key in ('time_s', 'electrode_id', 'region', 'stimulus', 'good_channel', 'good_trial')}
        result.update(x=d['lfp_uv'], event_id=d['event_id'] if 'event_id' in d.files else np.arange(d['lfp_uv'].shape[1]),
                      provenance=provenance)
    return result


def load_input(session, spec):
    source = resolve_path(spec['input'], session) if spec.get('input') else only_file(session, ['*.h5', '*.hdf5'])
    if source.suffix == '.npz':
        data = load_npz(source)
    else:
        data = load_h5(source, session, spec)
    x, t = data['x'], data['time_s']
    if len(t) < 500 or not np.allclose(np.diff(t), .001, atol=1e-9, rtol=0):
        raise ValueError('Require verified 1000 Hz LFPs; no implicit resampling')
    if t[0] > -1 or t[-1] < 2.25:
        raise ValueError('Need at least -1 to +2.25 s relative to verified stimulus onset')
    if len(np.unique(data['electrode_id'])) != len(data['electrode_id']):
        raise ValueError('Duplicate physical electrode IDs')
    channels = np.flatnonzero(data['good_channel'])
    if not len(channels):
        raise ValueError('No eligible channels')
    # Use the same events for representative-channel selection and both regions.
    good = np.all(data['good_trial'][channels], axis=0)
    if good.sum() < 3:
        raise ValueError('Fewer than three jointly accepted events across eligible channels')
    accepted = x[channels][:, good]
    if not np.isfinite(accepted).all() or np.any(np.ptp(accepted, axis=-1) == 0):
        raise ValueError('Accepted LFP contains nonfinite or constant traces; review upstream QC')
    counts = {str(label): int(np.sum(good & (data['stimulus'] == label))) for label in np.unique(data['stimulus'])}
    if any(n == 0 for n in counts.values()):
        raise ValueError('At least one selected stimulus has no jointly accepted events')
    groups = [np.flatnonzero(data['region'][channels] == r) for r in ('GC', 'PC')]
    if any(len(g) < 2 for g in groups):
        raise ValueError('Need at least two eligible channels per region for intra-region control')
    data.update(source=source, channels=channels, good=good, groups=groups, counts=counts)
    return data


def transform(data, profile, seed=0):
    """Two STFT passes select channels without retaining a full 64-channel phase tensor."""
    def channel_stft(index):
        channel = data['channels'][index]
        x = data['x'][channel, data['good']]
        f, t, z = stft(detrend(x, axis=-1), fs=1000, window='hann', nperseg=500,
                       noverlap=499 if profile == 'repository' else 495,
                       boundary='zeros', padded=True)
        fi = f < 100
        ti = (t >= 0) & (t < min(5., len(x[0])/1000))
        t += data['time_s'][0]
        return z[:, fi][:, :, ti], f[fi], t[ti]

    selected, powers = [], {}
    for group in data['groups']:
        mean = None
        for index in group:
            z, f, t = channel_stft(index)
            p = np.angle(z)
            if mean is None:
                mean = np.zeros_like(p)
            mean += p / len(group)
            powers[int(index)] = np.mean(abs(z)**2, axis=0)
        losses = []
        for index in group:
            z, _, _ = channel_stft(index)
            error = np.angle(z)-mean
            losses.append(np.sum(abs(error) if profile == 'repository' else error**2))
        selected.append(int(group[np.argmin(losses)]))
    rng = np.random.default_rng(seed)
    intra = np.array([rng.choice(group, 2, replace=False) for group in data['groups']])
    used = sorted(set(selected) | set(intra.ravel().tolist()))
    phases = []
    for index in used:
        z, f, t = channel_stft(index)
        phases.append(np.angle(z))
    phase = np.stack(phases, axis=1)
    local_selected = [used.index(i) for i in selected]
    local_intra = np.array([[used.index(int(i)) for i in pair] for pair in intra])
    return phase, np.stack([powers[i] for i in used]), f, t, np.array(used), local_selected, local_intra


def select_channels(phase, groups, profile):
    selected = []
    for group in groups:
        p = phase[:, group]
        error = p - p.mean(axis=1, keepdims=True)
        loss = np.abs(error) if profile == 'repository' else error**2
        selected.append(int(group[np.argmin(loss.sum(axis=(0, 2, 3)))]))
    return selected


def phase_locking(a, b):
    return np.abs(np.exp(1j*(a-b)).mean(axis=0))


def sustained_deviation(mean, bounds, time_s, minimum_s=.1):
    dev = (mean > bounds[1]).astype(np.int8) - (mean < bounds[0]).astype(np.int8)
    # Treat opposite signs separately and retain runs ending at the last sample.
    dt = np.median(np.diff(time_s))
    edges = np.r_[0, np.flatnonzero(np.diff(dev) != 0)+1, len(dev)]
    for start, stop in zip(edges[:-1], edges[1:]):
        if (stop-start)*dt < minimum_s-1e-9:
            dev[start:stop] = 0
    return dev


def calculate(phase, groups, f, t, profile, n_boot, seed, selected=None, intra=None):
    rng = np.random.default_rng(seed)
    if selected is None:
        selected = select_channels(phase, groups, profile)
    if intra is None:
        intra = [rng.choice(group, 2, replace=False) for group in groups]
    a, b = [phase[:, i] for i in selected]
    diff = np.exp(1j*(a-b))
    intra_diff = np.stack([np.exp(1j*(phase[:, pair[0]]-phase[:, pair[1]])) for pair in intra])
    bands = BANDS[profile]
    freq_inds = [(f >= low) & (f < high) for _, low, high in bands]
    if any(not mask.any() for mask in freq_inds):
        raise ValueError('Frequency grid misses a required band')
    boot = np.empty((n_boot, len(bands), len(t)))
    shuffle = np.empty_like(boot)
    within = np.empty((n_boot, 2, len(bands), len(t)))
    n = len(a)
    for k in range(n_boot):
        idx = rng.integers(n, size=n)
        value = np.abs(diff[idx].mean(axis=0))
        boot[k] = [value[mask].mean(axis=0) for mask in freq_inds]
        idx = rng.integers(n, size=n)
        value = np.abs(intra_diff[:, idx].mean(axis=1))
        within[k] = np.stack([value[:, mask].mean(axis=1) for mask in freq_inds], axis=1)
        # Upstream uses replacement sampling, not a permutation without replacement.
        idx = rng.integers(n, size=n)
        value = phase_locking(a[idx], b)
        shuffle[k] = [value[mask].mean(axis=0) for mask in freq_inds]
    baseline = (t > -.75) & (t < -.25)
    if baseline.sum() < 2:
        raise ValueError('Missing -750 to -250 ms baseline')
    bounds = np.stack([np.percentile(boot[:, band, baseline], [2.5, 97.5]) for band in range(len(bands))])
    dev = np.stack([sustained_deviation(mean, ci, t) for mean, ci in zip(boot.mean(axis=0), bounds)])
    return dict(coherence=abs(diff.mean(axis=0)), bootstrap_band=boot, shuffle_band=shuffle,
                intra_band=within, baseline_bounds=bounds, deviation=dev,
                selected_indices=np.array(selected), intra_indices=np.array(intra),
                frequency_hz=f, time_s=t, band=np.array([band[0] for band in bands]),
                band_limits_hz=np.array([band[1:] for band in bands]))


def save_figure(fig, stem):
    fig.tight_layout()
    for ext in ('png', 'pdf'):
        fig.savefig(stem.with_suffix('.'+ext), dpi=160)
    plt.close(fig)


def plot_session(result, powers, folder, title):
    t, f = result['time_s'], result['frequency_hz']
    fig, axes = plt.subplots(len(result['band']), 1, figsize=(10, 12), sharex=True)
    for i, ax in enumerate(axes):
        for key, label in [('bootstrap_band', 'GC–PC'), ('shuffle_band', 'Replacement shuffle')]:
            m, sd = result[key][:, i].mean(axis=0), result[key][:, i].std(axis=0)
            ax.plot(t, m, label=label)
            ax.fill_between(t, m-3*sd, m+3*sd, alpha=.18)
        for region in range(2):
            ax.plot(t, result['intra_band'][:, region, i].mean(axis=0), ls=':', alpha=.6,
                    label=('Within GC', 'Within PC')[region])
        for bound in result['baseline_bounds'][i]:
            ax.axhline(bound, color='red', alpha=.4, lw=.7)
        ax.fill_between(t, 0, 1, where=result['deviation'][i] != 0, color='gold', alpha=.15)
        ax.axvline(0, color='black', lw=.7)
        ax.set(ylabel=str(result['band'][i]), ylim=(0, 1))
    axes[0].legend(ncol=2, fontsize=8)
    axes[0].set_title(title+'\nPhase coherence; shading = ±3 bootstrap SD; gold = baseline deviation')
    axes[-1].set_xlabel('Time from declared alignment event (s)')
    save_figure(fig, folder/'phase_coherence')
    fig, axes = plt.subplots(1, 2, figsize=(12, 4), sharey=True)
    for ax, region, index in zip(axes, ('GC', 'PC'), result['selected_indices']):
        power = powers[index]
        sd = power.std(axis=-1, keepdims=True)
        z = np.divide(power-power.mean(axis=-1, keepdims=True), sd,
                      out=np.zeros_like(power), where=sd>0)
        image = ax.pcolormesh(t, f, z, shading='auto', cmap='RdBu_r', vmin=-3, vmax=3)
        ax.axvline(0, color='black', lw=.7)
        ax.set(title=region, xlabel='Time from declared alignment event (s)', ylabel='Frequency (Hz)')
        fig.colorbar(image, ax=ax, label='Temporal z-score of trial-mean STFT power')
    fig.suptitle(title+'; selected channels, pooled trials (descriptive power context)')
    save_figure(fig, folder/'power_context')


def digest_file(path):
    digest = hashlib.sha256()
    with path.open('rb') as handle:
        for block in iter(lambda: handle.read(1024*1024), b''):
            digest.update(block)
    return digest.hexdigest()


def resolve_manifest(list_path, config_path):
    sessions = read_session_list(list_path)
    cfg = json.loads(Path(config_path).read_text())
    overrides = {resolve_path(key, Path(config_path).resolve().parent): value
                 for key, value in cfg.get('sessions', {}).items()}
    unknown = set(overrides)-set(sessions)
    if unknown:
        raise ValueError(f'Config sessions absent from fin_inter_list: {sorted(map(str, unknown))}')
    manifest = []
    for session in sessions:
        spec = {**cfg.get('defaults', {}), **overrides.get(session, {})}
        manifest.append(dict(session=str(session), output_key=session.name+'-'+hashlib.sha256(str(session).encode()).hexdigest()[:10], spec=spec))
    return manifest


def aggregate(summaries, out):
    if not summaries:
        return
    first = summaries[0]
    for r in summaries[1:]:
        if r['deviation'].shape != first['deviation'].shape or not np.allclose(r['time_s'], first['time_s'], atol=1e-9, rtol=0):
            raise ValueError('Session time grids differ: per-session outputs saved; no invalid aggregate produced')
    dev = np.stack([r['deviation'] for r in summaries])
    means = np.stack([r['mean'] for r in summaries])
    np.savez_compressed(out/'aggregate.npz', deviation=dev, mean_coherence=means,
                        session=np.array([r['session'] for r in summaries]), time_s=first['time_s'], band=first['band'])
    fig, axes = plt.subplots(3, len(first['band']), figsize=(18, 9), squeeze=False)
    for i, name in enumerate(first['band']):
        t = first['time_s']
        dt = np.median(np.diff(t))
        edges = np.r_[t-dt/2, t[-1]+dt/2]
        mesh = axes[0, i].pcolormesh(edges, np.arange(len(summaries)+1)-.5, dev[:, i], shading='flat', cmap='coolwarm', vmin=-1, vmax=1)
        axes[0, i].set_yticks(np.arange(len(summaries)), labels=np.arange(1,len(summaries)+1))
        if i == len(first['band'])-1:
            bar = fig.colorbar(mesh, ax=axes[0, i], ticks=[-1,0,1])
            bar.ax.set_yticklabels(['Decrease','No change','Increase'])
        axes[0, i].set(title=name, ylabel='Session row (see aggregate.npz)')
        axes[1, i].plot(t, (dev[:, i] != 0).mean(axis=0), color='black')
        axes[1, i].set(ylabel='Fraction with baseline deviation', ylim=(-.03, 1.03))
        axes[2, i].plot(t, means[:, i].T, alpha=.3)
        axes[2, i].plot(t, means[:, i].mean(axis=0), color='black', lw=2)
        axes[2, i].set(ylabel='Phase coherence', xlabel='Time from alignment (s)', ylim=(0, 1))
        for ax in axes[:, i]:
            ax.axvline(0, color='black', lw=.6)
    fig.suptitle('All completed sessions; equal session weighting, descriptive (not animal-level inference)')
    save_figure(fig, out/'all_sessions')


def run(list_path, config_path, output, check=False, profile='repository', n_boot=500, seed=0):
    if n_boot < 2 or seed < 0:
        raise ValueError('Require at least two bootstrap draws and a nonnegative seed')
    manifest = resolve_manifest(list_path, config_path)
    output = Path(output).resolve()
    for row in manifest:
        session = Path(row['session'])
        if output == session or session in output.parents:
            raise ValueError('Output must be outside every source session directory')
    output.mkdir(parents=True, exist_ok=False)
    write_json(output/'manifest_resolved.json', manifest)
    # Export a plain absolute-path list accepted by the original scripts.
    (output/'fin_inter_list.txt').write_text(''.join(row['session']+'\n' for row in manifest))
    statuses, summaries = [], []
    for row in manifest:
        session, spec = Path(row['session']), row['spec']
        status = dict(session=str(session), output_key=row['output_key'], animal_id=spec.get('animal_id'), status='checking')
        try:
            data = load_input(session, spec)
            status.update(input=str(data['source']), accepted_trials=data['counts'],
                          eligible_electrodes=data['electrode_id'][data['channels']].tolist(), status='ready')
            if not check:
                folder = output/row['output_key']
                folder.mkdir()
                session_seed = int.from_bytes(hashlib.sha256((str(seed)+str(session)).encode()).digest()[:8], 'little')
                phase, powers, f, t, used, selected, intra = transform(data, profile, session_seed)
                result = calculate(phase, None, f, t, profile, n_boot, session_seed, selected, intra)
                ids = data['electrode_id'][data['channels']][used]
                result.update(trial_mean_stft_power=powers, phase_channel_electrode_id=ids, selected_electrode_id=ids[result['selected_indices']], intra_electrode_id=ids[result['intra_indices']],
                              accepted_event_id=data['event_id'][data['good']], stimulus=data['stimulus'][data['good']])
                np.savez_compressed(folder/'phase_coherence.npz', **result)
                plot_session(result, powers, folder, session.name)
                try:
                    revision = subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=Path(__file__).parent, text=True, stderr=subprocess.DEVNULL).strip()
                except (OSError, subprocess.CalledProcessError):
                    revision = 'unavailable'
                write_json(folder/'provenance.json', dict(source=str(data['source']), input_sha256=digest_file(data['source']),
                    session=str(session), animal_id=spec.get('animal_id'), profile=profile, bootstrap_samples=n_boot, seed=seed, session_seed=session_seed,
                    reference_sha=REFERENCE_SHA, implementation_revision=revision,
                    source_module_sha256=digest_file(Path(__file__)), preprocessing=data['provenance'],
                    settings=dict(window_samples=500, overlap_samples=499 if profile=='repository' else 495,
                                  max_frequency_exclusive=100, detrend='linear', boundary='zeros', padded=True,
                                  selection_time_s=[float(t[0]), float(t[-1])], pooling='all selected stimuli',
                                  qc='intersection across all eligible channels', bands=BANDS[profile]),
                    selected_electrode_id=result['selected_electrode_id'].tolist(), counts=data['counts']))
                summaries.append(dict(session=str(session), deviation=result['deviation'], mean=result['bootstrap_band'].mean(axis=0), time_s=t, band=result['band']))
                status.update(status='completed', output=str(folder))
        except Exception as exc:
            status.update(status='failed', error=f'{type(exc).__name__}: {exc}')
        statuses.append(status)
        write_json(output/'status.json', statuses)
        print(f"{session.name}: {status['status']} {status.get('error', '')}", flush=True)
    if not check:
        aggregate(summaries, output)
    return statuses


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('session_list')
    parser.add_argument('config')
    parser.add_argument('output')
    parser.add_argument('--check', action='store_true', help='Validate inputs and mapping without STFT/bootstrap')
    parser.add_argument('--profile', choices=tuple(BANDS), default='repository')
    parser.add_argument('--bootstrap', type=int, default=500)
    parser.add_argument('--seed', type=int, default=0)
    args = parser.parse_args()
    statuses = run(args.session_list, args.config, args.output, args.check, args.profile, args.bootstrap, args.seed)
    if any(row['status']=='failed' for row in statuses):
        raise SystemExit(1)


if __name__ == '__main__':
    main()
