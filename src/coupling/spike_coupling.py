"""Trial-differenced Spearman spike noise correlations (Mahmood-style)."""
import itertools
import json
import numpy as np
from .phase_coherence import only_file, resolve_path, flatten_ids
from .correlation import validate_trials, conditions, condition_indices, correlation_test, correct, progress


def load(session, spec):
    source = resolve_path(spec['input'], session) if spec.get('input') else only_file(session, ['*.h5', '*.hdf5'])
    if source.suffix == '.npz':
        with np.load(source, allow_pickle=False) as archive:
            data = {k: archive[k] for k in ('spikes', 'time_s', 'unit_id', 'region', 'good_unit',
                                           'trial_id', 'stimulus', 'stratum', 'trial_order', 'good_trial')}
        sources = [source]
    else:
        data, sources = load_h5(source, session, spec)
    x = data['spikes']
    if x.ndim != 3:
        raise ValueError('spikes must be trial x unit x sample')
    n, u, samples = x.shape
    validate_trials(data, n)
    if data['good_unit'].shape != (u,) or data['good_unit'].dtype.kind != 'b':
        raise ValueError('good_unit must be a boolean unit vector')
    if data['unit_id'].shape != (u,) or len(np.unique(data['unit_id'])) != u:
        raise ValueError('unit_id must contain unique unit identities')
    if data['region'].shape != (u,) or not set(data['region']).issubset({'GC', 'PC'}):
        raise ValueError('region must label each unit GC or PC')
    if any(not np.any(data['good_unit'] & (data['region'] == r)) for r in ('GC', 'PC')):
        raise ValueError('Need at least one eligible sorted unit in GC and PC')
    t = data['time_s']
    if t.shape != (samples,) or samples < 2 or not np.isfinite(t).all() or np.any(np.diff(t) <= 0) or not np.allclose(np.diff(t), np.diff(t)[0]):
        raise ValueError('time_s must be finite, increasing and uniformly sampled')
    accepted = x[data['good_trial']][:, data['good_unit']]
    if not np.isfinite(accepted).all() or np.any(accepted < 0) or np.any(accepted != np.floor(accepted)):
        raise ValueError('Accepted spikes must contain nonnegative finite integer counts')
    data['sources'] = sources
    return data


def load_h5(source, session, spec):
    import h5py
    try:
        import hdf5plugin  # noqa: F401 -- registers compression filters with HDF5
    except ImportError as error:
        raise ImportError("Install compression support: python -m pip install hdf5plugin") from error
    info_path = resolve_path(spec['info'], session) if spec.get('info') else only_file(session, ['*.info'])
    info = json.loads(info_path.read_text())
    keys = spec.get('region_keys', {'GC': 'gc', 'PC': 'pc'})
    if set(keys) != {'GC', 'PC'}:
        raise ValueError('region_keys must map GC and PC')
    layout = {r: flatten_ids(info['electrode_layout'][k]) for r, k in keys.items()}
    if set(layout['GC']) & set(layout['PC']):
        raise ValueError('Overlapping regional electrodes')
    fs, onset = float(spec['spike_fs_hz']), float(spec['onset_s'])
    if not np.isfinite([fs, onset]).all() or fs <= 0:
        raise ValueError('Invalid spike sampling rate/onset')
    mapping = spec['stimulus_by_din']
    if not mapping or len(set(mapping.values())) != len(mapping):
        raise ValueError('Require unique stimulus_by_din labels')
    arrays, labels, ids, strata = [], [], [], []
    with h5py.File(source, 'r') as h:
        units = np.asarray(h['unit_descriptor'])
        if not units.dtype.names or 'electrode_number' not in units.dtype.names:
            raise ValueError('unit_descriptor must have electrode_number')
        region = []
        for electrode in units['electrode_number']:
            hits = [r for r, electrodes in layout.items() if electrode in electrodes]
            region.append(hits[0] if len(hits) == 1 else 'excluded')
        if spec.get('single_units_only', True) and 'single_unit' not in units.dtype.names:
            raise ValueError('Missing single_unit flags; explicitly configure single_units_only=false after QC if appropriate')
        good = units['single_unit'].astype(bool) if spec.get('single_units_only', True) else np.ones(len(units), bool)
        selected = np.flatnonzero(np.isin(region, ['GC', 'PC']))
        if 'eligible_unit_ids' in spec:
            good &= np.isin(np.arange(len(units)), spec['eligible_unit_ids'])
        for din, label in sorted(mapping.items(), key=lambda item: int(item[0])):
            matches = [key for key in h['spike_trains'] if key.startswith('dig_in_') and int(key[7:]) == int(din)]
            if len(matches) != 1:
                raise ValueError(f'Missing or ambiguous spike_trains DIN {din}; generate sorted spike arrays first')
            group = h['spike_trains'][matches[0]]
            a = np.asarray(group['spike_array'])
            if a.ndim != 3 or a.shape[1] != len(units):
                raise ValueError('spike_array must be trial x unit_descriptor row x sample')
            if 'laser_durations' in group:
                laser = np.asarray(group['laser_durations'])
                lag = np.asarray(group['laser_onset_lag']) if 'laser_onset_lag' in group else np.zeros_like(laser)
                if laser.shape != (len(a),) or lag.shape != laser.shape or not np.isfinite(laser).all() or not np.isfinite(lag).all():
                    raise ValueError('Invalid per-trial laser metadata')
                strata.extend([f'laser:{d:g}:lag:{l:g}' if d else 'laser:off' for d, l in zip(laser, lag)])
            else:
                if info.get('laser_params', {}).get('dig_in_nums'):
                    raise ValueError('Optogenetic session lacks per-trial laser metadata; use reviewed NPZ')
                strata.extend(['laser:off']*len(a))
            arrays.append(a[:, selected])
            labels.extend([label]*len(a))
            ids.extend([f'DIN{int(din)}:trial{i}' for i in range(len(a))])
    x = np.concatenate(arrays)
    return dict(spikes=x, time_s=np.arange(x.shape[-1])/fs-onset,
                unit_id=selected, region=np.asarray(region)[selected], good_unit=good[selected],
                trial_id=np.array(ids), stimulus=np.array(labels), stratum=np.array(strata),
                trial_order=np.arange(len(x)), good_trial=np.ones(len(x), bool)), [source, info_path]


def windows(spec, time):
    result = [(str(k), float(v[0]), float(v[1])) for k, v in spec.get('windows_s', {'baseline': [-2, 0], 'evoked': [0, 2]}).items()]
    rolling = spec.get('rolling_s')
    if rolling:
        start, stop, width, step = map(float, rolling)
        if not np.isfinite([start, stop, width, step]).all() or width <= 0 or step <= 0 or stop-start < width:
            raise ValueError('rolling_s must be [start, stop, positive width, positive step]')
        result.extend((f'rolling_{i}', float(a), float(a+width)) for i, a in enumerate(np.arange(start, stop-width+1e-10, step)))
    if not result:
        raise ValueError('No analysis windows')
    for _, a, b in result:
        if not np.isfinite([a, b]).all() or a >= b or a < time[0]-1e-9 or b > time[-1]+np.diff(time)[0]+1e-9:
            raise ValueError('Window outside recorded sample interval or invalid bounds')
    return result


def analyze(data, spec, repeats, seed):
    rng = np.random.default_rng(seed)
    rows, nulls = [], {}
    units = np.flatnonzero(data['good_unit'])
    pairs = list(itertools.combinations(units, 2))
    mode = spec.get('detrend', 'first_difference')
    if mode not in ('first_difference', 'none'):
        raise ValueError('detrend must be first_difference or none')
    for label, a, b in windows(spec, data['time_s']):
        counts = data['spikes'][..., (data['time_s'] >= a-1e-10) & (data['time_s'] < b-1e-10)].sum(axis=-1, dtype=float)
        for stimulus, stratum in conditions(data):
            ix = condition_indices(data, stimulus, stratum)
            good = data['good_trial'][ix]
            values = counts[ix]
            if mode == 'first_difference':
                values = np.diff(values, axis=0)[good[1:] & good[:-1]]
            else:
                values = values[good]
            for i, j in progress(pairs, f'{stimulus}/{stratum} {label}: unit pairs'):
                stat, null = correlation_test(values[:, i], values[:, j], repeats, rng, spec.get('min_trials', 8))
                ri, rj = data['region'][i], data['region'][j]
                key = f'test_{len(rows)}'
                rows.append(dict(stimulus=stimulus, stratum=stratum, window=label, start_s=a, end_s=b,
                                 unit_a=str(data['unit_id'][i]), unit_b=str(data['unit_id'][j]),
                                 comparison=f'{ri}-{rj}' if ri == rj else 'GC-PC',
                                 detrend=mode, trials_total=len(ix), trials_accepted=int(good.sum()),
                                 null_key=key, **stat))
                nulls[key] = null
    correct(rows)
    return rows, nulls


def main():
    from .correlation_cli import main as run
    run('spike')


if __name__ == '__main__':
    main()
