"""Figure 2: simultaneous sorted spikes, alpha LFP and nasal pressure."""
import argparse
import json
import re
from pathlib import Path

import numpy as np
from scipy import signal
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from .events import write_csv
from .phase_coherence import only_file, flatten_ids
from .respiration_lfp import read_rows, fingerprint
from .respiration_spectra import load_inputs, FS


def load_sorted(session, events, raw_fs, time_s, h5=None, info=None):
    """Bin absolute blech_clust sorted-unit sample times on the LFP clock.

    Read original times rather than indexing spike_array by trial position: the
    latter can omit cutoff trials and uses a separate configured epoch duration.
    Preserve multiple spikes in a millisecond (counts, not binary occupancy).
    """
    import h5py
    try:
        import hdf5plugin  # noqa: F401
    except ImportError:
        pass  # Uncompressed HDF5 works without the optional filter package.
    session = Path(session).resolve()
    source = session / h5 if h5 else only_file(session, ['*.h5', '*.hdf5'])
    info_path = session / info if info else only_file(session, ['*.info'])
    if source.resolve().parent != session or info_path.resolve().parent != session:
        raise ValueError('Spike HDF5 and info must belong to the prepared session')
    layout = json.loads(info_path.read_text())['electrode_layout']
    regions = {r: flatten_ids(layout[r.lower()]) for r in ('GC', 'PC')}
    if set(regions['GC']) & set(regions['PC']):
        raise ValueError('Overlapping GC/PC electrode assignments')
    ids, labels, electrodes, arrays = [], [], [], []
    with h5py.File(source, 'r') as h:
        if 'sorted_units' not in h or 'unit_descriptor' not in h:
            raise ValueError('Run blech_clust spike sorting first: need sorted_units and unit_descriptor')
        desc = h['unit_descriptor'][:]
        if not {'electrode_number', 'single_unit'}.issubset(desc.dtype.names or ()):
            raise ValueError('Need electrode_number and single_unit in unit_descriptor')
        nodes = {}
        for name in h['sorted_units']:
            match = re.fullmatch(r'unit(\d+)', name)
            if match:
                uid = int(match[1])
                if uid in nodes:
                    raise ValueError('Duplicate sorted unit identity')
                nodes[uid] = name
        if set(nodes) != set(range(len(desc))):
            raise ValueError('sorted_units IDs do not match unit_descriptor rows')
        for uid, row in enumerate(desc):
            hits = [r for r, es in regions.items() if int(row['electrode_number']) in es]
            if not row['single_unit'] or len(hits) != 1:
                continue
            times = np.asarray(h['sorted_units'][nodes[uid]]['times']).reshape(-1)
            if not np.isfinite(times).all() or np.any(times < 0) or np.any(times != np.floor(times)):
                raise ValueError('Sorted spike times must be nonnegative raw sample indices')
            times = np.sort(times)
            epochs = []
            for event in events:
                # Match the rounded LFP/ADC center, not a separately rounded spike epoch.
                center = int(event['center_lfp_sample']) / FS * raw_fs
                edges = center + np.r_[time_s, time_s[-1] + 1/FS] * raw_fs
                indices = np.searchsorted(times, edges, side='left')
                epochs.append(np.diff(indices))
            arrays.append(epochs); ids.append(uid); labels.append(hits[0])
            electrodes.append(int(row['electrode_number']))
    if set(labels) != {'GC', 'PC'}:
        raise ValueError('Need at least one sorted single unit assigned to each of GC and PC')
    return (np.asarray(arrays, dtype=np.int32), np.array(ids), np.array(labels),
            np.array(electrodes), [source, info_path])


def alpha_filter(x):
    """Zero-phase 8–12 Hz, fourth-order Butterworth; filter full 7 s epochs."""
    return signal.sosfiltfilt(signal.butter(4, [8, 12], fs=FS, btype='bandpass', output='sos'), x, axis=-1)


def run(prepared, audit, output, h5=None, info=None, check=False):
    data, _, pressure, accepted, qc, meta, resp = load_inputs(prepared, audit)
    events = read_rows(Path(prepared)/'events.csv')
    t = data['time_s']
    spikes, units, regions, electrodes, sources = load_sorted(
        meta['source'], events, meta['raw_fs_hz'], t, h5, info)
    # Keep the same events in every row of a stimulus figure; choose deterministic
    # channels by QC coverage, never by alpha amplitude or apparent coupling.
    selections = []
    for stim in dict.fromkeys(data['stimulus']):
        candidates = accepted & (data['stimulus'] == stim)
        channels = []
        for region in ('GC', 'PC'):
            possible = np.flatnonzero(data['good_channel'] & (data['region'] == region))
            if not len(possible):
                raise ValueError('No eligible LFP channel in '+region)
            scores = [np.sum(candidates & data['good_trial'][c] & np.isfinite(data['lfp_uv'][c]).all(axis=1)) for c in possible]
            channels.append(int(possible[np.argmax(scores)]))
        keep = candidates.copy()
        for c in channels:
            keep &= data['good_trial'][c] & np.isfinite(data['lfp_uv'][c]).all(axis=1)
        selections.append((str(stim), channels, np.flatnonzero(keep)))
    if not any(len(js) for _, _, js in selections):
        raise ValueError('No simultaneous pressure/GC/PC LFP trials pass QC')
    if check:
        print(f'Ready: {len(units)} single units; '+', '.join(f'{s}: {len(js)} matched trials' for s, _, js in selections))
        return
    out = Path(output)
    out.mkdir(parents=True, exist_ok=False)
    write_csv(out/'pressure_trial_qc.csv', qc)
    write_csv(out/'units.csv', [dict(unit_id=int(u), region=r, electrode=int(e)) for u,r,e in zip(units,regions,electrodes)])
    archive = dict(time_s=t, unit_id=units, unit_region=regions)
    retention = []
    with np.load(Path(prepared)/'cleaned_lfp.npz') as z:
        channel_ids = z['electrode_id'] if 'electrode_id' in z else np.arange(len(data['region']))
    for number, (stim, channels, js) in enumerate(selections, 1):
        print(f'Figure 2 [{number}/{len(selections)}] {stim}: {len(js)} matched trials', flush=True)
        retention.append(dict(stimulus=stim, total_trials=int(np.sum(data['stimulus']==stim)), matched_trials=len(js),
                              GC_channel=str(channel_ids[channels[0]]), PC_channel=str(channel_ids[channels[1]])))
        if not len(js):
            continue
        j = js[0]
        alpha = np.stack([alpha_filter(data['lfp_uv'][c, js]) for c in channels])
        fig, axes = plt.subplots(6, 1, figsize=(11, 11), sharex=True, layout='constrained')
        axes[0].plot(t, pressure[j], color='black', lw=.8)
        axes[0].set(ylabel='Pressure\n(ADC counts)', title=f'A  {stim} | event {events[j]["event_id"]} | pressure polarity: {resp["sign"]}')
        for r, region in enumerate(('GC', 'PC')):
            axes[r+1].plot(t, alpha[r, 0], lw=.8, color=('tab:blue','tab:orange')[r])
            axes[r+1].set(ylabel=f'{region} alpha\n(µV)', title=f'B  8–12 Hz · electrode {channel_ids[channels[r]]}')
            us = np.flatnonzero(regions == region)
            axes[r+3].eventplot([t[spikes[u,j]>0] for u in us], colors=('tab:blue','tab:orange')[r], lineoffsets=np.arange(len(us)), linelengths=.8)
            axes[r+3].set(yticks=np.arange(len(us)), yticklabels=units[us], ylabel=f'{region} unit ID', title='C  Same example trial · sorted single units')
            # Non-overlapping 50 ms bins; equal unit mean, then matched-trial mean.
            counts = spikes[us][:,js].reshape(len(us),len(js),-1,50).sum(axis=-1)
            rate = counts.mean(axis=0)/.05
            bt = t.reshape(-1,50).mean(axis=1)
            axes[5].plot(bt, rate.mean(axis=0), label=f'{region}: {len(us)} units, {len(js)} trials')
            archive[f'{stim}__{region}_rate_hz'] = rate
        axes[5].set(ylabel='Spikes/s/unit', xlabel='Time from DIN command (s)', title='D  Matched-trial population firing rate · 50 ms bins')
        axes[5].legend(fontsize=8)
        for ax in axes:
            ax.axvline(0, color='.5', lw=.7, ls='--'); ax.set_xlim(-.5, 2.5)
            ax.spines[['top','right']].set_visible(False)
        axes[0].plot([2,2.2],[.08,.08],transform=axes[0].get_xaxis_transform(),color='black',lw=2)
        axes[0].text(2.1,.1,'200 ms',transform=axes[0].get_xaxis_transform(),ha='center',fontsize=8)
        fig.suptitle(Path(meta['source']).name+' · simultaneous recording examples')
        for ext in ('png','pdf'):
            fig.savefig(out/f'{number:02d}_figure2.{ext}', dpi=180)
        plt.close(fig)
        archive[f'{stim}__event_id'] = data['event_id'][js]
        archive[f'{stim}__pressure_counts'] = pressure[js]
        archive[f'{stim}__alpha_uv'] = alpha
        archive[f'{stim}__spike_counts'] = spikes[:,js]
    archive['rate_time_s'] = t.reshape(-1,50).mean(axis=1)
    np.savez_compressed(out/'plotted_data.npz', **archive)
    write_csv(out/'retention.csv',retention)
    inputs = sources + [Path(prepared)/n for n in ('cleaned_lfp.npz','events.csv','preprocessing.json')] + [Path(audit)/n for n in ('pressure_1khz.npz','events_respiration.csv','respiration_provenance.json')]
    (out/'provenance.json').write_text(json.dumps(dict(lfp=meta,respiration=resp,
        hashes={str(p.resolve()):fingerprint(p) for p in inputs},
        alpha_hz=[8,12], filter='Butterworth order 4, forward/backward SOS on full -2 to +5 s epochs; display -0.5 to +2.5 s',
        spikes='sorted_units/unitNNN/times: raw sample indices; unit_descriptor single_unit; 1 ms counts on rounded LFP clock',
        selection='Per stimulus, each region channel with most valid pressure-matched trials; first shared accepted event',
        limitations='Descriptive session validation, not coupling inference. Single-unit flags are not independent sorting/drift QC. Laser conditions are not separated. No histology or waveform-quality claim. DIN is not odor arrival. Unknown polarity does not identify inhalation.'),indent=2))
    return out


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('prepared'); parser.add_argument('audit'); parser.add_argument('output')
    parser.add_argument('--h5', help='HDF5 filename inside session (required if ambiguous)')
    parser.add_argument('--info', help='Region .info filename inside session (required if ambiguous)')
    parser.add_argument('--check', action='store_true')
    a = parser.parse_args()
    run(a.prepared,a.audit,a.output,a.h5,a.info,a.check)


if __name__ == '__main__':
    main()
