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
from .respiration_spectra import load_inputs, FS, PAIRS


def unit_regions(layout, prepared_regions=None):
    """Resolve anatomical names without guessing numeric CAR-group meanings."""
    if prepared_regions is not None:
        regions = {r: list(prepared_regions[r]) for r in ('GC', 'PC')}
    else:
        regions = {}
        for region in ('GC', 'PC'):
            matches = [key for key in layout if str(key).strip().upper() == region]
            if len(matches) != 1:
                raise ValueError(
                    f'Cannot identify {region} in electrode_layout keys {list(layout)}. '
                    'Use --unit-regions-from-prepared only when unit_descriptor electrode_number '
                    'uses the same electrode IDs as the prepared LFP mapping.')
            regions[region] = flatten_ids(layout[matches[0]])
    if any(not regions[r] for r in ('GC', 'PC')):
        raise ValueError('Empty GC/PC electrode mapping')
    if set(regions['GC']) & set(regions['PC']):
        raise ValueError('Overlapping GC/PC electrode assignments')
    return regions


def load_sorted(session, events, raw_fs, time_s, h5=None, info=None, prepared_regions=None):
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
    regions = unit_regions(layout, prepared_regions)
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


REGION_COLORS = ('#007C83', '#B56A3B')
RAW_COLOR = '#B8BEC4'
ENVELOPE_COLOR = '#9B4264'


def amplitude_envelope(x):
    """Descriptive analytic-signal magnitude; not an inhalation detector."""
    return np.abs(signal.hilbert(x, axis=-1))


def plot_column(axes, t, units, regions, spikes, pressure, alpha, channel_ids, event_id,
                broadband=None, raw_pressure=None):
    """Rasters first; all example traces use the first accepted event."""
    rates = {}
    for r, region in enumerate(('GC', 'PC')):
        color = REGION_COLORS[r]
        us = np.flatnonzero(regions == region)
        axes[r].eventplot([t[spikes[u,0]>0] for u in us], colors=color,
                          lineoffsets=np.arange(len(us)), linelengths=.8, linewidths=.8)
        # Sparse tick labels avoid overlap for large populations; all units remain plotted.
        ticks = np.arange(0, len(us), max(1, int(np.ceil(len(us)/12))))
        axes[r].set(yticks=ticks, yticklabels=units[us[ticks]], ylim=(-.7,len(us)-.3),
                    ylabel=f'{region} unit ID', title=f'A  {region} single units · event {event_id}')
        if broadband is not None:
            axes[r+2].plot(t,broadband[r,0],color=RAW_COLOR,lw=.65,alpha=.8,
                           label='Prepared LFP (1–100 Hz)',zorder=1)
        axes[r+2].plot(t,alpha[r,0],color=color,lw=1.1,label='Alpha (8–12 Hz)',zorder=2)
        axes[r+2].plot(t,amplitude_envelope(alpha[r,0]),color=ENVELOPE_COLOR,lw=1.25,
                       label='Alpha amplitude envelope',zorder=3)
        axes[r+2].set(ylabel=f'{region} LFP\n(µV)', title=f'B  LFP and alpha · electrode {channel_ids[r]}')
        axes[r+2].legend(loc='upper right',fontsize=7,frameon=True,facecolor='white',
                         edgecolor='none',framealpha=.85,ncol=3)
        counts = spikes[us].reshape(len(us),spikes.shape[1],-1,50).sum(axis=-1)
        rates[region] = counts.mean(axis=0)/.05
        axes[5].plot(t.reshape(-1,50).mean(axis=1),rates[region].mean(axis=0),color=color,
                     label=f'{region}: {len(us)} units, {spikes.shape[1]} trials')
    if raw_pressure is not None:
        raw_centered = raw_pressure[0]-np.mean(raw_pressure[0])
        axes[4].plot(t,raw_centered,color=RAW_COLOR,lw=.65,alpha=.8,
                     label='ADC (1 kHz, mean-centered)',zorder=1)
    axes[4].plot(t,pressure[0],color='#263238',lw=1.1,label='Pressure (0.5–20 Hz)',zorder=2)
    axes[4].plot(t,amplitude_envelope(pressure[0]),color=ENVELOPE_COLOR,lw=1.25,
                 label='Pressure amplitude envelope',zorder=3)
    axes[4].legend(loc='upper right',fontsize=7,frameon=True,facecolor='white',
                   edgecolor='none',framealpha=.85,ncol=2)
    axes[4].set(ylabel='Pressure\n(ADC counts)',title='C  Respiration · same example trial')
    axes[4].plot([2,2.2],[.08,.08],transform=axes[4].get_xaxis_transform(),color='black',lw=2)
    axes[4].text(2.1,.1,'200 ms',transform=axes[4].get_xaxis_transform(),ha='center',fontsize=8)
    axes[5].set(ylabel='Spikes/s/unit',xlabel='Time from DIN command (s)',
                title='D  Population firing rate · 50 ms bins')
    axes[5].legend(fontsize=8,frameon=False)
    for ax in axes:
        ax.axvline(0,color='#87929A',lw=.7,ls='--'); ax.set_xlim(-.5,2.5)
        ax.spines[['top','right']].set_visible(False)
        ax.tick_params(colors='#37434A',labelsize=9)
        for side in ('left','bottom'):
            ax.spines[side].set_color('#87929A'); ax.spines[side].set_linewidth(.6)
    return rates


def run(prepared, audit, output, h5=None, info=None, check=False, unit_regions_from_prepared=False):
    data, raw_pressure, pressure, accepted, qc, meta, resp = load_inputs(prepared, audit)
    events = read_rows(Path(prepared)/'events.csv')
    t = data['time_s']
    prepared_regions = None
    if unit_regions_from_prepared:
        with np.load(Path(prepared)/'cleaned_lfp.npz', allow_pickle=False) as z:
            ids = z['electrode_id']
        if ids.shape != data['region'].shape or len(np.unique(ids)) != len(ids):
            raise ValueError('Invalid prepared electrode IDs')
        prepared_regions = {r: ids[data['region'] == r].astype(int).tolist() for r in ('GC', 'PC')}
        print('Unit regions: using prepared LFP electrode IDs (including LFP-QC-excluded electrodes)', flush=True)
    spikes, units, regions, electrodes, sources = load_sorted(
        meta['source'], events, meta['raw_fs_hz'], t, h5, info, prepared_regions)
    # Keep the same events in every row of a stimulus figure; choose deterministic
    # channels by QC coverage, never by alpha amplitude or apparent coupling.
    selections = []
    for stim in dict.fromkeys(data['stimulus']):
        candidates = accepted & (data['stimulus'] == stim)
        paired = next((pair for pair in PAIRS.values() if stim in pair), (stim,))
        masks = [accepted & (data['stimulus'] == label) for label in paired]
        channels = []
        for region in ('GC', 'PC'):
            possible = np.flatnonzero(data['good_channel'] & (data['region'] == region))
            if not len(possible):
                raise ValueError('No eligible LFP channel in '+region)
            scores = []
            for c in possible:
                valid = data['good_trial'][c] & np.isfinite(data['lfp_uv'][c]).all(axis=1)
                counts = [int(np.sum(mask & valid)) for mask in masks]
                scores.append((min(counts),sum(counts)))
            channels.append(int(possible[max(range(len(possible)), key=lambda k: scores[k])]))
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
        fig, axes = plt.subplots(6, 1, figsize=(11, 13), sharex=True, layout='constrained',
                                 gridspec_kw={'height_ratios':[1.6,1.2,1,1,1,1]})
        rates = plot_column(axes,t,units,regions,spikes[:,js],pressure[js],alpha,
                            channel_ids[channels],events[j]['event_id'],
                            broadband=data['lfp_uv'][channels][:,js],raw_pressure=raw_pressure[js])
        for region, rate in rates.items():
            archive[f'{stim}__{region}_rate_hz'] = rate
        fig.suptitle(Path(meta['source']).name+f' · {stim} · simultaneous recording examples')
        for ext in ('png','pdf'):
            fig.savefig(out/f'{number:02d}_figure2.{ext}', dpi=180)
        plt.close(fig)
        archive[f'{stim}__electrode_id'] = channel_ids[channels]
        archive[f'{stim}__event_id'] = data['event_id'][js]
        archive[f'{stim}__pressure_counts'] = pressure[js]
        archive[f'{stim}__alpha_uv'] = alpha
        archive[f'{stim}__prepared_lfp_uv'] = data['lfp_uv'][channels][:,js]
        archive[f'{stim}__alpha_envelope_uv'] = amplitude_envelope(alpha)
        archive[f'{stim}__raw_pressure_counts'] = raw_pressure[js]
        archive[f'{stim}__pressure_envelope_counts'] = amplitude_envelope(pressure[js])
        archive[f'{stim}__spike_counts'] = spikes[:,js]
    for odor, pair in PAIRS.items():
        fig, axes = plt.subplots(6,2,figsize=(16,13),sharex=True,sharey='row',
                                 layout='constrained',gridspec_kw={'height_ratios':[1.6,1.2,1,1,1,1]})
        for col, (route,stim) in enumerate(zip(('ORTHONASAL','RETRONASAL'),pair)):
            if f'{stim}__event_id' not in archive:
                for ax in axes[:,col]:
                    ax.text(.5,.5,'No matched accepted trials',ha='center',transform=ax.transAxes)
                axes[0,col].set_title(f'{route} · {stim}')
                continue
            plot_column(axes[:,col],t,units,regions,archive[f'{stim}__spike_counts'],
                        archive[f'{stim}__pressure_counts'],archive[f'{stim}__alpha_uv'],
                        archive[f'{stim}__electrode_id'],archive[f'{stim}__event_id'][0],
                        broadband=archive[f'{stim}__prepared_lfp_uv'],
                        raw_pressure=archive[f'{stim}__raw_pressure_counts'])
            axes[0,col].set_title(f'{route} · {stim}\n'+axes[0,col].get_title(),fontweight='bold')
        fig.suptitle(f'{odor} · Orthonasal versus retronasal\n'+Path(meta['source']).name,
                     fontsize=15,fontweight='bold')
        for ext in ('png','pdf'):
            fig.savefig(out/f'{odor}_ortho_vs_retro.{ext}',dpi=180)
        plt.close(fig)
    archive['rate_time_s'] = t.reshape(-1,50).mean(axis=1)
    np.savez_compressed(out/'plotted_data.npz', **archive)
    write_csv(out/'retention.csv',retention)
    inputs = sources + [Path(prepared)/n for n in ('cleaned_lfp.npz','events.csv','preprocessing.json')] + [Path(audit)/n for n in ('pressure_1khz.npz','events_respiration.csv','respiration_provenance.json')]
    (out/'provenance.json').write_text(json.dumps(dict(lfp=meta,respiration=resp,
        hashes={str(p.resolve()):fingerprint(p) for p in inputs},
        unit_region_source="prepared LFP electrode IDs" if unit_regions_from_prepared else "session info electrode_layout",
        prepared_unit_regions=prepared_regions,
        display=dict(region_colors=dict(zip(('GC','PC'),REGION_COLORS)),raw_color=RAW_COLOR,
                     envelope_color=ENVELOPE_COLOR,
                     background_lfp='Prepared 1–100 Hz referenced LFP, not raw amplifier voltage',
                     background_pressure='Anti-aliased 1 kHz ADC, mean removed over full 7 s epoch for display',
                     envelope='abs(Hilbert(filtered signal)) over full epoch; pressure envelope is descriptive broadband amplitude, not inhalation timing'),
        alpha_hz=[8,12], filter='Butterworth order 4, forward/backward SOS on full -2 to +5 s epochs; display -0.5 to +2.5 s',
        spikes='sorted_units/unitNNN/times: raw sample indices; unit_descriptor single_unit; 1 ms counts on rounded LFP clock',
        selection='For matched odors, same electrodes across routes: maximize minimum route trial coverage, then total coverage; ties archive order. First shared accepted event per stimulus. Other stimuli maximize coverage.',
        limitations='Descriptive session validation, not coupling inference. Single-unit flags are not independent sorting/drift QC. Laser conditions are not separated. No histology or waveform-quality claim. DIN is not odor arrival. Unknown polarity does not identify inhalation.'),indent=2))
    return out


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('prepared'); parser.add_argument('audit'); parser.add_argument('output')
    parser.add_argument('--h5', help='HDF5 filename inside session (required if ambiguous)')
    parser.add_argument('--info', help='Region .info filename inside session (required if ambiguous)')
    parser.add_argument('--check', action='store_true')
    parser.add_argument('--unit-regions-from-prepared', action='store_true',
                        help='Use prepared electrode IDs/GC-PC labels for units; requires matching electrode numbering')
    a = parser.parse_args()
    run(a.prepared,a.audit,a.output,a.h5,a.info,a.check,a.unit_regions_from_prepared)


if __name__ == '__main__':
    main()
