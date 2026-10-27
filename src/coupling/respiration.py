"""Read-only TTL/pressure timing audit; physiology labels require verified polarity."""

from coupling.array_store import save_arrays
import argparse
import json
import math
from pathlib import Path
import numpy as np
from scipy.signal import butter, sosfiltfilt, resample_poly, find_peaks, welch
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from .events import sample_rate, read_binary, inventory, select_six, write_csv
from .prepare import check_time
from .qc import flat_fraction

FS = 1000
SETTINGS = dict(fs_hz=FS, bandpass_hz=[.5, 20], butterworth_order=4,
                minimum_cycle_hz=1., maximum_cycle_hz=15.,
                peak_prominence_scaled_mad=.5, search_after_ttl_s=2.,
                filter_edge_guard_s=5., raw_flat_ms=10., raw_bad_fraction=.01)


def downsample_counts(raw, fs):
    """Centered unsigned ADC counts; bounded raw-rate memory and FIR context."""
    if fs != int(fs):
        raise ValueError('Integer acquisition rate required')
    fs = int(fs)
    g = math.gcd(fs, FS)
    up, down = FS//g, fs//g
    result = np.empty(math.ceil(len(raw)*up/down), dtype=float)
    for start in range(0, len(raw), 30*fs):
        end = min(len(raw), start+30*fs)
        lo, hi = max(0, start-fs), min(len(raw), end+fs)
        x = np.asarray(raw[lo:hi], dtype=float)-32768.
        y = resample_poly(x, up, down)
        a, b = start*up//down, min(len(result), math.ceil(end*up/down))
        offset = (start-lo)*up//down
        result[a:b] = y[offset:offset+b-a]
    return result


def detect_cycles(filtered, minimum_hz=1., maximum_hz=15., prominence_mad=.5):
    """Negative-going crossing before a prominent trough, then positive crossing.

    Samples are explicitly pressure-phase candidates until polarity is verified.
    Local periods and trough prominence are screening heuristics, not probabilities.
    """
    if not np.isfinite([minimum_hz, maximum_hz, prominence_mad]).all() or not (0 < minimum_hz < maximum_hz <= 50) or prominence_mad <= 0:
        raise ValueError("Invalid cycle screening settings")
    x = np.asarray(filtered, dtype=float)
    if x.ndim != 1 or not np.isfinite(x).all():
        raise ValueError('Need finite 1D pressure')
    scale = 1.4826*np.median(np.abs(x-np.median(x)))
    if scale <= 0:
        return []
    neg = np.flatnonzero((x[:-1] >= 0) & (x[1:] < 0))+1
    pos = np.flatnonzero((x[:-1] <= 0) & (x[1:] > 0))+1
    peaks, props = find_peaks(-x, distance=math.ceil(FS/maximum_hz), prominence=prominence_mad*scale)
    rows = []
    seen = set()
    for peak, prominence in zip(peaks, props['prominences']):
        i = np.searchsorted(neg, peak, side='right')-1
        if i < 0 or i+1 >= len(neg):
            continue
        onset, next_onset = int(neg[i]), int(neg[i+1])
        j = np.searchsorted(pos, peak, side='right')
        if j == len(pos) or pos[j] >= next_onset or onset in seen:
            continue
        seen.add(onset)
        offset = int(pos[j])
        hz = FS/(next_onset-onset)
        reasons = []
        if not minimum_hz <= hz <= maximum_hz:
            reasons.append('cycle_rate_outside_configured_bounds')
        if min(offset-onset, next_onset-offset) < 10:
            reasons.append('phase_shorter_than_10ms')
        rows.append(dict(cycle_id=len(rows), negative_start_sample=onset,
                         negative_peak_sample=int(peak), positive_start_sample=offset,
                         next_negative_start_sample=next_onset, cycle_hz=float(hz),
                         prominence_scaled_mad=float(prominence/scale),
                         screen_pass=not reasons, reasons=';'.join(reasons)))
    return rows


def first_after(samples, onset, max_delay=2.):
    samples = np.asarray(samples, dtype=int)
    i = np.searchsorted(samples, onset, side='left')
    if i == len(samples) or samples[i]-onset > max_delay*FS:
        return None
    return int(samples[i])


def _latency(sample, ttl):
    return '' if sample is None else float((sample-ttl)/FS)


def plot_audit(out, unfiltered, filtered, events, cycles, sign):
    label = 'candidate inhalation' if sign != 'unknown' else 'candidate negative phase'
    pos_label = 'candidate exhalation' if sign != 'unknown' else 'candidate positive phase'
    fig, axes = plt.subplots(3, 2, figsize=(10, 9), layout='constrained')
    ax = axes[0, 0]
    # First in-bounds event, selected without optimizing for a compelling response.
    example = next((e for e in events if e['snippet_in_bounds']), None)
    if example:
        k = example['ttl_lfp_sample']
        idx = np.arange(k-2000, k+5000)
        t = (idx-k)/FS
        raw = unfiltered[idx]-np.median(unfiltered[idx])
        ax.plot(t, raw, color='.7', lw=.6, label='Anti-aliased pressure')
        filtered_ax = ax.twinx()
        filtered_ax.plot(t, filtered[idx], color='black', lw=.8, label='0.5–20 Hz filtered')
        filtered_ax.set_ylabel('Filtered ADC counts', fontsize=8)
        for c in cycles:
            for key, color in [('negative_start_sample','tab:blue'),('positive_start_sample','tab:orange')]:
                sample = c[key]
                if idx[0] <= sample <= idx[-1]:
                    filtered_ax.scatter((sample-k)/FS, filtered[sample], color=color, s=10)
        ax.axvline(0, color='red', ls='--', lw=.8)
        ax.set_title(f'A  First in-bounds trial: {example["stimulus"]}, event {example["event_id"]}')
        ax.legend(fontsize=7,loc='upper left')
        filtered_ax.legend(fontsize=7,loc='upper right')
    ax.set(xlabel='Time from TTL (s)', ylabel='Centered ADC counts')
    stimuli = ('BA','IA','S','oBA','oIA','BK')
    freq, psd = welch(unfiltered, fs=FS, nperseg=min(8000,len(unfiltered)))
    axes[0,1].semilogy(freq,psd,color='black',lw=.8)
    axes[0,1].axvspan(.5,20,color='tab:blue',alpha=.1)
    axes[0,1].set(xlim=(0,100),title='B  Anti-aliased pressure spectrum',xlabel='Frequency (Hz)',ylabel='ADC counts²/Hz')
    for axis, alignment, title in [(axes[1,0],'ttl_lfp_sample','C  TTL-aligned pressure'),
                                   (axes[1,1],'negative_start_sample','D  First negative-phase alignment')]:
        for si, stim in enumerate(stimuli):
            snippets = []
            for e in events:
                if e['stimulus'] != stim or e['timing_qc_reasons']:
                    continue
                center = e[alignment]
                if center == '':
                    continue
                idx = np.arange(int(center)-1000, int(center)+2001)
                if idx[0] < 0 or idx[-1] >= len(filtered):
                    continue
                y = filtered[idx]
                scale = np.std(y)
                if scale > 0:
                    snippets.append((y-np.mean(y))/scale)
            if snippets:
                axis.plot(np.arange(-1000,2001)/FS, np.mean(snippets,axis=0)+si*3,
                          lw=1, label=f'{stim} (n={len(snippets)})')
        axis.axvline(0,color='.4',ls='--',lw=.8)
        axis.set(title=title, xlabel='Time from alignment (s)', ylabel='Standardized pressure + offsets')
        axis.legend(fontsize=7, ncol=2)
    axis = axes[2,0]
    for i, stim in enumerate(stimuli):
        for key, offset, color in [('negative_latency_s',-.12,'tab:blue'),('positive_latency_s',.12,'tab:orange')]:
            vals = [e[key] for e in events if e['stimulus']==stim and not e['timing_qc_reasons'] and e[key]!='']
            jitter = np.linspace(-.05,.05,len(vals))
            axis.scatter(vals, i+offset+jitter, s=9, color=color, alpha=.5)
        rows = [e for e in events if e['stimulus']==stim]
        usable = sum(not e['timing_qc_reasons'] for e in rows)
        axis.text(2.03,i,f'{usable}/{len(rows)}',fontsize=7,va='center')
    axis.set(yticks=range(6),yticklabels=stimuli,xlim=(0,2.35),
             xlabel='Delay after TTL (s)',title='E  Candidate phase delays after TTL')
    axis.scatter([],[],color='tab:blue',label=label,s=10)
    axis.scatter([],[],color='tab:orange',label=pos_label,s=10)
    axis.legend(fontsize=7)
    rate_ax=axes[2,1]
    for passed,color in [(True,'tab:blue'),(False,'tab:orange')]:
        rates=[c['cycle_hz'] for c in cycles if c['screen_pass']==passed]
        rate_ax.hist(rates,bins=np.arange(0,31),alpha=.5,color=color,label='Screen pass' if passed else 'Flagged')
    rate_ax.set(xlim=(0,30),xlabel='Candidate cycle rate (Hz)',ylabel='Cycles',title='F  Cycle screening is not breath validation')
    rate_ax.legend(fontsize=8)
    fig.suptitle('Delivery/pressure timing audit · '+('polarity UNVERIFIED' if sign=='unknown' else 'user-verified polarity'),fontsize=12)
    for ext in ('png','pdf','svg'):
        fig.savefig(out/f'delivery_respiration.{ext}',dpi=300)
    plt.close(fig)


def audit(session, config_file, output, adc='board-ADC-00.dat', sign='unknown', polarity_note='', waveform_review_note=''):
    session, out = Path(session).resolve(), Path(output).resolve()
    if out.exists():
        raise FileExistsError(out)
    if out == session or session in out.parents:
        raise ValueError('Output must be outside recording')
    if sign not in ('unknown','negative','positive'):
        raise ValueError('Unknown inhalation sign')
    if sign != 'unknown' and not polarity_note.strip():
        raise ValueError('Verified polarity needs --polarity-note documenting independent evidence')
    config = json.loads(Path(config_file).read_text())
    fs = sample_rate(session/config.get('header','info.rhd'))
    if config.get('raw_fs_hz',fs) != fs:
        raise ValueError('Configured rate disagrees with RHD')
    records, inventory_rows = inventory(session,fs)
    selected = select_six({k:len(v[0]) for k,v in records.items()})
    labels = config['stimulus_by_din']
    from .power import STIMULI
    if set(labels) != set(map(str,selected)) or set(labels.values()) != set(STIMULI):
        raise ValueError('Stimulus labels must exactly cover selected DINs and expected stimuli')
    n = next(iter(records.values()))[2]
    tick = check_time(session/config.get('timestamps','time.dat'),n)
    adc_path = (session/adc).resolve()
    if adc_path.parent != session or not adc_path.name.startswith('board-ADC-'):
        raise ValueError('ADC must be a board-ADC file directly inside recording')
    raw = read_binary(adc_path,'<u2')
    if len(raw) != n:
        raise ValueError('ADC/DIN sample lengths differ')
    if np.ptp(raw[::max(1,int(fs/1000))]) == 0:
        raise ValueError('Constant ADC signal')
    print('Resampling ADC counts and filtering continuous pressure',flush=True)
    unfiltered = downsample_counts(raw,fs)
    filtered = sosfiltfilt(butter(4,[.5,20],fs=FS,btype='bandpass',output='sos'),unfiltered)
    if sign == 'positive':
        unfiltered = -unfiltered
        filtered = -filtered
    cycles = detect_cycles(filtered)
    passed = [c for c in cycles if c['screen_pass']]
    negative = np.array([c['negative_start_sample'] for c in passed],dtype=int)
    positive = np.array([c['positive_start_sample'] for c in passed],dtype=int)
    events = []
    for din in selected:
        for trial,(on,off) in enumerate(zip(*records[din][:2])):
            k = int(np.rint(on/fs*FS))
            bounds = k-2000>=0 and k+5000<=len(filtered)
            reasons = []
            rail, flat = None, None
            if not bounds:
                reasons.append('snippet_boundary')
            else:
                x = np.asarray(raw[max(0,int(on-2*fs)):min(n,int(on+5*fs))])
                rail = float(np.mean((x==0)|(x==65535)))
                flat = flat_fraction(x,math.ceil(.01*fs))
                if rail > .01:
                    reasons.append('raw_rail_fraction_gt_1pct')
                if flat > .01:
                    reasons.append('raw_flat_fraction_gt_1pct')
            if k < 5*FS or k+5000 > len(filtered)-5*FS:
                reasons.append('filter_recording_edge_guard')
            neg, pos = first_after(negative,k), first_after(positive,k)
            if neg is None:
                reasons.append('no_screened_negative_phase_within_2s')
            if pos is None:
                reasons.append('no_screened_positive_phase_within_2s')
            events.append(dict(din=din,stimulus=labels[str(din)],original_trial_index=trial,
                               ttl_raw_sample=int(on),offset_raw_sample=int(off),ttl_lfp_sample=k,
                               ttl_s=float(on/fs),pulse_ms=float((off-on)/fs*1000),
                               rounding_error_ms=float(k-on/fs*FS),snippet_in_bounds=bounds,
                               negative_start_sample='' if neg is None else neg,
                               positive_start_sample='' if pos is None else pos,
                               negative_latency_s=_latency(neg,k),positive_latency_s=_latency(pos,k),
                               inhalation_latency_s=_latency(neg,k) if sign!='unknown' and waveform_review_note.strip() else '',
                               exhalation_latency_s=_latency(pos,k) if sign!='unknown' and waveform_review_note.strip() else '',
                               raw_rail_fraction=rail,raw_flat_fraction=flat,
                               timing_qc_reasons=';'.join(reasons),polarity_verified=sign!='unknown',
                               review_flags=';'.join((["polarity_unverified"] if sign=='unknown' else [])+(["waveform_unreviewed"] if not waveform_review_note.strip() else [])),
                               usable_for_physiological_timing=not reasons and sign!='unknown' and bool(waveform_review_note.strip()),
                               waveform_reviewed=bool(waveform_review_note.strip())))
    events.sort(key=lambda e:(e['ttl_raw_sample'],e['din']))
    for i,e in enumerate(events):
        e['event_id']=i
    out.mkdir(parents=True)
    write_csv(out/'events_respiration.csv',events)
    if cycles:
        write_csv(out/'pressure_cycles.csv',cycles)
    write_csv(out/'din_inventory.csv',inventory_rows)
    save_arrays(out/'pressure_1khz.h5',unfiltered_counts=unfiltered.astype('float32'),
                        filtered_counts=filtered.astype('float32'),fs_hz=FS)
    freq,psd=welch(unfiltered,fs=FS,nperseg=min(8000,len(unfiltered)))
    band=psd[(freq>=.5)&(freq<=20)]
    spectral_flatness=float(np.exp(np.mean(np.log(np.maximum(band,np.finfo(float).tiny))))/np.mean(band)) if np.mean(band)>0 else None
    meta = dict(pressure_spectral_flatness_0p5_20hz=spectral_flatness,
                candidate_cycles=len(cycles),screen_passing_cycles=len(passed),session_id=config.get('session_id'),animal_id=config.get('animal_id'),
                source=str(session),settings=SETTINGS,raw_fs_hz=fs,first_timestamp=tick,
                adc_dtype='little-endian uint16',amplitude_units='centered ADC counts; not volts or pressure units',
                sign=sign,polarity_note=polarity_note,waveform_review_note=waveform_review_note,
                alignment='shared sample index and validated consecutive time.dat; hardware/catheter delay unmeasured',
                no_swallow_measurement=True,configuration=config,selected_dins=selected,total_trials=len(events),
                candidate_screen_pass_trials=sum(not e['timing_qc_reasons'] for e in events),
                physiological_timing_trials=sum(e['usable_for_physiological_timing'] for e in events),
                files=[dict(name=p.name,size_bytes=p.stat().st_size,mtime_ns=p.stat().st_mtime_ns)
                       for p in [adc_path,session/config.get('header','info.rhd'),
                                 session/config.get('timestamps','time.dat')]+sorted(session.glob('board-DIN-*.dat'))])
    (out/'respiration_provenance.json').write_text(json.dumps(meta,indent=2))
    plot_audit(out,unfiltered,filtered,events,cycles,sign)
    (out/'caption.txt').write_text('Single-session pressure timing audit. A: first in-bounds TTL trial, anti-aliased and filtered ADC counts on separate y axes, with candidate phase markers. B: anti-aliased pressure PSD. C/D: per-trial standardized pressure averaged by stimulus with vertical offsets, aligned to TTL or first screened negative-going phase start after TTL. Phase alignment sharpens pressure by construction and does not establish neural locking. E: first screened negative/positive phase starts within 2 s of TTL; counts are screen-passing/total trials. F: candidate cycle rates and screening flags. Filters are zero-phase; crossing timing is an operational candidate, not calibrated airflow onset. '+('Polarity unverified: no inhalation/exhalation inference. ' if sign=='unknown' else 'Physiology labels use user-verified polarity. ')+'No swallowing, odor arrival or neural coupling is measured. Automatic screening requires manual waveform review.\n')
    print(json.dumps({k:meta[k] for k in ['total_trials','candidate_screen_pass_trials','physiological_timing_trials']}))
    return out


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('session');p.add_argument('config');p.add_argument('output')
    p.add_argument('--adc',default='board-ADC-00.dat')
    p.add_argument('--inhalation-sign',choices=['unknown','negative','positive'],default='unknown')
    p.add_argument('--polarity-note',default='')
    p.add_argument('--waveform-review-note',default='',help='Independent review evidence that candidate phases correspond to breathing')
    a = p.parse_args()
    audit(a.session,a.config,a.output,a.adc,a.inhalation_sign,a.polarity_note,a.waveform_review_note)


if __name__ == '__main__':
    main()
