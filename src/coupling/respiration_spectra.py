"""DIN-aligned respiration power and channel-wise LFP/pressure coherence."""

from coupling.array_store import load_arrays, save_arrays
import argparse
import json
from pathlib import Path
import numpy as np
from scipy import signal
from scipy.integrate import trapezoid
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from .respiration_lfp import read_rows, fingerprint
from .events import write_csv

FS = 1000
PAIRS = {'BA': ('oBA', 'BA'), 'IA': ('oIA', 'IA')}


def spectra(x, y, fs=FS, window_s=1., hop_s=.1):
    """Average cross/auto spectra over trials before taking the coherence ratio.

    Input: matched trial x time arrays. No concatenation across trial boundaries.
    Each sliding window has one tapered Fourier coefficient per trial; at least
    three independent trials are required. Time resolution is window_s, not hop_s.
    """
    x, y = np.asarray(x, float), np.asarray(y, float)
    n, hop = int(round(fs*window_s)), int(round(fs*hop_s))
    if x.shape != y.shape or x.ndim != 2 or len(x) < 3:
        raise ValueError('Need >=3 matched trials with identical trial/time shapes')
    if n < 4 or hop < 1 or x.shape[1] < n or not np.isfinite(x).all() or not np.isfinite(y).all():
        raise ValueError('Invalid windows or nonfinite data')
    if np.any(x.std(axis=1) == 0) or np.any(y.std(axis=1) == 0):
        raise ValueError('Constant signal trial')
    taper = signal.windows.hann(n, sym=False)
    def transform(a):
        a = np.lib.stride_tricks.sliding_window_view(a, n, axis=-1)[:, ::hop]
        return np.fft.rfft((a-a.mean(axis=-1, keepdims=True))*taper, axis=-1)
    a, b = transform(x), transform(y)
    p, q = np.mean(abs(a)**2, axis=0), np.mean(abs(b)**2, axis=0)
    cross = np.mean(a*b.conj(), axis=0)
    coh = np.full_like(p, np.nan)
    np.divide(abs(cross)**2, p*q, out=coh, where=p*q>0)
    # One-sided PSD density, preserving physical input units squared / Hz.
    scale = np.full(n//2+1, 2/(fs*np.sum(taper**2)))
    scale[0] /= 2
    if n % 2 == 0: scale[-1] /= 2
    return np.fft.rfftfreq(n, 1/fs), (np.arange(p.shape[0])*hop+n/2)/fs, p*scale, q*scale, np.clip(coh, 0, 1)


def load_inputs(prepared, audit):
    prepared, audit = Path(prepared), Path(audit)
    meta = json.loads((prepared/'preprocessing.json').read_text())
    resp = json.loads((audit/'respiration_provenance.json').read_text())
    for key in ('raw_fs_hz', 'first_timestamp', 'configuration'):
        if meta[key] != resp[key]: raise ValueError('ADC/LFP mismatch: '+key)
    if Path(meta['source']).resolve() != Path(resp['source']).resolve():
        raise ValueError('ADC/LFP source mismatch')
    # Check source freshness, including ADC, clocks and DINs.
    for row in meta['files']+resp['files']:
        stat = (Path(meta['source'])/row['name']).stat()
        if (stat.st_size, stat.st_mtime_ns) != (row['size_bytes'], row['mtime_ns']):
            raise ValueError('Stale cache: '+row['name'])
    events, revents = read_rows(prepared/'events.csv'), read_rows(audit/'events_respiration.csv')
    if len(events) != len(revents): raise ValueError('Event count mismatch')
    for e, r in zip(events, revents):
        for key in ('event_id', 'din', 'stimulus', 'original_trial_index'):
            if e[key] != r[key]: raise ValueError('Event identity mismatch: '+key)
        if int(e['onset_sample']) != int(r['ttl_raw_sample']) or int(e['center_lfp_sample']) != int(r['ttl_lfp_sample']):
            raise ValueError('ADC/LFP onset mismatch')
    with load_arrays(prepared/'cleaned_lfp.h5', allow_pickle=False) as z:
        if json.loads(str(z['provenance_json'])) != meta:
            raise ValueError('Cleaned archive preprocessing provenance mismatch')
        data = {k:z[k] for k in ('lfp_uv','time_s','event_id','stimulus','region','good_channel','good_trial')}
    t, lfp = data['time_s'], data['lfp_uv']
    if lfp.shape[1:] != (len(events), len(t)) or not np.allclose(t, np.arange(-2000,5000)/FS):
        raise ValueError('Unexpected LFP trial/time axes')
    if not np.array_equal(data['event_id'], [int(e['event_id']) for e in events]) or not np.array_equal(data['stimulus'], [e['stimulus'] for e in events]):
        raise ValueError('Cleaned archive event order mismatch')
    if data['good_trial'].shape != lfp.shape[:2] or len(data['region']) != len(lfp) or len(data['good_channel']) != len(lfp):
        raise ValueError('Invalid QC dimensions')
    with load_arrays(audit/'pressure_1khz.h5', allow_pickle=False) as z:
        if float(z['fs_hz']) != FS: raise ValueError('Pressure sample rate mismatch')
        pressure, filtered = z['unfiltered_counts'].astype(float), z['filtered_counts'].astype(float)
    epochs = np.full((len(events), len(t)), np.nan)
    display = epochs.copy()
    accepted = np.zeros(len(events), bool)
    rows = []
    for j, (e,r) in enumerate(zip(events,revents)):
        k = int(e['center_lfp_sample'])
        # PSD/coherence do not require identifiable breaths; retain only signal QC reasons.
        reasons = [s for s in r['timing_qc_reasons'].split(';') if s and not s.startswith('no_screened_')]
        if k-2000 < 0 or k+5000 > len(pressure): reasons.append('boundary')
        if not reasons:
            epochs[j], display[j] = pressure[k-2000:k+5000], filtered[k-2000:k+5000]
            if not np.isfinite(epochs[j]).all() or np.std(epochs[j]) == 0: reasons.append('invalid_pressure')
        accepted[j] = not reasons
        rows.append(dict(event_id=e['event_id'], stimulus=e['stimulus'], accepted=bool(accepted[j]), reasons=';'.join(reasons)))
    return data, epochs, display, accepted, rows, meta, resp


def run(prepared, audit, output, window_s=1., hop_s=.1):
    data, pressure, display, accepted, rows, meta, resp = load_inputs(prepared, audit)
    out = Path(output)
    if out.exists(): raise FileExistsError(out)
    if not np.isfinite([window_s,hop_s]).all() or not 0 < hop_s <= window_s <= 2:
        raise ValueError('Need 0 < hop <= window <= 2 seconds')
    out.mkdir(parents=True)
    write_csv(out/'pressure_trial_qc.csv', rows)
    t = data['time_s']; channel_rows = []; results = {}; csv_rows = []
    for stim in dict.fromkeys(data['stimulus']):
        mask = accepted & (data['stimulus'] == stim)
        if mask.sum() >= 3:
            f, centers, p, _, _ = spectra(pressure[mask], pressure[mask], window_s=window_s, hop_s=hop_s)
            results[(str(stim), 'pressure')] = p.T
        for region in ('GC','PC'):
            cs, ps = [], []
            for c in np.flatnonzero(data['good_channel'] & (data['region'] == region)):
                keep = mask & data['good_trial'][c]
                valid = np.isfinite(data['lfp_uv'][c]).all(axis=-1) & (np.nanstd(data['lfp_uv'][c],axis=-1)>0)
                keep &= valid
                channel_rows.append(dict(stimulus=stim,region=region,channel_index=int(c),n_trials=int(keep.sum()),retained=bool(keep.sum()>=3)))
                if keep.sum() < 3: continue
                f, centers, _, p, coh = spectra(pressure[keep],data['lfp_uv'][c,keep],window_s=window_s,hop_s=hop_s)
                cs.append(coh.T); ps.append(p.T)
            if cs:
                # Average channel coherence, never the CAR-cancelled regional waveform.
                results[(str(stim),region)] = np.mean(cs,axis=0)
                results[(str(stim),region+'_power')] = np.mean(ps,axis=0)
    write_csv(out/'channel_trial_counts.csv',channel_rows)
    if not results: raise ValueError('No condition has >=3 valid trials; inspect QC CSVs')
    time = centers+t[0]; band = (f>=0)&(f<=12); post=(time>=0)&(time<=2)
    archive = dict(frequency_hz=f,time_s=time)
    for (stim,kind),arr in results.items():
        archive[stim+'__'+kind] = arr
        for fi in np.flatnonzero(band):
            csv_rows.append(dict(stimulus=stim,measure=kind,frequency_hz=float(f[fi]),post_0_2s_mean=float(np.mean(arr[fi,post]))))
    save_arrays(out/'spectra.h5',**archive)
    write_csv(out/'spectra_0_2s.csv',csv_rows)
    for odor, stimuli in PAIRS.items():
        fig, axes = plt.subplots(4,2,figsize=(11,12),layout='constrained',squeeze=False)
        for col, (route,stim) in enumerate(zip(('Ortho','Retro'),stimuli)):
            mask=accepted & (data['stimulus']==stim)
            ax=axes[0,col]
            if mask.any():
                j=np.flatnonzero(mask)[0]; z=display[j]; z=(z-z.mean())/z.std()
                ax.plot(t,z,color='0.35',lw=.9)
                y=-3.5;ax.plot([1.9,2.1],[y,y],color='black',lw=2);ax.text(2,y-.15,'200 ms',ha='center',va='top',fontsize=8)
                ax.set_ylim(min(-4.3,z.min()-.2),max(3,z.max()+.2))
            else: ax.text(.5,.5,'No accepted pressure trials',transform=ax.transAxes,ha='center')
            ax.axvline(0,color='black',lw=.8); ax.set(xlim=(-.5,2.5),ylabel='Pressure (z)',title=f'A  {route} · {stim} · example trial',xlabel='Time from DIN (s)')
            ax=axes[1,col]
            for kind,label,color in [('pressure','Respiration','black'),('GC_power','GC LFP','tab:blue'),('PC_power','PC LFP','tab:orange')]:
                arr=results.get((stim,kind))
                if arr is None:continue
                values=arr[band][:,post].mean(axis=1)
                # Explicit shape normalization; raw PSDs are exported.
                area=trapezoid(values,f[band])
                if area>0:ax.plot(f[band],values/area,label=label,color=color)
            ax.set(xlim=(0,12),xlabel='Frequency (Hz)',ylabel='Unit-area PSD',title='B  Power · window centers 0–2 s')
            if ax.has_data():ax.legend(fontsize=8)
            for row,region in enumerate(('GC','PC'),2):
                ax=axes[row,col]; arr=results.get((stim,region))
                if arr is not None:
                    mesh=ax.pcolormesh(time,f[band],arr[band],shading='auto',vmin=0,vmax=1,cmap='viridis')
                    fig.colorbar(mesh,ax=ax,label='Magnitude-squared coherence')
                else:ax.text(.5,.5,'Need ≥3 matched trials/channel',transform=ax.transAxes,ha='center')
                ax.axvline(0,color='white',lw=.8);ax.set(xlim=(-.5,2.5),ylim=(0,12),xlabel='Time from DIN (s)',ylabel='Frequency (Hz)',title=f'C  {region} LFP ↔ respiration')
        fig.suptitle(f'{odor}: single-session descriptive spectra | {window_s:g} s windows | pressure polarity: {resp["sign"]}')
        for ext in ('png','pdf'):fig.savefig(out/f'{odor}_respiration_lfp.{ext}',dpi=180)
        plt.close(fig)
    inputs=[Path(prepared)/n for n in ('preprocessing.json','events.csv','cleaned_lfp.h5')]+[Path(audit)/n for n in ('respiration_provenance.json','events_respiration.csv','pressure_1khz.h5')]
    (out/'provenance.json').write_text(json.dumps(dict(lfp=meta,respiration=resp,window_s=window_s,hop_s=hop_s,fs_hz=FS,
        hashes={str(p.resolve()):fingerprint(p) for p in inputs},minimum_trials=3,
        estimator='Hann-tapered trial-ensemble cross/auto spectra; magnitude-squared ratio then equal-channel mean',
        limitations='Descriptive single-session estimates, no significance inference. Stimulus-evoked common signals can contribute. Pressure PSD uses unfiltered anti-aliased ADC counts; display uses 0.5–20 Hz. LFP cache is 1–100 Hz: sub-1-Hz content attenuated. Pressure polarity does not change PSD/coherence; inhale labels require verification. DIN is not measured odor arrival. No confidence intervals or pooled-animal inference.'),indent=2))
    print(out)
    return out


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('prepared');parser.add_argument('audit');parser.add_argument('output')
    parser.add_argument('--window',type=float,default=1.);parser.add_argument('--hop',type=float,default=.1)
    a=parser.parse_args();run(a.prepared,a.audit,a.output,a.window,a.hop)


if __name__=='__main__':main()
