"""Optional paired DIN/pressure-event LFP plots from continuous cached LFPs."""

from coupling.array_store import load_arrays, save_arrays
import argparse
import csv
import hashlib
import json
from pathlib import Path
import subprocess
import numpy as np
from scipy.signal import butter, sosfiltfilt
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from .events import inventory, select_six, sample_rate, read_binary, write_csv
from .prepare import check_time
from .respiration import detect_cycles, first_after, FS
from .power import normalize_power, STIMULI, REGIONS, SETTINGS
from .qc import metrics

NAME = 'LFPs matched to brathing'
TIME = np.arange(-2000, 5000)/FS


def read_rows(path):
    with Path(path).open() as stream:
        return list(csv.DictReader(stream))


def fingerprint(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024*1024), b''):
            h.update(chunk)
    return h.hexdigest()


def event_key(landmark, sign):
    if landmark in ('negative', 'positive'):
        return landmark+'_start_sample'
    if landmark not in ('inspiration', 'expiration'):
        raise ValueError('Landmark must be inspiration, expiration, negative or positive')
    # Audit cache is sign-inverted when positive inspiration was specified.
    if sign == 'unknown':
        raise ValueError('Physiological landmarks require a polarity-labelled audit; use negative/positive candidates otherwise')
    return ('negative' if landmark == 'inspiration' else 'positive')+'_start_sample'


def select_event(cycles, key, ttl, delay, window):
    anchor = (ttl+delay)*FS
    samples = [c[key] for c in cycles if c['screen_pass']]
    return first_after(samples, anchor, window)


def bounds(center, length, guard=0):
    return center-2000 >= guard and center+5000 <= length-guard


def referenced_epoch(continuous, center, channels, members, method, source, raw_fs):
    """Recheck shifted raw/unreferenced trials and fixed reference donors."""
    x = np.asarray(continuous[:, center-2000:center+5000], dtype=float).copy()
    good = np.zeros(len(channels), dtype=bool)
    for c, ch in enumerate(channels):
        raw = read_binary(source/ch['filename'], '<i2')
        lo, hi = int(round((center/FS-2)*raw_fs)), int(round((center/FS+5)*raw_fs))
        good[c] = not metrics(np.asarray(raw[lo:hi], dtype=float)*.195, raw_fs)['bad'] and not metrics(x[c])['bad']
    ids = {int(c['electrode_id']): i for i, c in enumerate(channels)}
    if method == 'mean_car':
        for group, donors in members.items():
            donor_indices = [ids[int(i)] for i in donors]
            targets = [i for i, ch in enumerate(channels) if ch['reference_group'] == group]
            if len(donor_indices) < 3:
                good[targets] = False
                continue
            valid = good[donor_indices].all()
            x[targets] -= x[donor_indices].mean(axis=0)
            good[targets] &= valid
    elif method != 'none':
        raise ValueError('Unknown cached reference method')
    good &= np.array([not metrics(y)['bad'] for y in x])
    return x, good


def validate_inputs(prepared, audit, rules):
    meta = json.loads((prepared/'preprocessing.json').read_text())
    resp = json.loads((audit/'respiration_provenance.json').read_text())
    source = Path(meta['source']).resolve()
    if Path(resp['source']).resolve() != source:
        raise ValueError('LFP and ADC caches belong to different recordings')
    if meta['raw_fs_hz'] != resp['raw_fs_hz'] or meta['first_timestamp'] != resp['first_timestamp']:
        raise ValueError('ADC/LFP clock metadata differ')
    if meta['configuration'] != resp['configuration']:
        raise ValueError('ADC/LFP session configurations differ')
    if meta['filtering']['fs_hz'] != FS:
        raise ValueError('Continuous LFP must be sampled at 1000 Hz')
    # Reject stale source files rather than trusting cache names alone.
    for record in meta['files']+resp['files']:
        stat = (source/record['name']).stat()
        if (stat.st_size, stat.st_mtime_ns) != (record['size_bytes'], record['mtime_ns']):
            raise ValueError('Source changed since cache creation: '+record['name'])
    fs = sample_rate(source/meta['configuration'].get('header', 'info.rhd'))
    if fs != meta['raw_fs_hz']:
        raise ValueError('RHD sample rate changed')
    records, _ = inventory(source, fs)
    selected = select_six({k: len(v[0]) for k, v in records.items()})
    if selected != resp['selected_dins']:
        raise ValueError('DIN selection differs from cache')
    n = next(iter(records.values()))[2]
    if check_time(source/meta['configuration'].get('timestamps','time.dat'), n) != meta['first_timestamp']:
        raise ValueError('Timestamp origin changed')
    events = read_rows(prepared/'events.csv')
    audit_events = read_rows(audit/'events_respiration.csv')
    if len(events) != len(audit_events):
        raise ValueError('ADC/LFP event count mismatch')
    for a, b in zip(events, audit_events):
        for key in ('event_id','din','stimulus','original_trial_index'):
            if a[key] != b[key]:
                raise ValueError('ADC/LFP event identity mismatch')
        if int(a['onset_sample']) != int(b['ttl_raw_sample']):
            raise ValueError('ADC/LFP event timing mismatch')
        if int(a['onset_sample']) != records[int(a['din'])][0][int(a['original_trial_index'])]:
            raise ValueError('Cached event differs from raw DIN')
    if rules.get('expected_audit_sign', resp['sign']) != resp['sign']:
        raise ValueError('Alignment config expects a different audit polarity convention')
    if set(rules['stimuli']) != set(STIMULI):
        raise ValueError('Alignment rules must explicitly cover all six stimuli')
    for rule in rules['stimuli'].values():
        event_key(rule['landmark'], resp['sign'])
        if not rule['route'] or not np.isfinite([rule['delay_s'],rule['window_s']]).all() or rule['delay_s'] < 0 or rule['window_s'] <= 0:
            raise ValueError('Need route, nonnegative delay and positive bounded search window')
    return meta, resp, source, events, fs, n


def diagnostics(out, pressure, cycles, rows, provisional, sign):
    for stim in STIMULI:
        subset = [r for r in rows if r['stimulus'] == stim]
        retained = [r for r in subset if r['retained']]
        fig, axes = plt.subplots(3,1,figsize=(9,9),layout='constrained')
        example = next((r for r in subset if r['selected_sample'] != '' and bounds(r['selected_sample'], len(pressure))), None)
        if example:
            k = example['selected_sample']; idx = np.arange(k-2000,k+5000)
            axes[0].plot(TIME,pressure[idx],lw=.7,color='black')
            for key,color,label in [('negative_start_sample','tab:blue','negative phase'),('positive_start_sample','tab:orange','positive phase')]:
                ss = np.array([c[key] for c in cycles if idx[0]<=c[key]<=idx[-1]],dtype=int)
                if sign != 'unknown':
                    label = 'inspiration' if key == 'negative_start_sample' else 'expiration'
                axes[0].scatter((ss-k)/FS,pressure[ss],s=12,color=color,label=label)
            axes[0].axvline(example['din_s']-k/FS,color='red',ls='--',label='DIN')
            axes[0].axvline(0,color='green',label='selected event')
            axes[0].legend(fontsize=8)
        axes[0].set(title='First in-bounds candidate (may be excluded)',ylabel='Filtered ADC counts')
        if retained:
            data = np.stack([pressure[r['selected_sample']-2000:r['selected_sample']+5000] for r in retained])
            # Display only; trial amplitudes remain in counts in cached pressure.
            scale = np.std(data,axis=1,keepdims=True)
            data = (data-data.mean(axis=1,keepdims=True))/np.maximum(scale,1e-12)
            axes[1].imshow(data,aspect='auto',origin='lower',extent=(-2,5,-.5,len(data)-.5),cmap='RdBu_r',vmin=-3,vmax=3)
            axes[1].scatter([-r['latency_s'] for r in retained],range(len(retained)),s=10,c='lime',label='DIN offset')
            axes[1].legend(fontsize=8)
        else:
            axes[1].text(.5,.5,'No retained trials',transform=axes[1].transAxes,ha='center')
        axes[1].set(ylabel='Retained trial order',xlabel='Time from selected event (s)',title='Pressure heatmap (display z-score; ±3)')
        for keep, label in [(True,'retained'),(False,'excluded')]:
            values=[r['latency_s'] for r in subset if r['retained']==keep and r['latency_s']!='']
            if values: axes[2].hist(values,bins=15,alpha=.6,label=label)
        axes[2].set(xlabel='Selected event minus DIN (s)',ylabel='Trials',title=f'Retained {len(retained)} / {len(subset)}; missing events {sum(r["selected_sample"]=="" for r in subset)}')
        if axes[2].has_data(): axes[2].legend()
        fig.suptitle(f'{stim}: respiration-based proxy'+(' — PROVISIONAL' if provisional else ' — reviewed pressure landmark'))
        fig.savefig(out/f'{stim}_respiration.png',dpi=160);plt.close(fig)


def run(prepared, audit, rule_file, output=None, clim=3.):
    prepared, audit = Path(prepared).resolve(), Path(audit).resolve()
    rules = json.loads(Path(rule_file).read_text())
    meta, resp, source, events, fs, n = validate_inputs(prepared,audit,rules)
    out = Path(output).resolve() if output else source/NAME
    if out.exists(): raise FileExistsError(out)
    if out.name != NAME: raise ValueError('Output folder must be named exactly '+NAME)
    if not np.isfinite(clim) or clim <= 0: raise ValueError('Invalid color limit')
    continuous = np.load(prepared/'lfp_unreferenced.npy',mmap_mode='r')
    din_trials = np.load(prepared/'lfp_trials.npy',mmap_mode='r')
    channels = read_rows(prepared/'channels.csv')
    cq = read_rows(prepared/'channel_qc.csv')
    tq = read_rows(prepared/'trial_qc.csv')
    ids = {int(ch['electrode_id']):i for i,ch in enumerate(channels)}
    good = np.array([next(r['good_channel']=='True' for r in cq if int(r['electrode_id'])==int(ch['electrode_id'])) for ch in channels])
    original = np.zeros((len(channels),len(events)),bool)
    for r in tq: original[ids[int(r['electrode_id'])],int(r['event_id'])] = r['accepted']=='True'
    if continuous.shape != (len(channels),int(np.ceil(n/fs*FS))) or din_trials.shape != (len(channels),len(events),7000):
        raise ValueError('Unexpected continuous/trial cache shape')
    with load_arrays(audit/'pressure_1khz.h5',allow_pickle=False) as d:
        if float(d['fs_hz']) != FS: raise ValueError('Pressure cache rate mismatch')
        pressure0 = d['unfiltered_counts'].astype(float)
    if len(pressure0) != continuous.shape[1]: raise ValueError('ADC/LFP resampled lengths differ')
    detector = rules['detector']
    low, high = detector['bandpass_hz']
    if not 0 < low < high < FS/2: raise ValueError('Invalid pressure bandpass')
    pressure = sosfiltfilt(butter(4,[low,high],fs=FS,btype='bandpass',output='sos'),pressure0)
    cycles = detect_cycles(pressure,detector['minimum_hz'],detector['maximum_hz'],detector['prominence_mad'])
    adc_record = next(r for r in resp['files'] if r['name'].startswith('board-ADC-'))
    raw_adc = read_binary(source/adc_record['name'],'<u2')
    provisional = resp['sign']=='unknown' or not resp.get('waveform_review_note','').strip() or not rules.get('definition_verified',False)
    rows, pair_rows, sums, counts = [], [], {}, {}
    f,t,z,_ = normalize_power(np.sin(2*np.pi*8*TIME),TIME)
    crop = (t>=-1-1e-9)&(t<=1.75+1e-9); grid=t[crop]
    for j,event in enumerate(events):
        stim=event['stimulus'];rule=rules['stimuli'][stim];din=float(event['onset_s'])
        key=event_key(rule['landmark'],resp['sign'])
        center=select_event(cycles,key,din,rule['delay_s'],rule['window_s'])
        reasons=[]; pressure_flags=[]; pair=np.zeros(len(channels),bool)
        din_center=int(event['center_lfp_sample'])
        if center is None: reasons.append('no_screened_event_in_search_window')
        elif not bounds(center,len(pressure),5000) or not bounds(din_center,len(pressure),5000): reasons.append('epoch_or_filter_edge_boundary')
        else:
            # Screen the union, including the newly exposed shifted tail.
            a=int(min(din_center,center)/FS*fs-2*fs);b=int(max(din_center,center)/FS*fs+5*fs)
            raw=np.asarray(raw_adc[a:b],dtype=float)
            q=metrics(raw,fs)
            if q['flat_fraction']>.01 or not q['nonconstant']:pressure_flags.append('pressure_flat')
            if np.mean((raw==0)|(raw==65535))>.01:pressure_flags.append('pressure_clipping')
            if q['outlier_fraction']>.01:pressure_flags.append('pressure_raw_artifact')
            if metrics(pressure[min(din_center,center)-2000:max(din_center,center)+5000])['bad']:pressure_flags.append('pressure_filtered_artifact')
            reasons.extend(pressure_flags)
            if not reasons:
                shifted, valid = referenced_epoch(continuous,center,channels,meta['reference']['members'],meta['reference']['method'],source,fs)
                pair=valid & original[:,j] & good
                if not pair.any(): reasons.append('no_matched_lfp_channels')
                for c in np.flatnonzero(pair):
                    reference=np.asarray(din_trials[c,j],dtype=float)
                    for mode,x in [('DIN',reference),('respiration',shifted[c])]:
                        _,_,z,frequency_valid=normalize_power(x,TIME,baseline_reference=reference)
                        group=(mode,c,stim)
                        if group not in sums:
                            sums[group]=np.zeros((len(f),len(grid)));counts[group]=np.zeros(len(f),int)
                        sums[group]+=np.nan_to_num(z[:,crop]);counts[group]+=frequency_valid
        row=dict(event_id=int(event['event_id']),stimulus=stim,route=rule['route'],din=int(event['din']),
                 original_trial_index=int(event['original_trial_index']),din_s=din,
                 selected_sample='' if center is None else center,selected_event_s='' if center is None else center/FS,
                 latency_s='' if center is None else center/FS-din,delivery_delay_s=rule['delay_s'],
                 landmark=rule['landmark'],alignment_rule='first screened onset at/after DIN + delay',search_window_s=rule['window_s'],
                 provisional=provisional,quality_flags=';'.join(pressure_flags),exclusion_reason=';'.join(reasons),
                 retained=not reasons,matched_channels=int(pair.sum()))
        rows.append(row)
        for c,ch in enumerate(channels):
            pair_rows.append(dict(event_id=int(event['event_id']),electrode_id=ch['electrode_id'],region=ch['region'],
                                  retained=bool(pair[c]),reason='' if pair[c] else (';'.join(reasons) or ('original_lfp_qc' if not original[c,j] else 'shifted_raw_donor_or_post_car_qc'))))
    out.mkdir(parents=True)
    write_csv(out/'trial_alignment.csv',rows);write_csv(out/'paired_channel_trials.csv',pair_rows)
    if cycles:write_csv(out/'pressure_cycles.csv',cycles)
    panels=np.full((2,2,6,len(f),len(grid)),np.nan);n_channels=np.zeros((2,6,len(f)),int)
    valid_counts=[]
    for ri,region in enumerate(REGIONS):
        for si,stim in enumerate(STIMULI):
            for mi,mode in enumerate(('DIN','respiration')):
                means=[]
                for c,ch in enumerate(channels):
                    key=(mode,c,stim)
                    if ch['region']!=region or key not in sums:continue
                    denom=counts[key];mean=np.full_like(sums[key],np.nan)
                    np.divide(sums[key],denom[:,None],out=mean,where=denom[:,None]>0);means.append(mean)
                    valid_counts.append(dict(mode=mode,region=region,stimulus=stim,electrode_id=ch['electrode_id'],valid_trials_per_frequency=denom.tolist()))
                if means:
                    stack=np.stack(means);denom=np.isfinite(stack).sum(0)
                    np.divide(np.nansum(stack,axis=0),denom,out=panels[mi,ri,si],where=denom>0)
                    n_channels[ri,si]=denom[:,0]
            fig,axes=plt.subplots(1,2,figsize=(11,4),layout='constrained',sharey=True)
            for mi,ax in enumerate(axes):
                mesh=ax.pcolormesh(grid*1000,f,panels[mi,ri,si],cmap='turbo',vmin=-clim,vmax=clim,shading='auto')
                ax.axvline(0,color='white',lw=1)
                ax.set(title=('DIN','Respiratory proxy')[mi],xlabel='Time from alignment (ms)',ylabel='Hz',ylim=(0,70))
                if not np.isfinite(panels[mi,ri,si]).any():ax.text(.5,.5,'No retained data',transform=ax.transAxes,ha='center')
            fig.colorbar(mesh,ax=axes,label='Power z-score: shared DIN baseline')
            fig.suptitle(f'{region} / {stim} · matched channel/trial pairs'+(' · PROVISIONAL' if provisional else ''))
            fig.savefig(out/f'{region}_{stim}_paired_power.png',dpi=160);plt.close(fig)
    diagnostics(out,pressure,cycles,rows,provisional,resp['sign'])
    save_arrays(out/'paired_power.h5',z_power=panels,alignment=['DIN','respiration'],region=REGIONS,stimulus=STIMULI,
                        frequency_hz=f,time_s=grid,n_channels=n_channels)
    try: revision=subprocess.check_output(['git','rev-parse','HEAD'],text=True,cwd=Path(__file__).parent,stderr=subprocess.DEVNULL).strip()
    except (OSError,subprocess.CalledProcessError): revision='unknown'
    provenance=dict(cache_files=[dict(path=str(p),size_bytes=p.stat().st_size,mtime_ns=p.stat().st_mtime_ns) for p in [prepared/'lfp_unreferenced.npy',prepared/'lfp_trials.npy']], rules=rules,revision=revision,provisional=provisional,lfp_preprocessing=meta,respiration_audit=resp,
        inputs={str(p):fingerprint(p) for p in [prepared/'preprocessing.json',prepared/'events.csv',prepared/'channel_qc.csv',prepared/'trial_qc.csv',audit/'respiration_provenance.json',audit/'pressure_1khz.h5',Path(rule_file)]},
        power_settings=SETTINGS,color_limit=clim,channel_counts=valid_counts,
        baseline='Same DIN-relative -1.5 to -0.5 s STFT-center moments for both alignments, per channel/trial/frequency',
        source_clock='Shared raw sample index; validated time.dat; times are seconds from recording start, not wall clock',
        limitation='Pressure zero-crossing proxy, not measured odor arrival. Zero-phase filtering can move waveform crossings. Fixed 5 s edge guard; interior artifact propagation and physiological detection accuracy unvalidated.',
        retained_trials=sum(r['retained'] for r in rows),total_trials=len(rows))
    (out/'provenance.json').write_text(json.dumps(provenance,indent=2))
    (out/'caption.txt').write_text('Paired DIN and respiratory-event power on identical channel/trial pairs. Both use the same DIN-prestimulus baseline moments; 500 ms Hann window, 1 ms hop, equal trials within channel then equal channels. Pressure landmarks and latency distributions are conditional on screening and the bounded first-onset rule. Events during an ongoing phase align to its next onset. Alignment of the pressure trace sharpens it by construction, including for filtered noise. No odor arrival, swallowing, neural locking or coupling is established. '+('PROVISIONAL: definition, polarity or waveform verification incomplete.' if provisional else 'Pressure landmark reviewed; odor arrival remains unmeasured.'))
    print(out)
    return out


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('prepared_session');p.add_argument('respiration_audit');p.add_argument('alignment_config')
    p.add_argument('--output',help='New folder named exactly "'+NAME+'"; default is inside source recording')
    p.add_argument('--clim',type=float,default=3.)
    a=p.parse_args();run(a.prepared_session,a.respiration_audit,a.alignment_config,a.output,a.clim)


if __name__=='__main__':main()
