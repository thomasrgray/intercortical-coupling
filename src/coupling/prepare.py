"""Raw Intan DAT -> separate mapped/QC-masked LFP archive. Originals are read-only."""
import argparse
import csv
import json
import math
from pathlib import Path
import subprocess
import numpy as np
from scipy.signal import resample_poly, butter, sosfiltfilt, iirnotch, filtfilt, welch
from .events import sample_rate, read_binary, inventory, select_six, write_csv
from .qc import metrics, peer_flags


def check_time(path,n):
    a=read_binary(path,'<i4')
    if len(a)!=n:raise ValueError('Timestamp length differs from DIN length')
    for k in range(0,n-1,1_000_000):
        delta=np.diff(a[k:min(n,k+1_000_001)].astype(np.int64))
        if np.any(delta!=1):raise ValueError('Timestamp gaps/resets: split recording before processing')
    return int(a[0])


def load_channels(folder,config):
    if 'channels' in config:
        rows=config['channels']
    else:
        with (folder/config['electrode_layout']).open() as h:
            table=list(csv.DictReader(h))
        rows=[dict(filename=r['filename'],electrode_id=int(r['electrode_ind']),
                   region=r.get('original_CAR_group',r['CAR_group']).upper(),
                   reference_group=r['CAR_group']) for r in table]
    if not rows:raise ValueError('No configured amplifier channels')
    ids=[r['electrode_id'] for r in rows];names=[r['filename'] for r in rows]
    if len(ids)!=len(set(ids)) or len(names)!=len(set(names)):
        raise ValueError('Duplicate electrode ID or amplifier file')
    if set(r['region'] for r in rows)!={'GC','PC'}:raise ValueError('Explicit GC/PC mapping required')
    for r in rows:
        p=(folder/r['filename']).resolve()
        if p.parent!=folder.resolve() or not p.name.startswith('amp-'):
            raise ValueError('Amplifier files must be directly inside the source folder')
        if not p.is_file():raise FileNotFoundError(p)
        if not r.get('reference_group'):raise ValueError('Missing reference group')
    return rows


def resample_channel(raw,fs):
    """Bounded raw-rate memory; FIR context overlaps and output grid stays fixed."""
    if fs!=int(fs):raise ValueError('Integer raw sample rate required')
    fs=int(fs);g=math.gcd(fs,1000);up,down=1000//g,fs//g
    nout=(len(raw)*up+down-1)//down
    result=np.empty(nout,dtype=np.float64)
    # Whole seconds guarantee aligned polyphase output at every chunk boundary.
    for start in range(0,len(raw),30*fs):
        end=min(len(raw),start+30*fs);lo=max(0,start-fs);hi=min(len(raw),end+fs)
        x=np.asarray(raw[lo:hi],dtype=np.float64)*.195
        y=resample_poly(x,up,down)
        a=start*up//down;b=min(nout,(end*up+down-1)//down)
        offset=(start-lo)*up//down
        result[a:b]=y[offset:offset+b-a]
    return result


def filter_lfp(x,notch=True):
    # Resample_poly performs anti-alias filtering BEFORE decimation.
    y=sosfiltfilt(butter(4,[1,100],fs=1000,btype='bandpass',output='sos'),x)
    if notch:
        b,a=iirnotch(60,30,fs=1000);y=filtfilt(b,a,y)
    return y


def raw_quality(raw,fs):
    """Exact global integer histogram metrics and sustained-flat fraction."""
    hist=np.zeros(65536,dtype=np.int64);flat=0;carry=0;previous=None
    minimum=int(math.ceil(.01*fs))
    for k in range(0,len(raw),1_000_000):
        x=raw[k:k+1_000_000].astype(np.int32)
        hist+=np.bincount(x+32768,minlength=65536)
        edges=np.r_[0,np.flatnonzero(np.diff(x)!=0)+1,len(x)];runs=np.diff(edges)
        if previous is not None and x[0]==previous:runs[0]+=carry
        elif carry>=minimum:flat+=carry
        if len(runs)>1:flat+=int(runs[:-1][runs[:-1]>=minimum].sum())
        carry=int(runs[-1]);previous=x[-1]
    if carry>=minimum:flat+=carry
    values=np.arange(65536)-32768;n=len(raw)
    def median(v,h):
        order=np.argsort(v);cs=h[order].cumsum()
        return float((v[order[np.searchsorted(cs,(n-1)//2+1)]]+v[order[np.searchsorted(cs,n//2+1)]])/2)
    med=median(values,hist);dist=np.abs(values-med);mad=median(dist,hist)
    out=float(hist[dist>10*1.4826*mad].sum()/n)
    mean=float(np.sum(values*hist)/n);sd=float(np.sqrt(np.sum((values-mean)**2*hist)/n))
    return dict(median_uv=med*.195,sd_uv=sd*.195,flat_fraction=flat/n,outlier_fraction=out,
                bad=bool(mad==0 or flat/n>.01 or out>.01),rail_fraction=float((hist[0]+hist[-1])/n))


def prepare(folder,config_path,output):
    folder=Path(folder).resolve();out=Path(output).resolve()
    if out==folder or folder in out.parents:
        raise ValueError('Output must be outside the recording folder')
    if out.exists():raise FileExistsError(out)
    config=json.loads(Path(config_path).read_text())
    for key in ('animal_id','session_id','stimulus_by_din','reference_method'):
        if not config.get(key):raise ValueError(f'Missing configuration: {key}')
    method=config['reference_method']
    if method not in ('mean_car','none'):raise ValueError('reference_method must be mean_car or none')
    fs=sample_rate(folder/config.get('header','info.rhd'))
    if 'raw_fs_hz' in config and config['raw_fs_hz']!=fs:raise ValueError('Configured Fs disagrees with RHD')
    records,dins=inventory(folder,fs)
    selected=select_six({d:len(v[0]) for d,v in records.items()})
    labels=config['stimulus_by_din']
    from .power import STIMULI
    if set(labels)!=set(map(str,selected)) or set(labels.values())!=set(STIMULI):
        raise ValueError(f'Configured labels must cover exactly selected DINs {selected} and six stimulus names')
    n=next(iter(records.values()))[2]
    first_tick=check_time(folder/config.get('timestamps','time.dat'),n)
    channels=load_channels(folder,config)
    for row in channels:
        if (folder/row['filename']).stat().st_size!=2*n:raise ValueError('Amplifier/DIN length mismatch')
    events=[]
    for din in selected:
        for trial,(on,off) in enumerate(zip(*records[din][:2])):
            center=int(np.rint(on/fs*1000))
            events.append(dict(din=din,stimulus=labels[str(din)],original_trial_index=trial,
                onset_sample=int(on),offset_sample=int(off),onset_s=float(on/fs),center_lfp_sample=center,
                rounding_error_ms=float(center-on/fs*1000),
                in_bounds=bool(center-2000>=0 and center+5000<=math.ceil(n/fs*1000))))
    events.sort(key=lambda e:(e['onset_sample'],e['din']))
    for i,e in enumerate(events):e['event_id']=i
    out.mkdir(parents=True)
    write_csv(out/'din_inventory.csv',[dict(r,selected=r['din'] in selected) for r in dins])
    write_csv(out/'events.csv',events);write_csv(out/'channels.csv',channels)
    continuous=np.lib.format.open_memmap(out/'lfp_unreferenced.npy',mode='w+',dtype='float32',shape=(len(channels),math.ceil(n/fs*1000)))
    trials=np.lib.format.open_memmap(out/'lfp_trials.npy',mode='w+',dtype='float32',shape=(len(channels),len(events),7000))
    trials[:]=np.nan
    good_trial=np.zeros((len(channels),len(events)),bool);rows=[];trial_rows=[]
    for c,ch in enumerate(channels):
        print(f'Processing {ch["filename"]} ({c+1}/{len(channels)})',flush=True)
        raw=read_binary(folder/ch['filename'],'<i2');rq=raw_quality(raw,fs)
        y=filter_lfp(resample_channel(raw,fs),config.get('notch_60hz',True));continuous[c]=y
        q=metrics(y);q['bad']=q['bad'] or rq['bad']
        for j,e in enumerate(events):
            reasons=[]
            if not e['in_bounds']:reasons.append('boundary')
            else:
                k=e['center_lfp_sample'];trials[c,j]=y[k-2000:k+5000]
                a=int(e['onset_sample']-2*fs);b=int(e['onset_sample']+5*fs)
                rawq=metrics(np.asarray(raw[max(0,a):min(n,b)],dtype=float)*.195,fs)
                lq=metrics(trials[c,j])
                if rawq['bad']:reasons.append('raw_trial_qc')
                if lq['bad']:reasons.append('lfp_trial_qc')
            good_trial[c,j]=not reasons
            trial_rows.append(dict(electrode_id=ch['electrode_id'],event_id=j,initial_reasons=';'.join(reasons)))
        q['bad_trial_fraction']=float(1-good_trial[c].mean())
        rows.append(dict(**ch,**q,raw_flat_fraction=rq['flat_fraction'],raw_outlier_fraction=rq['outlier_fraction']))
    continuous.flush();trials.flush()
    reasons=peer_flags(rows,[r['region'] for r in channels])
    good=np.array([not x for x in reasons]);groups={}
    for c,ch in enumerate(channels):groups.setdefault(ch['reference_group'],[]).append(c)
    reference_members={}
    for group,indices in groups.items():
        if len({channels[c]['region'] for c in indices})!=1:raise ValueError('Reference group spans regions')
        members=[c for c in indices if good[c]];reference_members[group]=[channels[c]['electrode_id'] for c in members]
        if method=='mean_car' and len(members)<3:
            for c in indices:good[c]=False;reasons[c].append('fewer_than_3_reference_channels')
            continue
        if method=='mean_car':
            # Fixed membership; reject an event for the group if any reference donor fails.
            valid=good_trial[members].all(axis=0)
            for j in range(len(events)):
                ref=np.mean(np.asarray(trials[members,j]),axis=0)
                for c in indices:trials[c,j]-=ref
            good_trial[indices]&=valid
    for c,ch in enumerate(channels):
        for j in range(len(events)):
            if good_trial[c,j] and metrics(trials[c,j])['bad']:good_trial[c,j]=False
    for c,row in enumerate(rows):
        row.update(good_channel=bool(good[c]),reasons=';'.join(reasons[c]))
    for row in trial_rows:
        c=next(i for i,ch in enumerate(channels) if ch['electrode_id']==row['electrode_id']);j=row['event_id']
        row['accepted']=bool(good[c] and good_trial[c,j])
        row['final_reasons']=row['initial_reasons'] or ('' if row['accepted'] else ('channel_excluded' if not good[c] else 'reference_donor_or_post_car_qc'))
    write_csv(out/'channel_qc.csv',rows);write_csv(out/'trial_qc.csv',trial_rows)
    try:revision=subprocess.check_output(['git','rev-parse','HEAD'],cwd=Path(__file__).resolve().parent,text=True,stderr=subprocess.DEVNULL).strip()
    except (OSError,subprocess.CalledProcessError):revision='unknown'
    provenance=dict(source=str(folder),reference=dict(method=method,members=reference_members),
        filtering=dict(fs_hz=1000,bandpass_hz=[1,100],order=4,notch_60hz=config.get('notch_60hz',True),anti_alias='resample_poly'),
        qc_method='raw and LFP: finite/nonconstant, sustained flat >=10ms >1%, robust score >10 fraction >1%; leave-one-out same-region >3 scaled MAD; min 3 peers',
        event_alignment='TTL first-high sample, rounded to nearest 1ms; original indices in events.csv',
        raw_fs_hz=fs,first_timestamp=first_tick,revision=revision,configuration=config,
        files=[dict(name=p.name,size_bytes=p.stat().st_size,mtime_ns=p.stat().st_mtime_ns) for p in sorted(folder.glob('*.dat'))])
    (out/'preprocessing.json').write_text(json.dumps(provenance,indent=2))
    from .qc_plots import plot_qc
    plot_qc(out,rows,dins,selected,trials,good,good_trial)
    if not all(any(good[c] and channels[c]['region']==region for c in range(len(channels))) for region in ('GC','PC')):
        raise ValueError(f'No accepted channels in a region. Inspect {out}/channel_qc.csv; no cleaned archive written.')
    # Matched accepted epochs: measured reference effect, after identical filtering.
    # Channel-wise accumulation avoids another full trial-array copy.
    before_psd=np.full((len(channels),501),np.nan)
    after_psd=np.full_like(before_psd,np.nan)
    for c in np.flatnonzero(good):
        before=[];after=[]
        for j in np.flatnonzero(good_trial[c]):
            k=events[j]['center_lfp_sample']
            freq,p=welch(continuous[c,k-2000:k+5000],fs=1000,nperseg=1000)
            before.append(p)
            _,p=welch(trials[c,j],fs=1000,nperseg=1000);after.append(p)
        if before:
            before_psd[c]=np.mean(before,axis=0);after_psd[c]=np.mean(after,axis=0)
    np.savez_compressed(out/'cleaned_lfp.npz',lfp_uv=trials,time_s=np.arange(-2000,5000)/1000,
        stimulus=np.array([e['stimulus'] for e in events]),event_id=np.arange(len(events)),
        electrode_id=np.array([c['electrode_id'] for c in channels]),region=np.array([c['region'] for c in channels]),
        good_channel=good,good_trial=good_trial,session_id=config['session_id'],animal_id=config['animal_id'],
        provenance_json=json.dumps(provenance),reference_frequency_hz=np.arange(501),
        unreferenced_psd=before_psd,referenced_psd=after_psd)
    return out/'cleaned_lfp.npz'


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('session');p.add_argument('config');p.add_argument('output');p.add_argument('--plot-power',action='store_true')
    a=p.parse_args();clean=prepare(a.session,a.config,a.output)
    if a.plot_power:
        from .power import analyze
        analyze(clean,Path(a.output)/'power')
    print(clean)

if __name__=='__main__':main()
