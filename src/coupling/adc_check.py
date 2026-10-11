"""Read-only, fixed-segment ADC processing sensitivity diagnostic."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
from scipy.signal import butter, sosfiltfilt, welch
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from .events import read_binary, sample_rate
from .respiration import downsample_counts


def compare_segment(raw, fs):
    """Compare independent filter orderings; discard five seconds at either edge.

    Direct raw-rate filtering precedes integer decimation. The alternative uses
    the production polyphase anti-alias resampler before the same bandpass.
    """
    if fs < 1000 or fs % 1000:
        raise ValueError('Diagnostic requires raw rate to be an integer multiple of 1000 Hz')
    if len(raw) < 20*fs:
        raise ValueError('Need at least 20 seconds per segment')
    x = np.asarray(raw, dtype=float)-32768
    direct = sosfiltfilt(butter(4, [.5,20], fs=fs, btype='bandpass', output='sos'), x)[::int(fs/1000)]
    unfiltered = downsample_counts(raw, fs)
    after = sosfiltfilt(butter(4, [.5,20], fs=1000, btype='bandpass', output='sos'), unfiltered)
    keep = slice(5000,-5000)
    a,b,u = direct[keep],after[keep],unfiltered[keep]
    f,p = welch(u,fs=1000,nperseg=8000)
    band = p[(f>=.5)&(f<=20)]
    scale = np.std(a)
    metrics = dict(correlation=float(np.corrcoef(a,b)[0,1]) if scale>0 and np.std(b)>0 else None,
                   relative_rms_error=float(np.sqrt(np.mean((a-b)**2))/scale) if scale>0 else None,
                   spectral_flatness=float(np.exp(np.mean(np.log(np.maximum(band,np.finfo(float).tiny))))/np.mean(band)) if np.mean(band)>0 else None,
                   raw_min=int(np.min(raw)),raw_max=int(np.max(raw)),
                   rail_fraction=float(np.mean((raw==0)|(raw==65535))))
    return metrics,(u,a,b,f,p)


def audit(adc, header, output, seconds=30):
    adc,header,out = Path(adc).resolve(),Path(header).resolve(),Path(output).resolve()
    if out.exists():
        raise FileExistsError(out)
    if adc.parent == out or adc.parent in out.parents:
        raise ValueError('Output must be outside recording directory')
    fs = sample_rate(header)
    if not 20 <= seconds <= 120:
        raise ValueError('Segment duration must be 20–120 seconds')
    raw = read_binary(adc,'<u2')
    size = int(seconds*fs)
    if len(raw)<3*size:
        raise ValueError('Need three nonoverlapping segments; reduce --seconds (minimum 20)')
    starts = [0,(len(raw)-size)//2,len(raw)-size]
    rows=[]
    fig,axes=plt.subplots(3,3,figsize=(12,8),layout='constrained')
    for i,start in enumerate(starts):
        segment=np.array(raw[start:start+size])
        metrics,arrays=compare_segment(segment,fs)
        shuffled=segment.copy()
        np.random.default_rng(20261011+i).shuffle(shuffled)
        control,carrays=compare_segment(shuffled,fs)
        u,a,b,f,p=arrays
        row=dict(start_raw_sample=start,stop_raw_sample=start+size,start_s=start/fs,
                 sha256=hashlib.sha256(segment.tobytes()).hexdigest(),**metrics,
                 shuffled_control=control)
        rows.append(row)
        t=np.arange(min(5000,len(a)))/1000+start/fs+5
        axes[i,0].plot(t,u[:len(t)]-np.median(u),color='.7',lw=.5,label='Anti-aliased ADC')
        filtered_ax=axes[i,0].twinx()
        filtered_ax.plot(t,a[:len(t)],lw=.9,label='Filter at raw rate')
        filtered_ax.plot(t,b[:len(t)],ls='--',lw=.7,color='tab:orange',label='Resample then filter')
        limit=max(np.max(np.abs(a[:len(t)])),np.max(np.abs(carrays[1][:len(t)])),1)*1.1
        filtered_ax.set_ylim(-limit,limit)
        filtered_ax.set_ylabel('Filtered counts',fontsize=8)
        if i==0: filtered_ax.legend(fontsize=6,loc='upper right')
        axes[i,0].set(title=f'{["Early","Middle","Late"][i]} · r={metrics["correlation"]:.5f}' if metrics['correlation'] is not None else 'Constant segment',xlabel='Recording time (s)',ylabel='Centered counts')
        axes[i,1].semilogy(f,p,label='Recorded ADC')
        axes[i,1].semilogy(carrays[3],carrays[4],alpha=.7,label='Shuffled samples')
        axes[i,1].set(xlim=(.5,40),xlabel='Hz',ylabel='Counts²/Hz',title='Unfiltered, anti-aliased PSD')
        axes[i,2].plot(t,carrays[1][:len(t)],color='.4')
        axes[i,2].set_ylim(-limit,limit)
        axes[i,2].set(xlabel='Recording time (s)',ylabel='Centered counts',title='Shuffled raw samples, then filtered')
    axes[0,0].legend(fontsize=7);axes[0,1].legend(fontsize=7)
    fig.suptitle('ADC input audit · fixed segments · polarity and physiology unverified')
    out.mkdir(parents=True)
    for ext in ('png','pdf','svg'):
        fig.savefig(out/f'adc_segments.{ext}',dpi=250)
    plt.close(fig)
    meta=dict(adc=str(adc),header=str(header),dtype='<u2',raw_fs_hz=fs,total_samples=len(raw),
              source_size_bytes=adc.stat().st_size,source_mtime_ns=adc.stat().st_mtime_ns,
              segment_seconds=seconds,edge_discard_seconds=5,bandpass_hz=[.5,20],order=4,
              seed=20261011,segments=rows,
              interpretation='No breath validation, polarity inference, synchronization check or timing acceptance. Shuffle is a noise diagnostic, not a route-contrast surrogate.')
    (out/'adc_segments.json').write_text(json.dumps(meta,indent=2))
    print(json.dumps(rows,indent=2))
    return meta


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('adc');parser.add_argument('header');parser.add_argument('output')
    parser.add_argument('--seconds',type=int,default=30)
    args=parser.parse_args()
    audit(args.adc,args.header,args.output,args.seconds)

if __name__=='__main__':
    main()
