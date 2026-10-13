"""Figure 6: session-level QC, power, coherence and bivariate spectral Granger.

Electrode pairs and trials are technical observations, not biological replicates.
"""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import subprocess

import numpy as np
from scipy import signal, stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from .power import validate, STIMULI, REGIONS

WINDOWS = {'pre': (-1.75, -.25), 'post': (.25, 1.75)}
DEFAULT_PAIRS = [['BA', 'oBA'], ['IA', 'oIA']]


def mean_finite(x, axis=0):
    x = np.asarray(x)
    n = np.isfinite(x).sum(axis=axis)
    return np.divide(np.nansum(x, axis=axis), n,
                     out=np.full(np.shape(n), np.nan), where=n > 0)


def pooled_coherence(x, y, fs=1000, nperseg=500):
    """Pool auto/cross spectra over trials BEFORE taking their ratio."""
    if x.shape != y.shape or x.ndim != 2 or x.shape[-1] < nperseg:
        raise ValueError('Coherence requires matching trial × sample arrays')
    if not np.isfinite(x).all() or not np.isfinite(y).all():
        raise ValueError('Nonfinite accepted coherence data')
    kw = dict(fs=fs, window='hann', nperseg=nperseg, noverlap=nperseg//2,
              detrend='constant', axis=-1)
    f, xx = signal.welch(x, **kw)
    _, yy = signal.welch(y, **kw)
    _, xy = signal.csd(x, y, **kw)
    den = xx.mean(0) * yy.mean(0)
    c = np.divide(abs(xy.mean(0))**2, den, out=np.full_like(den, np.nan), where=den > 0)
    return f, np.clip(c, 0, 1)


def spectral_granger(a, sigma, f, fs):
    """Geweke bivariate spectrum, output [GC→PC, PC→GC, frequency].

    a[lag,target,source]; subtract source's conditional innovation contribution
    from target power, retaining instantaneous residual covariance in H Sigma H*.
    """
    phase = np.exp(-2j*np.pi*f[:, None]*np.arange(1, len(a)+1)/fs)
    h = np.linalg.inv(np.eye(2)[None] - np.einsum('fl,lij->fij', phase, a))
    s = h @ sigma @ h.conj().transpose(0, 2, 1)
    result = []
    for source, target in [(0, 1), (1, 0)]:
        total = s[:, target, target].real
        conditional = sigma[source, source] - sigma[source, target]**2/sigma[target, target]
        intrinsic = total - conditional * abs(h[:, target, source])**2
        if np.any(intrinsic <= 0) or np.any(total <= 0):
            raise ValueError('Nonpositive intrinsic spectral power')
        result.append(np.log(total/intrinsic))
    result = np.array(result)
    if not np.isfinite(result).all() or result.min() < -1e-8:
        raise ValueError('Invalid spectral Granger estimate')
    return np.maximum(result, 0)


def fit_granger(x, fs=250, orders=(2, 5, 10, 15, 20), whiteness_alpha=.01):
    """OLS ensemble VAR: no lag or residual diagnostic crosses trial boundaries.

    All BIC candidates use the same target samples (start at maximum order).
    Diagnostics are screening rules, not a proof of stationarity/model adequacy.
    """
    x = np.asarray(x, float)
    if x.ndim != 3 or x.shape[-1] != 2 or not np.isfinite(x).all():
        raise ValueError('VAR requires finite trial × sample × 2 data')
    orders = tuple(sorted(set(orders)))
    if not orders or orders[0] < 1 or x.shape[1] <= 3*orders[-1]:
        raise ValueError('Insufficient samples for requested VAR orders')
    scale = x.std(axis=(0, 1))
    if np.any(scale <= np.finfo(float).eps):
        raise ValueError('Constant VAR channel')
    half = x.shape[1]//2
    v1, v2 = x[:, :half].var(axis=1), x[:, half:].var(axis=1)
    variance_ratio = float(np.max(np.median(np.maximum(v1, v2)/np.maximum(np.minimum(v1, v2), 1e-15), axis=0)))
    mean_shift = float(np.max(np.median(abs(x[:, :half].mean(1)-x[:, half:].mean(1))/scale, axis=0)))
    x = (x-x.mean(1, keepdims=True))/scale
    start = orders[-1]
    y = x[:, start:].reshape(-1, 2)
    designs = [x[:, start-lag:x.shape[1]-lag] for lag in range(1, start+1)]
    full = np.concatenate(designs, axis=-1).reshape(len(y), -1)
    fits = []
    for p in orders:
        z = full[:, :2*p]
        if len(y) <= 10*z.shape[1]:
            continue
        coef, _, rank, sv = np.linalg.lstsq(z, y, rcond=None)
        e = y-z@coef
        sigma = e.T@e/len(e)
        sign, logdet = np.linalg.slogdet(sigma)
        if rank != z.shape[1] or sign <= 0:
            continue
        bic = float(logdet + np.log(len(y))*4*p/len(y))
        fits.append((bic, p, coef, sigma, e, float(sv[0]/sv[-1])))
    if not fits:
        raise ValueError('No full-rank adequately sampled VAR candidate')
    bic, p, coef, sigma, e, condition = min(fits, key=lambda item: item[0])
    a = coef.reshape(p, 2, 2).transpose(0, 2, 1)
    companion = np.zeros((2*p, 2*p))
    companion[:2] = np.concatenate(a, axis=1)
    companion[2:, :-2] = np.eye(2*(p-1))
    radius = float(max(abs(np.linalg.eigvals(companion))))
    residual = e.reshape(x.shape[0], x.shape[1]-start, 2)
    lags = min(max(2*p, 20), residual.shape[1]//3)
    # Multivariate residual portmanteau, asymptotic chi-square approximation.
    inv = np.linalg.inv(sigma)
    q = 0.
    max_corr = 0.
    for lag in range(1, lags+1):
        left, right = residual[:, lag:].reshape(-1, 2), residual[:, :-lag].reshape(-1, 2)
        cov = left.T@right/len(left)
        q += len(left)*np.trace(cov.T@inv@cov@inv)
        max_corr = max(max_corr, float(np.max(abs(cov)/np.sqrt(np.outer(np.diag(sigma), np.diag(sigma))))))
    white_p = float(stats.chi2.sf(q, max(1, 4*(lags-p))))
    reasons = []
    if radius >= 1: reasons.append('unstable')
    if condition > 1e8: reasons.append('ill_conditioned')
    if np.linalg.cond(sigma) > 1e8: reasons.append('singular_innovations')
    if white_p < whiteness_alpha: reasons.append('residual_autocorrelation')
    if variance_ratio > 4 or mean_shift > 1: reasons.append('nonstationarity_screen')
    diag = dict(order=p, bic=bic, spectral_radius=radius, design_condition=condition,
                residual_whiteness_p=white_p, residual_max_abs_correlation=max_corr,
                variance_ratio=variance_ratio, mean_shift=mean_shift,
                n_observations=len(y), passed=not reasons, reasons=';'.join(reasons),
                order_at_upper_bound=p == orders[-1])
    return a, sigma, diag


def save_csv(path, rows, fields=None):
    with Path(path).open('w', newline='') as h:
        writer = csv.DictWriter(h, fieldnames=fields or list(rows[0]))
        writer.writeheader(); writer.writerows(rows)


def save_fig(fig, out, name):
    fig.tight_layout()
    for ext in ('png', 'pdf'):
        fig.savefig(out/f'{name}.{ext}', dpi=160)
    plt.close(fig)


def matched_contrasts(values, stimuli, pairs):
    """Average within-electrode-pair ortho-minus-retro differences.

    Missing estimates must never change pair membership between routes.
    values: electrode pair × stimulus × ... × frequency.
    """
    contrasts, counts = [], []
    for retro, ortho in pairs:
        a, b = values[:, stimuli.index(retro)], values[:, stimuli.index(ortho)]
        valid = np.isfinite(a) & np.isfinite(b)
        contrasts.append(mean_finite(np.where(valid, b-a, np.nan)))
        counts.append(valid.sum(axis=0))
    return np.array(contrasts), np.array(counts)


def analyze(source, output, max_channels=5, min_trials=5, matched_pairs=None,
            orders=(2, 5, 10, 15, 20), remove_evoked=False):
    source, out = Path(source), Path(output)
    if out.exists(): raise FileExistsError(out)
    if max_channels < 1 or min_trials < 3: raise ValueError('Need >=1 channel and >=3 trials')
    if not orders or any(not isinstance(p, int) or p < 1 for p in orders):
        raise ValueError('VAR orders must be positive integers')
    pairs = DEFAULT_PAIRS if matched_pairs is None else matched_pairs
    if not pairs or any(len(p) != 2 or p[0] not in STIMULI or p[1] not in STIMULI or p[0] == p[1] for p in pairs):
        raise ValueError('Matched pairs must be explicit [retro, ortho] stimulus names')
    with np.load(source, allow_pickle=False) as d:
        provenance = validate(d)
        data = {k: d[k] for k in d.files}
    x, t = data['lfp_uv'], data['time_s']
    if not np.allclose(np.diff(t), .001, rtol=0, atol=1e-9):
        raise ValueError('Figure 6 requires 1000 Hz cleaned input')
    if t[0] > -1.75 or t[-1] < 1.75: raise ValueError('Need -1.75 to +1.75 s coverage')
    if 'event_id' not in data or len(np.unique(data['event_id'])) != x.shape[1] or data['event_id'].shape != (x.shape[1],):
        raise ValueError('Unique original event_id values are required')
    accepted = data['good_channel'][:, None] & data['good_trial']
    if not np.isfinite(x[accepted]).all(): raise ValueError('Nonfinite accepted data')
    selected = [sorted(np.flatnonzero(data['good_channel'] & (data['region'] == r)),
                       key=lambda c: data['electrode_id'][c])[:max_channels] for r in REGIONS]
    if any(not cs for cs in selected): raise ValueError('No accepted channels in a region')
    electrode_pairs = [(a, b) for a in selected[0] for b in selected[1]]
    f = np.arange(0, 81, 2.)
    gf = np.arange(0, 81, 1.)
    coh = np.full((len(electrode_pairs), 6, 2, len(f)), np.nan)
    shuffled = np.full_like(coh, np.nan)
    gc = np.full((len(electrode_pairs), 6, 2, 2, len(gf)), np.nan)
    power = np.full((x.shape[0], 6, 2, len(f)), np.nan)
    rows, diagnostics, retention = [], [], []
    masks = [(t >= lo-1e-9) & (t < hi-1e-9) for lo, hi in WINDOWS.values()]
    # Resample whole epochs before cropping, keeping filter edges away from windows.
    tx = t[0]+np.arange((len(t)+3)//4)/250
    gmasks = [(tx >= lo-1e-9) & (tx < hi-1e-9) for lo, hi in WINDOWS.values()]
    gx = {c: signal.resample_poly(x[c, accepted[c]], 1, 4, axis=-1) for cs in selected for c in cs}
    event_lookup = {c: {event: k for k, event in enumerate(np.flatnonzero(accepted[c]))} for c in gx}
    for c in range(len(x)):
        for s, stim in enumerate(STIMULI):
            ids = np.flatnonzero(accepted[c] & (data['stimulus'] == stim))
            retention.append(dict(electrode_id=str(data['electrode_id'][c]), region=str(data['region'][c]),
                stimulus=stim, good_channel=bool(data['good_channel'][c]),
                accepted_trials=len(ids), total_trials=int(sum(data['stimulus'] == stim)),
                selected_for_coupling=c in gx))
            if len(ids):
                for w, mask in enumerate(masks):
                    _, psd = signal.welch(x[c, ids][:, mask], fs=1000, nperseg=500, noverlap=250)
                    power[c, s, w] = psd.mean(0)[:len(f)]
    for pi, (a, b) in enumerate(electrode_pairs):
        print(f'Coupling pair {pi+1}/{len(electrode_pairs)}', flush=True)
        for si, stim in enumerate(STIMULI):
            ids = np.flatnonzero(accepted[a] & accepted[b] & (data['stimulus'] == stim))
            base = dict(gc_electrode=str(data['electrode_id'][a]), pc_electrode=str(data['electrode_id'][b]), stimulus=stim)
            for j in np.flatnonzero(data['stimulus'] == stim):
                rows.append(dict(**base, event_id=str(data['event_id'][j]), accepted=bool(j in ids)))
            for wi, (window, mask) in enumerate(zip(WINDOWS, masks)):
                diag = dict(**base, window=window, n_trials=len(ids), order='', bic='', spectral_radius='',
                    design_condition='', residual_whiteness_p='', residual_max_abs_correlation='',
                    variance_ratio='', mean_shift='', n_observations='', passed=False,
                    reasons='insufficient_trials', order_at_upper_bound='')
                if len(ids) >= min_trials:
                    xx, yy = x[a, ids][:, mask].copy(), x[b, ids][:, mask].copy()
                    if remove_evoked:
                        xx -= xx.mean(0); yy -= yy.mean(0)
                    _, cc = pooled_coherence(xx, yy)
                    coh[pi, si, wi] = cc[:len(f)]
                    _, cc = pooled_coherence(xx, np.roll(yy, 1, axis=0))
                    shuffled[pi, si, wi] = cc[:len(f)]
                    z = np.stack([gx[c][[event_lookup[c][j] for j in ids]][:, gmasks[wi]] for c in (a, b)], axis=-1)
                    if remove_evoked: z -= z.mean(0, keepdims=True)
                    try:
                        ar, sigma, fit = fit_granger(z, orders=orders)
                        diag.update(fit)
                        if fit['passed']:
                            gc[pi, si, wi] = spectral_granger(ar, sigma, gf, 250)
                    except (ValueError, np.linalg.LinAlgError) as exc:
                        diag.update(passed=False, reasons=str(exc))
                diagnostics.append(diag)
    cdelta, cn = matched_contrasts(coh, list(STIMULI), pairs)
    gdelta, gn = matched_contrasts(gc, list(STIMULI), pairs)
    out.mkdir(parents=True)
    save_csv(out/'retention.csv', retention)
    save_csv(out/'pair_event_qc.csv', rows)
    save_csv(out/'granger_diagnostics.csv', diagnostics)
    np.savez_compressed(out/'plotted_values.npz', frequency_hz=f, granger_frequency_hz=gf,
        power_uv2_per_hz=power, coherence=coh, trial_shift_coherence=shuffled, granger=gc,
        coherence_route_difference=cdelta, coherence_route_pair_counts=cn,
        granger_route_difference=gdelta, granger_route_pair_counts=gn,
        electrode_pairs=data['electrode_id'][np.array(electrode_pairs)], electrode_id=data['electrode_id'],
        stimulus=STIMULI, region=data['region'], windows=list(WINDOWS), directions=['GC_to_PC','PC_to_GC'],
        matched_pairs=np.array(pairs), event_id=data['event_id'], accepted=accepted,
        animal_id=data['animal_id'], session_id=data['session_id'])
    title = f"{data['session_id'].item()} | {data['animal_id'].item()}"
    if remove_evoked: title += ' | evoked-subtracted coupling'
    if 'synthetic' in json.dumps(provenance).lower(): title = 'SYNTHETIC VALIDATION — ' + title
    fig, axes = plt.subplots(2, 2, figsize=(13, 8))
    axes[0, 0].imshow(accepted, aspect='auto', interpolation='nearest', cmap='Greys', vmin=0, vmax=1)
    axes[0, 0].set(title='A | Accepted channel-events (black)', xlabel='Event position (IDs in NPZ/CSV)', ylabel='Channel position')
    for r in REGIONS:
        counts = [sum(row['accepted_trials'] for row in retention if row['region']==r and row['stimulus']==s) for s in STIMULI]
        axes[0, 1].plot(STIMULI, counts, 'o-', label=r)
    axes[0, 1].set(title='Retention by stimulus', ylabel='Accepted channel-events'); axes[0, 1].legend()
    diagnostic_values = {}
    freq = data.get('reference_frequency_hz')
    if freq is not None:
        diagnostic_values.update(reference_frequency_hz=freq, unreferenced_psd=data['unreferenced_psd'], referenced_psd=data['referenced_psd'])
        for r in REGIONS:
            cs = (data['region']==r) & data['good_channel']
            for key, style in [('unreferenced_psd','--'), ('referenced_psd','-')]:
                axes[1, 0].semilogy(freq, mean_finite(data[key][cs]), style, label=f'{r} {key}')
        axes[1, 0].set(xlim=(1,100), title='Measured reference effect (same accepted trials)', xlabel='Hz', ylabel='µV²/Hz')
        axes[1, 0].legend(fontsize=7)
    else:
        axes[1, 0].text(.03, .5, 'Reference PSD diagnostics unavailable.\nRe-run coupling-prepare to save them.\nNo before/after effect inferred.', transform=axes[1, 0].transAxes)
        axes[1, 0].set_axis_off()
    filtering = provenance['filtering']
    if isinstance(filtering, dict) and filtering.get('bandpass_hz') == [1, 100] and filtering.get('order') == 4:
        ff, hh = signal.sosfreqz(signal.butter(4,[1,100],fs=1000,btype='bandpass',output='sos'), worN=4096, fs=1000)
        response = abs(hh)**2
        if filtering.get('notch_60hz'):
            bb, aa = signal.iirnotch(60,30,fs=1000)
            _, hn = signal.freqz(bb,aa,worN=4096,fs=1000); response *= abs(hn)**2
        diagnostic_values.update(filter_frequency_hz=ff, filter_amplitude_gain=response)
        axes[1, 1].plot(ff,20*np.log10(np.maximum(response,1e-12)))
        axes[1, 1].set(xlim=(0,120),ylim=(-80,3),xlabel='Hz',ylabel='Amplitude gain (dB)',title='Declared zero-phase filter response; excludes resampler')
    else:
        axes[1, 1].text(.03,.5,'Filter response not reconstructed:\n'+str(filtering),wrap=True,transform=axes[1,1].transAxes)
        axes[1, 1].set_axis_off()
    np.savez_compressed(out/'reference_filter_diagnostics.npz', **diagnostic_values)
    fig.suptitle(title); save_fig(fig,out,'A_qc')
    fig, axes = plt.subplots(2,6,figsize=(18,6),sharex=True,sharey=True)
    for ri, r in enumerate(REGIONS):
        for si, s in enumerate(STIMULI):
            ax=axes[ri,si]
            for wi,w in enumerate(WINDOWS):
                yy=mean_finite(power[data['good_channel'] & (data['region']==r),si,wi])
                ax.semilogy(f[1:],yy[1:],label=w)
            ax.set(title=f'{r} / {s}',xlim=(1,80),xlabel='Hz')
        axes[ri,0].set_ylabel('PSD (µV²/Hz)')
    axes[0,0].legend(); fig.suptitle('B | Equal-channel power context — '+title); save_fig(fig,out,'B_power')
    for label, values, delta, freq, ylabel in [('C',coh,cdelta,f,'Magnitude-squared coherence'),('D',gc,gdelta,gf,'Spectral Granger (nats)')]:
        fig, axes=plt.subplots(2, max(6,len(pairs)),figsize=(18,7),squeeze=False)
        for si,s in enumerate(STIMULI):
            ax=axes[0,si]
            for wi,w in enumerate(WINDOWS):
                if label=='C':
                    ax.plot(freq,mean_finite(values[:,si,wi]),label=w)
                else:
                    for di,direction in enumerate(['GC→PC','PC→GC']):
                        ax.plot(freq,mean_finite(values[:,si,wi,di]),['-','--'][di],label=f'{w} {direction}')
            ax.set(title=s,xlim=(1,80),xlabel='Hz',ylabel=ylabel if si==0 else '')
            if label=='C':
                ax.plot(freq,mean_finite(shuffled[:,si,1]),':',color='gray',label='post trial-shift control')
                ax.set_ylim(0,1)
            if not np.isfinite(values[:,si]).any(): ax.text(.05,.5,'No valid estimates',transform=ax.transAxes)
        axes[0,0].legend(fontsize=7)
        for i,(retro,ortho) in enumerate(pairs):
            ax=axes[1,i]
            for wi,w in enumerate(WINDOWS):
                if label=='C': ax.plot(freq,delta[i,wi],label=w)
                else:
                    for di,direction in enumerate(['GC→PC','PC→GC']): ax.plot(freq,delta[i,wi,di],['-','--'][di],label=f'{w} {direction}')
            ax.axhline(0,color='gray',lw=.5);ax.set(title=f'{ortho} − {retro}',xlim=(1,80),xlabel='Hz',ylabel='Route difference')
            ax.legend(fontsize=7)
        for ax in axes[1,len(pairs):]: ax.set_axis_off()
        if label=='D' and len(pairs)<axes.shape[1]:
            ax=axes[1,-1];ax.set_axis_on()
            counts=[sum(bool(d['passed']) for d in diagnostics if d['stimulus']==s) for s in STIMULI]
            totals=[sum(d['stimulus']==s for d in diagnostics) for s in STIMULI]
            ax.bar(STIMULI,totals,color='.8',label='attempted');ax.bar(STIMULI,counts,label='passed')
            ax.set(title='VAR diagnostic retention',ylabel='Pair-window fits');ax.legend(fontsize=7)
        fig.suptitle(f'{label} | Descriptive session estimates — '+title);save_fig(fig,out,f'{label}_'+('coherence' if label=='C' else 'granger'))
    try: revision=subprocess.check_output(['git','rev-parse','HEAD'],cwd=Path(__file__).parent,text=True,stderr=subprocess.DEVNULL).strip()
    except (OSError,subprocess.CalledProcessError): revision='unknown'
    digest=hashlib.sha256()
    with source.open('rb') as h:
        for block in iter(lambda:h.read(1024*1024),b''): digest.update(block)
    settings=dict(windows=WINDOWS,max_channels=max_channels,min_trials=min_trials,matched_pairs=pairs,
        orders=list(orders),remove_evoked=remove_evoked,fs_hz=1000,var_fs_hz=250,
        coherence_segment_samples=500,coherence_overlap_samples=250,whiteness_alpha=.01,
        selection='first accepted physical IDs per region; fixed across stimuli',
        inference='descriptive session only; no biological replication or corrected significance',
        trial_shift='one cyclic trial shift within stimulus, not a significance threshold',
        granger='bivariate, not conditional on respiration or other electrodes')
    (out/'provenance.json').write_text(json.dumps(dict(input=str(source.resolve()),input_sha256=digest.hexdigest(),
        revision=revision,preprocessing=provenance,settings=settings,session_id=str(data['session_id'].item()),
        animal_id=str(data['animal_id'].item())),indent=2))
    return out


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('cleaned_npz');p.add_argument('output_directory')
    p.add_argument('--max-channels',type=int,default=5)
    p.add_argument('--min-trials',type=int,default=5)
    p.add_argument('--config',type=Path,help='JSON with matched_pairs and optional orders')
    p.add_argument('--remove-evoked',action='store_true',help='Sensitivity run: subtract within-stimulus trial mean waveform')
    args=p.parse_args()
    cfg=json.loads(args.config.read_text()) if args.config else {}
    analyze(args.cleaned_npz,args.output_directory,args.max_channels,args.min_trials,
            cfg.get('matched_pairs'),tuple(cfg.get('orders',[2,5,10,15,20])),args.remove_evoked)

if __name__=='__main__': main()
