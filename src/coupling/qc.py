"""Explicit channel/trial QC metrics. No spike sorting dependency."""
import numpy as np


def flat_fraction(x,min_run):
    if len(x)==0:return 1.
    edges=np.r_[0,np.flatnonzero(np.diff(x)!=0)+1,len(x)]
    lengths=np.diff(edges)
    return float(lengths[lengths>=min_run].sum()/len(x))


def metrics(x,fs=1000,flat_ms=10):
    x=np.asarray(x)
    finite=bool(np.isfinite(x).all())
    if not finite or not x.size:
        return dict(finite=finite,nonconstant=False,median_uv=None,sd_uv=None,
                    flat_fraction=1.,outlier_fraction=1.,bad=True)
    median=float(np.median(x));mad=float(np.median(np.abs(x-median)))
    fraction=float(np.mean(np.abs(x-median)>10*1.4826*mad)) if mad>0 else float(np.mean(x!=median))
    flat=flat_fraction(x,max(2,int(np.ceil(fs*flat_ms/1000))))
    nonconstant=bool(np.ptp(x)>0)
    return dict(finite=True,nonconstant=nonconstant,median_uv=median,sd_uv=float(np.std(x)),
                flat_fraction=flat,outlier_fraction=fraction,
                bad=bool(not nonconstant or mad==0 or flat>.01 or fraction>.01))


def peer_flags(rows,regions,min_peers=3):
    """Leave-one-out comparison with usable same-region peers; zero MAD is explicit.

    Review/insufficient peers is not silently accepted. Rows already failing
    waveform checks are excluded from the peer pool. One fixed comparison pass.
    """
    reasons=[[] for _ in rows]
    for i,row in enumerate(rows):
        if row['bad']:
            reasons[i].append('waveform_qc');continue
        peers=[r for j,r in enumerate(rows) if j!=i and regions[j]==regions[i] and not r['bad']]
        if len(peers)<min_peers:
            reasons[i].append('insufficient_same_region_peers');continue
        for key in ('median_uv','sd_uv','bad_trial_fraction'):
            vals=np.array([r[key] for r in peers]);center=np.median(vals)
            scale=1.4826*np.median(np.abs(vals-center));delta=abs(row[key]-center)
            if (scale==0 and delta>0) or (scale>0 and delta>3*scale):
                reasons[i].append('peer_'+key)
    return reasons
