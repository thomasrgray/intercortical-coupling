"""Reproducible synthetic Figure 6 smoke run, not study data."""

from coupling.array_store import save_arrays
import argparse
import json
from pathlib import Path
import numpy as np
from coupling.power import STIMULI
from coupling.lfp_coupling import analyze


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('output',type=Path)
    args=p.parse_args()
    args.output.mkdir(parents=True,exist_ok=False)
    rng=np.random.default_rng(20261013)
    n=8*6
    z=np.zeros((n,5200,2));noise=rng.normal(size=z.shape)
    for k in range(1,len(z[0])):
        z[:,k,0]=.97*z[:,k-1,0]+noise[:,k,0]
        z[:,k,1]=.95*z[:,k-1,1]+.15*z[:,k-1,0]+noise[:,k,1]
    # Generate directly on the canonical 1 kHz input grid.
    x=z[:,700:].transpose(2,0,1)
    t=np.arange(-2000,2500)/1000
    source=args.output/'synthetic_lfp.h5'
    save_arrays(source,lfp_uv=x,time_s=t,stimulus=np.repeat(STIMULI,8),
        electrode_id=[0,32],region=['GC','PC'],event_id=np.arange(n),
        good_channel=np.ones(2,bool),good_trial=np.ones((2,n),bool),
        animal_id='SYNTHETIC',session_id='SYNTHETIC_GC_TO_PC',
        provenance_json=json.dumps(dict(source='synthetic VAR seed 20261013',reference='none; synthetic',
            filtering='none; synthetic 1000 Hz VAR',qc_method='all finite synthetic trials',event_alignment='synthetic TTL')))
    analyze(source,args.output/'fig06_lfp_coupling')
    # Explicit model-order sensitivity: anti-alias filtering adds temporal dependence.
    # Preserve the default run and all failed fits; do not change diagnostic thresholds.
    analyze(source,args.output/'higher_order_sensitivity',orders=tuple(range(10,31)))

if __name__=='__main__':main()
