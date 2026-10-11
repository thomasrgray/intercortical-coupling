"""Read-only Intan event inventory; never guess between ambiguous count groups."""
from pathlib import Path
from collections import defaultdict
import argparse
import csv
import json
import re
import struct
import numpy as np


def sample_rate(header):
    with Path(header).open('rb') as h:
        b=h.read(12)
    if len(b)!=12:
        raise ValueError('Truncated RHD header')
    magic,major,minor,fs=struct.unpack('<Ihhf',b)
    if magic!=0xc6912702 or not np.isfinite(fs) or fs<1000:
        raise ValueError('Invalid RHD header/sample rate')
    return float(fs)


def read_binary(path,dtype):
    p=Path(path)
    if p.stat().st_size==0 or p.stat().st_size%np.dtype(dtype).itemsize:
        raise ValueError(f'Empty or truncated binary: {p}')
    return np.memmap(p,dtype=dtype,mode='r')


def pulses(path,chunk=1_000_000):
    """Return first-high onset and first-low offset (exclusive), sample coordinates.

    Non-binary files and pulses touching recording boundaries are invalid.
    Chunk boundaries do not add/delete edges. DIN bitfields are not supported.
    """
    a=read_binary(path,'<u2'); starts=[]; stops=[];prev=0
    if a[0]!=0 or a[-1]!=0:
        raise ValueError(f'Incomplete pulse at recording boundary: {path}')
    for k in range(0,len(a),chunk):
        x=np.asarray(a[k:k+chunk])
        if np.any(x>1):
            raise ValueError(f'Expected unpacked binary DIN values 0/1: {path}')
        delta=np.diff(x.astype(np.int8),prepend=prev)
        starts.extend((np.flatnonzero(delta==1)+k).tolist())
        stops.extend((np.flatnonzero(delta==-1)+k).tolist());prev=int(x[-1])
    if len(starts)!=len(stops):
        raise ValueError(f'Unpaired pulse edges: {path}')
    return np.asarray(starts,dtype=np.int64),np.asarray(stops,dtype=np.int64),len(a)


def select_six(counts):
    groups=defaultdict(list)
    for din,n in counts.items():
        if n>0:groups[n].append(din)
    candidates=[sorted(v) for v in groups.values() if len(v)==6]
    if len(candidates)!=1:
        raise ValueError(f'Need exactly one positive count group containing exactly six DINs; found {dict(groups)}. Check delivery logs; no trials were trimmed.')
    return candidates[0]


def inventory(folder,fs):
    records={};rows=[]
    for p in sorted(Path(folder).glob('board-DIN-*.dat')):
        m=re.fullmatch(r'board-DIN-(\d+)\.dat',p.name)
        if not m:continue
        din=int(m[1])
        if din in records:raise ValueError(f'Duplicate DIN {din}')
        on,off,n=pulses(p);records[din]=(on,off,n)
        widths=(off-on)/fs*1000
        rows.append(dict(din=din,filename=p.name,samples=n,trials=len(on),
                         median_pulse_ms=float(np.median(widths)) if len(on) else None))
    if not rows:raise ValueError('No DIN files found')
    if len({x['samples'] for x in rows})!=1:raise ValueError('DIN lengths differ')
    return records,rows


def write_csv(path,rows):
    if not rows:return
    with Path(path).open('w',newline='') as h:
        w=csv.DictWriter(h,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)


def scan(folder,output,header='info.rhd'):
    output=Path(output)
    if output.exists():raise FileExistsError(output)
    fs=sample_rate(Path(folder)/header);records,rows=inventory(folder,fs)
    output.mkdir(parents=True)
    write_csv(output/'din_inventory.csv',rows)
    selected=select_six({d:len(v[0]) for d,v in records.items()})
    (output/'selected_dins.json').write_text(json.dumps(dict(selected_dins=selected,fs_hz=fs),indent=2))
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,ax=plt.subplots(figsize=(7,4))
    ax.bar([str(r['din']) for r in rows],[r['trials'] for r in rows],
           color=['tab:blue' if r['din'] in selected else '0.7' for r in rows])
    ax.set(xlabel='DIN',ylabel='Complete pulses',title='Selected stimulus DINs (blue)')
    fig.tight_layout();fig.savefig(output/'din_counts.png',dpi=160);plt.close(fig)
    return selected,rows


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('session');p.add_argument('output')
    p.add_argument('--header',default='info.rhd');a=p.parse_args()
    selected,rows=scan(a.session,a.output,a.header)
    print(json.dumps(dict(selected_dins=selected,inventory=rows),indent=2))

if __name__=='__main__':main()
