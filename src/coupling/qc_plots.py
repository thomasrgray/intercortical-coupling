"""QC figures from saved metrics; no additional rejection decisions."""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def plot_qc(out,rows,dins,selected,trials,good,good_trial):
    fig,ax=plt.subplots(figsize=(7,4))
    ax.bar([str(r['din']) for r in dins],[r['trials'] for r in dins],color=['tab:blue' if r['din'] in selected else '0.7' for r in dins])
    ax.set(xlabel='DIN',ylabel='Complete pulses',title='Selected stimulus DINs (blue)')
    fig.tight_layout();fig.savefig(out/'din_counts.png',dpi=160);plt.close(fig)
    fig,axes=plt.subplots(3,1,figsize=(10,8),sharex=True)
    ids=[str(r['electrode_id']) for r in rows]
    for ax,key in zip(axes,['median_uv','sd_uv','bad_trial_fraction']):
        ax.bar(ids,[r[key] or 0 for r in rows],color=['tab:blue' if g else 'tab:red' for g in good]);ax.set_ylabel(key)
    axes[0].set_title('Channel QC: accepted blue / excluded red');axes[-1].set_xlabel('Physical electrode ID')
    axes[-1].tick_params(axis='x',rotation=90);fig.tight_layout();fig.savefig(out/'channel_qc.png',dpi=160);plt.close(fig)
    fig,ax=plt.subplots(figsize=(10,5));ax.imshow(good[:,None]&good_trial,aspect='auto',interpolation='nearest',vmin=0,vmax=1,cmap='gray_r')
    ax.set(xlabel='Original event index',ylabel='Array channel position',title='Accepted channel-trials (black)')
    fig.tight_layout();fig.savefig(out/'trial_retention.png',dpi=160);plt.close(fig)
    fig,axes=plt.subplots(2,1,figsize=(10,6),sharex=True)
    for ax,region in zip(axes,['GC','PC']):
        cs=[i for i,r in enumerate(rows) if r['region']==region and good[i] and good_trial[i].any()]
        if cs:
            c=cs[0];j=np.flatnonzero(good_trial[c])[0]
            ax.plot(np.arange(-2000,5000)/1000,trials[c,j],lw=.5)
            ax.set_title(f'{region}: electrode {rows[c]["electrode_id"]}, event {j} (first accepted)')
        else:ax.set_title(f'{region}: no accepted trace')
        ax.axvline(0,color='red');ax.set_ylabel('µV')
    axes[-1].set_xlabel('Time from TTL (s)');fig.tight_layout();fig.savefig(out/'representative_lfp.png',dpi=160);plt.close(fig)
