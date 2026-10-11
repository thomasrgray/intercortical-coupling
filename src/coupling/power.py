"""Baseline-normalized power from explicitly mapped, previously cleaned LFP trials.

This stage does not perform CAR or certify that an input is clean. The required
provenance and rejection masks make that upstream responsibility explicit.
"""
import argparse
import json
from pathlib import Path
import hashlib
import numpy as np
from scipy.signal import stft
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

STIMULI = ("BA", "IA", "S", "oBA", "oIA", "BK")
REGIONS = ("GC", "PC")
SETTINGS = dict(fs_hz=1000, window="hann", nperseg=500, noverlap=499,
                max_hz=70, baseline_s=[-1.5, -0.5], plot_s=[-1, 1.75],
                response_s=[0.25, 1.75], ddof=0, version=1,
                weighting="equal trials within channel, then equal channels",
                boundary=None, padded=False, detrend=False)

def normalize_power(x, time_s):
    """One accepted channel/trial, with baseline-only moments for each frequency."""
    x = np.asarray(x, dtype=float)
    time_s = np.asarray(time_s, dtype=float)
    if x.ndim != 1 or x.shape != time_s.shape or not np.isfinite(x).all():
        raise ValueError("STFT requires a finite one-dimensional accepted trial")
    if len(x) < 500 or not np.allclose(np.diff(time_s), .001, atol=1e-9, rtol=0):
        raise ValueError("Require uniformly sampled 1000 Hz input")
    f, t, z = stft(x, fs=1000, window="hann", nperseg=500, noverlap=499,
                   boundary=None, padded=False, detrend=False)
    t += time_s[0]
    keep = f <= 70
    f, p = f[keep], np.abs(z[keep])**2
    b = (t >= -1.5-1e-9) & (t < -.5-1e-9)
    if time_s[0] > -1.75+1e-9 or time_s[-1] < 1.999-1e-9 or b.sum() < 2:
        raise ValueError("Need input covering at least -1.75 to +1.999 seconds")
    mu, sd = p[:, b].mean(-1), p[:, b].std(-1)
    valid = np.isfinite(sd) & (sd > np.finfo(float).eps * np.maximum(1., np.abs(mu)))
    norm = np.full_like(p, np.nan)
    np.divide(p-mu[:, None], sd[:, None], out=norm, where=valid[:, None])
    return f, t, norm, valid

def validate(d):
    required = {"lfp_uv", "time_s", "stimulus", "electrode_id", "region",
                "good_channel", "good_trial", "session_id", "animal_id", "provenance_json"}
    if not required.issubset(d.files):
        raise ValueError(f"Missing fields: {sorted(required-set(d.files))}")
    x = d["lfp_uv"]
    if x.ndim != 3:
        raise ValueError("lfp_uv axes must be channel, event, sample")
    c, r, n = x.shape
    for key, shape in [("time_s", (n,)), ("stimulus", (r,)), ("electrode_id", (c,)),
                       ("region", (c,)), ("good_channel", (c,)), ("good_trial", (c,r))]:
        if d[key].shape != shape:
            raise ValueError(f"Wrong shape for {key}; expected {shape}")
    if len(np.unique(d["electrode_id"])) != c:
        raise ValueError("Duplicate physical electrode IDs")
    if set(d["stimulus"].tolist()) != set(STIMULI) or set(d["region"].tolist()) != set(REGIONS):
        raise ValueError("Require six named stimuli and explicit GC/PC region labels")
    if d["good_channel"].dtype != bool or d["good_trial"].dtype != bool:
        raise ValueError("QC masks must be boolean")
    p = json.loads(str(d["provenance_json"].item()))
    for key in ("source", "reference", "filtering", "qc_method", "event_alignment"):
        if not p.get(key):
            raise ValueError(f"Missing preprocessing provenance: {key}")
    if not np.isfinite(d["time_s"]).all():
        raise ValueError("Nonfinite time vector")
    return p

def plot_panel(out, name, f, t, z, clim):
    fig, ax = plt.subplots(figsize=(7,4))
    mesh = ax.pcolormesh(t*1000, f, z, shading="auto", cmap="turbo", vmin=-clim, vmax=clim)
    ax.axvline(0, color="white", lw=1)
    ax.set(xlim=(-1000,1750), ylim=(0,70), xlabel="Time from TTL onset (ms)", ylabel="Frequency (Hz)", title=name)
    fig.colorbar(mesh, ax=ax, label="Power z-score vs prestimulus baseline")
    fig.tight_layout(); fig.savefig(out/f"{name}_spectrogram.png", dpi=160); plt.close(fig)
    pre = (t >= -1-1e-9) & (t < -.5-1e-9)
    post = (t >= .25-1e-9) & (t <= 1.75+1e-9)
    # Pre ends at -0.5; post starts at +0.25 to keep 500 ms windows away from onset.
    pre_y, post_y = z[:,pre].mean(-1), z[:,post].mean(-1)
    fig, ax = plt.subplots(figsize=(6,4))
    ax.plot(f, pre_y, label="Pre (-1 to -0.5 s)")
    ax.plot(f, post_y, label="Post (+0.25 to +1.75 s)")
    ax.set(xlabel="Frequency (Hz)", ylabel="Baseline-normalized power (z)", title=name)
    ax.legend(); fig.tight_layout(); fig.savefig(out/f"{name}_pre_post.png", dpi=160); plt.close(fig)
    np.savetxt(out/f"{name}_pre_post.csv", np.c_[f,pre_y,post_y], delimiter=",",
               header="frequency_hz,pre_z,post_z", comments="")

def analyze(source, output, clim=3.):
    source, output = Path(source), Path(output)
    if not np.isfinite(clim) or clim <= 0:
        raise ValueError("Color limit must be positive and finite")
    # Never reuse an output folder: prevents stale plots or silent replacement.
    if output.exists():
        raise FileExistsError(output)
    with np.load(source, allow_pickle=False) as d:
        provenance = validate(d)
        rows, panels, counts = [], [], []
        grid = None
        for region in REGIONS:
            for stim in STIMULI:
                channel_means = []
                for c in np.flatnonzero((d["region"] == region) & d["good_channel"]):
                    total, count = None, None
                    trials = np.flatnonzero((d["stimulus"] == stim) & d["good_trial"][c])
                    for r in trials:
                        f,t,z,valid = normalize_power(d["lfp_uv"][c,r], d["time_s"])
                        crop = (t >= -1-1e-9) & (t <= 1.75+1e-9)
                        z = z[:,crop]; grid = (f,t[crop])
                        if total is None:
                            total, count = np.zeros_like(z), np.zeros(len(f), dtype=int)
                        total += np.nan_to_num(z); count += valid
                    if total is None:
                        continue
                    mean = np.full_like(total, np.nan)
                    np.divide(total, count[:,None], out=mean, where=count[:,None]>0)
                    channel_means.append(mean)
                    rows.append(dict(region=region, stimulus=stim, electrode_id=str(d["electrode_id"][c]),
                                     accepted_trials=len(trials), valid_trials_per_frequency=count.tolist()))
                if not channel_means:
                    raise ValueError(f"No accepted channels/trials for {region}/{stim}")
                means = np.stack(channel_means)
                denom = np.isfinite(means).sum(0)
                panel = np.full(means.shape[1:], np.nan)
                np.divide(np.nansum(means,axis=0),denom,out=panel,where=denom>0)
                panels.append(panel); counts.append(denom)
        f,t = grid
        output.mkdir(parents=True)
        arrays = np.stack(panels).reshape(2,6,len(f),len(t))
        np.savez_compressed(output/"session_power.npz", z_power=arrays, frequency_hz=f,time_s=t,
                            region=np.array(REGIONS),stimulus=np.array(STIMULI),
                            n_channels=np.stack(counts).reshape(arrays.shape),
                            session_id=d["session_id"],animal_id=d["animal_id"],
                            settings_json=json.dumps(SETTINGS,sort_keys=True))
        for i, region in enumerate(REGIONS):
            for j, stim in enumerate(STIMULI):
                plot_panel(output,f"{region}_{stim}",f,t,arrays[i,j],clim)
        digest = hashlib.sha256()
        with source.open("rb") as h:
            for block in iter(lambda:h.read(1024*1024),b""):
                digest.update(block)
        (output/"provenance.json").write_text(json.dumps(dict(input_sha256=digest.hexdigest(),
            preprocessing=provenance,settings=SETTINGS,color_limit=clim,channel_counts=rows),indent=2))

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("cleaned_npz");p.add_argument("output_directory");p.add_argument("--clim",type=float,default=3.)
    a=p.parse_args();analyze(a.cleaned_npz,a.output_directory,a.clim)

if __name__ == "__main__":
    main()
