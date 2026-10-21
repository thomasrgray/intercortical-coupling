"""Respiration, GC/PC LFP, and optional sorted-unit firing-rate spectra.

Uses existing coupling data contracts; LFP is never treated as spike activity.
"""
import argparse
import csv
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy import signal
from scipy.ndimage import gaussian_filter1d

from coupling.power import validate as validate_lfp
from coupling.respiration import FS, downsample_counts
from coupling.sniffing import trials_from_csv

ROUTES = {"BA": "Retro", "IA": "Retro", "oBA": "Ortho", "oIA": "Ortho"}
BAND = (1., 12.)
WINDOW = (-.5, 2.5)
REGIONS = ("GC", "PC")
COLORS = {"Respiration": "#777777", "GC firing rate": "#6375a1",
          "PC firing rate": "#b27c70", "GC LFP": "#6375a1", "PC LFP": "#b27c70"}


def validate_trial_table(trials):
    """Canonical chronological trial index and DIN:within-DIN identity."""
    if not trials:
        raise ValueError("No trials")
    ids = [int(r["trial"]) for r in trials]
    if ids != list(range(len(trials))):
        raise ValueError("trial_info_frame abs_trial_num must be 0..N-1 in chronological order")
    seen = {}
    keys = []
    for row in trials:
        din = int(row["din"])
        i = seen.get(din, 0)
        keys.append(f"DIN{din}:trial{i}")
        seen[din] = i + 1
    return keys


def welch_epoch(x):
    """One-sided density PSD; 1-kHz signal, 1-Hz bins, 50% overlap."""
    x = np.asarray(x, dtype=float)
    if x.ndim != 1 or len(x) < FS or not np.isfinite(x).all():
        raise ValueError("Require at least 1 s finite 1-D signal")
    return signal.welch(x, fs=FS, window="hann", nperseg=FS,
                        noverlap=FS//2, detrend="constant", scaling="density")


def epoch(x, time_s):
    t = np.asarray(time_s, dtype=float)
    if len(t) != len(x) or len(t) < 2 or not np.allclose(np.diff(t), 1/FS, atol=1e-8):
        raise ValueError("Signal time grid must have 1-kHz sampling")
    if t[0] > WINDOW[0]+1e-8 or t[-1] < WINDOW[1]-1/FS-1e-8:
        raise ValueError("Signal must cover -0.5 to +2.5 s relative to DIN onset")
    a = int(round((WINDOW[0]-t[0])*FS))
    b = a + int(round((WINDOW[1]-WINDOW[0])*FS))
    return np.asarray(x[a:b], dtype=float)


def summarize(psds, kind, label, route, source):
    if not psds:
        return None
    values = np.asarray(psds, dtype=float)
    if values.ndim != 2 or not np.isfinite(values).all():
        raise ValueError("Nonfinite power spectra")
    avg = values.mean(axis=0)
    sem = values.std(axis=0,ddof=1)/np.sqrt(len(values)) if len(values)>1 else np.zeros_like(avg)
    return dict(label=label, kind=kind, route=route, source=source,
                n_replicates=len(values), mean=avg, sem=sem)


def validate_source(provenance, adc):
    source = provenance.get("source", "")
    if not source or Path(source).name != Path(adc).parent.name:
        raise ValueError("Cleaned LFP source session does not match ADC parent session")
    if provenance.get("filtering",{}).get("fs_hz") != FS:
        raise ValueError("Cleaned LFP must be sampled at 1000 Hz")
    if not provenance.get("event_alignment","").startswith("TTL first-high"):
        raise ValueError("Need original DIN TTL-aligned cleaned LFP, not sniff-aligned epochs")


def respiration_spectra(adc, trials, raw_fs):
    raw = np.memmap(adc, mode="r", dtype="<u2")
    if int(raw_fs)!=raw_fs or raw_fs<FS:
        raise ValueError("raw_fs must be integer and at least 1000")
    if max(r["end_sample"] for r in trials) > len(raw):
        raise ValueError("Trial sample beyond ADC: recording mismatch or truncated file")
    centered = downsample_counts(raw, raw_fs)
    filtered = signal.sosfiltfilt(
        signal.butter(4,[.5,20],btype="bandpass",fs=FS,output="sos"), centered)
    before,after = -round(WINDOW[0]*FS), round(WINDOW[1]*FS)
    out = {"Retro":{}, "Ortho":{}}
    for index,trial in enumerate(trials):
        route=ROUTES.get(trial["stimulus"])
        if route is None:continue
        center=round(trial["start_sample"]*FS/raw_fs)
        if center-before<0 or center+after>len(filtered):continue
        x=filtered[center-before:center+after]
        sd=float(np.std(x))
        if not np.isfinite(sd) or sd<=0:continue
        freq, power=welch_epoch((x-x.mean())/sd)
        out[route][index]=power
    if not any(out.values()):
        raise ValueError("No valid in-bounds respiration trials")
    return freq,out


def lfp_spectra(lfp_npz, trials, adc):
    result={}
    with np.load(lfp_npz, allow_pickle=False) as d:
        provenance=validate_lfp(d)
        validate_source(provenance,adc)
        if "event_id" not in d.files or not np.array_equal(d["event_id"],np.arange(len(trials))):
            raise ValueError("Cleaned LFP event_id must map 1:1 to chronological trial frame")
        if not np.array_equal(d["stimulus"],np.asarray([t["stimulus"] for t in trials])):
            raise ValueError("LFP stimulus order differs from trial_info_frame")
        if d["lfp_uv"].shape[1] != len(trials):
            raise ValueError("LFP trial count differs from trial_info_frame")
        data,t=d["lfp_uv"],d["time_s"]
        for route in ("Retro","Ortho"):
            indices=[i for i,r in enumerate(trials) if ROUTES.get(r["stimulus"])==route]
            for region in REGIONS:
                channel_psds=[]
                for channel in np.flatnonzero(d["good_channel"] & (d["region"]==region)):
                    per_trial=[]
                    for i in indices:
                        if not d["good_trial"][channel,i]:continue
                        x=epoch(data[channel,i],t)
                        if not np.isfinite(x).all():
                            raise ValueError("Nonfinite accepted LFP trial")
                        freq,power=welch_epoch(x)
                        per_trial.append(power)
                    if per_trial:
                        channel_psds.append(np.mean(per_trial,axis=0))
                label=f"{region} LFP"
                result[(route,label)]=summarize(channel_psds,"lfp",label,route,str(lfp_npz))
    return freq,result


def spike_spectra(spikes, config, trials, keys, adc):
    from coupling.spike_coupling import load
    cfg=json.loads(Path(config).read_text())
    spec={**cfg.get("defaults",{}),"input":str(Path(spikes).resolve())}
    data=load(Path(adc).parent.resolve(),spec)
    if not np.allclose(np.diff(data["time_s"]),1/FS,atol=1e-8):
        raise ValueError("Require 1-kHz binned sorted spike arrays")
    ids=list(map(str,data["trial_id"]))
    if len(set(ids))!=len(ids):
        raise ValueError("Duplicate spike trial IDs")
    lookup={key:i for i,key in enumerate(ids)}
    if any(key not in keys for key in lookup):
        raise ValueError("Sorted spikes contain IDs absent from the trial frame")
    matches={i:lookup[key] for i,key in enumerate(keys) if key in lookup}
    if not matches:
        raise ValueError("No spike trials match DIN and within-DIN index")
    for i,j in matches.items():
        if data["stimulus"][j]!=trials[i]["stimulus"]:
            raise ValueError("Spike stimulus mismatch for "+keys[i])
    results={}
    for route in ("Retro","Ortho"):
        for region in REGIONS:
            unit_psds=[]
            for unit in np.flatnonzero(data["good_unit"] & (data["region"]==region)):
                per_trial=[]
                for i,j in matches.items():
                    if ROUTES.get(trials[i]["stimulus"])!=route or not data["good_trial"][j]:
                        continue
                    if data["stratum"][j]!="laser:off":
                        continue
                    counts=epoch(data["spikes"][j,unit],data["time_s"])
                    if not np.isfinite(counts).all() or np.any(counts<0):
                        raise ValueError("Invalid selected spike counts")
                    rate=gaussian_filter1d(counts*FS,sigma=25,mode="reflect")
                    freq,power=welch_epoch(rate)
                    per_trial.append(power)
                if per_trial:
                    unit_psds.append(np.mean(per_trial,axis=0))
            label=f"{region} firing rate"
            results[(route,label)]=summarize(unit_psds,"firing_rate",label,route,str(spikes))
    return freq,results


def normalized_psd(freq,power):
    mask=(freq>=BAND[0]) & (freq<=BAND[1])
    integral=np.trapz(power[mask],freq[mask])
    if not np.isfinite(integral) or integral<=0:
        return np.full_like(power,np.nan)
    return power/integral


def plot_spectra(freq,results,output):
    has_spikes=any(r["kind"]=="firing_rate" for r in results)
    fig,axes=plt.subplots(1,2,figsize=(10.5,4.2),constrained_layout=True)
    for ax,kind in zip(axes,("non_lfp","lfp")):
        for row in results:
            if kind=="lfp" and row["kind"]!="lfp":continue
            if kind=="non_lfp" and row["kind"]=="lfp":continue
            y=normalized_psd(freq,row["mean"])
            if not np.isfinite(y).any():continue
            ax.plot(freq,y,c=COLORS.get(row["label"]),lw=1.5,
                    ls="-" if row["route"]=="Ortho" else "--",
                    label=f'{row["route"]} {row["label"]}')
        ax.set(xlim=(0,12),xlabel="Frequency (Hz)",ylabel="Unit-area PSD (1–12 Hz)")
        ax.legend(frameon=False,fontsize=7)
    axes[0].set_title("B   Respiration + firing-rate power" if has_spikes else
                      "B   Respiration only: spikes unavailable",loc="left",weight="bold")
    axes[1].set_title("GC / PC LFP power (not firing rate)",loc="left",weight="bold")
    fig.suptitle("DIN-matched, single-session spectral shapes (descriptive)",fontsize=10)
    fig.savefig(output,dpi=300)
    plt.close(fig)


def run(adc,trial_info,lfp,output,spikes=None,spike_config=None,raw_fs=30000):
    adc,trial_info,lfp,output=map(Path,(adc,trial_info,lfp,output))
    if not lfp.exists():raise FileNotFoundError("Missing cleaned LFP archive: "+str(lfp))
    if output.exists():raise FileExistsError(output)
    trials=trials_from_csv(trial_info)
    keys=validate_trial_table(trials)
    freq,resp=respiration_spectra(adc,trials,raw_fs)
    results=[]
    for route,rows in resp.items():
        summary=summarize(list(rows.values()),"respiration","Respiration",route,str(adc))
        if summary is not None:results.append(summary)
    f_lfp,lfp_rows=lfp_spectra(lfp,trials,adc)
    if not np.array_equal(freq,f_lfp):raise ValueError("LFP frequency grid mismatch")
    results.extend(x for x in lfp_rows.values() if x is not None)
    if spikes:
        if not spike_config:raise ValueError("Spike config required for region mapping and QC")
        f_spikes,spike_rows=spike_spectra(spikes,spike_config,trials,keys,adc)
        if not np.array_equal(freq,f_spikes):raise ValueError("Spike frequency grid mismatch")
        results.extend(x for x in spike_rows.values() if x is not None)
        if not any(x["kind"]=="firing_rate" for x in results):
            raise ValueError("No eligible sorted units after QC and laser-off selection")
    if not results:raise ValueError("No eligible spectra")
    output.mkdir(parents=True)
    with (output/"spectra.csv").open("w",newline="") as h:
        fields=["route","signal","signal_type","frequency_hz","mean_psd","sem_psd",
                "normalized_shape","n_replicates","source"]
        w=csv.DictWriter(h,fieldnames=fields);w.writeheader()
        for row in results:
            norm=normalized_psd(freq,row["mean"])
            for i,hz in enumerate(freq):
                if hz>12:break
                w.writerow(dict(route=row["route"],signal=row["label"],signal_type=row["kind"],
                                frequency_hz=float(hz),mean_psd=float(row["mean"][i]),
                                sem_psd=float(row["sem"][i]),normalized_shape=float(norm[i]),
                                n_replicates=row["n_replicates"],source=row["source"]))
    for extension in ("pdf","png"):
        plot_spectra(freq,results,output/f"respiration_neural_spectra.{extension}")
    (output/"provenance.json").write_text(json.dumps(dict(adc=str(adc.resolve()),
        trial_info=str(trial_info.resolve()),lfp=str(lfp.resolve()),
        spikes=str(Path(spikes).resolve()) if spikes else None,
        spike_config=str(Path(spike_config).resolve()) if spike_config else None,
        raw_fs_hz=raw_fs,sample_rate_hz=FS,epoch_s=WINDOW,
        welch_nperseg=FS,welch_noverlap=FS//2,firing_rate_sigma_ms=25,
        lfp_input="QC-cleaned LFP; not spikes",
        spike_selection="eligible sorted units, laser:off, explicit DIN/trial matches",
        respiration="0.5-20 Hz filtered ADC, per-trial z-scored",
        normalization="Unit-area 1-12 Hz for plots only; unnormalized PSD in spectra.csv",
        limitations="Single-session descriptive, channels/units not biological replicates",
        n_trials=len(trials),n_spectra=len(results)),indent=2)+"\n")
    return results


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--adc",required=True)
    p.add_argument("--trial-info",required=True)
    p.add_argument("--lfp",required=True,help="cleaned_lfp.npz from coupling-prepare")
    p.add_argument("--spikes",help="Optional sorted spike NPZ or blech_clust HDF5")
    p.add_argument("--spike-config",help="Existing spike config JSON")
    p.add_argument("--raw-fs",type=int,default=30000)
    p.add_argument("--output",required=True)
    args=p.parse_args()
    result=run(args.adc,args.trial_info,args.lfp,args.output,args.spikes,
               args.spike_config,args.raw_fs)
    print(f"Saved {len(result)} route/signal spectra in {args.output}")


if __name__=="__main__":
    main()
