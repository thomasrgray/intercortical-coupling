"""Dewan-inspired sniff metrics from synchronized Intan ADC + trial_info_frame.

Retains BOTH peak-based Dewan timestamps and crossing-based Maier candidates.
All pressure landmarks are *candidates* until the experiment's polarity is validated.
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

from coupling.respiration import downsample_counts, detect_cycles, FS


def trials_from_csv(path):
    with Path(path).open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    required = {"start_taste", "end_taste", "taste", "dig_in_num_taste"}
    if not rows or not required.issubset(rows[0]):
        raise ValueError("Expected trial_info_frame columns: " + ", ".join(sorted(required)))
    trials = []
    for row in rows:
        try:
            start, end = int(float(row["start_taste"])), int(float(row["end_taste"]))
            din = int(float(row["dig_in_num_taste"]))
        except (TypeError, ValueError) as exc:
            raise ValueError("Invalid trial frame timing or DIN") from exc
        if end <= start:
            raise ValueError("Stimulus pulse must have positive duration")
        trials.append(dict(start_sample=start, end_sample=end,
                           stimulus=row["taste"], din=din,
                           trial=row.get("abs_trial_num", str(len(trials)))))
    return sorted(trials, key=lambda x: x["start_sample"])


def detect_peaks(x, fs=FS, polarity="positive"):
    """Dewan-style positive/negative peak detection; timestamps are in samples."""
    if polarity not in ("positive", "negative"):
        raise ValueError("polarity must be positive or negative")
    # Onset convention does not affect peak *selection*: positive vs negative is explicit.
    sign = 1 if polarity == "positive" else -1
    inh, _ = signal.find_peaks(sign*x, distance=round(.030*fs),
                               width=round(.010*fs), prominence=.1)
    exh, _ = signal.find_peaks(-sign*x, distance=round(.030*fs),
                               width=round(.010*fs), prominence=.1)
    # Suppress multiple candidate inhalation peaks lacking a separating exhalation.
    kept = []
    for idx in inh:
        if kept and not np.any((exh > kept[-1]) & (exh < idx)):
            if sign*x[idx] > sign*x[kept[-1]]:
                kept[-1] = idx
        else:
            kept.append(int(idx))
    return np.asarray(kept, dtype=int), np.asarray(exh, dtype=int)


def process(adc_path, trials_path, out_dir, raw_fs=30000, polarity="unknown",
            before_s=2., after_s=2., bin_ms=100, step_ms=50):
    adc_path, trials_path, out = Path(adc_path), Path(trials_path), Path(out_dir)
    if out.exists():
        raise FileExistsError("Refusing to overwrite existing analysis: " + str(out))
    if raw_fs <= 0 or before_s <= 0 or after_s <= 0 or bin_ms <= 0 or step_ms <= 0:
        raise ValueError("All sampling/window arguments must be positive")
    if polarity not in ("unknown", "positive", "negative"):
        raise ValueError("polarity must be unknown, positive or negative")
    if adc_path.stat().st_size % 2:
        raise ValueError("ADC uint16 binary has an odd byte count")
    raw = np.memmap(adc_path, dtype="<u2", mode="r")
    trials = trials_from_csv(trials_path)
    if max(t["end_sample"] for t in trials) > len(raw):
        raise ValueError("Trial start/end exceed ADC length: verify recording and sampling clock")
    # Intan RHD USB ADC storage, centered counts; voltage conversion depends on board.
    centered = downsample_counts(raw, raw_fs)
    sos = signal.butter(4, [.5, 20.], btype="bandpass", fs=FS, output="sos")
    filtered = signal.sosfiltfilt(sos, centered)
    # Dewan performs a further high-pass at 1 Hz and z-scores each trial.
    hp = signal.butter(2, 1., btype="highpass", fs=FS, output="sos")
    dewan_filtered = signal.sosfiltfilt(hp, filtered)
    cycles = detect_cycles(filtered)
    good_crossings = np.array([c["negative_start_sample"] for c in cycles if c["screen_pass"]], dtype=int)
    positive_crossings = np.array([c["positive_start_sample"] for c in cycles if c["screen_pass"]], dtype=int)
    out.mkdir(parents=True)
    summary, bins = [], []
    plot_data = {}
    n_before, n_after = round(before_s*FS), round(after_s*FS)
    for trial in trials:
        start = round(trial["start_sample"]*FS/raw_fs)
        end = round(trial["end_sample"]*FS/raw_fs)
        row = dict(**trial, event_sample_1khz=start,
                   pulse_ms=(trial["end_sample"]-trial["start_sample"])*1000/raw_fs,
                   polarity=polarity, status="ready")
        lo, hi = start-n_before, start+n_after
        if lo < 0 or hi > len(filtered):
            row["status"] = "out_of_bounds"
            summary.append(row)
            continue
        x = dewan_filtered[lo:hi]
        z = (x-np.mean(x))/np.std(x) if np.std(x) > 0 else np.zeros_like(x)
        # Unknown polarity: emit peak candidates but NEVER designate inhalations.
        peak_a, ex_a = detect_peaks(z, polarity="positive")
        peak_b, ex_b = detect_peaks(z, polarity="negative")
        local = (lo, hi)
        for tag, peaks, exh in [("positive_peak", peak_a, ex_a), ("negative_peak", peak_b, ex_b)]:
            rel_ms = (peaks+lo-start)*1000/FS
            pre = rel_ms[(rel_ms >= -before_s*1000) & (rel_ms < 0)]
            post = rel_ms[(rel_ms >= 0) & (rel_ms < after_s*1000)]
            row[tag+"_pre_count"] = int(len(pre))
            row[tag+"_post_count"] = int(len(post))
            for j in range(3):
                row[tag+f"_latency_{j+1}_ms"] = float(post[j]) if j<len(post) else ""
            # Dewan-like peak width, measured at half prominence, NOT full physiological inspiration duration.
            widths = signal.peak_widths(z if tag == "positive_peak" else -z, peaks, rel_height=.5)[0] if len(peaks) else []
            row[tag+"_mean_half_prominence_width_ms"] = float(np.mean(widths)) if len(widths) else ""
            for center_ms in np.arange(-before_s*1000+bin_ms/2, after_s*1000-bin_ms/2+1, step_ms):
                count = int(((rel_ms >= center_ms-bin_ms/2) & (rel_ms < center_ms+bin_ms/2)).sum())
                bins.append(dict(trial=trial["trial"], stimulus=trial["stimulus"],
                                 din=trial["din"], peak_polarity=tag,
                                 bin_center_ms=float(center_ms), count=count,
                                 rate_hz=count*1000/bin_ms))
        for tag, arr in [("negative_crossing",good_crossings),
                         ("positive_crossing",positive_crossings)]:
            rel = (arr[(arr>=lo)&(arr<hi)]-start)*1000/FS
            post = rel[rel>=0]
            row[tag+"_first_after_ms"] = float(post[0]) if len(post) else ""
            row[tag+"_pre_count"] = int(((rel<0)&(rel>=-before_s*1000)).sum())
            row[tag+"_post_count"] = int(((rel>=0)&(rel<after_s*1000)).sum())
        if polarity == "positive":
            row["verified_inhalation_peak_after_ms"] = row["positive_peak_latency_1_ms"]
            row["verified_inhalation_onset_after_ms"] = row["positive_crossing_first_after_ms"]
        elif polarity == "negative":
            row["verified_inhalation_peak_after_ms"] = row["negative_peak_latency_1_ms"]
            row["verified_inhalation_onset_after_ms"] = row["negative_crossing_first_after_ms"]
        else:
            row["verified_inhalation_peak_after_ms"] = ""
            row["verified_inhalation_onset_after_ms"] = ""
        summary.append(row)
        if trial["stimulus"] not in plot_data:
            plot_data[trial["stimulus"]] = []
        plot_data[trial["stimulus"]].append(z)
    def save_csv(path, rows):
        if not rows:
            return
        columns = list(dict.fromkeys(k for row in rows for k in row))
        with path.open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=columns)
            writer.writeheader()
            writer.writerows(rows)
    save_csv(out/"trials.csv", summary)
    save_csv(out/"binned_sniff_counts.csv", bins)
    save_csv(out/"cycles.csv", cycles)
    t = np.arange(-n_before,n_after)/FS
    fig, axes = plt.subplots(len(plot_data),1,figsize=(11,max(3,2.2*len(plot_data))),sharex=True,squeeze=False)
    for ax,(stim,traces) in zip(axes[:,0],sorted(plot_data.items())):
        stack = np.asarray(traces)
        mean = stack.mean(axis=0)
        sem = stack.std(axis=0,ddof=1)/np.sqrt(len(stack)) if len(stack)>1 else np.zeros_like(mean)
        ax.plot(t,mean,color="black",lw=1)
        ax.fill_between(t,mean-sem,mean+sem,alpha=.2,color="gray")
        ax.axvline(0,color="red",ls="--",lw=.8)
        ax.set(ylabel="Z",title=f"{stim} (n={len(stack)})")
    if plot_data:
        axes[-1,0].set_xlabel("Seconds from DIN rising edge")
        fig.suptitle("Dewan-style pressure traces | " + ("POLARITY NOT VERIFIED" if polarity=="unknown" else "polarity: "+polarity))
        fig.tight_layout()
        fig.savefig(out/"respiration_by_stimulus.png",dpi=180)
        fig.savefig(out/"respiration_by_stimulus.pdf")
    plt.close(fig)
    routes={"oBA":"ortho","oIA":"ortho","BK":"ortho_blank","BA":"retro","IA":"retro","S":"taste"}
    groups={}
    for row in summary:
        if row["status"]!="ready":continue
        group=(row["stimulus"],routes.get(row["stimulus"],"unmapped"))
        groups.setdefault(group,[]).append(row)
    grouped=[]
    for (stim,route),rows in sorted(groups.items()):
        item=dict(stimulus=stim,route=route,n_trials=len(rows),
                  polarity=polarity,mean_pulse_ms=float(np.mean([r["pulse_ms"] for r in rows])))
        for key in ["positive_peak_post_count","negative_peak_post_count",
                    "negative_crossing_first_after_ms","positive_crossing_first_after_ms"]:
            vals=[float(r[key]) for r in rows if r.get(key,"")!=""]
            item["mean_"+key]=float(np.mean(vals)) if vals else ""
        grouped.append(item)
    save_csv(out/"session_summary.csv",grouped)
    (out/"provenance.json").write_text(json.dumps(dict(adc=str(adc_path.resolve()),
        trial_info=str(trials_path.resolve()),raw_fs_hz=raw_fs,output_fs_hz=FS,
        filter_hz=[.5,20],dewan_highpass_hz=1,polarity=polarity,
        n_trials=len(summary),n_ready=sum(r["status"]=="ready" for r in summary),
        caveat="Crossings are candidates; polarity must be independently verified. Peak widths are not inspiration duration."),indent=2))
    return summary


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--adc",required=True)
    p.add_argument("--trial-info",required=True)
    p.add_argument("--output",required=True)
    p.add_argument("--raw-fs",type=int,default=30000)
    p.add_argument("--polarity",choices=["unknown","positive","negative"],default="unknown")
    p.add_argument("--before",type=float,default=2.)
    p.add_argument("--after",type=float,default=2.)
    args=p.parse_args()
    rows=process(args.adc,args.trial_info,args.output,args.raw_fs,args.polarity,args.before,args.after)
    print(f"Processed {len(rows)} trials; outputs: {args.output}")


if __name__=="__main__":
    main()
