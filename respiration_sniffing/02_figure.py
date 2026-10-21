"""Reference-inspired A-D respiration figure (NOT neural coherence).
Run: python respiration_sniffing/02_figure.py --adc PATH --trial-info PATH --output results/figure.pdf
"""
import argparse
from pathlib import Path
import csv
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import signal
from coupling.respiration import downsample_counts, FS
from coupling.sniffing import trials_from_csv, detect_peaks

ROUTES = {"oBA": "Orthonasal", "oIA": "Orthonasal", "BA": "Retronasal", "IA": "Retronasal"}
PALETTE = {"Orthonasal": "#817fae", "Retronasal": "#656565"}


def scale_bar(ax, seconds=.2, label="200 ms"):
    """Data-coordinate scale bar, not dependent on figure DPI or panel size."""
    left, right = ax.get_xlim()
    if right-left <= seconds*1.3:
        raise ValueError("Waveform panel too short for 200-ms scale bar")
    bottom, top = ax.get_ylim()
    x1 = right - .075*(right-left)
    y = bottom + .11*(top-bottom)
    ax.plot([x1-seconds, x1], [y, y], color="black", lw=2.2, solid_capstyle="butt", zorder=10)
    ax.text(x1-seconds/2, y+.035*(top-bottom), label, ha="center", va="bottom", fontsize=8)
    return (x1-seconds, x1)


def get_snippets(adc, trials, raw_fs, pre=.5, post=2.5):
    raw = np.memmap(adc, mode="r", dtype="<u2")
    filtered = signal.sosfiltfilt(
        signal.butter(4, [.5, 20], btype="bandpass", fs=FS, output="sos"),
        downsample_counts(raw, raw_fs))
    n0, n1 = round(pre*FS), round(post*FS)
    snippets = {}
    for trial in trials:
        stim = trial["stimulus"]
        if stim not in ROUTES:
            continue
        center = round(trial["start_sample"]*FS/raw_fs)
        if center < n0 or center+n1 > len(filtered):
            continue
        x = filtered[center-n0:center+n1].copy()
        sd = x.std()
        if sd == 0:
            continue
        z = (x-x.mean())/sd
        snippets.setdefault(stim, []).append(z)
    if not snippets:
        raise ValueError("No in-bounds valid ortho/retro trials")
    return snippets, np.arange(-n0, n1)/FS


def route_arrays(snippets):
    groups = {}
    for odor, trials in snippets.items():
        groups.setdefault(ROUTES[odor], []).extend(trials)
    return {route: np.asarray(rows) for route, rows in groups.items()}


def figure(snippets, times, output, peak_polarity="unknown", spectra_csv=None):
    """Figure anatomy mirrors reference; C is pressure POWER, not neural coherence."""
    grouped = route_arrays(snippets)
    if not grouped:
        raise ValueError("No route groups")
    fig, axs = plt.subplots(2, 2, figsize=(9.4, 6.9), constrained_layout=True)
    ax = axs[0, 0]
    # A: one illustrative trial per route, vertically offset, with candidate peaks.
    for i, route in enumerate(("Orthonasal", "Retronasal")):
        if route not in grouped:
            continue
        y = grouped[route][0]
        offset = (1-i)*4
        ax.plot(times, y+offset, lw=.9, c=PALETTE[route])
        peak_sign = "positive" if peak_polarity == "unknown" else peak_polarity
        candidates, _ = detect_peaks(y, polarity=peak_sign)
        good = candidates[(times[candidates] >= -.5) & (times[candidates] <= 2.5)]
        marker_top = offset + 2.5
        ax.vlines(times[good], marker_top, marker_top+.28, color=PALETTE[route], lw=.7)
        ax.text(.02, .95-i*.14, route, transform=ax.transAxes, va="top", fontsize=8,
                color=PALETTE[route])
    ax.axvline(0, lw=.8, color="black")
    ax.set(xlim=(-.5, 2.5), xlabel="Time from DIN (s)", ylabel="Standardized pressure + offset")
    ax.set_yticks([])
    scale_bar(ax)
    ax.set_title("A   Respiration traces / candidate peaks", loc="left", weight="bold")

    # B: Optional actual sorted-spike/respiration spectra imported from stage 03.
    # If sorted units aren't available, label LFP curves as LFP, NEVER firing rate.
    ax = axs[1, 0]
    if spectra_csv is not None:
        with Path(spectra_csv).open(newline="") as source:
            rows = list(csv.DictReader(source))
        has_spiking = any(row["signal_type"] == "firing_rate" for row in rows)
        allowed = {"respiration", "firing_rate"} if has_spiking else {"respiration", "lfp"}
        grouped_psd = {}
        for row in rows:
            if row["signal_type"] not in allowed:
                continue
            key = (row["route"], row["signal"])
            grouped_psd.setdefault(key, []).append((float(row["frequency_hz"]), float(row["normalized_shape"])))
        if not grouped_psd:
            raise ValueError("No usable respiration or neural PSD rows in " + str(spectra_csv))
        palette = {"Respiration":"#666666", "GC firing rate":"#6375a1",
                   "PC firing rate":"#b27c70", "GC LFP":"#6375a1","PC LFP":"#b27c70"}
        for (route, label), values in sorted(grouped_psd.items()):
            values.sort()
            ax.plot([x for x, _ in values], [y for _,y in values],
                    label=f"{route} {label}", color=palette.get(label, ".3"),
                    ls="-" if route == "Ortho" else "--", lw=1.3)
        ax.set(xlim=(0, 12), xlabel="Frequency (Hz)", ylabel="Unit-area PSD (1–12 Hz)")
        ax.set_title("B   Respiration / firing rate" if has_spiking else
                     "B   Respiration / LFP (not firing)", loc="left", weight="bold")
    else:
        for route, stack in grouped.items():
            psds = []
            for x in stack:
                f, p = signal.welch(x, fs=FS, nperseg=min(1024, len(x)), noverlap=512)
                psds.append(p)
            power = np.mean(psds, axis=0)
            ax.plot(f, power, color=PALETTE[route], label=f"{route} (n={len(stack)})", lw=1.5)
        ax.set(xlim=(0, 12), xlabel="Frequency (Hz)", ylabel="Pressure power (Z²/Hz)")
        ax.set_title("B   Respiration power spectrum", loc="left", weight="bold")
    ax.legend(frameon=False, fontsize=7)

    # C: pressure time-frequency power, explicitly NOT respiration-spike coherence.
    ax = axs[0, 1]
    all_stfts = []
    for stack in grouped.values():
        for x in stack:
            f, ts, z = signal.stft(x, fs=FS, nperseg=256, noverlap=224, boundary=None)
            all_stfts.append(np.abs(z)**2)
    # Use all available trials in descriptive pooled power figure.
    pwr = np.mean(all_stfts, axis=0)
    mask = f <= 12
    img = ax.pcolormesh(times[0]+ts, f[mask], 10*np.log10(np.maximum(pwr[mask], 1e-12)),
                        cmap="viridis", shading="auto", rasterized=True)
    ax.axvline(0, color="black", lw=1)
    ax.set(xlim=(-.5, 2.5), ylim=(0, 12), xlabel="Time from DIN (s)", ylabel="Frequency (Hz)")
    ax.set_title("C   Respiration power (not coherence)", loc="left", weight="bold")
    fig.colorbar(img, ax=ax, label="Power (dB)")

    # D: route-specific, *descriptive* pressure-band time courses.
    ax = axs[1, 1]
    for route, stack in grouped.items():
        curves = []
        kernel = np.ones(201)/201
        for x in stack:
            # Mean-squared standardized pressure, smoothed over 200-ms windows.
            curves.append(np.convolve(x*x, kernel, mode="same"))
        curves = np.asarray(curves)
        avg = curves.mean(axis=0)
        sem = curves.std(axis=0, ddof=1)/np.sqrt(len(curves)) if len(curves)>1 else np.zeros_like(avg)
        ax.plot(times, avg, color=PALETTE[route], lw=1.4, label=f"{route}")
        ax.fill_between(times, avg-sem, avg+sem, color=PALETTE[route], alpha=.18)
    ax.axvline(0, color="black", lw=.8)
    ax.set(xlim=(-.5, 2.5), xlabel="Time from DIN (s)", ylabel="Smoothed Z² (± SEM)")
    ax.set_title("D   Pressure-energy by route", loc="left", weight="bold")
    ax.legend(frameon=False, fontsize=8)
    note = "Peaks are candidates; inhalation polarity not independently confirmed" if peak_polarity=="unknown" else (
        f"Candidate peak polarity: {peak_polarity}; physiological verification required")
    fig.suptitle("Stimulus-aligned respiration | "+note, fontsize=10)
    output = Path(output)
    if output.suffix.lower() not in (".pdf",".png",".svg"):
        raise ValueError("Output must be .pdf, .png or .svg")
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=300)
    plt.close(fig)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--adc", required=True)
    p.add_argument("--trial-info", required=True)
    p.add_argument("--output", required=True)
    p.add_argument("--raw-fs", type=int, default=30000)
    p.add_argument("--peak-polarity", choices=["unknown","positive","negative"], default="unknown")
    p.add_argument("--spectra", default=None, help="Optional spectra.csv from 03_neural_spectra.py for panel B")
    args=p.parse_args()
    snippets, t=get_snippets(args.adc, trials_from_csv(args.trial_info), args.raw_fs)
    figure(snippets,t,args.output,args.peak_polarity,args.spectra)
    print(args.output)


if __name__ == "__main__":
    main()
