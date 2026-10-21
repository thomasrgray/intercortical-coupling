# Respiration / sniffing processing (Dewan-inspired)

Independent workflow alongside `lfp_power` and `lfp_phase_coherence`.

Adapted conceptually from [Dewan Lab sniffing](https://github.com/OlfactoryBehaviorLab/sniffing), especially their filtered/z-scored traces, peak counts, first-three-peak latencies, half-prominence durations, and moving-window sniff counts. This implementation reads **Intan RHD unsigned board-ADC-00.dat** and **blech_clust trial_info_frame.csv**, not DewanH5. It reuses `coupling.respiration.downsample_counts` and `detect_cycles` for negative/positive zero-crossing candidates.

## Example, from the repository root

```bash
python -m pip install -e .
python respiration_sniffing/01_process.py \
  --adc "/media/thomas/Data/TG22/Feb2_OrthoRetro_TG22_260202_170014/board-ADC-00.dat" \
  --trial-info "/path/to/trial_info_frame.csv" \
  --raw-fs 30000 \
  --output "results/sniff_TG22_Feb2_$(date +%Y%m%d_%H%M%S)"
```

Or `python -m coupling.sniffing --adc ... --trial-info ... --output ...`.

**Before physiological interpretation**: default `--polarity unknown` produces positive/negative peak and zero-crossing *candidates*, not inhalation labels. Inspect the raw waveform and sensor wiring to establish physiological polarity. Only then use `--polarity positive` or `--polarity negative` to populate verified inhalation timing columns. In the present implementation positive/negative specify the sign of the **peak**; crossing labels are independent; do not assume that a positive peak implies a positive-going start crossing without validating polarity and waveform morphology.

## Outputs

* `trials.csv` — stimulus, DIN number, absolute event sample, pulse width, peak-based counts, first three peak latencies, half-prominence widths, candidate zero-crossing timing, QC status.
* `binned_sniff_counts.csv` — overlapping 100-ms windows at 50-ms increments, per trial and candidate peak sign.
* `cycles.csv` — candidate negative-to-positive pressure cycles from existing crossing detector, with screening flags.
* `session_summary.csv` — per-stimulus trial counts, nominal route labels, mean pulse widths and candidate timing.
* `respiration_by_stimulus.pdf/png` — mean ± SEM z-scored waveform per stimulus (DIN-onset aligned).
* `provenance.json` — source filenames, sampling rate, filters, polarity and eligible trial counts.

## Methodological differences and safeguards

* Dewan code uses an initial Chebyshev bandpass 0.01–100 Hz plus a 1-Hz high-pass and per-trial z-score. We use the repository's existing 0.5–20 Hz fourth-order Butterworth offline filter, then a 1-Hz high-pass and trial z-scoring for **Dewan-comparable peak detection**. The Maier-style candidates remain based on the 0.5–20-Hz filtered signal. Exact numeric equivalence is **not** claimed.
* Positive and negative peak candidates are kept separately. The existing `detect_cycles` uses negative-going crossings before troughs; this is a *candidate* onset detector until pressure polarity is verified.
* Reported peak widths are SciPy half-prominence widths. They are **not** full inspiration duration. The standalone Dewan code's midpoint measure is related but not identical.
* DIN timestamps are in **raw Intan samples**, not milliseconds. We convert samples to the downsampled 1-kHz timebase. The example provided had 180 trials across six repeated stimulus DINs; pulse lengths can differ sharply between ortho and retro.
* Do not treat a Bpod final-valve event or DIN onset as physical odor arrival. PID recording and separate calibration are required for that.
* The input CSV and ADC must originate from the exact same session/clock. The code rejects out-of-bounds events but cannot prove that files belong together. Inspect the QC plot and provenance before comparing groups.
* Summary rows are descriptive. Across-animal statistical tests must use animals as independent units, not channels, sniffs, or repeated trials.

No existing LFP workflows are modified by this addition.

## Reference-style A–D figure with 200-ms scale bar

The attached literature panel motivated an optional paper-style 2×2 figure. Generate it with:

```bash
python respiration_sniffing/02_figure.py \\
  --adc /path/to/board-ADC-00.dat \\
  --trial-info /path/to/trial_info_frame.csv \\
  --raw-fs 30000 \\
  --output results/respiration_reference_style.pdf
```

Panels: **A** shows a single example orthonasal/retronasal pressure trace each, with candidate peak ticks and an exact **200 ms horizontal scale bar** in data coordinates; **B** plots route-specific Welch respiration power spectra; **C** shows the pressure time–frequency *power* spectrogram; **D** overlays route-specific smoothed standardized pressure energy ± trial SEM. These plots match the reference's compact A–D layout but **do not calculate neural coherence**: the provided recording has pressure ADC and DIN timestamps, not neural spike times. Coherence requires synchronized valid spiking/LFP events and an explicitly defined estimator and controls. This figure deliberately labels C as power and D as energy so it cannot be mistaken for a respiration–spiking coherence analysis. Route traces are illustrative; use animal-level replication for inference.

## Respiration–LFP and sorted-spike firing-rate spectra

Use `03_neural_spectra.py` to compare DIN-matched nasal respiration, QC-cleaned GC/PC LFP and (optionally) true sorted-unit GC/PC firing-rate spectra. The reference-style `02_figure.py` accepts `--spectra results/RUN/spectra.csv` to populate panel B without relabeling LFP as spiking. Detailed methods, paths, requirements and limitations: [neural_spectra.md](neural_spectra.md).
