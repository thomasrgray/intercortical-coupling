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

## Respiration power and LFP–respiration coherence (reference panels A–C)

Use the **same session/configuration** for `coupling-prepare` and
`coupling-respiration-qc` first (see the repository's respiration audit docs).
Then run:

```bash
python -m coupling.respiration_spectra \
  /path/to/prepared_session /path/to/respiration_audit \
  results/respiration_lfp_spectra --window 1 --hop 0.1
```

The output directory must be new. The workflow reads `cleaned_lfp.npz`,
`events.csv`, `preprocessing.json`, and the audit's `pressure_1khz.npz`,
`events_respiration.csv`, `respiration_provenance.json`. It validates clock,
source, original DIN/trial identity, onset sample, archive order and source-file
freshness before combining signals. No raw recordings are distributed in this
repository; run this on local session files.

`BA_respiration_lfp.{png,pdf}` and `IA_respiration_lfp.{png,pdf}` compare matching
odors across ortho/retro columns:

- **A:** first accepted pressure trial, DIN at zero, 200 ms scale bar.
- **B:** respiration and GC/PC LFP power spectra, unit-area normalized over
  0–12 Hz for shape comparisons. Raw PSD density remains in `spectra.npz`.
- **C:** GC and PC LFP–respiration magnitude-squared coherence heatmaps, with a
  shared 0–1 scale. These are pressure–LFP estimates, distinct from GC–PC coherence.

All six stimuli are processed separately and exported to `spectra.npz` and
`spectra_0_2s.csv` when eligible, including sucrose and blank controls. Plots keep
BA and IA separate; they do not pool distinct odors into a route average.
CSV summaries average window centers from 0 to 2 s; windows at zero straddle
stimulus onset. Pressure QC exclusions and accepted trial counts per channel
are exported. LFP channel/trial masks are applied, and >=3 matched trials per
channel/condition are required. Three is a computational minimum, not a claim
of reliable estimation. Unequal retained trial counts can bias route contrasts;
inspect counts before interpreting differences.

The estimator uses Hann-tapered Fourier transforms of each trial in sliding
windows, averages cross/auto spectra **across trials before** taking
`abs(mean(Sxy))**2 / (mean(Sxx)*mean(Syy))`, then averages channel coherence within
region. It never concatenates trials or averages common-average-referenced
regional LFP waveforms (which can cancel). A 1 s window gives 1 Hz frequency
bins and 1 s temporal resolution; the 0.1 s hop does not improve that resolution.

Pressure power uses the anti-aliased unfiltered ADC cache (counts²/Hz); only the
example waveform uses the audit's 0.5–20 Hz filter. LFP power is µV²/Hz and
inherits the prepared 1–100 Hz filter, so sub-1-Hz content is attenuated. Unknown
polarity permits spectra but no inhale/exhale labels. Missing screened cycles
alone do not exclude a trial; signal/boundary QC failures do. Pressure validity
still needs waveform review. Coherence includes stimulus-locked common signals
and does not establish causal coupling. Single-session descriptive estimates
have no inferential error bands; no spike trains are needed. DIN zero is a
command timestamp, not measured odor arrival. Provenance records inputs,
hashes, filters and estimation settings.
