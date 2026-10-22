# Respiration / sniffing processing (Dewan-inspired)

Process nasal pressure and generate respiration power, GC/PC LFP power, and LFP–respiration coherence for one session. This folder sits alongside `lfp_power` and `lfp_phase_coherence`.

## Start here: copy and paste the whole block

Run this in Bash with your analysis Python environment active. It uses the TG22 Feb2 session from the completed LFP–respiration spectra run, and includes the standalone sniffing analysis in the same sequence. The standalone command also completed on TG23 Feb1 in the supplied session log; the combined block below has not itself been rerun on the recordings.

The repository is assumed to be at `$HOME/intercortical-coupling`. For this TG22 example, the session needs raw `amp-*.dat` files, stimulus DIN files, `board-ADC-00.dat`, and `trial_info_frame.csv`. The sniffing step uses the CSV; LFP preparation and the respiration audit use the raw session and JSON configuration.

```bash
(
  set -e

  cd "$HOME/intercortical-coupling"
  git switch master
  git pull --ff-only
  python -m pip install -e .

  # Session settings: keep the ADC, trial CSV, and configuration matched.
  SESSION="/media/thomas/Data/TG22/Feb2_OrthoRetro_TG22_260202_170014"
  CONFIG="$PWD/configs/tg22_feb2.example.json"
  TRIAL_INFO="$SESSION/trial_info_frame.csv"
  RUN="$PWD/results/$(basename "$SESSION")_respiration_lfp_$(date +%Y%m%d_%H%M%S)"

  # Check required inputs before starting the expensive LFP preparation.
  test -d "$SESSION"
  test -f "$CONFIG"
  test -f "$SESSION/board-ADC-00.dat"
  test -f "$TRIAL_INFO"
  mkdir -p "$PWD/results"
  mkdir "$RUN"

  # 1. Standalone sniffing: pressure traces, peak/cycle candidates, and tables.
  python respiration_sniffing/01_process.py \
    --adc "$SESSION/board-ADC-00.dat" \
    --trial-info "$TRIAL_INFO" \
    --raw-fs 30000 \
    --polarity unknown \
    --output "$RUN/sniffing"

  # 2. Prepare cleaned GC/PC LFPs from the raw amplifier recordings.
  python -m coupling.prepare \
    "$SESSION" \
    "$CONFIG" \
    "$RUN/prepared"

  # 3. Prepare and audit respiration using the SAME session and configuration.
  python -m coupling.respiration \
    "$SESSION" \
    "$CONFIG" \
    "$RUN/respiration_audit" \
    --adc board-ADC-00.dat \
    --inhalation-sign unknown

  # 4. Combine the prepared signals into power and coherence figures.
  python -m coupling.respiration_spectra \
    "$RUN/prepared" \
    "$RUN/respiration_audit" \
    "$RUN/spectra" \
    --window 1 \
    --hop 0.1

  printf '\nCompleted. Outputs saved to:\n%s\n' "$RUN"
  printf '\nSniffing plots/tables: %s\nLFP–respiration figures: %s\n' \
    "$RUN/sniffing" "$RUN/spectra"
  if command -v xdg-open >/dev/null 2>&1; then
    xdg-open "$RUN/spectra" || true
  fi
)
```

The parentheses keep settings local to this run; `set -e` stops the sequence if a processing command fails. Each run gets a timestamped directory, and each stage writes to a new subdirectory. No intermediate path needs to be copied into the next command.

**For another session:** change `SESSION`, `CONFIG`, and (if named differently) `TRIAL_INFO` in the settings section. Use that session's verified channel/region and DIN/stimulus mapping; do not reuse the TG22 configuration unchanged for TG23 or another recording. Update `--raw-fs` and the JSON sampling rate together if acquisition was not at 30 kHz. For **sniffing only**, run through step 1 and omit steps 2–4; raw amplifier files and the LFP configuration are not needed by step 1.

## What a completed run looks like

The supplied TG22 spectra run processed all 64 amplifier files, wrote `prepared/cleaned_lfp.npz`, reported 180 respiration trials with 180 candidate-screen passes, and finished the spectra directory. The separate TG23 sniffing run printed `Processed 180 trials`.

An audit result of `"physiological_timing_trials": 0` is expected with `--inhalation-sign unknown`: the run can produce power and coherence without verified inhale/exhale labels. It does **not** establish physiological breathing validity. Review the pressure waveform and polarity before interpreting candidate timings.

| Directory inside the run | Main results |
| --- | --- |
| `sniffing/` | `respiration_by_stimulus.png/pdf`, peak/cycle tables, and provenance |
| `prepared/` | `cleaned_lfp.npz`, `events.csv`, and `preprocessing.json` |
| `respiration_audit/` | `pressure_1khz.npz`, `events_respiration.csv`, and `respiration_provenance.json` |
| `spectra/` | `BA_respiration_lfp.png/pdf`, `IA_respiration_lfp.png/pdf`, `spectra.npz`, and `spectra_0_2s.csv` |

For each matched odor, the BA and IA figures contain:

- **A:** first accepted pressure trial, DIN at zero, 200 ms scale bar.
- **B:** respiration and GC/PC LFP power spectra, unit-area normalized over 0–12 Hz for shape comparisons. Raw PSD density remains in `spectra.npz`.
- **C:** GC and PC LFP–respiration magnitude-squared coherence heatmaps with a shared 0–1 scale. These estimate pressure–LFP coherence, distinct from GC–PC coherence.

## Troubleshooting the command block

| Where it stops | What to check |
| --- | --- |
| Repository update | Resolve any Git error before rerunning; the block does not discard local work. |
| Input checks before step 1 | Confirm the session directory, JSON, ADC, and exact trial CSV filename. |
| LFP preparation | Confirm the raw amplifier/DIN files and session-specific channel/stimulus configuration. |
| Spectra input validation | Use `prepared/` and `respiration_audit/` from the same run and session. Both preceding stages must finish. |
| Output directory already exists | Rerun with a fresh timestamp; the spectra output directory must be new. |

## Standalone sniffing outputs and interpretation

Adapted conceptually from [Dewan Lab sniffing](https://github.com/OlfactoryBehaviorLab/sniffing), especially their filtered/z-scored traces, peak counts, first-three-peak latencies, half-prominence durations, and moving-window sniff counts. This implementation reads **Intan RHD unsigned board-ADC-00.dat** and **blech_clust trial_info_frame.csv**, not DewanH5. It reuses `coupling.respiration.downsample_counts` and `detect_cycles` for negative/positive zero-crossing candidates.

`python -m coupling.sniffing` is equivalent to the step-1 wrapper when given the same arguments.

**Before physiological interpretation**: default `--polarity unknown` produces positive/negative peak and zero-crossing *candidates*, not inhalation labels. Inspect the raw waveform and sensor wiring to establish physiological polarity. Only then use `--polarity positive` or `--polarity negative` to populate verified inhalation timing columns. In the present implementation positive/negative specify the sign of the **peak**; crossing labels are independent; do not assume that a positive peak implies a positive-going start crossing without validating polarity and waveform morphology.

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

## Spectra inputs and estimation details

The spectra stage reads the prepared and audit files listed above. It validates clock, source, original DIN/trial identity, onset sample, archive order and source-file freshness before combining signals. No raw recordings are distributed in this repository; run this on local session files.

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
