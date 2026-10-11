# Raw recording -> LFP power

## Run in PyCharm terminal

```bash
python -m pip install -e .
python -m unittest discover -s tests -v
coupling-scan-dins /path/to/session /path/to/results/din_scan
coupling-prepare /path/to/session configs/tg22_feb2.example.json /path/to/results/lfp_run --plot-power
```

The TG22 example explicitly maps C-000…031 to GC and C-032…063 to PC;
verify this mapping for each session. It uses conventional arithmetic mean CAR
within each region, not blech_clust's median/MAD-normalized CAR. This is a
configurable methodological choice, not a claim of equivalence. Use
`reference_method: none` for a deliberately unreferenced sensitivity run.
A `channels` list is sufficient: no blech_clust output is required. Alternatively,
`electrode_layout` can name a CSV containing filename, electrode_ind, CAR_group,
and optionally original_CAR_group; original groups supply anatomical labels.

## DIN selection

Scan every board-DIN-*.dat as unpacked uint16 binary 0/1. Count complete high
pulses with first-high onset/first-low exclusive offset. Select the unique positive
pulse-count group containing EXACTLY SIX DINs. Do not assume the largest/smallest
count or lowest six IDs. A seven-channel tie, two six-channel candidates, a missing
stimulus trial, nonbinary values, or boundary-truncated pulse stops processing.
`coupling-scan-dins` saves the inventory before selection errors so it is reviewable.
Labels must exactly cover the selected DINs. These rules cannot identify stimuli
by count alone; the explicit setup mapping supplies their meanings.

Actual uploaded example: DIN 0–5 each 30 pulses; DIN 6–7 each 90. The scanner selects
0–5. Pulse medians: 15 ms for 0/1/3; about 2 s for 2/4/5; 10/20 ms for 6/7.
No condition is truncated to equalize counts. Event IDs and boundary exclusions
are retained in events.csv. Onset rounding to 1 kHz is recorded per trial.

## Processing and QC

- Read header sample rate, validate config rate, all amplifier/DIN lengths and
  consecutive int32 timestamps. Timestamp resets/gaps must be resolved first.
- Evaluate exact raw integer histogram median/SD/outlier fractions and sustained
  identical-value runs. Inspect raw trial QC before filtering can conceal dropouts.
- Convert Intan signed int16 counts using 0.195 µV/count. Resample with anti-alias
  polyphase FIR, in 30-second chunks with 1-second overlapping FIR context.
- Filter the continuous 1 kHz trace: fourth-order Butterworth 1–100 Hz, zero phase;
  optional 60 Hz Q=30 notch (on by default). The raw-rate memory footprint is
  bounded; downsampled continuous/trial arrays use disk-backed NPY files.
- Extract -2 to +5 s snippets. Retain original event identities; mark boundary
  snippets rejected rather than shortening other conditions.
- Raw and unreferenced LFP channels/trials: finite, nonconstant, nonzero MAD;
  sustained-flat fraction >1% fails, robust-outlier fraction >1% fails using
  abs(x-median)/(1.4826*MAD)>10. Sustained flat = identical run >=10 ms (provisional).
- Same-region leave-one-out peer comparison of LFP median voltage, SD and initial
  bad-trial fraction: >3 scaled MAD fails. Require >=3 usable peers; insufficient
  peers exclude with a review reason. With zero peer MAD, equal values pass and
  unequal values fail; no arbitrary epsilon rescues them. This conservative rule
  can be harsh; review the metrics before adopting thresholds for the manuscript.
- Freeze good reference members. Mean CAR requires >=3 accepted donors per group
  and never mixes regions. The implementation subtracts the fixed group's mean
  from trial voltages AFTER identical linear filtering; this is algebraically
  equivalent to referencing continuously before those filters for fixed members.
- If any donor fails a trial, reject that trial for the entire reference group.
  Recheck referenced trial waveforms. No artifact replacement/interpolation.

Filters can spread artifacts in time; this first version checks raw and filtered
trial windows but does not yet dilate artifact intervals by a validated filter
settling time. Review boundary/nearby-artifact sensitivity before final inference.
The initial artifact-trial fractions drive peer QC; final post-CAR retention is
saved separately. Channel QC is not silently re-optimized after CAR.

## Outputs

DIN inventory/count plot; channels.csv; events.csv; channel_qc.csv and metric plot;
trial_qc.csv and retention plot; representative accepted traces; preprocessing.json;
lfp_unreferenced.npy; lfp_trials.npy; cleaned_lfp.npz. Original recordings are never
opened for writing. Archives duplicate some disk-backed arrays for convenient
plotting; reserve several GB per full session. NPZ plotting loads the trial tensor
into RAM once; approximately 323 MB at 64 × 180 × 7000 float32 samples, plus overhead.
Source file sizes/mtimes and source path are logged, not full raw-data content hashes.

With --plot-power: 12 region/stimulus spectrograms, 12 pre/post normalized frequency
plots and CSVs, plus baseline-normalized session arrays/provenance. Parameters:
1 kHz, Hann 500 samples, overlap 499, <=70 Hz, turbo, -1000…1750 ms. Baseline centers
[-1.5,-0.5) s; normalize each trial/channel/frequency before averaging. These remain
provisional analysis windows. Power is not a coupling estimate.

The nasal ADC and VDD streams are not required or consumed in this TTL-aligned
LFP stage. Breathing alignment is a separate upcoming component.
