# Delivery/respiration audit: first single-session increment

## Question

Does a difference in TTL-to-pressure-phase timing provide an alternative
explanation for ortho/retro response timing? Competing explanations include
respiratory sampling, catheter/amplifier delay, oral stimulation, valve timing,
weak/noisy pressure, and genuinely different neural processing.

This increment audits the pressure record. It does not realign LFP/spikes,
measure swallowing/odor arrival, or establish neural response reliability.

## Copy-paste first run

From the repository root in PyCharm's terminal:

```bash
python -m pip install -e .
python -m unittest discover -s tests -v
SESSION="/media/thomas/Data/TG22/Feb2_OrthoRetro_TG22_260202_170014"
RESP_RUN="$PWD/results/TG22_Feb2_respiration_$(date +%Y%m%d_%H%M%S)"
coupling-respiration-qc "$SESSION" configs/tg22_feb2.example.json "$RESP_RUN"
```

Defaults deliberately keep polarity **unknown**. This produces a diagnostic
figure and candidate pressure phases, with empty inhale/exhale latency fields
and zero trials marked usable for physiological timing. This is a successful
audit, not a failure to detect breathing. No waveform or polarity is guessed.

If this recording used the same differential wiring as the Maier setup, the published polarity convention is positive-going inspiration. After confirming the wiring, use a NEW run and document that basis:

```bash
INHALE_SIGN="positive"
POLARITY_NOTE="Same Honeywell 24PCAFA6G and Grass P55 differential wiring as Maier et al. 2012; lab setup confirmation; paper defines inspiration as positive-going."
read -r -p "Evidence from manual waveform/breath review: " WAVEFORM_NOTE
RESP_VERIFIED="$PWD/results/TG22_Feb2_respiration_reviewed_$(date +%Y%m%d_%H%M%S)"
coupling-respiration-qc "$SESSION" configs/tg22_feb2.example.json "$RESP_VERIFIED" --inhalation-sign "$INHALE_SIGN" --polarity-note "$POLARITY_NOTE" --waveform-review-note "$WAVEFORM_NOTE"
```

Use the positive setting only if the differential wiring/orientation matches the
published setup. The separate waveform review note must describe actual inspection
of this recording, not the literature convention. An attractive alignment plot is
not independent verification. The tool records user review; it cannot certify
that the note is correct or replace trial-by-trial waveform inspection.

## Inputs and output identity

Inputs: `board-ADC-00.dat`, all `board-DIN-*.dat`, `info.rhd`, `time.dat`, and the
same session config used for raw preparation. Intan board ADC DAT is little-endian
uint16 [2]; amplifier DAT scaling must not be applied to this stream. We retain
centered ADC counts, avoiding invented voltage/pressure calibration. Board-mode
voltage scaling and the sensor/amplifier transfer function remain unverified.
Use `--adc board-ADC-01.dat` only if records identify that as the pressure channel.

All ADC/DIN sample counts must agree, timestamps must be consecutive, and config
sample rate/labels must agree. The unique six equal-count DIN selection is reused.
Event IDs use the same chronological ordering as raw preparation, with original
DIN and within-DIN trial indices preserved. Matching sample lengths/indexes does
not measure physical catheter, electronics or valve delay.

Outputs in a new folder outside the recording:

- `delivery_respiration.png`, `.pdf`, `.svg` and `caption.txt`.
- `events_respiration.csv`: all trials, independent first negative/positive phase
  starts after TTL, latency/rounding, raw QC, polarity/review status and timing
  eligibility. No missing event is imputed or shifted to zero.
- `pressure_cycles.csv`: candidate zero crossings, troughs, periods, prominence
  and screening reasons. Failed candidates remain in the inventory.
- `pressure_1khz.npz`, `din_inventory.csv`, `respiration_provenance.json`.

## Figure and controls

A: first in-bounds event chosen chronologically, not for effect size; anti-aliased
and filtered pressure have **separate y axes** and candidate phase markers.
B: pressure spectrum, useful for identifying tonal contamination or broadly
noise-like signals. C/D: per-trial standardized pressure averaged by condition,
with vertical offsets, under TTL and first negative-phase alignment.
E: first negative/positive phase delays (screen-passing/total counts).
F: candidate cycle rates with rejected candidates displayed.

Alignment to a phase detected from the same pressure trace necessarily sharpens
that trace, even for noise. It is not evidence of respiratory neural locking.
S remains sucrose; BK remains orthonasal blank. Neither replaces intraoral water.

## Provisional implementation choices

Read-only ADC counts are centered and anti-alias resampled to 1 kHz in bounded
30 s raw chunks with 1 s FIR context. Continuous fourth-order Butterworth 0.5–20 Hz
is applied forward/backward. There is no 60 Hz notch in this pressure stage.

Candidate negative peaks require prominence ≥0.5 scaled global MAD and ≥67 ms
spacing. A candidate cycle pairs a preceding negative-going zero crossing with
its trough and a following positive crossing before the next negative crossing.
Cycles outside 1–15 Hz or with either phase shorter than 10 ms are flagged.
These are explicit engineering screens, not literature-validated confidence
probabilities. Candidate rates should not be reported as validated sniff rates.

Each TTL is associated independently with the first screened negative and first
screened positive phase start within 2 s. If TTL occurs during an ongoing phase,
this is the **next phase start**, not the onset of current airflow. Missing starts
remain missing. The 2 s search interval censors later events; distributions apply
only to screened events in that interval and are not estimates of all breaths.

Trial screening uses raw −2 to +5 s windows: >1% rail samples or >1% sustained
flat samples (identical runs ≥10 ms) flags the trial. Recording-edge filter guard
is provisionally 5 s. Filter transient dilation around interior artifacts has
not been validated. These checks are insufficient to establish signal quality:
manual waveform, PSD, rate and polarity review are required. Zero-phase filtering
can spread/shift features; threshold crossings are operational pressure markers,
not precisely measured airflow onset.

## Real TG22 first audit

Uploaded ADC/header/timestamps and all eight DIN streams were inspected. The
recording has 102,814,080 samples at 30 kHz (3427.136 s); all lengths agree and
selected DINs yield 180 original trials. In unknown-polarity mode all 180 have
screened phase candidates, but **zero are approved for physiological timing**.

The first 0.5–20 Hz spectrum has spectral flatness 0.9994 (near 1 indicates a flat spectrum), and candidate cycles are often fast (median 11.1 Hz across all candidates). These values are descriptive diagnostics, not validated rejection thresholds. Tom confirms this session used the same Honeywell pressure sensor and Grass differential preamplifier chain as the Maier laboratory setup. Maier et al. (2012) specify Honeywell 24PCAFA6G and Grass P55 at 100×, and define inhalation as a positive-going pressure voltage deflection (upward zero crossing); exhalation is negative-going [3]. If the session's differential input wiring matches that setup, the expected inhalation polarity is therefore positive. The ADC header itself does not store sensor model or differential wiring, so this is lab-provided setup provenance plus the published polarity convention, not a fact recovered from the DAT header.

The flat full-session spectrum and similar shuffled-control segment spectra still raise a concern about signal quality. The named sensor/amplifier model is no longer an open question; remaining candidates include ADC input assignment/differential wiring, catheter patency, amplifier configuration, or a weak/noisy signal. The segment test shows resampling is not the cause. Before interpreting candidate latency, inspect the fixed segments and annotate recognizable breaths across conditions. Detector cutoffs were not relaxed.

## Literature and next gate

Samuelsen & Fontanini [1] used respiration to distinguish intraoral odor timing
from early oral responses. This motivates respiratory review for retro as well
as ortho; it does not prescribe this detector or numerical thresholds.

1. Samuelsen CL, Fontanini A (2017). *Processing of Intraoral Olfactory and Gustatory
   Signals in the Gustatory Cortex of Awake Rats.*
   [doi:10.1523/JNEUROSCI.1926-16.2016](https://doi.org/10.1523/JNEUROSCI.1926-16.2016).
2. Intan Technologies. *RHD2000 Data File Formats.*
   [Official format documentation](https://www.intantech.com/files/Intan_RHD2000_data_file_formats.pdf).
3. Maier JX, Wachowiak M, Katz DB (2012). *Chemosensory convergence on primary
   olfactory cortex.* The paper specifies the Honeywell 24PCAFA6G, Grass P55 at
   100×, and inspiratory positive-going voltage deflection.
   [doi:10.1523/JNEUROSCI.3540-12.2012](https://doi.org/10.1523/JNEUROSCI.3540-12.2012).

Next: use the published positive-going inspiration convention only if the
recording differential wiring matches the lab-confirmed setup; manually annotate
a fixed trial sample across all six stimuli, then test detector agreement with
those annotations. Only afterward add inhale/exhale-aligned LFP or pPC spike reliability.
Full sorted spike arrays, unit-quality information and verified unit-region mapping
are needed for the latter and were not supplied in this increment.

## Independent ADC subsample diagnostic

`python -m coupling.adc_check ADC HEADER OUTPUT [--seconds 30]` reads only
fixed early/middle/late segments (20–120 s each). It applies the same fourth-order
0.5–20 Hz zero-phase filter directly at the raw rate, then decimates to 1 kHz;
the comparison path uses the production polyphase resampler before filtering.
Five seconds at each segment edge are discarded. This is a numerical sensitivity
check, not a replacement timing detector. Raw-rate filtering provides anti-alias
attenuation before the diagnostic decimation; no bare stride sampling is used
as an unfiltered spectrum. The current diagnostic requires integer multiples of
1 kHz. Output folders must be new and outside the ADC source directory.

A fixed-seed permutation of raw samples retains the amplitude distribution while
removing temporal ordering. It demonstrates that smooth cycles can arise from
filtering noise. It is not a physiological reference signal or a valid surrogate
for route-dependent neural coupling. Plot raw/filtered counts on separate axes;
filtered recorded/control panels use matching limits per row. Do not compare
amplitudes across different axes. JSON retains source metadata, exact segment
sample offsets, SHA256 fingerprints, and both processing-path metrics.

TG22 uploaded-input check: RHD v1.5, board mode 0, 30 kHz, only ADC-00 enabled;
its custom name is also ADC-00 and notes are empty. The user confirms the physical Honeywell/Grass P55 chain; the header itself does
not encode that sensor or wiring. The published matching setup defines inspiration
as positive-going. Uint16 little-endian agrees with Intan's
per-channel ADC format [2]; changing interpretation to produce a plausible trace
would be unjustified. Counts are centered only for display, not calibrated pressure.
For the fixed 30 s segments, processing-path correlations were 0.99923, 0.99950,
0.99919; RMS differences were 4.00%, 3.20%, 4.08% of direct-filter SD. No rail
samples occurred. Segment spectral flatness was 0.875, 0.895, 0.867, versus shuffled
controls 0.890, 0.891, 0.867. Short-segment Welch estimates average fewer windows
than the full-session PSD, so these values are not directly comparable to the
full-session 0.9994. No flatness cutoff or equivalence test was applied.

This reduces concern that resampling alone created the observed signal. The confirmed
sensor/amplifier model narrows the hardware question, but DAT metadata do not establish
that ADC-00 receives the signal with the expected polarity, nor that its current
waveform has usable breath SNR. Review independently identified breaths
before interpreting any candidate latency. Neural coupling,
respiratory neural locking, swallowing and odor arrival remain unresolved.
