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

After independent polarity verification, replace the unknown label in a NEW run:

```bash
read -r -p "Verified inhalation deflection (negative or positive): " INHALE_SIGN
read -r -p "Independent evidence for polarity: " POLARITY_NOTE
read -r -p "Evidence from manual waveform/breath review: " WAVEFORM_NOTE
RESP_VERIFIED="$PWD/results/TG22_Feb2_respiration_reviewed_$(date +%Y%m%d_%H%M%S)"
coupling-respiration-qc "$SESSION" configs/tg22_feb2.example.json "$RESP_VERIFIED" --inhalation-sign "$INHALE_SIGN" --polarity-note "$POLARITY_NOTE" --waveform-review-note "$WAVEFORM_NOTE"
```

Only supply these notes after verifying them against acquisition documentation,
calibration or independently observed breaths. An attractive alignment plot is
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

The first 0.5–20 Hz spectrum has spectral flatness 0.9994 (near 1 indicates a flat spectrum), and candidate cycles are often fast (median 11.1 Hz across all candidates). These values are descriptive diagnostics, not validated rejection thresholds. This
raises a concern about noise, weak catheter signal or ADC identity, and prevents
biological timing interpretation. Successful automatic detection does not resolve
that concern. Confirm the ADC assignment, nasal catheter patency/pressure chain,
polarity and observable breaths before aligning neural data to these candidates.
The detector cutoffs were not relaxed to make this session pass.

## Literature and next gate

Samuelsen & Fontanini [1] used respiration to distinguish intraoral odor timing
from early oral responses. This motivates respiratory review for retro as well
as ortho; it does not prescribe this detector or numerical thresholds.

1. Samuelsen CL, Fontanini A (2017). *Processing of Intraoral Olfactory and Gustatory
   Signals in the Gustatory Cortex of Awake Rats.*
   [doi:10.1523/JNEUROSCI.1926-16.2016](https://doi.org/10.1523/JNEUROSCI.1926-16.2016).
2. Intan Technologies. *RHD2000 Data File Formats.*
   [Official format documentation](https://www.intantech.com/files/Intan_RHD2000_data_file_formats.pdf).

Next: resolve ADC identity/signal quality/polarity, manually evaluate a fixed
sample of trials across all six stimuli, then test detector agreement with those
annotations. Only afterward add inhale/exhale-aligned LFP or pPC spike reliability.
Full sorted spike arrays, unit-quality information and verified unit-region mapping
are needed for the latter and were not supplied in this increment.
