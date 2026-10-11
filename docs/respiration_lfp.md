# Optional respiration-aligned LFP comparison

This extension reuses `lfp_unreferenced.npy`, the original referenced DIN trials,
channel masks and fixed CAR donors from a completed `coupling-prepare` run. It
reuses the anti-aliased continuous ADC cache from `coupling-respiration-qc`.
Raw amplifier filtering is not repeated. Newly exposed raw amplifier segments
are still read for shifted-epoch quality control. Existing DIN outputs are unchanged.

## Method and interpretation

Maier, Wachowiak & Katz (2012), *Chemosensory convergence on primary olfactory
cortex*, **Measurement and analysis of respiration**, explicitly define inhalation
onsets as **upward zero crossings of the 0.5–20 Hz filtered respiration signal**.
Their Honeywell 24PCAFA6G / Grass P55 wiring gives positive voltage during
inspiration and negative voltage during expiration.
[Primary methods](https://pmc.ncbi.nlm.nih.gov/articles/PMC3525522/),
[DOI](https://doi.org/10.1523/JNEUROSCI.3540-12.2012).
A pressure trough is not this onset. Positive pressure *voltage* is not necessarily
positive physical pressure. Verify the actual differential wiring for this session.

For orthonasal odor, the example uses the first inspiratory onset after DIN (plus
any separately documented calibrated delivery delay). For intraoral odor it uses
expiration as a **hypothesis about retronasal sampling**, not as a measured arrival
or as a claim that Maier prescribed this first-expiration rule. Maier (2017)
discusses exhalation-driven intraoral olfactory stimulation:
[DOI](https://doi.org/10.1152/jn.00802.2016).
Sucrose uses an explicitly labelled intraoral-taste route and an expiratory
sampling reference only; it is not an odor-arrival estimate. BK is orthonasal blank.
Every stimulus rule can be changed independently, including inspiration for all
stimuli when testing that particular hypothesis.

The detector operates on the audit's sign-normalized cache. For a positive-sign
audit, the cached voltage was inverted, so a downward crossing in that cache is
an upward crossing of recorded voltage (inspiration). The verified example uses
physiological names and enforces the expected audit sign. The unknown-polarity
example uses recorded positive/negative candidate phases and enforces unknown
polarity. This prevents accidentally reusing a candidate rule with an inverted cache.

Anti-alias resampling comes from the audit (polyphase FIR with chunk context),
then configurable continuous fourth-order zero-phase Butterworth filtering and
prominent-trough/zero-crossing cycle screening. Rate range 1–15 Hz, minimum phase
10 ms, prominence 0.5 scaled MAD and search window 2 s are **provisional engineering
screens**, not rat physiology certification or thresholds specified by Maier.
Cycle CSV retains irregular-cycle reasons. Filtering has no causal group delay,
but can move crossings and spread transients. No timing correction is invented.
The 5 s recording-edge guard is conservative and unvalidated; interior artifact
spread remains a limitation. Full-session spectral diagnostics remain in the audit.

The rule is first screened onset at/after DIN + delay, within the configured
window. A phase already underway is not backdated; it waits for the next onset.
Missing events and incomplete epochs are excluded, never replaced with DIN.
Timing is quantized to the 1 kHz grid; source DIN timestamps remain unrounded.
ADC counts remain counts; no unverified volts/pressure conversion is applied.
`delay_s=0` means **no calibrated delay applied**, not proof of instantaneous arrival.
Nasal pressure alone cannot measure odor arrival, swallowing or catheter delay.

**TG22 caution:** the existing audit documents noise-like ADC spectra and a
similar shuffled-noise control. Smooth filtered oscillations and a sharp
phase-aligned pressure average do not demonstrate breathing. Until independent
waveform review establishes recognizable breaths, outputs stay PROVISIONAL.
This extension does not repair bad ADC data or relax thresholds to retain trials.

## Run locally

Install from the feature branch or after merging it:

```bash
python -m pip install -e .
SESSION="/media/thomas/Data/TG22/Feb2_OrthoRetro_TG22_260202_170014"
read -r -e -p "Completed coupling-prepare folder: " PREPARED
read -r -e -p "Completed unknown-polarity respiration audit folder: " RESP_AUDIT
python -m coupling.respiration_lfp "$PREPARED" "$RESP_AUDIT" configs/respiration_alignment.example.json
```

The two caches must refer to the same unchanged source, config, clock and event
identities. Stale file size/mtime, timestamp gaps, differing sampling rates,
ambiguous DIN selection, or changed event mappings fail explicitly. If needed,
create the audit first using the **same session config as preparation**:

```bash
RESP_AUDIT="$PWD/results/TG22_respiration_$(date +%Y%m%d_%H%M%S)"
python -m coupling.respiration "$SESSION" configs/tg22_feb2.example.json "$RESP_AUDIT"
```

Output defaults to the new session subfolder:

```text
/media/thomas/Data/TG22/Feb2_OrthoRetro_TG22_260202_170014/LFPs matched to brathing/
```

This explicit new subfolder is the exception requested to the older audit/raw
preparation convention of writing outside the recording. Source files remain
read-only. Existing output folders are never overwritten. For another comparison,
use `--output "$PWD/results/review2/LFPs matched to brathing"`.

After independently verifying polarity and reviewing actual breaths, create a
new audit with `--inhalation-sign positive --polarity-note "actual evidence"
--waveform-review-note "actual review evidence"`, then use
`configs/respiration_alignment.maier.json`. Do not copy placeholder notes as evidence.
The named sensor model alone does not establish usable waveform quality.

## Matched power and outputs

For each channel and event the same accepted pairs contribute to DIN and shifted
power. Shifted raw and filtered epochs are rechecked; CAR uses the original fixed
donors and rejects an event when a donor fails. Original channel/trial exclusions
remain enforced. GC and PC can retain different pairs; this is not a cross-region
matched-trial statistical test.

Both alignments use **the same DIN-relative baseline moments**, STFT centers
−1.5 to −0.5 s, independently per channel/trial/frequency. This prevents a shifted
baseline from including post-DIN responses. The existing Hann 500 ms window and
1 ms hop are reused; the latter is not 1 ms temporal resolution. Zero baseline
variance remains missing. Trials average within channel, then channels equally;
frequency-level counts are recorded. Shared display limits default to ±3 z.

- `GC_*_paired_power.png`, `PC_*_paired_power.png`: 12 side-by-side comparisons.
- `*_respiration.png`: six example traces with phases/DIN/selected event,
  event-aligned pressure heatmaps and DIN offsets, and latency/retention summaries.
- `trial_alignment.csv`: every original trial, route, timestamps, delay, rule,
  bounded search, provisional flag and exclusions. `retained` means technical
  screening only, not physiological validation. Missing latencies stay empty.
- `paired_channel_trials.csv`: auditable per-channel/event retention and reasons.
- `pressure_cycles.csv`: screened and rejected cycle candidates (if any).
- `paired_power.npz`: numerical paired panels and frequency-level channel counts.
- `provenance.json`, `caption.txt`: configuration, sources, review notes, cached
  processing provenance, source stats, small-input hashes, revision and limitations.

Empty conditions produce clearly labelled empty panels rather than silently
substituting another alignment. Pressure heatmaps standardize amplitude only for
display; sharpening after self-alignment is expected even for filtered noise.

## Validation

23 unit tests pass, including a synthetic raw ADC/DIN + cached-LFP end-to-end run
that produces all 18 figures, and tests for known phase timing, sign handling,
no-breath/search bounds, resampling continuity, clock mismatch, shifted-tail
artifacts, group-wide CAR donor failures, original trial exclusions, boundaries and shared-baseline arithmetic.
No complete biological amplifier session was available in this workspace;
real-data validation and manual breath annotations remain outstanding.
