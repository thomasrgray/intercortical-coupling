# Validation status

Thirteen unittest cases pass (including two saved-panel comparison tests). They cover exact DIN onset/offset detection across
chunk boundaries; unique-six selection and ambiguous/unequal count rejection;
chunked versus full polyphase resampling and anti-alias rejection; flat/outlier
and zero-MAD QC; raw histogram metrics; same-region peer handling; STFT grid and
baseline normalization; zero baseline variance; sample-rate validation; and power
pipeline rejection masks, channel counts, numerical output and overwrite refusal.

All eight real uploaded DIN streams were scanned: DINs 0–5 have 30 trials each;
DINs 6–7 have 90 pulses each. Only 0–5 were selected.

A separate end-to-end synthetic run used 16 amplifier channels (8 per region),
3 kHz raw sampling, timestamp/header files and eight DIN streams. Raw processing,
QC, mean CAR, cleaned archive creation and power plotting completed. It produced
12 spectrograms, 12 pre/post plots and 12 plotted-value CSVs. Representative channel
QC and spectrogram images were visually inspected. Synthetic output is not study data.

Initial validation used only two biological amplifier channels. Subsequently,
the complete TG22 Feb2 session executed locally and generated power/QC outputs
(23 GC and 15 PC channels retained). Biological validation remains incomplete:
review exclusions, reference membership, artifact propagation and sensitivity
before manuscript inference. Continuous filtering can spread transients; artifact
interval dilation has not been validated. Coherence, Granger, respiration alignment
and grand averages remain unimplemented. The provisional power baseline is
[-1.5,-0.5) seconds. See `raw_processing.md` for current assumptions and commands.

## Subsequent single-session execution

PR #2 added two comparison tests, bringing the suite to 13 passing tests.
User-generated TG22 Feb2 full-session power outputs and QC tables were inspected:
23 GC and 15 PC channels retained; GC IA/oIA comparison exports rendered. This
supersedes the earlier statement that no full biological session had executed,
but does not complete biological validation. Histological mapping, exclusions,
reference sensitivity and filter artifact spreading remain to be reviewed.
The README includes executable run/QC/comparison examples and primary references;
the cited physiological studies do not prescribe our provisional QC thresholds.
