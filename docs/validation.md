# Validation status

Eleven unittest cases pass. They cover exact DIN onset/offset detection across
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

No complete biological session has passed this pipeline: only two raw amplifier
channels were supplied. Run on a complete local recording and inspect exclusion
reasons, reference membership, artifact propagation and sensitivity before using
results in the manuscript. Continuous filtering can spread transients; artifact
interval dilation has not been validated. Coherence, Granger, respiration alignment
and grand averages remain unimplemented. The provisional power baseline is
[-1.5,-0.5) seconds. See `raw_processing.md` for current assumptions and commands.
