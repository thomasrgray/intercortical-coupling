# Validation of this starter

Five unittest cases pass: baseline mean/SD, STFT grid, later-response independence
of baseline, zero baseline variance, finite/sample-rate validation, plus a pipeline
case checking rejection masks, channel counts, numerical output and overwrite refusal.
The synthetic demonstration generated 12 spectrograms, 12 pre/post plots and 12
plotted-value CSVs. A representative spectrogram was visually inspected.

No full biological session has passed this pipeline. No raw preparation, CAR,
coherence, Granger or grand-average implementation is claimed in this release.
The input NPZ must already be prepared and QC-masked. The provisional baseline
is [-1.5,-0.5) seconds, explicitly documented for future revision.
