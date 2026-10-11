Current implementation: see `raw_processing.md`. Remaining steps below describe scientific validation and future coupling analyses.

# LFP build order and acceptance gates

1. **Figure 6B panel renderer — implemented.** Validate baseline normalization,
   STFT parameters, frequency/time axes, rejection masks, shared color limits,
   finite-input failure, zero-baseline handling. Synthetic run only so far.
2. **Raw preparation and Figure 6A — initial implementation available.** Inventory every layout electrode;
   validate sample rates/lengths and timestamp continuity; cross-check DIN events.
   Define reference groups explicitly, pre-screen bad channels, anti-alias before
   resampling, filter continuous data with adequate context, preserve originals.
   Conventional CAR and blech_clust's median/MAD-normalized reference are different
   choices: record which is used and run reference sensitivity comparisons.
3. **Requested QC — initial implementation available.** Require finite/nonconstant signals. Flag robust outliers using
   abs(x-median)/(1.4826*MAD)>10 and fail fractions >1%. Handle zero MAD explicitly.
   Flat fraction uses sustained identical-value runs: provisional minimum 10 ms,
   not isolated adjacent equal samples. Retain both metrics. Compute raw and LFP QC
   separately; apply sustained-flat detection before filters can hide dropouts.
   Trial flags use the same checks in each event window. Compare median voltage,
   SD and artifact-trial fraction with same-region peers (>3 scaled MADs fails).
   Insufficient peers and zero peer MAD must produce explicit review states,
   not an epsilon trick or silent acceptance. Freeze masks before spectra.
4. **Session power on real data.** Generate all six conditions for each region;
   inspect excluded-channel examples and accepted-trial counts. Power is not
   evidence for coupling. Show absolute pre/post PSD later as a separate product.
5. **Figure 6C coherence.** Prespecify estimator, epoch duration and frequency
   resolution; calculate per-trial cross/auto spectra using multiple segments or
   tapers (a single Fourier estimate yields trivial coherence). Control unequal
   epoch lengths/trial counts. Compare matched BA/oBA and IA/oIA. Save channel-pair
   identity and trial-shuffled nulls; summarize session then animal. Consider
   imaginary coherence as a sensitivity analysis, not a universal correction.
6. **Figure 6D Granger.** Use time-domain cleaned data, not power spectrograms.
   Prespecify multivariate model, model-order selection, stationarity checks,
   stability, residual whiteness and held-out validation. Check reference,
   filtering/downsampling and unobserved/common-drive sensitivity. Report GC→PC
   and PC→GC separately. No Granger implementation is claimed yet.
7. **grand_average_LFP.py — planned.** Read only normalized session outputs;
   reject duplicate sessions or incompatible alignment/reference/filter/normalization
   settings or axes. Average sessions within animal, then animals by default.
   Bootstrap animals for uncertainty; optional equal-session descriptive summary
   must be labeled. Never pool trials as independent animals or rewrite recordings.

Open methods decisions: reference groups/method; final baseline/pre/post windows;
minimum trials/peers; treatment of S/BK controls and opto trials; how nasal
inspiration and retro neural-state alignment enter sensitivity analyses.
