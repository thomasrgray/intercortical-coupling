Continue developing https://github.com/thomasrgray/intercortical-coupling for my manuscript "Differential inter-cortical coupling in orthonasal and retronasal processing."

Read README.md, docs/raw_processing.md, docs/validation.md and the current open PR/branch before editing. Avoid re-auditing unrelated legacy repositories. Preserve existing code, raw data and uncommitted user changes. Work in small tested increments; communicate briefly. I use PyCharm to pull, run, commit and push.

The repo is organized by eight manuscript figures. LFP is the first focus. New code provides raw Intan DAT processing, DIN discovery, channel/trial QC, explicit electrode/region mapping, conventional within-region mean CAR, anti-alias resampling to 1 kHz, continuous filtering, trial extraction, and power/QC plots. Check current implementation/status rather than assuming a PR was merged.

Inputs: amp-*.dat signed int16, info.rhd, time.dat, board-DIN-*.dat unpacked uint16. Config explicitly assigns filenames/electrode IDs to GC/PC and reference groups, and DINs to stimuli. Usually 0–31 GC, 32–63 PC; never infer a half split after selecting channels. No spike sorting is needed for LFP; include electrodes without units. Optional ADC-00 is nasal pressure, not used by current TTL-aligned LFP processing. Never rewrite source recordings or existing blech_clust HDF5 files.

Select the UNIQUE positive pulse-count group containing EXACTLY six DINs, then require the configured labels to match. Stop on ambiguity or mismatched counts; never truncate trials. Confirmed real example: DINs 0–5 each 30 pulses; 6–7 each 90. Mapping 0 BA, 1 S, 2 oIA, 3 IA, 4 oBA, 5 BK. Matched odor contrasts BA/oBA and IA/oIA; S and BK are separate controls.

STFT power: Fs 1000, Hann 500 samples, overlap 499, <=70 Hz. Z-score per trial/channel/frequency using baseline-only moments, then average trials and channels. Turbo; plot -1000…1750 ms. Current provisional baseline [-1.5,-0.5) s; snippets -2…+5 s. Save pre/post curves and normalized session arrays for future animal-level grand averages.

QC: finite/nonconstant, flat fraction >1%, robust-outlier fraction >1% with abs(x-median)/(1.4826*MAD)>10; same-region peer deviations >3 scaled MAD for median, SD, artifact-trial fraction. Flat means sustained identical runs >=10 ms, not isolated repeated ADC values. Review zero-MAD/insufficient-peer handling, fixed CAR donors, trial exclusions and filter artifact spreading.

Critical upstream finding: blech_clust CAR modifies HDF5 /raw but legacy LFP extraction reads original amplifier DAT, bypassing it. Our mean CAR differs from blech_clust's median/MAD-normalized reference. Record provenance; never double-reference unknown inputs.

Next priority: run the new processing on my COMPLETE local 64-channel session, inspect QC/exclusions and power plots, resolve any failures without loosening thresholds blindly. Only two raw amplifier channels were supplied to the previous chat, so complete biological CAR/peer QC has NOT been validated. DIN scanning was validated on all eight real uploaded streams. Coherence, Granger, respiration alignment, HMMs, spike coupling and grand averages remain future work. Ask for local command output only when necessary; do not request re-uploading an entire recording when local execution suffices.
