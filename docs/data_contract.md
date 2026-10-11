# Data contract

## Original sources (read only)

Per-channel Intan `amp-*.dat`, `.info`, `.params`, electrode layout CSV,
`info.rhd`, timestamp stream, stimulus events and breathing ADC. A processed
blech_clust HDF5 is an optional separate source, never assumed equivalent to DAT.
Record exact pipeline version, reference method/groups and filtering history.
The audited blech_clust CAR modifies HDF5 `/raw`; the audited legacy LFP extractor
reads amplifier DAT files. Do not infer CAR merely from `use_auto_CAR=true`.

## Canonical tables for the planned preparation stage

- sessions: animal_id, session_id, source_path, raw_fs_hz, sample_count, source revision.
- channels: array_position, electrode_id, filename, region, reference_group, inclusion reason.
- events: event_id, DIN, stimulus, route, onset_sample, offset_sample,
  original_trial_index, within_stimulus_index, alignment source, opto condition.
- channel_qc: median_uV, SD_uV, finite status, nonconstant status,
  adjacent-repeat fraction, sustained-flat fraction, robust-outlier fraction,
  artifact-trial fraction, same-region peer deviations, exclusion reasons.
- trial_qc: event_id, electrode_id, accepted flag and reasons.

Never guess half/half regions or derive the LFP electrode set from sorted units.
Never shorten all conditions to the minimum trial count. Exclusions preserve IDs.

## Implemented cleaned NPZ input (no pickled objects)

| Key | Shape/type | Meaning |
|---|---|---|
| lfp_uv | channel × event × sample, float | Cleaned time-domain LFP |
| time_s | sample, float | Relative to the declared alignment event, at 1 kHz |
| electrode_id | channel | Unique physical electrode IDs |
| region | channel, Unicode | GC or PC |
| stimulus | event, Unicode | BA, IA, S, oBA, oIA, BK |
| good_channel | channel, bool | Explicit upstream QC acceptance |
| good_trial | channel × event, bool | Channel-specific trial acceptance |
| session_id / animal_id | scalar Unicode | Stable biological identifiers |
| provenance_json | scalar Unicode | JSON with source, reference, filtering, qc_method, event_alignment |

The first implementation is TTL-aligned; retain at least -2 to +2 seconds
(prefer the existing -2 to +5 s snippets). A future adapter must also save the
canonical event table alongside this archive. Different trial counts are allowed.
No voltage interpolation across rejected artifacts is performed by the power stage.
Rejected samples may be nonfinite; accepted trials may not.

## Normalization and outputs

Per-trial/channel/frequency power = abs(STFT)**2. Baseline STFT center times are
[-1.5,-0.5) s; this is an explicit provisional default, not a user-confirmed
biological choice. All 500 ms windows used for baseline lie before TTL onset.
Use baseline mean and population SD, then average trials per channel, then channels.
500 samples at 1 kHz gives 2 Hz frequency-bin spacing; 1 ms hops are not independent
observations. Plot -1000 to +1750 ms using turbo with a common fixed color scale.
Pre/post curves are normalized power, not absolute PSD, and are labeled accordingly.

`session_power.npz`: region × stimulus × frequency × time z_power; matching
n_channels at each bin, named coordinates, animal/session IDs and settings.
`provenance.json`: source SHA256, settings, reference/filter declaration and
accepted/valid trial counts per electrode/frequency. Archive these with plots.
