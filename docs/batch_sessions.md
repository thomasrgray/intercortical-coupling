# Multiple-session LFP processing

`coupling-batch` discovers sessions, validates each folder against its own metadata,
and executes preparation plus Figure 6 A–D separately for each valid session.
It never pools electrodes/trials from different sessions into biological replicates.
Raw folders are read-only. This is serial batch processing to keep memory bounded;
`.params` CPU settings are not used.

## Tom's TG22 recordings

After checking out this change and installing `python -m pip install -e .`:

```bash
# Find folders containing .info files below the animal folder.
coupling-batch init /media/thomas/Data/TG22 configs/my_tg22_sessions.json

# Inspect/edit the manifest, then scan without filtering amplifier data.
coupling-batch check configs/my_tg22_sessions.json results/tg22_check_01

# Review results/tg22_check_01/preflight.json and summary.csv, then run.
coupling-batch run configs/my_tg22_sessions.json results/tg22_batch_01
```

For several animals, use `/media/thomas/Data` as the discovery root and edit the
manifest to retain the intended session folders. Discovery uses `.info` files,
not a guessed date pattern; `.bak` files are ignored. TG-number animal IDs are
inferred only when the session name/parent agree unambiguously. Review all IDs;
other naming schemes require filling `animal_id` in the manifest.

The [single-session example manifest](../configs/tg22_batch.example.json) points
to the supplied February 2 path and works as a starter without discovery.
Add another session object for each folder. Relative paths are resolved against
the manifest's directory; generated paths are absolute. IDs and paths must be
unique. Session IDs become output folder names and cannot contain path separators.

Both `check` and `run` require a **new output directory**. The `run` command repeats
the checks rather than trusting an older report. Every session is checked before
processing begins. Blocked sessions are skipped, valid sessions are processed,
and a runtime failure is recorded while subsequent sessions continue. Exit code
1 means at least one session was blocked/failed; inspect `summary.csv` for details.
There is no implicit resume or overwrite: to retry, make a manifest containing
only failed sessions and choose a new output directory. Completed outputs remain
available even when the batch exits with code 1.

## Required files and how they are interpreted

| Source | Use |
|---|---|
| Exactly one `*.info` JSON | Six stimulus names/DINs/trial counts; explicit anatomical electrode layout; acquisition settings; laser declarations |
| Exactly one `*electrode_layout*.csv` | Amplifier filenames, physical electrode IDs, reference groups; anatomical mapping cross-check |
| `info.rhd` | Authoritative Intan sample rate; must agree with declared rates |
| `time.dat` | Full int32 timestamp stream; must match sample count and increment by one |
| All mapped `amp-<port>-<number>.dat` | Int16 raw voltages; lengths must match DINs |
| All canonical `board-DIN-<number>.dat` | Uint16 unpacked 0/1 pulses; full measured inventory |
| Optional `*.params` JSON | Sampling-rate cross-check only; spike filters, waveform thresholds, clustering and auto-CAR settings are ignored |

The current batch contract supports the six labels **BA, IA, S, oBA, oIA, BK**.
It does not silently rename other stimulus labels or choose a mapping based on
trial count alone. Each session supplies its own mapping from `.info`.

Use original recording filenames. Uploaded copy suffixes such as `(1)` on DINs
are not accepted as canonical acquisition filenames. Multiple `.info`, CSV or
`.params` files are ambiguous: explicitly select the intended filename using
`info`, `layout`, or `params` in the session object. `header` and `timestamps` can
likewise select alternate names, directly inside that recording folder. Backup
metadata is never promoted automatically.

## Checks before processing

- Metadata JSON structure, safe/unique session IDs and required animal IDs.
- RHD validity; sample-rate agreement with `.info` and `.params`.
- Exactly matching `.info`/CSV electrode sets; unique IDs/filenames; explicit GC/PC
  assignment; no reference group spanning anatomical regions.
- All mapped amplifier files present and sample lengths consistent. Unmapped
  amplifier files are listed and ignored; channels are not automatically dropped
  just because they are missing or have zero predicted spike clusters.
- At least four mapped channels in each region, because existing QC requires
  three same-region peers per channel; at least three donors per reference group
  before mean-CAR. Upstream QC can
  still reject donors and make a group unusable during preparation.
- Scan every DIN for complete 0/1 pulses and consistent recording lengths. Select
  exactly one group of six equal positive trial counts, then require its identities
  and measured counts to agree with stimulus metadata. Extra DIN metadata counts
  are also checked. No trial truncation repairs a mismatch.
- Scan the full timestamp stream for gaps/resets and count boundary events that
  cannot provide the full epoch. Such events retain IDs but fail epoch QC.
- Block declared laser/opto sessions: the current figure does not stratify opto
  conditions. Extra DINs alone are not assumed to represent laser activation.
- Output location outside raw folders; rough working-array disk estimate and a
  20% disk margin before each session. This is not a guaranteed peak RAM/disk
  bound: filesystem compression and plot/CSV size vary.

Settings in `defaults` can be overridden for one session using its `options`
object: `reference_method`, `notch_60hz`, `max_channels`, `min_trials`, `orders`,
and `matched_pairs`. Stimulus/anatomical mappings remain metadata-derived and
cannot be replaced by these numerical options. The reference groups are taken
from CSV `CAR_group`; anatomy is taken from `.info` and checked against
`original_CAR_group` where available. A `.params` auto-CAR flag does not establish
that the original DAT files have already been referenced.

## Output layout

Each batch saves `manifest_resolved.json`, incremental `preflight.json`, and
`summary.csv`. The summary includes accepted channel-event counts and the number
of Granger models passing diagnostics; “completed” means the workflow completed,
not that every model passed or a biological effect exists.

Each session subfolder contains:

- `config.json`: resolved input mapping, options, metadata hashes and acquisition
  settings. This is usable with the existing `coupling-prepare` command.
- `prepared/`: cleaned LFP archive, upstream QC, event tables and preprocessing
  provenance.
- `fig06_lfp_coupling/`: A–D PNG/PDF panels, numerical estimates, model diagnostics
  and provenance.
- `run.log`, `status.json`, and on failure `failure.txt`.

Sessions remain separate, even within one animal. This change implements batch
processing, **not cross-session statistics or an animal-level group figure**.
See [Figure 6 methods](../figures/fig06_lfp_coupling/README.md) for estimator details.

## What the supplied February 2 examples establish

The uploaded `.info` and CSV agree on 64 electrodes, GC IDs 0–31 and PC IDs 32–63.
The mapping is DIN 0=BA, 1=S, 2=oIA, 3=IA, 4=oBA, 5=BK, with 30 trials each;
DINs 6 and 7 are declared as 90 trials each. The example `.params` contains
300–3000 Hz spike filters, which must not be used for this LFP analysis.

Read-only checks on the attached time, DIN-04 and amplifier-061 files found
102,814,080 samples per stream, uninterrupted timestamps starting at zero, and
30 DIN-04 pulses with a median width of 2 seconds. At the mutually declared
30 kHz rate this is 3427.136 seconds (57 min 7.136 s). The header rate cannot be
independently checked from these attachments because `info.rhd` was not included.
The single amplifier passed the global raw flat/outlier screen; this does not
establish channel retention after peer, epoch and reference QC.

No extra metadata is needed to implement the batch mechanics. Before scientific
interpretation, confirm that the per-session `.info` files are current, that
DINs 6/7 are not undeclared opto conditions, and that BA/oBA and IA/oIA identify
the intended matched odors. Full local folders need the RHD header, all mapped
amplifiers and DINs. You do not need to upload all recordings to run this locally.
The supplied partial files alone cannot support a complete session analysis.
