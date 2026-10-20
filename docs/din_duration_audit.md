# Compare DIN pulse durations across sessions

Run from the repository root after `python -m pip install -e .`:

```bash
coupling-audit-dins configs/fin_inter_list.txt "results/din_duration_$(date +%Y%m%d_%H%M%S)"
```

Equivalent: `python scripts/audit_din_durations.py SESSION_LIST OUTPUT`.

The list uses one session directory per line; relative paths resolve against the list.
Every board-DIN-*.dat file is inspected as little-endian unpacked uint16 0/1.
Each session's sample rate comes from info.rhd. No HDF5 plugins, LFP extraction,
electrode mapping, or six-equal-count requirement is needed.

## Outputs

- `din_duration_summary.csv`: stimulus labels from .info; complete/partial counts;
  min, 5th percentile, median, 95th percentile and max complete-pulse widths.
- `pulse_widths.csv`: every observed high segment, onset, exclusive offset, observed
  width, and censoring flags. Inspect for mixed pulse widths and outliers.
- `duration_comparison_01.png/.pdf`: annotated session-by-DIN median duration map.
  Additional pages cover lists longer than 15 sessions. Blank cells lack complete
  pulse medians. The logarithmic color scale is shared across pages. Row indices
  distinguish duplicate session names; CSVs retain absolute session paths.
- `audit.json`: session/metadata paths and duration-check settings.

## Interpretation

`--split-ms 100` is a configurable long/short heuristic, intended to separate
roughly 15 ms retro from roughly 2000 ms ortho pulses in this task. Labels oBA/oIA/BK
expect long pulses, BA/IA/S expect short pulses; other labels are unassigned.
`REVIEW` flags disagreement with this heuristic, not an automatic relabeling decision.
A median can conceal a mixed train: inspect min/max and the per-pulse CSV too.

Duration cannot identify BA versus IA, establish which valve was connected, measure
odor arrival, or establish physical delivery duration. Verify against Bpod wiring
and task logs. Extra cue lines remain visible and are never assigned identities.

Boundary-high segments are censored and EXCLUDED from complete-width statistics.
The blocked TG23 Feb4 recording can therefore be inspected without weakening the
production LFP pipeline's strict pulse checks. Raw files and metadata are read-only.
Malformed metadata leaves labels unavailable but does not prevent width measurement.
Nonbinary files, missing headers, duplicate physical DINs and missing declared DINs
produce error rows while other files/sessions continue. Unequal file lengths warn.
Exit status 1 means error rows exist; usable plots and tables are still saved.
Warnings and REVIEW flags alone do not cause a nonzero exit. Output folders must be new.

## Validation

`python -m unittest discover -s tests -p 'test_din_audit.py' -v`

Tests cover chunk boundaries, censored pulses, idle/nonbinary/truncated files,
route heuristics, malformed metadata, duplicate/missing DINs, read-only behavior,
and end-to-end CSV/PNG/PDF export.

## Single all-session summary

Every audit now also writes `all_sessions_summary.png`, `.pdf`, and `.svg`.
Panel A shows all sessions × DINs, with median pulse duration and stimulus labels.
Panel B shows complete-pulse counts and explicitly annotated boundary-partial counts.
An asterisk flags duration/label disagreement. Sessions with failed reads remain visible;
missing measurements are gray rather than represented as zero. The figure grows with
session count and is not paginated. The earlier paginated plots remain available for
very large lists. No averaging across different stimulus identities or sessions is done.
