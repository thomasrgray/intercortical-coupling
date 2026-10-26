# Archived Figure 2 display experiments

The active Figure 2 implementation is restored exactly from commit
`e27e8c2159e52d852ca9d61589a8cdc63ec87bcd`: blue/orange rasters first, alpha-only
LFP traces, filtered pressure, and firing rates; includes matched-odor comparisons.
The session-region mapping fix is retained in that version.

Archived source snapshots (not installed entry points):

- `20261026_overlays/figure2.py`: muted palette, background traces and envelopes;
  source commit `31e6972822955e1b663c3028a4bc83e17b9599e6`.
- `20261026_context_strips/figure2.py`: separate background strips and larger export;
  source commit `f98d21604aa6470c10ee37a6e5c7815c0d44a3a3`.

All earlier versions, tests and documentation also remain in Git history. These
snapshots require their corresponding package revision to reproduce exactly;
they are archival sources, not standalone scripts.

Actual recording images are generated on the acquisition/analysis workstation
and are not included in GitHub. After regenerating the preferred plots into a new
`figure2_preferred_*` directory, run:

```bash
python scripts/archive_figure2_outputs.py "$RUN"
```

This moves the original `fig02_recording_validation` and the timestamped
`figure2_comparison_*`, `figure2_overlays_*`, and `figure2_readable_*` directories
into `$RUN/archive/figure2_<timestamp>/`, including images, PDFs, NPZ, CSV and
provenance. It writes a manifest listing original/archive locations and files.
No files are deleted. Prepared LFP, respiration caches, preferred outputs and
unrelated analyses are left in place. Re-running with nothing eligible is a no-op.
