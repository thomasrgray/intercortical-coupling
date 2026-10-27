# All-session LFP power and phase-coupling summary

This reads completed results. It does not recalculate LFPs, STFTs or bootstraps.

From the repository root, after `python -m pip install -e .`:

```bash
coupling-summarize-lfp \
  --phase-run "$LFP_RUN/phase" \
  --power-root "$LFP_RUN" \
  --output "figures/fig06_lfp_coupling/summary_$(date +%Y%m%d_%H%M%S)"
```

`LFP_RUN` is the existing raw-to-phase run directory from the batch instructions.
If opening a new terminal, set it to that actual directory first. The phase run
must contain status.json and its completed session output folders. Power discovery
recursively finds session_power.h5 under one run root. Do not use the whole results/
directory: duplicate runs of the same animal/session are rejected.

You can supply only --phase-run or only --power-root. The script entry point is
`python lfp_phase_coherence/04_summary.py` with the same arguments.

## Figures

- `all_sessions_power`: one row per power session; GC ortho/retro and PC ortho/retro
  spectrograms on the same baseline-z scale.
- `all_sessions_phase`: one row per completed matched phase session, one column
  per frequency band; blue ortho and orange retro traces.
- `animal_mean_power`: GC/PC ortho, retro and ortho-minus-retro spectrograms,
  averaging sessions within each animal first and then weighting animals equally.
- `animal_mean_phase`: thin per-animal traces, thick equal-animal mean, and the
  signed ortho-minus-retro contrast. No significance stars or confidence intervals.

All figures are PNG, PDF and editable-text SVG, saved in a new output directory.
Source arrays, per-bin power animal counts, session/animal IDs and source hashes
are exported alongside the plots. Power and phase cohorts are listed independently
in summary_manifest.json: a session with valid power may have failed phase QC.
Failed or missing matched phase conditions are recorded as exclusions, not zeros.

## Comparisons and interpretation

Default `--comparison route` uses equal BA/IA means. For power, both odors must
contribute at each bin; for phase it reads the existing within-odor coherence
route summaries rather than pooling phase vectors. `--comparison BA` or
`--comparison IA` generates matched-odor summaries separately in a new directory.
Sucrose and blank controls are not included in these matched-odor route summaries.

Repeated sessions do not increase an animal's weight in the mean. Thin traces
show animal means, not bootstrap draws. These are descriptive summaries, not
animal-level hypothesis tests. Channels and trials are not independent animals.
Review reference/filtering/QC provenance before treating sessions as comparable.

Time/frequency grids, alignment labels, band definitions, phase profile/settings
and power normalization/settings must agree within their respective summaries.
No interpolation or silent mixing of settings is performed. Power NaNs remain
missing and animal contributions can vary by bin (counts exported). Missing animal
IDs, duplicate power sessions, invalid arrays and incompatible grids fail explicitly.
All figures are staged so a failed export leaves no completed output directory.

Tests: `python -m unittest discover -s tests -p 'test_lfp_summary.py' -v`.
