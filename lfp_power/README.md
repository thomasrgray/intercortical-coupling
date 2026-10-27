# LFP plotting and power spectra: Julia workflow

This is the default LFP power workflow, ported from
[thomasrgray/julia_blech_clust_testing](https://github.com/thomasrgray/julia_blech_clust_testing/tree/1a4a4eb6346b00226e907228f4a2c3f658223907).
It is bundled in this repository; no second checkout or PYTHONPATH setup is needed.

```bash
cd "$HOME/intercortical-coupling"
python -m pip install -e '.[lfp,phase,spikes]'
DATA_LIST="$PWD/configs/fin_inter_list.txt"
RUN="$PWD/results/julia_lfp_$(date +%Y%m%d_%H%M%S)"
coupling-lfp-pipeline --data-list "$DATA_LIST" --output-dir "$RUN" --jobs 1
```

The runner executes all five stages, streams progress, saves logs and `run.json`,
and stops after a failed stage. Use `--jobs 2` (or more) if memory permits;
STFT intermediates can be large. Individual commands use `COUPLING_LFP_JOBS`
(default 1). The source STFT progress bars remain enabled.

| Stage | Individual command | Outputs |
| --- | --- | --- |
| Channel QC and voltage plots | `coupling-lfp-validate --data-list "$DATA_LIST"` | Each session's `LFP_channel_check/` |
| Artifact removal | `coupling-lfp-clean --data-list "$DATA_LIST"` | Each session's `cleaned_LFP_artifacts/` |
| Spectrograms and pre/post spectra | `coupling-power --data-list "$DATA_LIST" --output-dir "$RUN/session_reports"` | Each session's `LFP_analyses/`; validation JSON in the run |
| Grand averages | `coupling-lfp-grand-average --data-list "$DATA_LIST" --output-dir "$RUN/grand_average"` | Regional PNG/SVG spectra and spectrograms |
| Animal-level comparisons | `coupling-compare-power --data-list "$DATA_LIST" --output-dir "$RUN/statistics"` | GC/PC and within-route comparisons, HDF5 values, statistical audit |

Numbered Python wrappers in this folder expose the same stages. `02_power.py`
also runs the new session plotting stage. The older `01_prepare.py` remains a
raw-preparation compatibility command; it is not part of this five-stage pipeline.

For one session, use `--dataset /absolute/path/to/session` and
`--stages validate clean plot average`. Statistics require at least two complete
animals for each comparison; insufficient groups are reported as untested and
return nonzero. A completed descriptive plot is not a significant statistical result.
Review channel QC before interpretation. To pause after QC, run
`--stages validate`, then use a new run directory with `--stages clean plot average compare`.

## What is retained from the source

- HDF5 `/Parsed_LFP` arrays; `.info` region/taste metadata; `.params` epoch timing.
- Source channel QC, approved-channel voltage plots, MAD trial removal and exclusions.
- Detrended Hann STFT: 500 ms windows, 1 ms hops, frequencies below 40 Hz.
- Per-trial/channel/frequency baseline z-scoring (−2 to −0.25 s), median
  spectrograms across trials/channels, and mean pre/post spectra (post 0–1.75 s;
  spectra displayed through 35 Hz). These are normalized power responses,
  not absolute PSD in µV²/Hz.
- Session-weighted descriptive grand averages; animal-weighted inferential comparisons.
- Confirmed legacy labels: `waterIA`/`waterBA` → retro `IA`/`BA`; accompanying
  plain `IA`/`BA` → `oIA`/`oBA`. Canonical labels are unchanged; collisions fail.
- Source plot layouts, colors, PNG/SVG exports and statistical correction.

See [statistical methods and limitations](../docs/LFP_spectrogram_statistics.md)
and [port provenance](../src/coupling/julia_lfp/PROVENANCE.md).
`--exclude-taste SESSION_BASENAME:TASTE` and `--mad-threshold` are supported by
the runner; `--animal-map`, `--permutations`, and `--seed` configure comparisons.

## Archive previous results and start again

```bash
cd "$HOME/intercortical-coupling"
# First inspect the proposed moves (no changes):
python scripts/archive_results.py --data-list configs/fin_inter_list.txt
# Move old outputs into a timestamped archive:
python scripts/archive_results.py --data-list configs/fin_inter_list.txt --apply
# Recompute QC, clean trials, spectra, plots, averages, and comparisons:
coupling-lfp-pipeline \
  --data-list configs/fin_inter_list.txt \
  --output-dir "results/julia_lfp_$(date +%Y%m%d_%H%M%S)" \
  --jobs 1
```

Archives default to the sibling `../intercortical-coupling_archives/TIMESTAMP/`.
The script moves all repository `results/`, generated figure artifacts under
`figures/`, and the three session output directories listed above. This includes
previous respiration/spike/phase outputs inside `results/`, as requested for a
clean results tree. It preserves source code and READMEs. Add repeatable
`--output-dir /other/generated/results` for outputs saved elsewhere; these cannot
be discovered reliably from the session list. `--archive-root` changes the destination.
Raw recordings, HDF5 files, `.info`, `.params`, and trial/DIN CSV metadata are
not moved. Paths containing raw files, source code, or symlinks are refused.

Each archive has a `manifest.json` with original and archived paths and per-move
status. To restore, first move any newer output aside, then move each `archived`
path back to its `source`; do not overwrite newer results. A failed/interrupted
archive retains completed moves and its manifest for recovery. Moves across
filesystems require enough destination space for the copied outputs.

This resets **derived analysis outputs**, not the original LFP extraction:
existing `/Parsed_LFP` and `/Parsed_LFP_channels` are reused. Missing arrays or
channel metadata trigger the upstream extractor, which writes derived nodes to
HDF5 and requires raw amplifier files and `trial_info_frame.csv`. The archive
script never deletes HDF5 nodes. Verify epoch timing against the acquisition
settings; the source workflow assumes 1 kHz LFPs.

## Earlier workflows

`coupling-power` and `coupling-compare-power` now accept session arguments above,
not the previous positional HDF5 arguments. Those older interfaces remain available
as `coupling-power-legacy` and `coupling-compare-power-legacy`. Existing Python
modules are retained for respiration, coherence, Figure 2 and legacy Figure 6
compatibility; they are not used by the new power pipeline. The source migration
does not replace phase-coherence or respiration analyses. Old `--plot-power`,
`coupling-lfp-figure`, `coupling-export-figure6`, and power summaries still use
the legacy HDF5 format; do not feed the new STFT pickle to them.

The [HDF5 workflow guide](../docs/hdf5_workflow.md) documents the repository-wide
archive migration and the optional `phase` stage using Julia-cleaned trial IDs
and channel QC. Single-array NPY caches remain supported.
