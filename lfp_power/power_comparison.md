# Two-panel sample power figure

Use saved `power/session_power.h5`; no raw processing is repeated.

```bash
python -m pip install -e .
coupling-compare-power-legacy "$RUN/power/session_power.h5" results/sample_gc_ia --left GC:IA --right GC:oIA
```

For an additional descriptive difference panel:

```bash
coupling-compare-power-legacy "$RUN/power/session_power.h5" results/sample_gc_ia_difference --left GC:IA --right GC:oIA --difference --difference-clim 1
```

Compare regions with `--left GC:S --right PC:S`, or the other matched odor pair
with `--left PC:BA --right PC:oBA`. Labels are case-sensitive. PC is displayed as
pPC. Each output directory must be new. If `$RUN` is no longer set, replace the
input argument with the full path to your existing session_power.h5.

Outputs: 600 dpi PNG, PDF and SVG (editable text, rasterized heatmaps), numerical
left/right/difference arrays with contributor counts, provenance JSON and caption.
The two source panels share turbo and symmetric ±3 z limits; `--clim` changes both
together. The optional difference uses a diverging map and independent symmetric
limits. White areas are missing values, not zero responses. Limits may clip values;
consult saved arrays and document any display-limit changes.

Difference is LEFT minus RIGHT of the already normalized session means. It is
not a paired-trial test, percent power change, coherence or a significant contrast.
Even cross-region z differences compare standardized responses, not absolute
voltage power. One animal/session remains one biological replicate; electrodes
are not independent animals. Counts can vary by frequency/time after exclusions.
The 500 ms window smooths timing; a 1 ms grid does not give 1 ms temporal resolution.
TTL timing and route-specific delivery/respiration remain potential confounds.

Confirm anatomy and review upstream QC/CAR before manuscript use. This command
only renders existing outputs; it does not change normalization or QC thresholds.
