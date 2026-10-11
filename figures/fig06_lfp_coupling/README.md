# LFP coupling — first implementation target

## Panels

A: channel/trial QC and retention, reference/filter diagnostics. B: GC/PC × six stimuli power context. C: GC-PC magnitude-squared coherence and route contrasts for matched odors. D: directional spectral Granger, diagnostics and matched-odor route contrasts.

## Required inputs

Explicitly mapped cleaned voltages, QC masks, original event IDs, reference provenance, session/animal IDs.

## Interpretation and validation

Power is contextual, not a coupling measure. Shared referencing, respiration, volume conduction and trial-locked responses can alter apparent coupling. Treat directionality as predictive dependence, not causal intervention.

## Output contract

Save panel PNG/PDF, plotted values, counts/exclusions, configuration and source revision under `results/<run>/fig06_lfp_coupling/`. Panel letters are provisional. No synthetic output may be labeled as study data.
