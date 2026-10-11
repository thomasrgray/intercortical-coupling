# Simultaneous recording validation

## Panels

A: GC/PC placement. B: example spike waveforms and unit QC. C: simultaneous spikes, pressure and LFP traces. D: animals/sessions/units/electrodes retained.

## Required inputs

Electrode layout, histology, unit metadata, voltage, respiration and masks.

## Interpretation and validation

An electrode can have valid LFP without a sorted unit. Report animal/session denominators separately.

## Output contract

Save panel PNG/PDF, plotted values, counts/exclusions, configuration and source revision under `results/<run>/fig02_recording_validation/`. Panel letters are provisional. No synthetic output may be labeled as study data.
