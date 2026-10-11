# Delivery and alignment

## Panels

A: delivery schematic and TTL/physical onset calibration. B: nasal-catheter pressure and detected inspiration around ortho trials. C: retro delivery and available swallow proxy. D: alignment uncertainty and event-count QC.

## Required inputs

Nasal catheter pressure, TTLs, valve/PID calibration, event metadata; direct swallow measurement if available.

## Interpretation and validation

HMM state transitions are neural timing proxies, not direct evidence that swallowing occurred. Keep TTL, inhalation and HMM alignment as distinct analyses.

## Output contract

Save panel PNG/PDF, plotted values, counts/exclusions, configuration and source revision under `results/<run>/fig01_delivery_alignment/`. Panel letters are provisional. No synthetic output may be labeled as study data.
