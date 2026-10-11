# PC route comparison

## Panels

A: BA vs oBA and IA vs oIA PSTHs. B: response-similarity distribution. C: held-out similarity after allowed retro time shift. D: route-specific population dynamics.

## Required inputs

Matched odor pairs, PC spike tensors, trial identities, alignment uncertainty.

## Interpretation and validation

S and BK are separate controls, not a matched odor pair. Define equivalence margins; nonsignificance is not sameness. Fit shifts on training trials, evaluate on held-out trials. Dynamics can use population rather than within-unit comparisons.

## Output contract

Save panel PNG/PDF, plotted values, counts/exclusions, configuration and source revision under `results/<run>/fig04_pc_route_dynamics/`. Panel letters are provisional. No synthetic output may be labeled as study data.
