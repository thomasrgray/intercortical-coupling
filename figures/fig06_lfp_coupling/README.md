# Figure 6 — Mahmood 2023 phase coherence

The preferred figure workflow is now [`coupling-phase-figure`](../../docs/mahmood_phase_coherence.md).
It computes GC–PC phase coherence with representative electrodes, replacement
bootstrap/shuffle controls, baseline deviations and all-session figures.
Each stimulus has its own band curves and GC/PC power context. Matched-odor
comparisons (oBA vs BA; oIA vs IA) accompany equal-odor ortho/retro summaries.
See the guide for event-specific time labels and route weighting.

Start with [`fin_inter_list.template.txt`](../../configs/fin_inter_list.template.txt)
and [`phase_coherence.example.json`](../../configs/phase_coherence.example.json).
The methods, file layout, command examples and flowchart are in the workflow guide.

The former A–D magnitude-squared-coherence/Granger contract remains available in
[`legacy_A_D.md`](legacy_A_D.md) for reproducing old outputs.
