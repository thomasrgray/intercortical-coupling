# Differential inter-cortical coupling in orthonasal and retronasal processing

Figure-first manuscript repository. GC = gustatory cortex; PC = posterior piriform cortex.
Each figure has a panel contract; reusable processing lives in `src/coupling`.
No legacy scripts or recordings are changed. No study outcomes are assumed.

## Start here

1. Read `docs/figure_plan.md` and the panel contracts in `figures/`.
2. Read `docs/data_contract.md` before preparing any session.
3. Follow `docs/lfp_roadmap.md`: Fig 6A quality control, Fig 6B power context,
   then Fig 6C coherence and Fig 6D Granger.

## What runs now

The cleaned-array power stage computes 1000 Hz Hann STFTs (500 samples,
499 overlap, 0–70 Hz); baseline z-scores are computed per channel/trial/frequency
BEFORE trial and channel means. Output: 12 GC/PC × stimulus spectrograms,
12 pre/post frequency plots, their CSVs, session arrays, counts and provenance.
An all-invalid frequency is retained as NaN, not converted into a zero response.

```bash
python -m pip install -e .
python -m unittest discover -s tests -v
python scripts/demo_power.py
coupling-power results/demo_input.npz results/demo_power
```

The demo is SYNTHETIC; it is a rendering/numerical check, not a biological result.
An existing output directory is rejected; choose a new directory for each run.

## Raw recording processing — available now

```bash
coupling-scan-dins /path/to/recording /path/to/results/din_scan
coupling-prepare /path/to/recording configs/tg22_feb2.example.json /path/to/results/session_run --plot-power
```

Edit the example's explicit amplifier/region/reference-group mapping and stimulus
labels for your recording. All configured amplifier files must exist. The output
must be a NEW directory OUTSIDE the recording directory.
Read `docs/raw_processing.md` for defaults, QC and limitations.

## Still planned

Respiration alignment, coherence, Granger, spike analyses, HMMs, GLMs and grand
averages. Raw preparation and automatic QC are now implemented, but have not yet
been validated end-to-end on a complete biological session. No study outcomes
are assumed from successful synthetic tests.

## Directory organization

| Directory | Purpose |
|---|---|
| `figures/fig01_delivery_alignment/` | Delivery and alignment validation |
| `figures/fig02_recording_validation/` | Simultaneous GC/PC, breathing, representative LFPs |
| `figures/fig03_pc_reliability/` | PC response fractions, latency, locking, split-half |
| `figures/fig04_pc_route_dynamics/` | Matched odors, temporal shifts, population dynamics |
| `figures/fig05_gc_dynamics/` | GC responses and three-phase hypothesis |
| `figures/fig06_lfp_coupling/` | QC, power, coherence, directional analyses — FIRST |
| `figures/fig07_transition_coupling/` | Coupling of inferred population transitions |
| `figures/fig08_spike_coupling/` | Noise correlations, cross-correlations, Poisson GLMs |
| `src/coupling/` | Shared tested computations |
| `configs/` | Session/analysis configuration |
| `scripts/` | Small entry points and synthetic demo |
| `tests/` | Numerical and input-contract tests |
| `docs/` | Methods decisions, data contracts, build order |

`data/` and `results/` are ignored by Git. Keep actual recordings outside this repo.
See `docs/github_setup.md` for creating the private remote repository.
