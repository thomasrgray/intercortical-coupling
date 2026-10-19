# Adapted spike noise and transition correlations

**The default is now `--profile repository`.** For the closer reproduction of the
reference setup and aggregation, follow [the repository-profile guide](mahmood_repository_parity.md).
The analysis behavior described below is preserved with `--profile adapted`; input
file contracts apply to both profiles.

This adds runnable Figure 8 noise correlations and Figure 7 transition correlations
alongside the existing phase-coherence workflow. It is an independent implementation
adapted to GC–PC and the six stimulus conditions, not a claim of exact paper replication.
Source recordings and fitted models are opened read-only. Outputs require a new directory.

## Reference code and explicit differences

Reference repository revision: `9c54acb57e227dc9d13597750165b3ebc7be68e9`.
The implementation was checked against executable scripts rather than the broader
feature claims in the reference READMEs:

- [Spike setup](https://github.com/abuzarmahmood/Mahmood_2023_BLA_GC/blob/9c54acb57e227dc9d13597750165b3ebc7be68e9/spike_noise_corrs/inter_region_noise_corrs_setup.py):
  sum 2000–4000 ms samples, first-difference along chronological trials separately
  by taste, z-score, Spearman correlation, 10,000 trial permutations; inter- and
  intra-region pairs. We omit the z-score because it cannot change Spearman ranks
  for nonconstant vectors. Constant vectors are explicitly excluded.
- [Transition engine](https://github.com/abuzarmahmood/Mahmood_2023_BLA_GC/blob/9c54acb57e227dc9d13597750165b3ebc7be68e9/transition_corrs/all_tastes/inter_region/run_inter_corrs.py):
  integer-bin posterior mode, Spearman correlation, shuffles and absolute timing
  error. Its `calc_mse` actually computes mean **absolute** error; we name it MAE.
- [All-to-all transition helper](https://github.com/abuzarmahmood/Mahmood_2023_BLA_GC/blob/9c54acb57e227dc9d13597750165b3ebc7be68e9/transition_corrs/all_tastes/corr_test.py):
  compares every cross-region transition combination. We use this instead of
  assuming that GC transition 1 and PC transition 1 represent the same state.

The source transition runner pools trials. Here, tests and permutations are separate
for every **stimulus × stratum** (laser duration/onset, or user-defined nuisance block).
This avoids interpreting different condition means as trial-to-trial coordination.
The spike default window is 0–2 s after the declared onset, plus a -2–0 s baseline.
Optional rolling windows extend the core calculation; they are not an exact port
of the reference rolling script.

We report two-sided finite-sample permutation p-values `(extreme + 1)/(B + 1)` for
rho, a lower-tail p-value for MAE, reference-style percentiles, and Bonferroni p-values
across all finite tests **within each session**, separately for rho and MAE.
There is no correction across sessions; cross-session output is descriptive.

**Exchangeability caveat:** neighboring first differences share a trial, so their
permutation null is a reference-code control, not a guarantee of calibrated inference
under serial dependence. Check trial drift/serial dependence and compare
`"detrend": "none"`; do not interpret small permutation p-values alone as definitive
coupling. Stimulus/laser stratification does not remove shared breathing or movement.

## Run spike noise correlations first

```bash
python -m pip install -e '.[spikes]'
cp configs/spike_noise.example.json configs/spike_noise.local.json
```

Edit `preprocessing` to document actual verification; do not merely delete the warning.
Confirm the spike-array sample rate and onset from the spike extraction configuration;
these are **not** the raw Intan sample rate or automatically the LFP settings.
The example DIN map follows the confirmed TG22 mapping, but must be verified per session.
Use `sessions` overrides, keyed by absolute session directory, when settings differ.
`animal_id` can be provided per session and is retained in the combined CSV.

```bash
RUN_STAMP=$(date +%Y%m%d_%H%M%S)
coupling-spike-noise configs/fin_inter_list.txt configs/spike_noise.local.json \
  "results/spike_check_${RUN_STAMP}" --profile adapted --check &&
coupling-spike-noise configs/fin_inter_list.txt configs/spike_noise.local.json \
  "results/spike_noise_${RUN_STAMP}" --profile adapted --shuffles 10000 --seed 0
```

`--check` checks every listed session, writes `preflight.json`, and exits unsuccessfully
if any input is unusable. It does not run correlations. Pair progress appears during
analysis; `status.json` tracks completed sessions and records runtime failures.
Optional `"rolling_s": [-1, 2, 0.5, 0.25]` adds 500 ms windows stepped by 250 ms.
Start without rolling windows to keep the number of tests and output size manageable.

### Blech-clust HDF5 inputs

Expected nodes:

- `/spike_trains/dig_in_N/spike_array`: **trial × unit × sample** integer spike counts.
- `/unit_descriptor`: structured rows with `electrode_number` and, by default,
  `single_unit`. Array unit columns must correspond to these rows.
- Session `.info` with `electrode_layout.gc` and `.pc`; nested electrode lists work.
- Per-DIN `laser_durations` and optional `laser_onset_lag` arrays where relevant.

No anatomy is inferred from unit count or row order. Units outside the two configured
regions are omitted; non-single units are excluded by default. `eligible_unit_ids`
can further restrict reviewed units using original descriptor-row IDs. Explicitly
setting `single_units_only=false` includes multiunits and must be justified in QC.
An optogenetic session declared in `.info` without per-trial laser durations fails.
Arrays in each DIN must retain chronological order; generated IDs are
`DIN0:trial0`, etc. Different DINs are never differenced together.

HDF5 loading starts with all parsed trials accepted. For artifact masks or additional
nuisance strata, export the NPZ contract below. Nonfinite accepted data fails; silent
neurons produce a `constant` result rather than a fabricated zero correlation.

### Reviewed spike NPZ alternative

Set `"input": "reviewed_spikes.npz"`. Required arrays (no Python objects/pickles):

| Key | Shape / interpretation |
|---|---|
| `spikes` | trial × unit × sample; nonnegative integer counts |
| `time_s` | sample; regular grid relative to declared event |
| `unit_id`, `region`, `good_unit` | unit; stable IDs, `GC`/`PC`, boolean inclusion |
| `trial_id` | trial; unique strings |
| `stimulus`, `stratum` | trial; nonempty strings; e.g. `BA`, `laser:off` |
| `trial_order` | trial; unique finite chronology index |
| `good_trial` | trial; boolean joint acceptance mask |

First differences use adjacent trials within each condition's ordered sequence;
if either endpoint is rejected, the difference is rejected. The code never bridges
a masked trial. A default minimum of eight **differences** therefore requires at
least nine accepted consecutive within-condition trials. Raw mode needs eight trials.
Use strata for acquisition blocks when differences should not cross block boundaries.

## Transition correlations: fit inputs are a separate prerequisite

These commands do **not** fit changepoints, choose state number, or establish model
convergence. Independently fit the GC and PC populations first. The original analysis
also loads existing fits. Do not align both regions to a shared inferred transition,
and do not use correlation strength to choose the number of states.

Export one NPZ per region with `tau` shaped **posterior draw × trial × transition**
in nonnegative **bin coordinates relative to fit start**. Export `trial_id`, `stimulus`,
`stratum`, `trial_order`, and `good_trial` as above. Both regions must contain exactly
the same trial IDs; differing row order is resolved by IDs, never by positional pairing.
Rejected trials remain present with masks. Per-trial stimulus/stratum/order must agree.

`metadata_json` is a scalar JSON string. Example for GC (PC must use its own model ID):

```json
{
  "region": "GC",
  "model_id": "TG22_Feb2_GC_4state_fit1",
  "independently_fitted": true,
  "fit_qc_passed": true,
  "fit_qc_notes": "Describe actual convergence and posterior predictive checks here",
  "bin_width_ms": 50,
  "fit_start_ms": 0,
  "fit_end_ms": 2000,
  "alignment_event": "DIN"
}
```

Metadata assertions record your review; the loader cannot independently verify
convergence. Fit intervals/alignment must match; bin widths may differ because times
are converted to milliseconds. Input tau must be finite and ordered per draw/trial.
Different transition counts across regions are supported.

Export from already loaded fitting outputs without loading an untrusted pickle:

```python
np.savez_compressed(
    "transition_posteriors/GC.npz", tau=tau_samples,
    trial_id=trial_ids, stimulus=stimulus_labels, stratum=stratum_labels,
    trial_order=chronological_indices, good_trial=joint_qc_mask,
    metadata_json=json.dumps(metadata),
)
```

Use the fit's actual trial metadata, not a guessed flattened taste order.
If a fit stores `draw × stimulus × trial × transition`, flatten the stimulus/trial
axes and the associated labels/IDs in the same order before exporting.

```bash
cp configs/transition_corr.example.json configs/transition_corr.local.json
# Edit preprocessing and confirm posterior file paths and metadata first.
RUN_STAMP=$(date +%Y%m%d_%H%M%S)
coupling-transition-corr configs/fin_inter_list.txt configs/transition_corr.local.json \
  "results/transition_check_${RUN_STAMP}" --profile adapted --check &&
coupling-transition-corr configs/fin_inter_list.txt configs/transition_corr.local.json \
  "results/transition_corr_${RUN_STAMP}" --profile adapted --shuffles 10000 --seed 0
```

Posterior mode uses integer truncation of nonnegative bin samples, matching the
reference; ties choose the smaller bin. Timing variance uses the original samples.
Independent posterior draws from each regional fit give a rho sensitivity interval;
this is **not** an animal-level confidence interval or a replacement for model checks.

## Outputs and interpretation

Each session writes under `fig08_spike_coupling/` or `fig07_transition_coupling/`:

- `correlations.csv`: all tested pairs, counts, exclusions/status, rho, null summaries,
  permutation p-values and correction. Spike rows distinguish baseline/evoked/rolling
  and GC–PC/GC–GC/PC–PC; transition rows include MAE and uncertainty.
- `nulls.npz`: shuffle samples keyed by `null_key`; transition MAE and sampled posterior
  rho have `_mae` and `_posterior_rho` suffixes when enough trials exist.
- `summary.png` and `.pdf`: session plots; not a complete manuscript panel set.
- `retention.json`: trial labels and joint masks, plus unit inclusion or fit metadata.
- `provenance.json`: resolved configuration, seeds, input hashes, code/reference revisions.
- `route_contrasts.csv`, when configured matched pairs are estimable: ortho minus retro
  rho for the same pair/stratum/window. This is a descriptive contrast, not a route test.
- Transitions additionally save `posterior_summary.npz` with matched trial IDs, masks,
  modes and variances for inspecting uncertain trials.

The run root contains `all_sessions.csv` with session and animal labels. Units/pairs
are not independent animals, and units are not tracked across sessions. Use a separate
animal-aware analysis for population claims; unequal route trial counts are reported,
not automatically corrected or subsampled. Bonferroni can make 10,000 shuffles too
coarse for very large families; increase shuffles rather than interpreting the minimum
p-value as zero.

Not implemented here: Bayesian changepoint fitting, within-region split-population
fits/correlations, uncertainty-ranked moving-window analysis, cross-correlograms,
Poisson GLMs, respiration nuisance regression, or animal-level route inference.
These remain distinct follow-up analyses; the existing figure contracts are broader
than this contribution. Synthetic tests verify computation, not biological validity.
