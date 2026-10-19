# Run the reference repository profile

`coupling-spike-noise` and `coupling-transition-corr` now default to
`--profile repository`. This follows the executable core setup/aggregation scripts
at Mahmood_2023_BLA_GC revision `9c54acb57e227dc9d13597750165b3ebc7be68e9`.
It reproduces that repository revision, which includes later edits, rather than
asserting exact reproduction of every analysis in the 2023 publication.

The older implementation remains available with `--profile adapted`. Its baseline,
route-separated transition analyses and pairwise permutation Bonferroni outputs
are not the reference profile's primary statistics. Existing result directories
are not modified; run into a fresh output directory.

## Spike setup and aggregation

| Step | Reference repository profile |
|---|---|
| Input | Sorted spike arrays, not electrode-level detected timestamps |
| Unit inclusion | All mapped sorted units, matching the source's `get_region_units()` use; an explicit `eligible_unit_ids` list / NPZ mask still applies |
| Region order | Larger population first, used as the shuffled population |
| Window | 0–2 s after stimulus onset; equivalent to source samples 2000:4000 for your 1 ms, 2 s pre-onset arrays |
| Preprocessing | Count spikes, first differences across chronological trials within each stimulus, then trial-axis z-score |
| Statistic | Spearman rho and SciPy's two-sided asymptotic Spearman p-value |
| Pair shuffle | 10,000 permutations of the first population; reset MT19937 seed 0 per pair and condition, matching `np.random.seed(0)` |
| Pair flag | `p_spearman < 0.05`, **not** a Bonferroni-corrected permutation p-value |
| Population statistic | Fraction of pairs meeting that flag, per stimulus and session |
| Population control | Same fraction for each repeat across synchronized pair shuffles |
| Reference interval | 5th–95th percentiles (`inverted_cdf`), after the source's tiny ±0.0005 jitter |
| Additional control | GC–GC and PC–PC fractions; same method, clearly labeled |

This is a reproduction of the reference significance-fraction analysis, not a
change of multiplicity correction chosen to make the current data significant.
The old independent-per-pair shuffle archive cannot supply the synchronized
fraction null: **rerun spike analysis**. We retain the pair rho/p/percentile values
and both raw and jittered fraction nulls for audit.

The repository profile overrides `windows_s`, `rolling_s`, `detrend`,
`single_units_only`, and `min_trials` in old configs. Resolved settings are saved in
preflight/provenance. It uses evoked only, no rolling windows, first differences
and z-score, all mapped units, and at least three differences for finite inference.
Reviewed NPZ good-unit/trial masks remain effective. Invalid/constant pairs count
as non-significant in the fraction denominator, as the source's `NaN < .05`
comparison does; both total and valid pair counts are exported.

The source aggregates inter-region fractions by stimulus and intra-region fractions
by region across stimuli. Our main table keeps stimulus-specific values for both;
the intra-region pooled table also reports the source-style across-stimulus summary.
No flavor/mixture-specific renaming or hard-coded animal exclusions are applied.

### Feb2 commands

Use your already verified `configs/spike_noise.feb2.local.json` and single-session list.
No change to their timing or electrode keys is needed.

```bash
python -m pip install -e '.[spikes]'
RUN_STAMP=$(date +%Y%m%d_%H%M%S)
coupling-spike-noise \
  configs/fin_inter_list.spike_feb2.local.txt \
  configs/spike_noise.feb2.local.json \
  "results/spike_repository_check_${RUN_STAMP}" \
  --profile repository --check &&
coupling-spike-noise \
  configs/fin_inter_list.spike_feb2.local.txt \
  configs/spike_noise.feb2.local.json \
  "results/spike_repository_${RUN_STAMP}" \
  --profile repository --shuffles 10000 --seed 0
```

Primary session outputs under `fig08_spike_coupling/`:

- `significant_fractions.pdf` / `.png` / `.svg`: observed fractions vs 5–95% shuffle ranges.
- `significant_fractions.csv`: stimulus, region comparison, denominators, fraction and bounds.
- `fraction_nulls.npz`: repeat-wise fractions, with raw and jittered versions.
- `intra_region_pooled.csv`: across-stimulus within-region fraction summaries.
- `correlations.csv`: rho, `p_spearman`, percentile and pair-level flags. No blanket pairwise Bonferroni column in this profile.
- `nulls.npz`: rho and `_p_spearman` arrays for each `null_key`.
- `summary.pdf`: secondary descriptive rho scatterplot; not the primary reference statistic.

The root `all_session_significant_fractions.csv` and matching figures combine
completed sessions without treating unit pairs as independent animals.
`route_contrasts.csv` remains an optional descriptive rho contrast, not a reference
repository hypothesis test.

## Transition setup and aggregation

Independent four-state models (three transitions per region) are required, with
matching fit intervals, alignment and bin widths. The input NPZ contract is unchanged:
see [posterior export instructions](mahmood_spike_transition.md#transition-correlations-fit-inputs-are-a-separate-prerequisite).
Model fitting itself is upstream, as in the source scripts.

The source uses two related analyses; this profile runs both:

1. **Setup runner:** corresponding transition indices, integer-bin posterior modes,
   Spearman rho, mean absolute timing error, 10,000 trial shuffles, and percentile ranks.
   It also orders trials by summed regional posterior variance and computes moving
   correlations with step 5. For numerical parity, we retain the source's 29-element
   slice in its nominal 30-trial window, and export `n=29` and the source rank endpoint.
2. **All-to-all aggregation:** nine cross-region transition pairs over pooled stimuli,
   5,000 shuffles by default, and flags at/above the 90th percentile of the shuffle null.
   Across sessions, count flags per transition pair, use an upper-tail binomial null
   with probability 0.10, and Bonferroni-correct across transition comparisons.

`--shuffles` controls the all-to-all stage; the matched setup stage retains 10,000.
Seeds make our transition and jitter outputs reproducible. The source's parallel
transition and jitter generation is not fully seeded, so individual random draws
are not bit-for-bit reproductions. We use the complete binomial survival function,
including the all-successes outcome omitted by one source PMF slicing expression.

```bash
RUN_STAMP=$(date +%Y%m%d_%H%M%S)
coupling-transition-corr configs/fin_inter_list.txt \
  configs/transition_corr.local.json "results/transition_repo_check_${RUN_STAMP}" \
  --profile repository --check &&
coupling-transition-corr configs/fin_inter_list.txt \
  configs/transition_corr.local.json "results/transition_repo_${RUN_STAMP}" \
  --profile repository --seed 0
```

Outputs include `matched_transitions.csv`, `matched_nulls.npz`,
`uncertainty_ranked_windows.csv` (when >=30 accepted trials), all-to-all
`correlations.csv` / `nulls.npz`, and run-level `transition_population.csv` plus
heatmaps. The source's misleading `mse` name is replaced by `mae_ms`; raw bin timing
errors are converted to milliseconds, while ranks/percentiles are unchanged.

## Deliberate differences and limits

- GC/PC labels replace BLA/GC; six stimulus labels are preserved without hard-coded
  taste renaming. Transition stimuli are pooled in this profile, so shared condition
  effects can contribute to correlation. Use `--profile adapted` for separated results.
- Laser strata remain separate, and joint masks / explicit trial IDs prevent accidental
  cross-condition or mismatched-trial comparisons. The reference has no such general guard.
- Region inclusion and sampling metadata remain explicit. Source recording HDF5 files
  are read-only; portable CSV/NPZ outputs replace in-place PyTables writes and unsafe
  model-pickle loading.
- Cross-session binomial calculations reproduce the source's session-level approach.
  Repeated sessions from the same animal are not independent biological replicates;
  these outputs do not establish animal-level inference.
- Trial differencing can introduce serial dependence. Matching the source does not
  eliminate that inferential limitation or validate spike sorting / model convergence.
- This ports the core spike setup/fraction aggregation and inter-region transition
  setup/all-to-all aggregation. It does **not** reproduce every repository script:
  separate split-population model fits, rolling-noise variants, PCA/firing-rate drift
  plots, mixture/child normalization and ancillary regressions remain outside this profile.

Source links:

- [Spike setup](https://github.com/abuzarmahmood/Mahmood_2023_BLA_GC/blob/9c54acb57e227dc9d13597750165b3ebc7be68e9/spike_noise_corrs/inter_region_noise_corrs_setup.py)
- [Spike aggregation](https://github.com/abuzarmahmood/Mahmood_2023_BLA_GC/blob/9c54acb57e227dc9d13597750165b3ebc7be68e9/spike_noise_corrs/inter_region_noise_corrs_aggregate.py)
- [Transition setup](https://github.com/abuzarmahmood/Mahmood_2023_BLA_GC/blob/9c54acb57e227dc9d13597750165b3ebc7be68e9/transition_corrs/all_tastes/inter_region/run_inter_corrs.py)
- [Transition aggregation](https://github.com/abuzarmahmood/Mahmood_2023_BLA_GC/blob/9c54acb57e227dc9d13597750165b3ebc7be68e9/transition_corrs/all_tastes/inter_region/inter_corrs_aggregate.py)
