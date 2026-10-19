"""Portable reproduction of the reference repository's core statistics.

Source algorithms: spike_noise_corrs/{setup,aggregate}, transition_corrs/all_tastes/
inter_region/{run_inter_corrs,inter_corrs_aggregate}; revision is pinned in the driver.
This module intentionally does not load pickles or modify recording HDF5 files.
"""
import itertools
import numpy as np
from scipy import stats
from .correlation import conditions, condition_indices, progress


def spearman_p(rho, n):
    """SciPy spearmanr's two-sided asymptotic t approximation, vectorized."""
    r = np.clip(np.asarray(rho, float), -1, 1)
    with np.errstate(divide='ignore', invalid='ignore'):
        t = r*np.sqrt((n-2)/((r+1)*(1-r)).clip(0))
    return 2*stats.t.sf(abs(t), n-2)


def rank_test(x, y, repeats, seed=0):
    """Reset MT19937 per pair, as np.random.seed(0) in reference calc_corrs.

Permuting average ranks gives the same Spearman statistic as reranking each
permuted vector. Integer/half-integer rank dot products resolve percentile ties
without floating point differences between BLAS and scipy.corrcoef.
"""
    x, y = np.asarray(x), np.asarray(y)
    if x.shape != y.shape or x.ndim != 1 or repeats < 1:
        raise ValueError('Require matched vectors and positive shuffle count')
    result = dict(n=len(x), rho=np.nan, p_spearman=np.nan, percentile=np.nan,
                  null_mean=np.nan, null_low=np.nan, null_high=np.nan, status='ok')
    empty = np.full(repeats, np.nan)
    if len(x) < 3:
        result['status'] = 'too_few_trials'
        return result, empty, empty.copy()
    if not np.isfinite(x).all() or not np.isfinite(y).all() or not np.ptp(x) or not np.ptp(y):
        result['status'] = 'constant_or_invalid'
        return result, empty, empty.copy()
    a, b = stats.rankdata(x), stats.rankdata(y)
    a -= a.mean()
    b -= b.mean()
    denominator = np.linalg.norm(a)*np.linalg.norm(b)
    score = a @ b
    rng = np.random.RandomState(seed)
    scores = np.empty(repeats)
    for start in range(0, repeats, 256):
        count = min(256, repeats-start)
        scores[start:start+count] = np.array([rng.permutation(a) for _ in range(count)]) @ b
    rho = float(np.clip(score/denominator, -1, 1))
    null = np.clip(scores/denominator, -1, 1)
    result.update(rho=rho, p_spearman=float(spearman_p(rho, len(x))),
                  percentile=float(stats.percentileofscore(scores, score)),
                  null_mean=float(null.mean()), null_low=float(np.percentile(null, 2.5)),
                  null_high=float(np.percentile(null, 97.5)))
    return result, null, spearman_p(null, len(x))


def spike_analyze(data, spec, repeats, seed):
    # Reference sorts populations by size before choosing the shuffled side.
    groups = [np.flatnonzero(data['good_unit'] & (data['region'] == r)) for r in ('GC', 'PC')]
    groups = [groups[i] for i in np.argsort([len(g) for g in groups])[::-1]]
    pairs = list(itertools.product(*groups))
    pairs += [pair for g in groups for pair in itertools.combinations(g, 2)]
    t = data['time_s']
    if not np.allclose(np.diff(t), .001) or t[0] > 0 or t[-1]+.001 < 2-1e-9:
        raise ValueError('Repository spike profile requires 1 ms bins covering 0–2 s')
    counts = data['spikes'][..., (t >= -1e-10) & (t < 2-1e-10)].sum(-1, dtype=float)
    rows, nulls = [], {}
    for stimulus, stratum in conditions(data):
        ix = condition_indices(data, stimulus, stratum)
        good = data['good_trial'][ix]
        differences = np.diff(counts[ix], axis=0)[good[:-1] & good[1:]]
        with np.errstate(divide='ignore', invalid='ignore'):
            values = stats.zscore(differences, axis=0) if len(differences) else differences
        for i, j in progress(pairs, f'{stimulus}/{stratum}: repository spike pairs'):
            result, null, null_p = rank_test(values[:, i], values[:, j], repeats, seed)
            key = f'test_{len(rows)}'
            nulls[key], nulls[key+'_p_spearman'] = null, null_p
            ri, rj = data['region'][i], data['region'][j]
            rows.append(dict(stimulus=stimulus, stratum=stratum, window='evoked', start_s=0., end_s=2.,
                unit_a=str(data['unit_id'][i]), unit_b=str(data['unit_id'][j]),
                comparison=f'{ri}-{rj}' if ri == rj else 'GC-PC',
                detrend='first_difference_zscore', trials_total=len(ix), trials_accepted=int(good.sum()),
                null_key=key, **result, significant_spearman=bool(result['p_spearman'] < .05)))
    return rows, nulls


def spike_fractions(rows, nulls, seed=0):
    """Source aggregation: fraction p<.05 by taste and repeat, 5–95% interval.

All eligible pairs form the denominator, as in the reference comparison of
NaN p-values with .05. Invalid-pair counts and statuses remain explicit.
"""
    groups = sorted({(r['stimulus'], r['stratum'], r['comparison']) for r in rows})
    summaries, distributions = [], {}
    jitter_rng = np.random.RandomState(seed)
    for stimulus, stratum, comparison in groups:
        selected = [r for r in rows if (r['stimulus'], r['stratum'], r['comparison']) == (stimulus, stratum, comparison)]
        values = np.array([nulls[r['null_key']+'_p_spearman'] < .05 for r in selected]).mean(0)
        actual = float(np.mean([r['p_spearman'] < .05 for r in selected]))
        # The source adds +/-0.0005 jitter before extracting inverted-CDF limits.
        jittered = values+(jitter_rng.random_sample(len(values))-.5)*.001
        low, high = np.percentile(jittered, [5, 95], method='inverted_cdf')
        key = f'fraction_{len(summaries)}'
        distributions[key] = values
        distributions[key+'_jittered'] = jittered
        summaries.append(dict(stimulus=stimulus, stratum=stratum, comparison=comparison,
            n_pairs=len(selected), n_valid_pairs=sum(r['status']=='ok' for r in selected),
            fraction_significant=actual, shuffle_mean=float(values.mean()),
            shuffle_low_5=max(0., float(low)), shuffle_high_95=float(high),
            above_shuffle_95=bool(actual >= high), null_key=key))
    return summaries, distributions


def intra_region_pooled(rows):
    """Source intra-region aggregation pools pair/taste flags within each region."""
    output = []
    for stratum, comparison in sorted({(r['stratum'], r['comparison']) for r in rows if r['comparison'] != 'GC-PC'}):
        selected = [r for r in rows if (r['stratum'], r['comparison']) == (stratum, comparison)]
        output.append(dict(stratum=stratum, comparison=comparison, pair_stimulus_tests=len(selected),
                           valid_tests=sum(r['status']=='ok' for r in selected),
                           fraction_significant=float(np.mean([r['p_spearman'] < .05 for r in selected]))))
    return output


def transition_analyze(data, spec, repeats, seed):
    """Pooled-stimulus all-to-all analysis from inter_corrs_aggregate.py.

Keep laser strata separate as a declared experimental adaptation. Four-state
fits (three transitions) and matching fit bin widths follow the source selection.
"""
    gc, pc = data['GC'], data['PC']
    if gc['mode_ms'].shape[1] != 3 or pc['mode_ms'].shape[1] != 3:
        raise ValueError('Repository transition profile requires independently fitted four-state models (3 transitions)')
    if gc['metadata']['bin_width_ms'] != pc['metadata']['bin_width_ms']:
        raise ValueError('Repository profile requires matching bin widths')
    rows, nulls = [], {}
    for stratum in sorted(set(data['stratum'])):
        ix = np.flatnonzero(data['good_trial'] & (data['stratum'] == stratum))
        # Source flattens stimuli before trials; chronology is retained within each.
        ix = np.array(sorted(ix, key=lambda k: (str(data['stimulus'][k]), data['trial_order'][k])), dtype=int)
        for a, b in progress(list(itertools.product(range(3), repeat=2)), f'{stratum}: pooled transition pairs'):
            x, y = gc['mode_ms'][ix, a], pc['mode_ms'][ix, b]
            # Source helper shuffles the SECOND region for all-to-all correlations.
            result, null, null_p = rank_test(y, x, repeats, seed+len(rows))
            key = f'test_{len(rows)}'
            nulls[key], nulls[key+'_p_spearman'] = null, null_p
            rows.append(dict(stimulus='ALL_POOLED', stratum=str(stratum), gc_transition=a, pc_transition=b,
                             comparison='GC-PC', null_key=key, **result,
                             above_shuffle_90=bool(result['percentile'] >= 90)))
    return rows, nulls


def transition_setup(data, seed=0, repeats=10000):
    """Matched-transition runner: rho, MAE, uncertainty-ranked moving windows."""
    rows, moving, nulls = [], [], {}
    gc, pc = data['GC'], data['PC']
    for stratum in sorted(set(data['stratum'])):
        ix = np.flatnonzero(data['good_trial'] & (data['stratum'] == stratum))
        ix = np.array(sorted(ix, key=lambda k: (str(data['stimulus'][k]), data['trial_order'][k])), dtype=int)
        for a in range(3):
            x, y = gc['mode_ms'][ix, a], pc['mode_ms'][ix, a]
            result, null, _ = rank_test(x, y, repeats, seed+len(rows))
            key = f'matched_{len(rows)}'
            rng = np.random.RandomState(seed+len(rows))
            mae = float(np.mean(abs(x-y))) if len(ix) else np.nan
            mae_null = np.array([np.mean(abs(rng.permutation(x)-y)) for _ in range(repeats)]) if len(ix) else np.full(repeats, np.nan)
            nulls[key], nulls[key+'_mae_ms'] = null, mae_null
            rows.append(dict(stratum=str(stratum), transition=a, null_key=key, **result,
                             mae_ms=mae, mae_percentile=float(stats.percentileofscore(mae_null, mae))))
            uncertainty = gc['variance_ms2'][ix, a]+pc['variance_ms2'][ix, a]
            order = np.argsort(uncertainty)
            # Reproduce the source's [start:start+29] slices for window_size=30.
            for start in range(0, len(ix)-30+1, 5):
                chosen = order[start:start+29]
                xx, yy = x[chosen], y[chosen]
                rho, p = stats.spearmanr(xx, yy) if np.ptp(xx) and np.ptp(yy) else (np.nan, np.nan)
                moving.append(dict(stratum=str(stratum), transition=a, uncertainty_rank_start=start,
                                   source_max_t=start+29, n=len(chosen), rho=float(rho), p_spearman=float(p)))
    return rows, moving, nulls


def transition_population(rows):
    """Source count of sessions above the 90th percentile, binomial upper tail."""
    result = []
    groups = sorted({(r['stratum'], r['gc_transition'], r['pc_transition']) for r in rows})
    for stratum, a, b in groups:
        selected = [r for r in rows if (r['stratum'], r['gc_transition'], r['pc_transition']) == (stratum, a, b)]
        valid = [r for r in selected if r['status'] == 'ok']
        n = len(valid)
        k = sum(r['above_shuffle_90'] for r in valid)
        p = float(stats.binom.sf(k-1, n, .1)) if n else np.nan
        result.append(dict(stratum=stratum, gc_transition=a, pc_transition=b, sessions_total=len(selected),
                           sessions_valid=n, sessions_above_90=k, fraction=k/n if n else np.nan,
                           p_binomial=p, p_binomial_bonferroni=min(1., p*len(groups)) if n else np.nan))
    return result
