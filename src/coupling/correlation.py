"""Shared rank/permutation statistics for independently implemented Mahmood analyses."""
import csv
import sys
import time
import numpy as np
from scipy.stats import rankdata, percentileofscore


def correlation_test(x, y, repeats, rng, min_trials=8):
    x, y = np.asarray(x, float), np.asarray(y, float)
    if x.shape != y.shape or x.ndim != 1:
        raise ValueError('Correlation inputs must be matched trial vectors')
    if repeats < 1:
        raise ValueError('repeats must be positive')
    result = dict(n=len(x), rho=np.nan, p_perm=np.nan, percentile=np.nan,
                  null_mean=np.nan, null_low=np.nan, null_high=np.nan, status='ok')
    empty = np.full(repeats, np.nan)
    if len(x) < min_trials:
        result['status'] = 'too_few_trials'
        return result, empty
    if not np.isfinite(x).all() or not np.isfinite(y).all():
        raise ValueError('Nonfinite accepted observations')
    a, b = rankdata(x), rankdata(y)
    a -= a.mean()
    b -= b.mean()
    norm = np.linalg.norm(a)*np.linalg.norm(b)
    if norm == 0:
        result['status'] = 'constant'
        return result, empty
    rho = float(np.dot(a, b)/norm)
    # Bounded temporary memory; average ranks preserve ties exactly.
    null = np.empty(repeats)
    for start in range(0, repeats, 256):
        count = min(256, repeats-start)
        permutations = np.array([rng.permutation(a) for _ in range(count)])
        null[start:start+count] = permutations @ b / norm
    result.update(rho=rho, p_perm=float((1+np.sum(abs(null) >= abs(rho)-1e-12))/(repeats+1)),
                  percentile=float(percentileofscore(null, rho)), null_mean=float(null.mean()),
                  null_low=float(np.quantile(null, .025)), null_high=float(np.quantile(null, .975)))
    return result, null


def correct(rows):
    """Bonferroni across all finite tests in a session, separately for each metric."""
    for field in ('p_perm', 'p_mae'):
        valid = [r for r in rows if np.isfinite(r.get(field, np.nan))]
        for r in rows:
            r[field+'_bonferroni'] = min(1., r[field]*len(valid)) if np.isfinite(r.get(field, np.nan)) else np.nan


def write_csv(path, rows):
    if not rows:
        raise ValueError('No output rows')
    with path.open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def validate_trials(data, n):
    for key in ('trial_id', 'stimulus', 'stratum', 'trial_order', 'good_trial'):
        if np.asarray(data[key]).shape != (n,):
            raise ValueError(f'{key} must have one value per trial')
    if len(np.unique(data['trial_id'])) != n:
        raise ValueError('Duplicate trial IDs')
    if data['good_trial'].dtype.kind != 'b':
        raise ValueError('good_trial must be boolean')
    if not np.isfinite(data['trial_order']).all() or len(np.unique(data['trial_order'])) != n:
        raise ValueError('trial_order must be finite and unique')
    if any(np.any(np.asarray(data[k]).astype(str) == '') for k in ('trial_id', 'stimulus', 'stratum')):
        raise ValueError('Empty trial labels')


def conditions(data):
    return sorted(set(zip(data['stimulus'].astype(str), data['stratum'].astype(str))))


def condition_indices(data, stimulus, stratum):
    ix = np.flatnonzero((data['stimulus'] == stimulus) & (data['stratum'] == stratum))
    return ix[np.argsort(data['trial_order'][ix])]


def progress(items, label):
    """Bounded-rate pair progress, including elapsed time and an ETA."""
    total = len(items)
    start = last = time.monotonic()
    print(f'{label}: 0/{total}', file=sys.stderr, flush=True)
    for done, item in enumerate(items, 1):
        yield item
        now = time.monotonic()
        if done == total or now-last >= 2:
            elapsed = now-start
            eta = elapsed/done*(total-done)
            filled = int(20*done/total)
            print(f"{label} [{'#'*filled}{'-'*(20-filled)}] {done}/{total} "
                  f'elapsed {elapsed:.0f}s | ETA {eta:.0f}s', file=sys.stderr, flush=True)
            last = now
