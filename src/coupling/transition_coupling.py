"""Correlation of independently fitted GC/PC changepoint posterior modes."""
import json
import itertools
import numpy as np
from scipy.stats import mode, spearmanr
from .phase_coherence import resolve_path
from .correlation import validate_trials, conditions, condition_indices, correlation_test, correct, progress


def load(session, spec):
    sources = []
    regions = []
    for region in ('GC', 'PC'):
        if region not in spec.get('posterior_files', {}):
            raise ValueError(f'Missing posterior_files.{region}: fit each region independently and export tau posterior NPZ first')
        path = resolve_path(spec['posterior_files'][region], session)
        with np.load(path, allow_pickle=False) as archive:
            data = {k: archive[k] for k in ('tau', 'trial_id', 'stimulus', 'stratum', 'trial_order', 'good_trial')}
            metadata = json.loads(str(archive['metadata_json'].item()))
        tau = data['tau']
        if tau.ndim != 3 or min(tau.shape) < 1 or tau.shape[0] < 2:
            raise ValueError('tau must be posterior draw x trial x transition; at least two draws')
        validate_trials(data, tau.shape[1])
        if not np.isfinite(tau).all() or np.any(tau < 0) or np.any(np.diff(tau, axis=-1) < 0):
            raise ValueError('tau must be finite, nonnegative, ordered bin coordinates')
        if metadata.get('region') != region or metadata.get('independently_fitted') is not True or metadata.get('fit_qc_passed') is not True:
            raise ValueError('Metadata must declare region, independently_fitted=true and fit_qc_passed=true')
        if not metadata.get('model_id') or not metadata.get('fit_qc_notes'):
            raise ValueError('Record model_id and fit_qc_notes (convergence and posterior predictive assessment)')
        for key in ('bin_width_ms', 'fit_start_ms', 'fit_end_ms'):
            if not np.isfinite(float(metadata[key])):
                raise ValueError(f'Invalid {key}')
        width = float(metadata['bin_width_ms'])
        start, end = float(metadata['fit_start_ms']), float(metadata['fit_end_ms'])
        if width <= 0 or start >= end or np.any(tau*width+start > end):
            raise ValueError('Posterior coordinates exceed declared fit interval')
        data['mode_ms'] = mode(tau.astype(np.int64), axis=0, keepdims=False).mode*width+start
        data['tau_ms'] = tau*width+start
        data['variance_ms2'] = np.var(data['tau_ms'], axis=0)
        data['metadata'] = metadata
        regions.append(data)
        sources.append(path)
    gc, pc = regions
    if gc['metadata']['model_id'] == pc['metadata']['model_id']:
        raise ValueError('GC and PC model_id must identify separate fits')
    if set(gc['trial_id']) != set(pc['trial_id']):
        raise ValueError('GC/PC trial ID sets differ; export the same trial universe with explicit masks')
    lookup = {v: i for i, v in enumerate(pc['trial_id'])}
    order = np.array([lookup[v] for v in gc['trial_id']])
    for key in ('trial_id', 'stimulus', 'stratum', 'trial_order', 'good_trial', 'mode_ms', 'variance_ms2'):
        pc[key] = pc[key][order]
    pc['tau_ms'] = pc['tau_ms'][:, order]
    for key in ('stimulus', 'stratum', 'trial_order'):
        if not np.array_equal(gc[key], pc[key]):
            raise ValueError(f'GC/PC {key} mismatch after trial matching')
    for key in ('fit_start_ms', 'fit_end_ms', 'alignment_event'):
        if key not in gc['metadata'] or gc['metadata'][key] != pc['metadata'].get(key):
            raise ValueError(f'GC/PC {key} must match')
    return dict(GC=gc, PC=pc, sources=sources, **{k: gc[k] for k in ('trial_id', 'stimulus', 'stratum', 'trial_order')},
                good_trial=gc['good_trial'] & pc['good_trial'])


def analyze(data, spec, repeats, seed):
    rng = np.random.default_rng(seed)
    gc, pc = data['GC'], data['PC']
    rows, nulls = [], {}
    pairs = list(itertools.product(range(gc['mode_ms'].shape[1]), range(pc['mode_ms'].shape[1])))
    posterior_draws = int(spec.get('posterior_draws', 200))
    if posterior_draws < 1:
        raise ValueError('posterior_draws must be positive')
    for stimulus, stratum in conditions(data):
        ix = condition_indices(data, stimulus, stratum)
        ix = ix[data['good_trial'][ix]]
        for a, b in progress(pairs, f'{stimulus}/{stratum}: transition pairs'):
            x, y = gc['mode_ms'][ix, a], pc['mode_ms'][ix, b]
            stat, null = correlation_test(x, y, repeats, rng, spec.get('min_trials', 8))
            key = f'test_{len(rows)}'
            extra = dict(mae_ms=np.nan, p_mae=np.nan, gc_posterior_variance_ms2=np.nan,
                         pc_posterior_variance_ms2=np.nan, posterior_rho_low=np.nan, posterior_rho_high=np.nan)
            nulls[key] = null
            if len(ix):
                extra.update(mae_ms=float(np.mean(abs(x-y))),
                             gc_posterior_variance_ms2=float(gc['variance_ms2'][ix, a].mean()),
                             pc_posterior_variance_ms2=float(pc['variance_ms2'][ix, b].mean()))
            if len(ix) >= spec.get('min_trials', 8):
                mae_null = np.array([np.mean(abs(rng.permutation(x)-y)) for _ in range(repeats)])
                extra['p_mae'] = float((1+np.sum(mae_null <= extra['mae_ms']+1e-12))/(repeats+1))
                nulls[key+'_mae'] = mae_null
                draw_rho = []
                for _ in range(posterior_draws):
                    gx = gc['tau_ms'][rng.integers(gc['tau_ms'].shape[0]), ix, a]
                    py = pc['tau_ms'][rng.integers(pc['tau_ms'].shape[0]), ix, b]
                    if np.ptp(gx) and np.ptp(py):
                        draw_rho.append(float(spearmanr(gx, py).statistic))
                if draw_rho:
                    extra['posterior_rho_low'], extra['posterior_rho_high'] = map(float, np.quantile(draw_rho, [.025, .975]))
                nulls[key+'_posterior_rho'] = np.array(draw_rho)
            rows.append(dict(stimulus=stimulus, stratum=stratum, gc_transition=a, pc_transition=b,
                             comparison='GC-PC', null_key=key, **stat, **extra))
    correct(rows)
    return rows, nulls


def main():
    from .correlation_cli import main as run
    run('transition')


if __name__ == '__main__':
    main()
