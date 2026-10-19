"""Parity checks against small direct translations of reference formulas."""
import contextlib
import io
import json
from pathlib import Path
import tempfile
import unittest
import numpy as np
from scipy import stats
from coupling import mahmood_repository as repo
from coupling.correlation_cli import run
from test_spike_transition import trial_data, save_posteriors
from coupling.transition_coupling import load as load_transitions


class RepositoryTests(unittest.TestCase):
    def test_rank_test_matches_literal_source_with_ties(self):
        x = np.array([0, 1, 2, 4, 5, 1, 3, 3, 7, 0, 4, 8.])
        y = np.array([1, 2, 4, 3, 2, 2, 4, 5, 6, 2, 0, 7.])
        x, y = stats.zscore(x), stats.zscore(y)
        rng = np.random.RandomState(0)
        direct = np.array([stats.spearmanr(rng.permutation(x), y) for _ in range(100)])
        actual, null, null_p = repo.rank_test(x, y, 100)
        np.testing.assert_allclose(null, direct[:, 0], atol=1e-14)
        np.testing.assert_allclose(null_p, direct[:, 1], atol=1e-14)
        self.assertAlmostEqual(actual['rho'], stats.spearmanr(x, y).statistic)
        self.assertAlmostEqual(actual['p_spearman'], stats.spearmanr(x, y).pvalue)
        self.assertAlmostEqual(actual['percentile'], stats.percentileofscore(direct[:, 0], actual['rho']))

    def test_synchronized_shuffles_and_fraction_denominator(self):
        x = np.arange(20)
        first, null, ps = repo.rank_test(x, x, 100)
        second, null2, ps2 = repo.rank_test(x*2, x*3, 100)
        np.testing.assert_equal(null, null2)
        np.testing.assert_equal(ps, ps2)
        rows = [dict(stimulus='BA', stratum='off', comparison='GC-PC', null_key='a', **first),
                dict(stimulus='BA', stratum='off', comparison='GC-PC', null_key='b', **second),
                dict(stimulus='BA', stratum='off', comparison='GC-PC', null_key='c', p_spearman=np.nan, status='constant_or_invalid')]
        values = dict(a_p_spearman=ps, b_p_spearman=ps2, c_p_spearman=np.full(100, np.nan))
        summary, nulls = repo.spike_fractions(rows, values)
        self.assertEqual(summary[0]['n_valid_pairs'], 2)
        self.assertAlmostEqual(summary[0]['fraction_significant'], 2/3)
        np.testing.assert_equal(nulls['fraction_0'], (ps<.05)*2/3)
        rng = np.random.RandomState(0)
        jitter = (ps<.05)*2/3+(rng.random_sample(100)-.5)*.001
        np.testing.assert_allclose(nulls['fraction_0_jittered'], jitter)
        self.assertAlmostEqual(summary[0]['shuffle_high_95'], np.percentile(jitter,95,method='inverted_cdf'))

    def test_spike_engine_matches_source_counts_diff_zscore(self):
        rng = np.random.default_rng(3)
        x = rng.poisson(.005, (24, 4, 4000))
        data = dict(spikes=x, time_s=np.arange(4000)/1000-2, unit_id=np.arange(4),
                    region=np.array(['GC','PC','PC','PC']), good_unit=np.ones(4,bool), **trial_data())
        with contextlib.redirect_stderr(io.StringIO()):
            rows, nulls = repo.spike_analyze(data, {}, 99, 0)
        row = rows[0]
        # Larger region PC is array1: this matters for synchronized shuffles.
        self.assertEqual((row['unit_a'], row['unit_b']), ('1','0'))
        self.assertEqual(len(rows), 6)
        self.assertEqual(repo.intra_region_pooled(rows)[0]['pair_stimulus_tests'], 3)
        counts = x[...,2000:4000].sum(-1)
        z = stats.zscore(np.diff(counts,axis=0),axis=0)
        self.assertAlmostEqual(row['rho'], stats.spearmanr(z[:,1],z[:,0]).statistic)
        rng = np.random.RandomState(0)
        expected = [stats.spearmanr(rng.permutation(z[:,1]),z[:,0]).statistic for _ in range(99)]
        np.testing.assert_allclose(nulls[row['null_key']], expected, atol=1e-14)

    def test_pooled_transition_and_population(self):
        with tempfile.TemporaryDirectory() as temp:
            folder = Path(temp)
            spec = save_posteriors(folder)
            for path in spec['posterior_files'].values():
                with np.load(path) as z:
                    values = {k:z[k] for k in z.files}
                # Three transitions/four states, with explicit two-stimulus labels.
                values['tau'] = np.concatenate([values['tau'], np.full((8,24,1),39.9)],axis=-1)
                values['stimulus'] = np.array(['BA' if int(v[5:])<12 else 'oBA' for v in values['trial_id']])
                np.savez(path, **values)
            data = load_transitions(folder, spec)
            with contextlib.redirect_stderr(io.StringIO()):
                rows, _ = repo.transition_analyze(data, {}, 50, 0)
            self.assertEqual(len(rows),9)
            self.assertEqual({r['n'] for r in rows},{24})
            self.assertEqual({r['stimulus'] for r in rows},{'ALL_POOLED'})
            pop = repo.transition_population([dict(r,session='a') for r in rows])
            self.assertEqual(pop[0]['sessions_above_90'],1)
            self.assertAlmostEqual(pop[0]['p_binomial'], .1)
            self.assertAlmostEqual(pop[0]['p_binomial_bonferroni'], .9)
            spec['preprocessing'] = 'synthetic independent models'
            (folder/'list.txt').write_text(str(folder)+'\n')
            (folder/'config.json').write_text(json.dumps({'defaults':spec}))
            with contextlib.redirect_stderr(io.StringIO()),contextlib.redirect_stdout(io.StringIO()):
                run('transition', folder/'list.txt', folder/'config.json', folder/'output', 20)
            self.assertTrue((folder/'output/transition_population.pdf').exists())
            self.assertTrue((folder/'output'/folder.name/'fig07_transition_coupling/matched_transitions.csv').exists())


    def test_source_moving_window_off_by_one_is_explicit(self):
        n=40
        x=np.tile(np.arange(n)[:,None],(1,3)).astype(float)
        data=dict(GC={'mode_ms':x,'variance_ms2':x}, PC={'mode_ms':x,'variance_ms2':x},
                  stimulus=np.array(['BA']*n), stratum=np.array(['off']*n),
                  good_trial=np.ones(n,bool),trial_order=np.arange(n))
        rows,moving,nulls=repo.transition_setup(data, repeats=20)
        self.assertEqual(len(rows),3)
        self.assertEqual(len(moving),9)
        self.assertEqual({r['n'] for r in moving},{29})
        self.assertEqual(moving[-1]['source_max_t'],39)

    def test_repository_batch_outputs_and_profile_override(self):
        with tempfile.TemporaryDirectory() as temp:
            folder=Path(temp);session=folder/'session';session.mkdir()
            rng=np.random.default_rng(2)
            np.savez(session/'spikes.npz', spikes=rng.poisson(.01,(24,2,4000)),
                     time_s=np.arange(4000)/1000-2,unit_id=np.arange(2),region=np.array(['GC','PC']),
                     good_unit=np.ones(2,bool),**trial_data())
            (folder/'list.txt').write_text(str(session)+'\n')
            (folder/'config.json').write_text(json.dumps({'defaults': {'input':'spikes.npz','preprocessing':'synthetic',
                 'windows_s':{'baseline':[-2,0]},'detrend':'none','single_units_only':True}}))
            with contextlib.redirect_stderr(io.StringIO()),contextlib.redirect_stdout(io.StringIO()):
                run('spike', folder/'list.txt', folder/'config.json', folder/'out', 20)
            out=folder/'out/session/fig08_spike_coupling'
            self.assertTrue((out/'significant_fractions.csv').exists())
            self.assertTrue((out/'significant_fractions.pdf').exists())
            self.assertTrue((folder/'out/all_session_significant_fractions.csv').exists())
            prov=json.loads((out/'provenance.json').read_text())
            self.assertEqual(prov['analysis_profile'],'repository')
            self.assertEqual(prov['config']['windows_s'],{'evoked':[0,2]})
            self.assertFalse(prov['config']['single_units_only'])


if __name__=='__main__':
    unittest.main()
