import contextlib
import io
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest
import numpy as np
from scipy.stats import spearmanr
from coupling.correlation import correlation_test
from coupling import spike_coupling as spike, transition_coupling as transition
from coupling.correlation_cli import run, route_contrasts


def trial_data(n=24):
    return dict(trial_id=np.array([f'trial{i}' for i in range(n)]), stimulus=np.array(['BA']*n),
                stratum=np.array(['off']*n), trial_order=np.arange(n), good_trial=np.ones(n, bool))


def spike_data():
    rng = np.random.default_rng(12)
    x = rng.poisson(1, (24, 4, 40))
    x[:, 2] = x[:, 0]
    return dict(spikes=x, time_s=np.arange(40)/10-2, unit_id=np.arange(4),
                region=np.array(['GC', 'GC', 'PC', 'PC']), good_unit=np.ones(4, bool), **trial_data())


def save_posteriors(folder, mismatch=False):
    rng = np.random.default_rng(3)
    tau = np.sort(rng.uniform(1, 40, (24, 2)), axis=-1)
    result = {}
    for region in ('GC', 'PC'):
        order = np.arange(24) if region == 'GC' else np.arange(24)[::-1]
        metadata = dict(region=region, model_id=region+'_fit', independently_fitted=True,
                        fit_qc_passed=True, fit_qc_notes='Synthetic fixture', bin_width_ms=50,
                        fit_start_ms=0, fit_end_ms=2000, alignment_event='DIN')
        data = {k: v[order] for k, v in trial_data().items()}
        if mismatch and region == 'PC':
            data['trial_id'][0] = 'other'
        path = folder/f'{region}.npz'
        np.savez(path, tau=np.tile(tau[order], (8, 1, 1)), metadata_json=json.dumps(metadata), **data)
        result[region] = str(path)
    return {'posterior_files': result}


class CorrelationTests(unittest.TestCase):
    def test_rank_ties_and_reproducibility(self):
        x = np.array([1, 1, 2, 3, 5, 5, 8, 9])
        y = np.array([2, 3, 1, 3, 4, 6, 5, 8])
        a, null = correlation_test(x, y, 99, np.random.default_rng(0))
        b, other = correlation_test(x, y, 99, np.random.default_rng(0))
        self.assertAlmostEqual(a['rho'], spearmanr(x, y).statistic)
        np.testing.assert_array_equal(null, other)
        self.assertGreaterEqual(a['p_perm'], .01)
        self.assertEqual(a, b)

    def test_degenerate(self):
        self.assertEqual(correlation_test(np.ones(12), np.arange(12), 9, np.random.default_rng())[0]['status'], 'constant')
        self.assertEqual(correlation_test(np.arange(3), np.arange(3), 9, np.random.default_rng())[0]['status'], 'too_few_trials')

    def test_spike_difference_pairs_and_masks(self):
        data = spike_data()
        data['good_trial'][4] = False
        with contextlib.redirect_stderr(io.StringIO()):
            rows, _ = spike.analyze(data, {'windows_s': {'evoked': [0, 2]}}, 99, 1)
        self.assertEqual(len(rows), 6)  # 4 inter + 1 GC-GC + 1 PC-PC
        self.assertEqual(rows[0]['n'], 21)  # do not bridge rejected trial
        same = next(r for r in rows if r['unit_a'] == '0' and r['unit_b'] == '2')
        self.assertAlmostEqual(same['rho'], 1)
        self.assertLess(same['p_perm'], .05)
        counts = data['spikes'][:, :, 20:].sum(-1)
        keep = data['good_trial'][1:] & data['good_trial'][:-1]
        self.assertAlmostEqual(rows[0]['rho'], spearmanr(np.diff(counts[:, 0])[keep], np.diff(counts[:, 1])[keep]).statistic)

    def test_strata_never_mixed(self):
        data = spike_data()
        data['stratum'] = np.array(['off', 'on']*12)
        with contextlib.redirect_stderr(io.StringIO()):
            rows, _ = spike.analyze(data, {'windows_s': {'evoked': [0, 2]}}, 9, 1)
        self.assertEqual(len(rows), 12)
        self.assertEqual({r['n'] for r in rows}, {11})

    def test_transition_reordering_modes_uncertainty(self):
        with tempfile.TemporaryDirectory() as temp:
            folder = Path(temp)
            spec = save_posteriors(folder)
            data = transition.load(folder, spec)
            np.testing.assert_equal(data['GC']['mode_ms'], data['PC']['mode_ms'])
            with contextlib.redirect_stderr(io.StringIO()):
                rows, null = transition.analyze(data, {**spec, 'posterior_draws': 10}, 99, 0)
            self.assertEqual(len(rows), 4)
            self.assertAlmostEqual(rows[0]['rho'], 1)
            self.assertEqual(rows[0]['mae_ms'], 0)
            self.assertAlmostEqual(rows[0]['gc_posterior_variance_ms2'], 0)
            self.assertIn('test_0_posterior_rho', null)
            bad = save_posteriors(folder, mismatch=True)
            with self.assertRaisesRegex(ValueError, 'trial ID sets'):
                transition.load(folder, bad)

    @unittest.skipUnless(importlib.util.find_spec('h5py') and importlib.util.find_spec('hdf5plugin'), 'install .[spikes] for HDF5 tests')
    def test_h5_region_mapping_and_opto(self):
        import h5py
        import hdf5plugin
        with tempfile.TemporaryDirectory() as temp:
            folder = Path(temp)
            info = {'electrode_layout': {'gc': [[8]], 'pc': [[3]]}, 'laser_params': {'dig_in_nums': [7]}}
            (folder/'session.info').write_text(json.dumps(info))
            path = folder/'session.h5'
            units = np.array([(3, 1), (8, 1), (8, 0)], dtype=[('electrode_number', 'i4'), ('single_unit', 'i4')])
            with h5py.File(path, 'w') as h:
                h['unit_descriptor'] = units
                group = h.create_group('spike_trains/dig_in_0')
                group.create_dataset('spike_array', data=np.zeros((10, 3, 400)), **hdf5plugin.Blosc())
                group['laser_durations'] = np.array([0, 500]*5)
                group['laser_onset_lag'] = np.zeros(10)
            spec = {'spike_fs_hz': 100, 'onset_s': 2, 'stimulus_by_din': {'0': 'BA'}}
            data = spike.load(folder, spec)
            self.assertEqual(data['region'].tolist(), ['PC', 'GC', 'GC'])
            self.assertEqual(data['good_unit'].tolist(), [True, True, False])
            self.assertEqual(len(set(data['stratum'])), 2)
            with h5py.File(path, 'a') as h:
                del h['spike_trains/dig_in_0/laser_durations']
            with self.assertRaisesRegex(ValueError, 'Optogenetic'):
                spike.load(folder, spec)

    def test_batch_both_end_to_end_readonly(self):
        with tempfile.TemporaryDirectory() as temp:
            folder = Path(temp)
            session = folder/'session'
            session.mkdir()
            source = session/'spikes.npz'
            np.savez(source, **spike_data())
            before = source.read_bytes()
            sessions = folder/'sessions.txt'
            sessions.write_text(str(session)+'\n')
            for kind in ('spike', 'transition'):
                spec = {'input': 'spikes.npz', 'preprocessing': 'Synthetic fixture, verified trial order'} if kind == 'spike' else {**save_posteriors(session), 'preprocessing': 'Synthetic independently fitted example', 'posterior_draws': 5}
                cfg = folder/f'{kind}.json'
                cfg.write_text(json.dumps({'defaults': spec}))
                with contextlib.redirect_stderr(io.StringIO()), contextlib.redirect_stdout(io.StringIO()):
                    run(kind, sessions, cfg, folder/(kind+'_check'), 9, 0, True)
                    run(kind, sessions, cfg, folder/kind, 9, 0)
                self.assertEqual(json.loads((folder/kind/'status.json').read_text())['state'], 'complete')
                self.assertTrue((folder/kind/'all_sessions.csv').exists())
                self.assertEqual(len(list((folder/kind).rglob('summary.png'))), 1)
            self.assertEqual(source.read_bytes(), before)

    def test_route_pairing_and_invalid_inputs(self):
        row = dict(stimulus='BA', stratum='off', window='evoked', unit_a='1', unit_b='2',
                   comparison='GC-PC', status='ok', rho=.1, n=20)
        rows = [row, {**row, 'stimulus': 'oBA', 'rho': .4},
                {**row, 'stimulus': 'oBA', 'stratum': 'on', 'rho': .9}]
        contrasts = route_contrasts(rows, {'BA': ['BA', 'oBA']}, 'spike')
        self.assertEqual(len(contrasts), 1)
        self.assertAlmostEqual(contrasts[0]['delta_ortho_minus_retro'], .3)
        with tempfile.TemporaryDirectory() as temp:
            folder = Path(temp)
            data = spike_data()
            data['spikes'] = data['spikes'].astype(float)
            data['spikes'][0, 0, 0] = np.nan
            np.savez(folder/'input.npz', **data)
            with self.assertRaisesRegex(ValueError, 'Accepted spikes'):
                spike.load(folder, {'input': 'input.npz'})
            with self.assertRaisesRegex(ValueError, 'outside recorded'):
                spike.windows({'windows_s': {'bad': [-4, 0]}}, data['time_s'])

    def test_preflight_collects_missing_posteriors(self):
        with tempfile.TemporaryDirectory() as temp:
            folder = Path(temp)
            (folder/'sessions.txt').write_text(str(folder)+'\n')
            (folder/'config.json').write_text(json.dumps({'defaults': {'preprocessing': 'reviewed'}}))
            with self.assertRaisesRegex(ValueError, 'Preflight failed'):
                run('transition', folder/'sessions.txt', folder/'config.json', folder/'out', check=True)
            report = json.loads((folder/'out/preflight.json').read_text())[0]
            self.assertIn('fit each region independently', report['error'])


if __name__ == '__main__':
    unittest.main()
