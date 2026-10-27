
from coupling.array_store import load_arrays, save_arrays
import hashlib
import json
from pathlib import Path
import tempfile
import unittest
try:
    import h5py
except ImportError:
    h5py = None
import numpy as np
from scipy.signal import stft, detrend
from coupling.phase_coherence import (phase_locking, select_channels, calculate,
    sustained_deviation, read_session_list, load_input, run, transform, equal_odor_result, alignment_label, analyze_conditions)


class PhaseTests(unittest.TestCase):
    def test_locked_and_independent(self):
        rng = np.random.default_rng(34)
        a = rng.uniform(-np.pi, np.pi, (2000, 2, 3))
        np.testing.assert_allclose(phase_locking(a, a+.7), 1)
        self.assertLess(phase_locking(a, rng.uniform(-np.pi, np.pi, a.shape)).max(), .07)

    def test_reference_selection_and_bootstrap(self):
        rng = np.random.default_rng(3)
        phase = rng.normal(size=(8, 4, 50, 100))
        groups = [np.array([0, 1]), np.array([2, 3])]
        expected = []
        for group in groups:
            p = phase[:, group]
            expected.append(group[np.argmin(abs(p-p.mean(axis=1)[:, None]).sum(axis=(0,2,3)))])
        self.assertEqual(select_channels(phase, groups, 'repository'), expected)
        f, t = np.arange(0, 100, 2), np.arange(100)*.01-1
        r = calculate(phase, groups, f, t, 'repository', 4, 2)
        again = calculate(phase, groups, f, t, 'repository', 4, 2)
        np.testing.assert_array_equal(r['bootstrap_band'], again['bootstrap_band'])
        np.testing.assert_allclose(r['coherence'], phase_locking(phase[:, expected[0]], phase[:, expected[1]]))
        # Independent literal implementation of upstream matched-pair resampling.
        rng = np.random.default_rng(2)
        for g in groups:
            rng.choice(g, 2, replace=False)
        idx = rng.integers(8, size=8)
        full = abs(np.exp(1j*(phase[idx, expected[0]]-phase[idx, expected[1]])).mean(axis=0))
        np.testing.assert_allclose(r['bootstrap_band'][0, 0], full[(f>=3)&(f<8)].mean(axis=0))

    def test_equal_odor_coherence_does_not_cancel_opposite_phase_lags(self):
        # Perfectly locked odors with opposite GC-PC lags should each give 1.
        f, t = np.arange(0, 100, 2), np.arange(100)*.01-1
        results = []
        for lag, n in [(0, 3), (np.pi, 9)]:
            phase = np.zeros((n, 4, len(f), len(t)))
            phase[:, 2:] = lag
            results.append(calculate(phase, None, f, t, 'repository', 4, n,
                                     [0, 2], np.array([[0, 1], [2, 3]])))
        route = equal_odor_result(results)
        np.testing.assert_allclose(route['coherence'], 1)
        np.testing.assert_allclose(route['bootstrap_band'], 1)
        # Different trial counts do not change the equal-odor weighting.
        results[1]['coherence'][:] = 0
        np.testing.assert_allclose(equal_odor_result(results)['coherence'], .5)

    def test_alignment_is_not_inferred_from_onset_offset(self):
        self.assertIn('declared parsed-trial event', alignment_label({'onset_s': 2}))
        self.assertEqual(alignment_label({'event_alignment': 'stimulus DIN onset'}),
                         'Time from stimulus DIN onset (s)')
        self.assertIn('inhalation', alignment_label({'event_alignment': 'first inhalation'}))

    def test_signed_terminal_duration(self):
        t = np.arange(60)*.005
        y = np.r_[np.zeros(10), np.ones(20)*2, np.ones(10)*-2, np.ones(20)*2]
        d = sustained_deviation(y, [-1,1], t)
        self.assertTrue(np.all(d[10:30] == 1))
        self.assertTrue(np.all(d[30:40] == 0))
        self.assertTrue(np.all(d[40:] == 1))

    def test_arrays_adapter_preserves_qc_and_time(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            session = root/'raw'; session.mkdir()
            prepared = root/'prepared'; prepared.mkdir()
            rng = np.random.default_rng(1)
            x = rng.normal(size=(4, 6, 5000))
            good = np.ones((4,6), dtype=bool); good[0,0] = False
            # Two trials per condition: rejecting one still preserves each condition.
            x = np.repeat(x, 2, axis=1)
            good = np.repeat(good, 2, axis=1);good[:,1] = True
            provenance = {k:'test' for k in ('source','reference','filtering','qc_method','event_alignment')}
            save_arrays(prepared/'cleaned_lfp.h5', lfp_uv=x, time_s=np.arange(5000)/1000-2,
                     stimulus=np.repeat(['BA','IA','S','oBA','oIA','BK'],2),
                     electrode_id=[1,2,3,4], region=['GC','GC','PC','PC'],
                     good_channel=np.ones(4,dtype=bool), good_trial=good,
                     session_id='test', animal_id='test', provenance_json=json.dumps(provenance))
            data = load_input(session, {'input':'../prepared/cleaned_lfp.h5'})
            self.assertEqual(data['good'].sum(), 11)
            self.assertEqual(data['counts']['BA'], 1)
            self.assertEqual(data['time_s'][0], -2)

    def fixture(self, root):
        session = root/'data'/'session one'
        session.mkdir(parents=True)
        (session/'session.info').write_text(json.dumps({'electrode_layout': {'gc': [[10,12]], 'pc': [[40,45]]}}))
        rng = np.random.default_rng(9)
        t = np.arange(5000)/1000
        with h5py.File(session/'session.h5', 'w') as h:
            h['Parsed_LFP_channels'] = [40,10,45,12]  # Deliberately non-sorted physical mapping.
            for din in [0,1]:
                x = rng.normal(size=(4,3,5000)) + np.sin(2*np.pi*8*t)[None,None,:]
                h[f'Parsed_LFP/dig_in_{din}_LFPs'] = x
        spec = {'fs_hz': 1000, 'onset_s': 2, 'region_keys': {'GC':'gc','PC':'pc'},
                'stimulus_by_din': {'0':'A','1':'B'}, 'eligible_electrode_ids': [10,12,40,45],
                'preprocessing': 'SYNTHETIC VALIDATION ONLY'}
        return session, spec

    @unittest.skipIf(h5py is None, 'Install .[phase] for HDF5 integration tests')
    def test_file_contract_and_end_to_end(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            session, spec = self.fixture(root)
            listing = root/'fin_inter_list.txt'
            listing.write_text('# comment\n\ndata/session one\n')
            self.assertEqual(read_session_list(listing), [session])
            cfg = root/'config.json';cfg.write_text(json.dumps({'defaults':spec}))
            data = load_input(session, spec)
            np.testing.assert_array_equal(data['groups'][0], [1,3])
            before = hashlib.sha256((session/'session.h5').read_bytes()).hexdigest()
            rows = run(listing, cfg, root/'output', profile='paper', n_boot=3)
            self.assertEqual(rows[0]['status'], 'completed', rows)
            self.assertTrue((root/'output'/'all_sessions.pdf').exists())
            folder = Path(rows[0]['output'])
            conditions = json.loads((folder/'conditions.json').read_text())
            self.assertEqual(conditions['ortho']['status'], 'skipped')
            for label in ('A', 'B'):
                target = folder/conditions[label]['directory']
                self.assertTrue((target/'power_context.pdf').exists())
                with load_arrays(target/'phase_coherence.h5') as result:
                    self.assertEqual(set(result['stimulus']), {label})
                    self.assertEqual(len(result['accepted_event_id']), 3)
            self.assertEqual(before, hashlib.sha256((session/'session.h5').read_bytes()).hexdigest())
            with self.assertRaises(FileExistsError):
                run(listing, cfg, root/'output', check=True)
            with self.assertRaises(ValueError):
                run(listing, cfg, session/'output', check=True)
            (session/'extra.h5').touch()
            with self.assertRaises(ValueError):
                load_input(session, spec)
            # Explicit paths resolve an otherwise ambiguous source directory.
            self.assertEqual(load_input(session, {**spec, 'input':'session.h5'})['counts'], {'A':3,'B':3})

    @unittest.skipIf(h5py is None, 'Install .[phase] for HDF5 integration tests')
    def test_stft_matches_loader_and_mapping(self):
        with tempfile.TemporaryDirectory() as directory:
            session, spec = self.fixture(Path(directory))
            data = load_input(session, spec)
            phase, powers, f, t, used, selected, intra, by_stim = transform(data, 'paper', 42, by_stimulus=True)
            for local, original in enumerate(used):
                x = data['x'][data['channels'][original], data['good']]
                ef, et, z = stft(detrend(x, axis=-1), fs=1000, window='hann', nperseg=500, noverlap=495)
                z = z[:, ef<100][:, :, et<5]
                np.testing.assert_allclose(phase[:,local], np.angle(z))
                np.testing.assert_allclose(powers[local], (abs(z)**2).mean(axis=0))
                for label in ('A', 'B'):
                    mask = data['stimulus'][data['good']] == label
                    np.testing.assert_allclose(by_stim[label][local], (abs(z[mask])**2).mean(axis=0))
            np.testing.assert_allclose(t, et[et<5]-2)

    @unittest.skipIf(h5py is None, 'Install .[phase] for HDF5 integration tests')
    def test_bad_mapping_and_qc(self):
        with tempfile.TemporaryDirectory() as directory:
            session, spec = self.fixture(Path(directory))
            with self.assertRaises(ValueError):
                load_input(session, {**spec, 'fs_hz':30000})
            with self.assertRaises(ValueError):
                load_input(session, {**spec, 'eligible_electrode_ids':[10,40]})
            with self.assertRaises(ValueError):
                load_input(session, {**spec, 'preprocessing':'REPLACE_ME'})
            listing = Path(directory)/'list.txt'
            listing.write_text(str(session)+'\n'+str(session)+'\n')
            with self.assertRaises(ValueError):
                read_session_list(listing)


if __name__ == '__main__':
    unittest.main()
