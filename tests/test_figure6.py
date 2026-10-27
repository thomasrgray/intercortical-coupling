
from coupling.array_store import load_arrays, save_arrays
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
import numpy as np
import matplotlib.pyplot as plt
from coupling.figure6 import export
from coupling.phase_coherence import calculate


class Figure6Tests(unittest.TestCase):
    def fixture(self, root):
        source = root/'session'; source.mkdir()
        metadata = {}
        t, f = np.linspace(-1, 1, 11), np.array([4., 8., 12.])
        for label, value in [('oBA', .7), ('BA', .3)]:
            folder = source/label; folder.mkdir()
            save_arrays(folder/'phase_coherence.h5', time_s=t, frequency_hz=f,
                     band=['theta'], band_limits_hz=[[3, 8]], selected_electrode_id=[1, 33],
                     alignment_label='Time from stimulus DIN onset (s)',
                     bootstrap_band=np.full((4, 1, len(t)), value),
                     coherence=np.full((len(f), len(t)), value))
            metadata[label] = {'status': 'completed', 'directory': label}
        (source/'conditions.json').write_text(json.dumps(metadata))
        return source

    def test_real_exports_and_source_values(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp); source = self.fixture(root)
            out = export(root/'panels', phase_session=source)
            for ext in ('pdf', 'svg', 'png'):
                self.assertGreater((out/f'phase_BA_band_0.{ext}').stat().st_size, 100)
            with load_arrays(out/'phase_BA_values.h5') as values:
                np.testing.assert_allclose(values['ortho_coherence']-values['retro_coherence'], .4)
            manifest = json.loads((out/'manifest.json').read_text())
            self.assertEqual(manifest['phase']['omitted_comparisons'], ['IA', 'route'])
            self.assertEqual(len(manifest['sources']), 3)
            self.assertFalse(plt.get_fignums())
            with self.assertRaises(FileExistsError):
                export(out, phase_session=source)

    def test_grid_mismatch_leaves_no_export(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp); source = self.fixture(root)
            path = source/'BA/phase_coherence.h5'
            with load_arrays(path) as d:
                values = dict(d)
            values['time_s'] += .01
            save_arrays(path, **values)
            with self.assertRaisesRegex(ValueError, 'different time_s'):
                export(root/'panels', phase_session=source)
            self.assertFalse((root/'panels').exists())
            self.assertFalse(list(root.glob('.figure6-*')))

    def test_render_failure_is_atomic_and_closes_figures(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp); source = self.fixture(root)
            with patch('matplotlib.figure.Figure.savefig', side_effect=OSError('disk full')):
                with self.assertRaises(OSError):
                    export(root/'panels', phase_session=source)
            self.assertFalse((root/'panels').exists())
            self.assertFalse(plt.get_fignums())

    def test_condition_path_cannot_escape(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp); source = self.fixture(root)
            (source/'conditions.json').write_text(json.dumps({'BA': {'status': 'completed', 'directory': '..'}}))
            with self.assertRaisesRegex(ValueError, 'escapes'):
                export(root/'panels', phase_session=source)

    def test_power_labels_nan_and_dimensions(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            z = np.ones((2, 4, 3, 11)); z[0, 0, 0] = np.nan
            path = root/'session_power.h5'
            save_arrays(path, z_power=z, time_s=np.linspace(-1, 1, 11), frequency_hz=[4,8,12],
                     region=['GC', 'PC'], stimulus=['BA','oBA','IA','oIA'])
            (root/'provenance.json').write_text(json.dumps({'preprocessing': {'event_alignment': 'inhalation'}}))
            with patch('coupling.figure6.heatmap') as draw:
                out = export(root/'panels', power=path)
                self.assertEqual(draw.call_count, 8)
                self.assertIn('inhalation', draw.call_args.args[7])
            with load_arrays(out/'power_values.h5') as data:
                self.assertTrue(np.isnan(data['z_power'][0,0,0]).all())
            save_arrays(path, z_power=z[:, :, :, :-1], time_s=np.linspace(-1,1,11),
                     frequency_hz=[4,8,12], region=['GC','PC'], stimulus=['BA','oBA','IA','oIA'])
            with self.assertRaisesRegex(ValueError, 'Power axes'):
                export(root/'invalid', power=path)

    def test_failed_session_does_not_block_next_preflight(self):
        from coupling.phase_coherence import run
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root/'bad').mkdir(); (root/'good').mkdir()
            listing = root/'sessions.txt'
            listing.write_text(str(root/'bad') + '\n' + str(root/'good') + '\n')
            config = root/'config.json'; config.write_text('{}')
            ready = {'source': root/'good/input.h5', 'counts': {'BA': 3},
                     'electrode_id': np.array([1, 33]), 'channels': np.array([0, 1])}
            with patch('coupling.phase_coherence.load_input', side_effect=[ValueError('bad mapping'), ready]):
                rows = run(listing, config, root/'check', check=True)
            self.assertEqual([r['status'] for r in rows], ['failed', 'ready'])
            saved = json.loads((root/'check/status.json').read_text())
            self.assertIn('bad mapping', saved[0]['error'])

    def test_invalid_bootstrap_fails_before_allocation(self):
        with self.assertRaisesRegex(ValueError, 'two bootstrap'):
            calculate(np.zeros((3,4,2,10)), None, np.arange(2), np.arange(10), 'repository', 0, 0)


if __name__ == '__main__':
    unittest.main()
