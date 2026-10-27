
from coupling.array_store import load_arrays, save_arrays
import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch
import numpy as np
from coupling.compare import compare, select_panel


class ComparisonTests(unittest.TestCase):
    def test_label_selection_and_difference_preserve_missing_values(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            source = root/'input.h5'
            # Reordered labels catch accidental positional indexing.
            z = np.zeros((2, 2, 2, 3))
            z[1, 1] = 4
            z[1, 0] = 1
            z[1, 0, 0, 0] = np.nan
            save_arrays(source, z_power=z, n_channels=np.ones_like(z, dtype=int)*7,
                     region=['PC', 'GC'], stimulus=['oIA', 'IA'],
                     frequency_hz=[2, 10], time_s=[-1, 0, 1.75],
                     session_id='test', animal_id='rat',
                     settings_json=json.dumps({'baseline_s': [-1.5, -.5]}))
            with patch('matplotlib.figure.Figure.savefig'):
                out = compare(source, root/'out', difference=True)
            with load_arrays(out/'comparison_values.h5') as data:
                self.assertTrue(np.isnan(data['difference_z'][0, 0]))
                np.testing.assert_array_equal(data['difference_z'][1], [3, 3, 3])
                np.testing.assert_array_equal(data['left_z'], np.full((2, 3), 4))
            with self.assertRaises(FileExistsError):
                compare(source, out)
            with self.assertRaises(ValueError):
                compare(source, root/'bad', left='GC:missing')
            self.assertFalse((root/'bad').exists())
            with self.assertRaises(ValueError):
                compare(source, root/'limit', clim=0)

    def test_selector_validation(self):
        with self.assertRaises(ValueError):
            select_panel({}, 'GC_IA')


if __name__ == '__main__':
    unittest.main()
