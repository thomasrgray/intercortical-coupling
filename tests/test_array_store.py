import ast
from pathlib import Path
import h5py
import numpy as np
import pytest
from coupling.array_store import load_arrays, save_arrays, is_array_store


def test_roundtrip_types_shapes_and_native_hdf5(tmp_path):
    values = dict(real=np.array([np.nan, np.inf, -2.], dtype=np.float32),
                  complex=np.array([1+2j]), integer=np.arange(6).reshape(2, 3),
                  boolean=np.array([[True, False]]), scalar=np.asarray(4.),
                  label=np.asarray('GC μV'), labels=np.array(['BA', 'oIA']),
                  empty=np.array([], dtype='U5'), empty_numeric=np.empty((0, 2)))
    path = tmp_path/'arrays.h5'
    save_arrays(path, **values)
    assert is_array_store(path)
    with h5py.File(path, 'r') as h:
        assert h['integer'].shape == (2, 3)
    with load_arrays(path) as actual:
        assert set(actual) == set(values)
        for name, value in values.items():
            assert actual[name].dtype == value.dtype
            np.testing.assert_array_equal(actual[name], value)


def test_never_overwrites_recordings_and_failed_write_is_atomic(tmp_path):
    raw = tmp_path/'recording.h5'
    with h5py.File(raw, 'w') as h:
        h.create_group('Parsed_LFP')
    before = raw.read_bytes()
    with pytest.raises(ValueError, match='Refusing'):
        save_arrays(raw, value=[1])
    assert raw.read_bytes() == before
    path = tmp_path/'derived.h5'
    save_arrays(path, value=[1])
    before = path.read_bytes()
    with pytest.raises(ValueError, match='Unsupported dtype'):
        save_arrays(path, objects=[{'bad': 1}])
    assert path.read_bytes() == before
    with pytest.raises(ValueError, match='HDF5'):
        load_arrays(tmp_path/'old.npz')


def test_analysis_code_does_not_read_or_write_zip_archives():
    import coupling
    root = Path(coupling.__file__).parent
    for path in root.rglob('*.py'):
        for node in ast.walk(ast.parse(path.read_text())):
            if isinstance(node, ast.Call) and isinstance(node.func, ast.Attribute):
                assert node.func.attr not in {'savez', 'savez_compressed'}, path
                if isinstance(node.func.value, ast.Name) and node.func.value.id == 'np' and node.func.attr == 'load':
                    # Only explicitly .npy caches and the Julia cleaned NPY loader remain.
                    assert '.npy' in ast.unparse(node) or path.name in {'lfp_common.py', 'session.py'}, path
