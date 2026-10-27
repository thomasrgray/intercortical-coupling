"""Versioned HDF5 storage for derived arrays, separate from recording HDF5.

No pickle or ZIP-array archives are read. Strings round-trip as NumPy Unicode
arrays; numeric, complex and boolean dtypes/shapes are preserved. Writes are
atomic and never replace an unrecognized recording HDF5 file.
"""
from collections.abc import Mapping
import os
from pathlib import Path
import tempfile
import h5py
import numpy as np

SCHEMA = 'intercortical-coupling.arrays'
VERSION = 1


def is_array_store(path):
    path = Path(path)
    if path.suffix.lower() not in {'.h5', '.hdf5'}:
        raise ValueError('Expected HDF5 (.h5/.hdf5); rerun earlier stages to regenerate old outputs')
    with h5py.File(path, 'r') as handle:
        return handle.attrs.get('schema') == SCHEMA


class ArrayStore(Mapping):
    """Read-only array mapping. Open in a context manager to release the file."""
    def __init__(self, path):
        if not is_array_store(path):
            raise ValueError(f'{path} is not a derived array store; use the session HDF5 loader')
        self._handle = h5py.File(path, 'r')
        if self._handle.attrs.get('schema_version') != VERSION:
            self._handle.close()
            raise ValueError('Unsupported array-store schema version')
        self.files = list(self._handle.keys())

    def __getitem__(self, key):
        dataset = self._handle[key]
        if 'numpy_unicode_dtype' in dataset.attrs:
            return np.asarray(dataset.asstr()[()], dtype=dataset.attrs['numpy_unicode_dtype'])
        return np.asarray(dataset[()])

    def __iter__(self):
        return iter(self.files)

    def __len__(self):
        return len(self.files)

    def close(self):
        self._handle.close()

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.close()


def load_arrays(path, *, allow_pickle=False):
    if allow_pickle:
        raise ValueError('Pickled objects are not supported')
    return ArrayStore(path)


def save_arrays(path, **arrays):
    path = Path(path)
    if path.suffix.lower() not in {'.h5', '.hdf5'}:
        raise ValueError('Derived array stores require .h5 or .hdf5')
    if path.is_symlink():
        raise ValueError('Refusing to overwrite a symlink')
    if path.exists() and not is_array_store(path):
        raise ValueError(f'Refusing to replace recording/non-array HDF5: {path}')
    # Validate before creating anything; never serialize Python objects.
    values = {}
    for key, value in arrays.items():
        if not key or '/' in key or key in {'.', '..'}:
            raise ValueError(f'Invalid array key: {key!r}')
        a = np.asarray(value)
        if a.dtype.kind not in 'biufcUS':
            raise ValueError(f'Unsupported dtype for {key}: {a.dtype}')
        values[key] = a
    fd, temporary = tempfile.mkstemp(prefix='.'+path.name, suffix='.tmp', dir=path.parent)
    os.close(fd)
    try:
        with h5py.File(temporary, 'w') as handle:
            handle.attrs.update(schema=SCHEMA, schema_version=VERSION)
            for key, a in values.items():
                options = {'compression': 'gzip', 'shuffle': True} if a.ndim and a.size and a.dtype.kind != 'U' else {}
                if a.dtype.kind == 'U':
                    dataset = handle.create_dataset(key, data=a.astype(object), dtype=h5py.string_dtype('utf-8'))
                    dataset.attrs['numpy_unicode_dtype'] = a.dtype.str
                else:
                    handle.create_dataset(key, data=a, **options)
            handle.flush()
        os.replace(temporary, path)
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)
