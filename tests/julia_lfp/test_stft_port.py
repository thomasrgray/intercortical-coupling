import numpy as np
from scipy.signal import detrend, stft
from coupling.julia_lfp.ephys_data import ephys_data
from coupling.julia_lfp.lfp_common import compute_stft


def test_packaged_stft_matches_source_formula(monkeypatch):
    monkeypatch.setenv('COUPLING_LFP_JOBS', '1')
    rng = np.random.default_rng(14)
    traces = rng.normal(size=(2, 2, 3000))
    dat = object.__new__(ephys_data)
    dat.lfp_array_list = [traces.copy()]
    dat.check_file_type = lambda: True
    compute_stft(dat, {'fs': 1000, 'pre_stim_s': 1, 'total_s': 3})
    f, t, z = stft(detrend(traces, axis=-1), fs=1000, window='hann', nperseg=500, noverlap=499)
    keep_f, keep_t = f < 40, (t >= 0) & (t < 3)
    expected = np.abs(z[..., keep_f, :][..., keep_t])**2
    np.testing.assert_array_equal(dat.freq_vec, f[keep_f])
    np.testing.assert_array_equal(dat.time_vec, t[keep_t])
    np.testing.assert_allclose(dat.amplitude_array_list[0], expected, rtol=1e-12, atol=1e-14)
