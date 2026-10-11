import tempfile
from pathlib import Path
import unittest
import numpy as np
from coupling.events import pulses,select_six
from coupling.prepare import resample_channel,raw_quality
from coupling.qc import metrics,peer_flags
from scipy.signal import resample_poly

class EventTests(unittest.TestCase):
    def test_chunk_boundaries_and_exact_onsets(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'din.dat';x=np.zeros(40,dtype='<u2');x[6:15]=1;x[21:27]=1;x.tofile(p)
            on,off,n=pulses(p,chunk=7)
            np.testing.assert_array_equal(on,[6,21]);np.testing.assert_array_equal(off,[15,27])
            x[0]=1;x.tofile(p)
            with self.assertRaises(ValueError):pulses(p)
    def test_selection_not_smallest_or_largest(self):
        self.assertEqual(select_six({**dict.fromkeys(range(6),30),6:90,7:90,8:0}),list(range(6)))
        with self.assertRaises(ValueError):select_six(dict.fromkeys(range(8),30))
        with self.assertRaises(ValueError):select_six({**dict.fromkeys(range(6),30),**dict.fromkeys(range(6,12),40)})
        with self.assertRaises(ValueError):select_six({**dict.fromkeys(range(5),30),5:29})

class SignalTests(unittest.TestCase):
    def test_resample_chunks_equal_full_and_alias_rejection(self):
        fs=3000;t=np.arange(fs*32)/fs
        x=(1000*np.sin(2*np.pi*20*t)+1000*np.sin(2*np.pi*980*t)).astype('<i2')
        y=resample_channel(x,fs);ref=resample_poly(x.astype(float)*.195,1,3)
        np.testing.assert_allclose(y,ref,atol=1e-9)
        expected=195*np.sin(2*np.pi*20*np.arange(len(y))/1000)
        self.assertLess(np.std(y[100:-100]-expected[100:-100]),1)
    def test_qc_flat_outlier_and_zero_mad(self):
        self.assertTrue(metrics(np.zeros(1000))['bad'])
        x=np.random.default_rng(5).normal(size=10000);x[:200]=100
        self.assertGreater(metrics(x)['outlier_fraction'],.01)
        self.assertTrue(metrics(x)['bad'])
    def test_exact_raw_histogram_and_flat_runs(self):
        x=np.random.default_rng(5).integers(-100,101,1_000_200,dtype=np.int16)
        x[999950:1000050]=10
        q=raw_quality(x,1000)
        self.assertAlmostEqual(q['median_uv'],np.median(x)*.195)
        self.assertAlmostEqual(q['sd_uv'],np.std(x)*.195)
        self.assertGreaterEqual(q['flat_fraction'],100/len(x))
    def test_peer_review_and_zero_scale(self):
        rows=[dict(bad=False,median_uv=0,sd_uv=1,bad_trial_fraction=0) for _ in range(5)]
        rows[-1]['sd_uv']=2
        flags=peer_flags(rows,['GC']*5)
        self.assertIn('peer_sd_uv',flags[-1]);self.assertEqual(flags[0],[])
        self.assertIn('insufficient_same_region_peers',peer_flags(rows[:2],['GC']*2)[0])

if __name__=='__main__':unittest.main()
