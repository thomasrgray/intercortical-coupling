import unittest
import numpy as np
from coupling.adc_check import compare_segment

class ADCCheckTests(unittest.TestCase):
    def test_known_signal_and_alias_rejection(self):
        fs=3000
        t=np.arange(30*fs)/fs
        raw=np.rint(37000+500*np.sin(2*np.pi*4*t)+500*np.sin(2*np.pi*704*t)).astype('<u2')
        metrics,(_,a,b,f,p)=compare_segment(raw,fs)
        expected=500*np.sin(2*np.pi*4*np.arange(5000,25000)/1000)
        self.assertGreater(np.corrcoef(a,expected)[0,1],.999)
        self.assertGreater(metrics['correlation'],.999)
        self.assertLess(metrics['relative_rms_error'],.01)
        self.assertEqual(f[np.argmax(p)],4)

    def test_constant_and_invalid_input(self):
        metrics,_=compare_segment(np.full(90000,37000,dtype='<u2'),3000)
        self.assertIsNone(metrics['spectral_flatness'])
        with self.assertRaises(ValueError): compare_segment(np.zeros(100),3000)
        with self.assertRaises(ValueError): compare_segment(np.zeros(90000),2500)
