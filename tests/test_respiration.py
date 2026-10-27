
from coupling.array_store import load_arrays
import hashlib
import json
import struct
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch
import numpy as np
from scipy.signal import resample_poly
from coupling.respiration import detect_cycles, first_after, downsample_counts, audit


class RespirationTests(unittest.TestCase):
    def test_known_phase_timing_and_rate(self):
        t = np.arange(10000)/1000
        cycles = detect_cycles(-np.sin(2*np.pi*4*t))
        passed = [c for c in cycles if c['screen_pass']]
        self.assertGreater(len(passed), 30)
        for c in passed:
            self.assertAlmostEqual(c['cycle_hz'],4,delta=.03)
            self.assertAlmostEqual((c['positive_start_sample']-c['negative_start_sample'])/1000,.125,delta=.002)
        self.assertEqual(detect_cycles(np.zeros(1000)),[])
        with self.assertRaises(ValueError):
            detect_cycles([0,np.nan,1])

    def test_no_imputed_phase_and_search_limit(self):
        self.assertIsNone(first_after([],100))
        self.assertIsNone(first_after([3101],1000))
        self.assertEqual(first_after([100,3100],1100),3100)
        self.assertEqual(first_after([100,200],100),100)
        slow = -np.sin(2*np.pi*.5*np.arange(15000)/1000)
        self.assertFalse(any(c['screen_pass'] for c in detect_cycles(slow)))

    def test_adc_unsigned_chunk_grid(self):
        fs = 3000
        t = np.arange(33*fs)/fs
        raw = np.rint(37000+500*np.sin(2*np.pi*4*t)).astype('<u2')
        got = downsample_counts(raw,fs)
        expected = resample_poly(raw.astype(float)-32768,1,3)
        np.testing.assert_allclose(got,expected,atol=1e-8)

    def test_audit_polarity_identity_and_read_only(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder);session = root/'session';session.mkdir()
            fs = 3000;n=30*fs;t=np.arange(n)/fs
            (session/'info.rhd').write_bytes(struct.pack('<Ihhf',0xc6912702,2,0,fs))
            np.arange(n,dtype='<i4').tofile(session/'time.dat')
            raw=np.rint(37000-600*np.sin(2*np.pi*4*t)).astype('<u2')
            raw.tofile(session/'board-ADC-00.dat')
            for din in range(8):
                x=np.zeros(n,dtype='<u2')
                for center in ([12*fs] if din<6 else [10*fs,12*fs,14*fs]):
                    x[center:center+30]=1
                x.tofile(session/f'board-DIN-{din:02d}.dat')
            config=root/'config.json'
            config.write_text(json.dumps(dict(raw_fs_hz=fs,session_id='test',animal_id='rat',
                stimulus_by_din=dict(zip(map(str,range(6)),['BA','S','oIA','IA','oBA','BK'])))))
            before=hashlib.sha256((session/'board-ADC-00.dat').read_bytes()).hexdigest()
            with patch('coupling.respiration.plot_audit'):
                out=audit(session,config,root/'out')
                verified=audit(session,config,root/'verified',sign='negative',polarity_note='Synthetic waveform ground truth',waveform_review_note='Synthetic known 4 Hz sinusoid')
                inverted=audit(session,config,root/'inverted',sign='positive',polarity_note='Synthetic inverted polarity check')
            import csv
            with (out/'events_respiration.csv').open() as stream:
                events=list(csv.DictReader(stream))
            self.assertEqual([int(e['event_id']) for e in events],list(range(6)))
            self.assertTrue(all(e['inhalation_latency_s']=='' for e in events))
            meta=json.loads((out/'respiration_provenance.json').read_text())
            self.assertEqual(meta['physiological_timing_trials'],0)
            self.assertEqual(meta['candidate_screen_pass_trials'],6)
            self.assertEqual(json.loads((verified/'respiration_provenance.json').read_text())['physiological_timing_trials'],6)
            with load_arrays(verified/'pressure_1khz.h5') as a,load_arrays(inverted/'pressure_1khz.h5') as b:
                np.testing.assert_allclose(a['filtered_counts'],-b['filtered_counts'])
            self.assertEqual(json.loads((inverted/'respiration_provenance.json').read_text())['physiological_timing_trials'],0)
            self.assertEqual(before,hashlib.sha256((session/'board-ADC-00.dat').read_bytes()).hexdigest())
            with self.assertRaises(FileExistsError):audit(session,config,out)
            with self.assertRaises(ValueError):audit(session,config,session/'bad')
            with self.assertRaises(ValueError):audit(session,config,root/'nopolarity',sign='negative')
            raw[:-1].tofile(session/'board-ADC-00.dat')
            with self.assertRaises(ValueError):audit(session,config,root/'mismatch')


if __name__ == '__main__':
    unittest.main()
