import tempfile
import unittest
from pathlib import Path
import numpy as np
from coupling.sniffing import detect_peaks, trials_from_csv, process


class SniffingTests(unittest.TestCase):
    def test_peak_detection(self):
        t=np.arange(3000)/1000.
        x=np.sin(2*np.pi*5*t)
        positive, negative=detect_peaks(x)
        self.assertGreaterEqual(len(positive),12)
        self.assertGreaterEqual(len(negative),12)
        self.assertTrue(np.all(x[positive] > 0))

    def test_csv_parser_and_reject_mismatch(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d)
            frame=root/"trials.csv"
            frame.write_text("start_taste,end_taste,taste,dig_in_num_taste,abs_trial_num\n3000,3060,oBA,4,0\n")
            rows=trials_from_csv(frame)
            self.assertEqual(rows[0]["stimulus"],"oBA")
            adc=root/"adc.dat"
            np.zeros(2000,dtype="<u2").tofile(adc)
            with self.assertRaises(ValueError):
                process(adc,frame,root/"out",raw_fs=1000)
