import json
import tempfile
import unittest
from pathlib import Path
import h5py
import numpy as np
from coupling.events import write_csv
from coupling.figure2 import load_sorted, alpha_filter, run, unit_regions, amplitude_envelope


class Figure2Tests(unittest.TestCase):
    def test_region_key_variants_and_explicit_prepared_mapping(self):
        self.assertEqual(unit_regions({'GC': [[1]], 'PC': [[2]]}), {'GC':[1], 'PC':[2]})
        self.assertEqual(unit_regions({' gc ': [1], 'Pc': [2]}), {'GC':[1], 'PC':[2]})
        with self.assertRaisesRegex(ValueError, 'unit-regions-from-prepared'):
            unit_regions({'0': [1], '1': [2]})
        with self.assertRaisesRegex(ValueError, 'Cannot identify'):
            unit_regions({'GC':[1], 'gc':[2], 'PC':[3]})
        self.assertEqual(unit_regions({'0':[1,2]}, {'GC':[1], 'PC':[2]}), {'GC':[1], 'PC':[2]})

    def test_sorted_units_with_prepared_mapping(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d); ev,t=self.fixture(root)
            (root/'session.info').write_text(json.dumps({'electrode_layout':{'0':[[1,2]]}}))
            x,ids,regions,_,_=load_sorted(root,ev,30000,t,prepared_regions={'GC':[1], 'PC':[2]})
            np.testing.assert_array_equal(regions,['GC','PC'])
            self.assertEqual(x[0,0,2000],2)

    def fixture(self, root):
        (root/'session.info').write_text(json.dumps({'electrode_layout':{'gc':[[1]],'pc':[[2]]}}))
        with h5py.File(root/'session.h5','w') as h:
            h['unit_descriptor'] = np.array([(1,1),(2,1),(1,0)],dtype=[('electrode_number','i4'),('single_unit','i4')])
            for u in range(3):
                # Two counts in the first millisecond after zero, plus a pre-onset spike.
                h[f'sorted_units/unit{u:03d}/times'] = [299970,300000,300010,600000]
        events=[dict(event_id=0,din=0,stimulus='BA',original_trial_index=0,onset_sample=300001,center_lfp_sample=10000)]
        return events, np.arange(-2000,5000)/1000

    def test_absolute_clock_counts_and_units(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d); ev,t=self.fixture(root)
            x,ids,regions,_,_=load_sorted(root,ev,30000,t)
            self.assertEqual(x.shape,(2,1,7000))
            np.testing.assert_array_equal(ids,[0,1])
            np.testing.assert_array_equal(regions,['GC','PC'])
            self.assertEqual(x[0,0,2000],2)
            self.assertEqual(x[0,0,1999],1)

    def test_reject_ambiguous_h5_and_region_overlap(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d); ev,t=self.fixture(root)
            (root/'extra.h5').touch()
            with self.assertRaisesRegex(ValueError,'expected one'):
                load_sorted(root,ev,30000,t)
            (root/'session.info').write_text(json.dumps({'electrode_layout':{'gc':[1],'pc':[1]}}))
            with self.assertRaisesRegex(ValueError,'Overlapping'):
                load_sorted(root,ev,30000,t,h5='session.h5')

    def test_envelope_preserves_sinusoidal_amplitude(self):
        t=np.arange(7000)/1000
        np.testing.assert_allclose(amplitude_envelope(3*np.sin(2*np.pi*10*t)),3,atol=1e-10)

    def test_alpha_passband(self):
        t=np.arange(7000)/1000
        x=np.sin(2*np.pi*10*t)+np.sin(2*np.pi*30*t)
        y=alpha_filter(x)[1000:-1000]; tm=t[1000:-1000]
        self.assertGreater(abs(np.mean(y*np.exp(-2j*np.pi*10*tm))),.45)
        self.assertLess(abs(np.mean(y*np.exp(-2j*np.pi*30*tm))),.01)

    def test_real_cache_end_to_end_and_clock_rejection(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d); ev,t=self.fixture(root)
            prep=root/'prep'; audit=root/'audit';prep.mkdir();audit.mkdir()
            meta=dict(source=d,raw_fs_hz=30000,first_timestamp=0,configuration={},files=[])
            (prep/'preprocessing.json').write_text(json.dumps(meta))
            (audit/'respiration_provenance.json').write_text(json.dumps(dict(meta,sign='unknown')))
            write_csv(prep/'events.csv',ev)
            re=[dict(event_id=0,din=0,stimulus='BA',original_trial_index=0,ttl_raw_sample=300001,ttl_lfp_sample=10000,timing_qc_reasons='')]
            write_csv(audit/'events_respiration.csv',re)
            pressure=np.sin(2*np.pi*4*np.arange(20000)/1000)
            np.savez(audit/'pressure_1khz.npz',fs_hz=1000,unfiltered_counts=pressure,filtered_counts=pressure)
            lfp=np.stack([np.sin(2*np.pi*10*t),np.cos(2*np.pi*10*t)])[:,None,:]
            np.savez(prep/'cleaned_lfp.npz',provenance_json=json.dumps(meta),lfp_uv=lfp,time_s=t,event_id=[0],stimulus=['BA'],region=['GC','PC'],good_channel=[True,True],good_trial=np.ones((2,1),bool),electrode_id=[1,2])
            run(prep,audit,root/'check',check=True)
            self.assertFalse((root/'check').exists())
            out=run(prep,audit,root/'out')
            self.assertTrue((out/'01_figure2.pdf').exists())
            self.assertTrue((out/'BA_ortho_vs_retro.pdf').exists())
            self.assertTrue((out/'IA_ortho_vs_retro.png').exists())
            with np.load(out/'plotted_data.npz') as z:
                self.assertEqual(z['BA__spike_counts'][0,0,2000],2)
                self.assertEqual(z['BA__prepared_lfp_uv'].shape,(2,1,7000))
                self.assertEqual(z['BA__raw_pressure_counts'].shape,(1,7000))
                np.testing.assert_allclose(z['BA__alpha_envelope_uv'],amplitude_envelope(z['BA__alpha_uv']))
                self.assertAlmostEqual(z['BA__GC_rate_hz'][0,40],40)
            re[0]['ttl_raw_sample']=300002
            write_csv(audit/'events_respiration.csv',re)
            with self.assertRaisesRegex(ValueError,'onset mismatch'):
                run(prep,audit,root/'bad')
