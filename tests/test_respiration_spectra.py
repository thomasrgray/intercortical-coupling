import json
import tempfile
import unittest
from pathlib import Path
import numpy as np
from coupling.events import write_csv
from coupling.respiration_spectra import spectra, run, load_inputs


class RespirationSpectraTests(unittest.TestCase):
    def test_shared_frequency_and_independent_control(self):
        rng=np.random.default_rng(7);t=np.arange(7000)/1000
        phases=rng.uniform(0,2*np.pi,(80,1))
        x=np.sin(2*np.pi*5*t+phases)+rng.normal(0,.2,(80,len(t)))
        y=np.sin(2*np.pi*5*t+phases+.5)+rng.normal(0,.2,x.shape)
        f,_,p,_,c=spectra(x,y)
        _,_,_,_,null=spectra(x,rng.normal(size=x.shape))
        k=np.argmin(abs(f-5))
        self.assertEqual(f[np.argmax(p.mean(axis=0))],5)
        self.assertGreater(c[:,k].mean(),.95)
        self.assertLess(null[:,k].mean(),.1)
        self.assertTrue(np.all((c>=0)&(c<=1)))

    def test_reject_single_trial(self):
        with self.assertRaises(ValueError):spectra(np.ones((1,7000)),np.ones((1,7000)))

    def test_full_workflow_and_alignment_rejection(self):
        rng=np.random.default_rng(4)
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);prep=root/'prep';audit=root/'audit';prep.mkdir();audit.mkdir()
            meta=dict(source=d,raw_fs_hz=1000,first_timestamp=0,configuration={},files=[])
            (prep/'preprocessing.json').write_text(json.dumps(meta))
            (audit/'respiration_provenance.json').write_text(json.dumps(dict(meta,sign='unknown')))
            stimuli=['oBA']*4+['BA']*4+['oIA']*4+['IA']*4
            centers=np.arange(16)*8000+8000;t=np.arange(centers[-1]+8000)/1000
            pressure=np.sin(2*np.pi*5*t)+rng.normal(0,.1,len(t))
            np.savez(audit/'pressure_1khz.npz',fs_hz=1000,unfiltered_counts=pressure,filtered_counts=pressure)
            ev=[dict(event_id=j,din=j//4,stimulus=s,original_trial_index=j%4,onset_sample=int(k),center_lfp_sample=int(k)) for j,(s,k) in enumerate(zip(stimuli,centers))]
            re=[dict(event_id=e['event_id'],din=e['din'],stimulus=e['stimulus'],original_trial_index=e['original_trial_index'],ttl_raw_sample=e['onset_sample'],ttl_lfp_sample=e['center_lfp_sample'],timing_qc_reasons='') for e in ev]
            write_csv(prep/'events.csv',ev);write_csv(audit/'events_respiration.csv',re)
            epochs=np.stack([pressure[k-2000:k+5000] for k in centers])
            lfp=np.stack([epochs+rng.normal(0,.2,epochs.shape) for _ in range(2)])
            np.savez(prep/'cleaned_lfp.npz',provenance_json=json.dumps(meta),lfp_uv=lfp,time_s=np.arange(-2000,5000)/1000,event_id=np.arange(16),stimulus=stimuli,region=['GC','PC'],good_channel=[True,True],good_trial=np.ones((2,16),bool))
            out=run(prep,audit,root/'out')
            self.assertTrue((out/'BA_respiration_lfp.png').exists())
            with np.load(out/'spectra.npz') as z:self.assertEqual(z['oBA__GC'].shape[0],501)
            re[0]['ttl_raw_sample']+=1;write_csv(audit/'events_respiration.csv',re)
            with self.assertRaisesRegex(ValueError,'onset mismatch'):load_inputs(prep,audit)
