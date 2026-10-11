import json
import struct
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch
import numpy as np
from scipy.signal import stft
from coupling.events import write_csv
from coupling.respiration import audit, detect_cycles
from coupling.respiration_lfp import run, event_key, select_event, bounds, validate_inputs, referenced_epoch
from coupling.power import normalize_power


class AlignmentTests(unittest.TestCase):
    def test_polarity_missing_events_and_boundaries(self):
        t=np.arange(10000)/1000
        cycles=detect_cycles(-np.sin(2*np.pi*4*t))
        self.assertEqual(event_key('inspiration','positive'),'negative_start_sample')
        self.assertEqual(event_key('expiration','negative'),'positive_start_sample')
        with self.assertRaises(ValueError):event_key('inspiration','unknown')
        self.assertIsNone(select_event([], 'negative_start_sample', 2.,0.,1.))
        k=select_event(cycles,'negative_start_sample',2.1,0.,.2)
        self.assertAlmostEqual(k/1000,2.25,delta=.002)
        self.assertIsNone(select_event(cycles,'negative_start_sample',2.1,0.,.1))
        self.assertFalse(bounds(1999,10000));self.assertTrue(bounds(2000,10000))
        self.assertFalse(bounds(5001,10000))

    def test_shared_baseline_not_renormalized(self):
        t=np.arange(-2000,5000)/1000
        rng=np.random.default_rng(42);x=rng.normal(size=len(t));y=4*x
        f,tt,z,valid=normalize_power(y,t,baseline_reference=x)
        _,_,a=stft(x,fs=1000,window='hann',nperseg=500,noverlap=499,boundary=None,padded=False,detrend=False)
        p=np.abs(a[:len(f)])**2;b=(tt>=-1.5-1e-9)&(tt<-.5-1e-9)
        expected=(16*p-p[:,b].mean(1)[:,None])/p[:,b].std(1)[:,None]
        np.testing.assert_allclose(z,expected,atol=1e-10)

    def test_shifted_car_donor_failure_is_group_wide(self):
        with tempfile.TemporaryDirectory() as td:
            source=Path(td);rng=np.random.default_rng(12)
            channels=[dict(filename=f'amp-{i}.dat',electrode_id=i,region='GC' if i<3 else 'PC',reference_group='GC' if i<3 else 'PC') for i in range(6)]
            continuous=rng.normal(size=(6,10000))
            for ch in channels:
                raw=rng.integers(-1000,1000,10000,dtype=np.int16)
                if ch['electrode_id']==0:raw[7900:8000]=0
                raw.tofile(source/ch['filename'])
            x,good=referenced_epoch(continuous,3000,channels,{'GC':[0,1,2],'PC':[3,4,5]},'mean_car',source,1000)
            self.assertFalse(good[:3].any());self.assertTrue(good[3:].all())
            np.testing.assert_allclose(x[3:].sum(0),0,atol=1e-14)

    def test_cached_session_end_to_end(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td);source=root/'raw';source.mkdir();prep=root/'prepared';prep.mkdir()
            fs=3000;n=30*fs;t=np.arange(n)/fs
            (source/'info.rhd').write_bytes(struct.pack('<Ihhf',0xc6912702,2,0,fs))
            np.arange(n,dtype='<i4').tofile(source/'time.dat')
            np.rint(37000+600*np.sin(2*np.pi*4*t)).astype('<u2').tofile(source/'board-ADC-00.dat')
            labels=['BA','S','oIA','IA','oBA','BK']
            for din in range(8):
                x=np.zeros(n,dtype='<u2')
                for sec in ([12.1] if din<6 else [10.,12.,14.]):
                    start=round(sec*fs);x[start:start+30]=1
                x.tofile(source/f'board-DIN-{din:02d}.dat')
            config=dict(raw_fs_hz=fs,session_id='synthetic',animal_id='synthetic',stimulus_by_din=dict(zip(map(str,range(6)),labels)))
            cfg=root/'session.json';cfg.write_text(json.dumps(config))
            with patch('coupling.respiration.plot_audit'):
                audit(source,cfg,root/'audit')
            resp=json.loads((root/'audit/respiration_provenance.json').read_text())
            channels=[dict(filename=f'amp-{r}.dat',electrode_id=i,region=r,reference_group=r) for i,r in enumerate(['GC','PC'])]
            rng=np.random.default_rng(13)
            continuous=rng.normal(size=(2,30000)).astype('float32')
            for ch in channels:rng.integers(-1000,1000,n,dtype=np.int16).tofile(source/ch['filename'])
            np.save(prep/'lfp_unreferenced.npy',continuous)
            trials=np.stack([continuous[:,10100:17100]]*6,axis=1)
            np.save(prep/'lfp_trials.npy',trials)
            events=[dict(event_id=i,din=i,stimulus=labels[i],original_trial_index=0,onset_sample=36300,onset_s=12.1,center_lfp_sample=12100) for i in range(6)]
            write_csv(prep/'events.csv',events);write_csv(prep/'channels.csv',channels)
            write_csv(prep/'channel_qc.csv',[dict(electrode_id=i,good_channel=True) for i in range(2)])
            write_csv(prep/'trial_qc.csv',[dict(electrode_id=i,event_id=j,accepted=not(i==0 and j==1)) for i in range(2) for j in range(6)])
            meta=dict(source=str(source),raw_fs_hz=fs,first_timestamp=0,configuration=config,filtering=dict(fs_hz=1000),reference=dict(method='none',members={'GC':[0],'PC':[1]}),files=resp['files'])
            (prep/'preprocessing.json').write_text(json.dumps(meta))
            rules=dict(definition_verified=False,detector=dict(bandpass_hz=[.5,20],minimum_hz=1.,maximum_hz=15.,prominence_mad=.5),stimuli={s:dict(route='synthetic',landmark='positive',delay_s=0.,window_s=2.) for s in labels})
            path=root/'rules.json';path.write_text(json.dumps(rules))
            out=run(prep,root/'audit',path)
            self.assertEqual(out.name,'LFPs matched to brathing')
            result=json.loads((out/'provenance.json').read_text())
            self.assertTrue(result['provisional']);self.assertEqual(result['retained_trials'],6)
            self.assertEqual(len(list(out.glob('*paired_power.png'))),12)
            with np.load(out/'paired_power.npz') as d:
                self.assertTrue(np.isnan(d['z_power'][:,0,2]).all()) # S excluded in GC only
                self.assertEqual(d['z_power'].shape[:3],(2,2,6))
            with self.assertRaises(FileExistsError):run(prep,root/'audit',path)
            # Real shifted-tail artifact must fail new QC despite old DIN acceptance.
            raw=np.fromfile(source/'amp-GC.dat',dtype='<i2')
            raw[51310:51850]=0;raw.tofile(source/'amp-GC.dat')
            _,good=referenced_epoch(continuous,12300,channels,meta['reference']['members'],'none',source,fs)
            self.assertFalse(good[0])
            resp_path=root/'audit/respiration_provenance.json';resp['first_timestamp']=1;resp_path.write_text(json.dumps(resp))
            with self.assertRaisesRegex(ValueError,'clock'):validate_inputs(prep,root/'audit',rules)


if __name__=='__main__':unittest.main()
