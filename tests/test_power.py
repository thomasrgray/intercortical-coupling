import unittest
import numpy as np
from coupling.power import normalize_power

class PowerTests(unittest.TestCase):
    def setUp(self):
        self.t=np.arange(-2000,2500)/1000
        self.x=np.random.default_rng(4).normal(size=len(self.t))
    def test_baseline_and_axes(self):
        f,t,z,valid=normalize_power(self.x,self.t)
        b=(t>=-1.5-1e-9)&(t<-.5-1e-9)
        np.testing.assert_allclose(z[valid][:,b].mean(-1),0,atol=1e-12)
        np.testing.assert_allclose(z[valid][:,b].std(-1),1,atol=1e-12)
        np.testing.assert_allclose(np.diff(t),.001,atol=1e-12)
        np.testing.assert_array_equal(f,np.arange(0,71,2))
    def test_baseline_unaffected_by_later_response(self):
        f,t,a,v=normalize_power(self.x,self.t)
        changed=self.x.copy();changed[self.t>0]*=10
        _,_,b,_=normalize_power(changed,self.t)
        np.testing.assert_allclose(a[:,t<-.5],b[:,t<-.5])
    def test_constant_is_invalid(self):
        _,_,z,valid=normalize_power(np.zeros_like(self.x),self.t)
        self.assertFalse(valid.any());self.assertTrue(np.isnan(z).all())
    def test_bad_input_rejected(self):
        x=self.x.copy();x[0]=np.nan
        with self.assertRaises(ValueError): normalize_power(x,self.t)
        with self.assertRaises(ValueError): normalize_power(self.x,self.t*2)


class PipelineTests(unittest.TestCase):
    def test_masks_and_equal_channel_weighting(self):
        import tempfile,json
        from pathlib import Path
        from unittest.mock import patch
        from coupling.power import analyze,STIMULI
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp); t=np.arange(-2000,2500)/1000
            rng=np.random.default_rng(55)
            x=rng.normal(size=(3,7,len(t)))
            x[2]=np.nan  # excluded channel must never enter STFT
            x[:,6]=np.nan # excluded event must never enter STFT
            mask=np.ones((3,7),bool);mask[:,6]=False
            np.savez(p/'in.npz',lfp_uv=x,time_s=t,stimulus=np.array([*STIMULI,'BA']),
                electrode_id=[12,59,13],region=['GC','PC','GC'],good_channel=[True,True,False],
                good_trial=mask,animal_id='test',session_id='test',provenance_json=json.dumps(
                dict(source='test',reference='test',filtering='test',qc_method='test',event_alignment='TTL')))
            with patch('coupling.power.plot_panel'):
                analyze(p/'in.npz',p/'out')
            with np.load(p/'out/session_power.npz') as out:
                self.assertEqual(out['z_power'].shape[:3],(2,6,36))
                self.assertTrue(np.isfinite(out['z_power']).all())
                self.assertTrue((out['n_channels']==1).all())
                f,tt,z,_=normalize_power(x[0,0],t)
                keep=(tt>=-1-1e-9)&(tt<=1.75+1e-9)
                np.testing.assert_allclose(out['z_power'][0,0],z[:,keep])
            with self.assertRaises(FileExistsError):analyze(p/'in.npz',p/'out')

if __name__=="__main__":unittest.main()
