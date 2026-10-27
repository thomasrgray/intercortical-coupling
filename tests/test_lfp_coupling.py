
from coupling.array_store import load_arrays, save_arrays
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
import numpy as np
from coupling.lfp_coupling import pooled_coherence, spectral_granger, fit_granger, matched_contrasts, analyze
from coupling.power import STIMULI


def simulate(seed=1, trials=20, samples=1000, coupling=.6):
    rng=np.random.default_rng(seed)
    x=np.zeros((trials,samples+300,2))
    noise=rng.normal(size=x.shape)
    for t in range(1,x.shape[1]):
        x[:,t,0]=.6*x[:,t-1,0]+noise[:,t,0]
        x[:,t,1]=.4*x[:,t-1,1]+coupling*x[:,t-1,0]+noise[:,t,1]
    return x[:,300:]


class CouplingTests(unittest.TestCase):
    def test_known_direction_and_channel_swap(self):
        x=simulate()
        a,s,d=fit_granger(x,orders=(1,2,5))
        self.assertTrue(d['passed'],d)
        g=spectral_granger(a,s,np.arange(126),250)
        self.assertGreater(g[0].mean(),10*g[1].mean())
        a,s,d=fit_granger(x[:,:,::-1],orders=(1,2,5))
        np.testing.assert_allclose(spectral_granger(a,s,np.arange(126),250),g[::-1],atol=1e-10)

    def test_instantaneous_covariance_is_not_lagged_direction(self):
        g=spectral_granger(np.zeros((1,2,2)),np.array([[1,.8],[.8,1]]),np.arange(126),250)
        np.testing.assert_allclose(g,0,atol=1e-12)

    def test_null_and_boundary_offsets(self):
        x=simulate(coupling=0)
        offsets=np.random.default_rng(8).normal(size=(20,1,2))*1000
        a,s,d=fit_granger(x+offsets,orders=(1,2,5))
        g=spectral_granger(a,s,np.arange(126),250)
        self.assertLess(g.max(),.01)
        a2,s2,_=fit_granger(x,orders=(1,2,5))
        np.testing.assert_allclose(g,spectral_granger(a2,s2,np.arange(126),250),atol=1e-9)

    def test_coherence(self):
        rng=np.random.default_rng(12)
        x=rng.normal(size=(30,1500));y=rng.normal(size=x.shape)
        _,same=pooled_coherence(x,x);_,null=pooled_coherence(x,y)
        np.testing.assert_allclose(same,1,atol=1e-12)
        self.assertLess(null.mean(),.03)

    def test_route_difference_uses_common_pairs(self):
        x=np.array([[[1.],[3.]],[[100.],[np.nan]]])
        delta,n=matched_contrasts(x,['BA','oBA'],[['BA','oBA']])
        self.assertEqual(delta.item(),2);self.assertEqual(n.item(),1)

    def test_underfit_residuals_are_rejected(self):
        rng=np.random.default_rng(19)
        x=rng.normal(size=(20,1000,2))
        for k in range(10,1000): x[:,k] += .7*x[:,k-10]
        _,_,d=fit_granger(x,orders=(1,))
        self.assertFalse(d['passed'])
        self.assertIn('residual_autocorrelation',d['reasons'])

    def test_constant_rejected(self):
        with self.assertRaises(ValueError):fit_granger(np.zeros((8,1000,2)))

    def test_pipeline_masks_missing_condition_and_ids(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp);t=np.arange(-2000,2000)/1000
            rng=np.random.default_rng(10)
            x=rng.normal(size=(3,18,len(t)))
            good=np.ones((3,18),bool);good[1,:3]=False
            x[1,:3]=np.nan;x[2]=np.nan
            save_arrays(p/'input.h5',lfp_uv=x,time_s=t,electrode_id=[2,7,99],region=['GC','PC','GC'],
                good_channel=np.array([True,True,False]),good_trial=good,stimulus=np.repeat(STIMULI,3),
                event_id=np.arange(18)+100,animal_id='SYNTHETIC',session_id='SYNTHETIC',
                provenance_json=json.dumps(dict(source='synthetic',reference='none',filtering='none',qc_method='test',event_alignment='TTL')))
            with patch('coupling.lfp_coupling.save_fig'):
                analyze(p/'input.h5',p/'out',min_trials=3,orders=(1,2))
            import matplotlib.pyplot as plt
            plt.close('all')
            with load_arrays(p/'out/plotted_values.h5') as d:
                self.assertTrue(np.isnan(d['coherence'][:,0]).all())
                self.assertTrue(np.isnan(d['coherence_route_difference'][0]).all())
                self.assertEqual(d['electrode_pairs'].tolist(),[[2,7]])
                self.assertEqual(d['event_id'][0],100)
                self.assertTrue(np.isfinite(d['coherence'][:,1:]).all())
            with self.assertRaises(FileExistsError):analyze(p/'input.h5',p/'out')

if __name__=='__main__':unittest.main()
