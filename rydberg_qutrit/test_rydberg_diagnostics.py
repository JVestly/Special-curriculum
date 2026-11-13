"""Focused tests for extended Rydberg noise diagnostics."""
import unittest
import numpy as np
from rydberg_diagnostics import detuning_quadrature,state_diagnostics,outcomes,integrate_batch,_initial_kets

class NoiseDiagnosticsTests(unittest.TestCase):
    def test_detuning_quadrature_recovers_gaussian(self):
        t2=3.0
        delta,w=detuning_quadrature(t2,order=9)
        self.assertAlmostEqual(float(np.sum(w)),1.0,places=13)
        self.assertAlmostEqual(float(np.sum(w*delta)),0.0,places=13)
        self.assertAlmostEqual(float(np.sum(w*delta**2)),2/t2**2,places=13)
        t=0.4
        coh=np.sum(w*np.exp(1j*delta*t))
        self.assertAlmostEqual(float(coh.real),np.exp(-(t/t2)**2),places=11)

    def test_initial_00_remains_00(self):
        k=_initial_kets(('00',))
        t,series=integrate_batch(np.array([np.outer(k[0],k[0].conj())]),grid_points=3)
        out=outcomes(series[-1,0]);self.assertAlmostEqual(out['computational'],1.0,places=10)
        self.assertAlmostEqual(out['any_loss'],0.0,places=10)
        self.assertAlmostEqual(out['rydberg'],0.0,places=10)

    def test_trace_and_state_symmetry(self):
        t,res=state_diagnostics(labels=('01','10','rr'),points=4)
        for k,values in res.items():
            self.assertTrue(np.allclose(values['trace'],1,atol=1e-9),msg=k)
            self.assertTrue(np.all(values['any_loss']>=-1e-10))
            self.assertTrue(np.all(values['any_loss']<=1+1e-10))
        self.assertAlmostEqual(res['01']['any_loss'][-1],res['10']['any_loss'][-1],places=9)
        self.assertGreater(res['rr']['any_loss'][-1],res['01']['any_loss'][-1])

if __name__=='__main__':unittest.main()
