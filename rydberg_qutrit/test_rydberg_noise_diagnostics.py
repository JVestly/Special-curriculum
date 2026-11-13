"""Focused tests for the state-resolved Rydberg noise extension."""
import unittest
import numpy as np
from rydberg_noise_diagnostics import detuning_samples, state_diagnostics

class RydbergNoiseDiagnosticsTests(unittest.TestCase):
    def test_gaussian_t2star_quadrature(self):
        samples=detuning_samples(3.0,order=9)
        delta=np.array([x for x,_ in samples])
        weights=np.array([w for _,w in samples])
        self.assertAlmostEqual(weights.sum(),1.,places=12)
        self.assertAlmostEqual(np.dot(weights,delta),0.,places=12)
        self.assertAlmostEqual(np.dot(weights,delta**2),2./9.,places=12)
        t=.4
        self.assertAlmostEqual(
            np.real(np.dot(weights,np.exp(1j*delta*t))),
            np.exp(-(t/3.)**2),places=10
        )

    def test_state_dependent_loss(self):
        results=state_diagnostics(("00","01","10","11","rr"))
        self.assertAlmostEqual(results["00"]["any_L"],0.,places=9)
        self.assertAlmostEqual(results["01"]["any_L"],results["10"]["any_L"],places=8)
        self.assertGreater(results["rr"]["any_L"],results["11"]["any_L"])
        for values in results.values():
            self.assertAlmostEqual(values["trace"],1.,places=7)
            self.assertAlmostEqual(
                values["computational"]+values["residual_rydberg"]+values["any_L"],
                1.,places=7
            )

if __name__=="__main__":
    unittest.main()
