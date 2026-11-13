"""Checks for the reduced, adiabatically eliminated intermediate-state model."""
import unittest
import numpy as np
from rydberg_scattering import (
    ScatteringParameters, effective_local_jump, integrate_channel,
    analyze, averaged_channel
)
from rydberg_full_channel import construct_channel

class ScatteringTests(unittest.TestCase):
    def test_jump_coefficients_from_published_laser_parameters(self):
        p=ScatteringParameters()
        k=effective_local_jump(0,p)
        self.assertEqual(k.shape,(4,4))
        self.assertAlmostEqual(k[1,1].real,
            np.sqrt(1/p.intermediate_lifetime_us)*p.omega_blue_mhz/(2*p.intermediate_detuning_mhz))
        self.assertAlmostEqual(k[1,2].real,
            np.sqrt(1/p.intermediate_lifetime_us)*p.omega_red_mhz/(2*p.intermediate_detuning_mhz))

    def test_dark_state_cancelation(self):
        p=ScatteringParameters()
        v=np.array([0,-p.omega_red_mhz,p.omega_blue_mhz,0],complex)
        v/=np.linalg.norm(v)
        self.assertLess(np.linalg.norm(effective_local_jump(0,p)@v),1e-14)

    def test_no_scattering_reproduces_prior_channel(self):
        a=integrate_channel(include_scattering=False)
        b=construct_channel()
        self.assertLess(np.max(abs(a-b)),2e-7)

    def test_scattering_model_cptp_and_choi(self):
        data=analyze(integrate_channel())
        self.assertLess(abs(data['normalization']-1),1e-8)
        self.assertLess(data['trace_error'],2e-8)
        self.assertLess(data['reconstruction_error'],2e-8)
        self.assertGreater(data['min_choi_eig'],-1e-8)
        self.assertGreater(data['nonidentity'],0)

    def test_scattering_scale_zero_reproduces_no_scattering(self):
        a=integrate_channel(scattering=ScatteringParameters(scattering_scale=0))
        b=integrate_channel(include_scattering=False)
        self.assertLess(np.max(abs(a-b)),1e-11)

    def test_branching_target_changes_channel(self):
        r1=analyze(integrate_channel(scattering=ScatteringParameters(target='1')))
        rL=analyze(integrate_channel(scattering=ScatteringParameters(target='L')))
        self.assertGreater(rL['outcomes']['any_loss'],r1['outcomes']['any_loss'])

if __name__=='__main__':
    unittest.main()
