"""Smoke/physics checks for the reduced channel. Execute with unittest."""
import unittest
import numpy as np
from rydberg_qutrit import Evered2023, rydberg_hamiltonian_numpy
from rydberg_full_channel import (all_decay_operators,construct_channel,
    choi_from_channel,kraus_from_choi,channel_report,COMP)

class Tests(unittest.TestCase):
    def test_4level_hamiltonian(self):
        H=rydberg_hamiltonian_numpy([0,2],omega_r_mhz=4.6,
            c6_mhz_um6=Evered2023().inferred_c6_mhz_um6,include_loss_level=True)
        self.assertEqual(H.shape,(16,16))
        np.testing.assert_allclose(H,H.conj().T)
    def test_branches(self):
        ops=all_decay_operators()
        self.assertEqual(len(ops),6)
        G=sum(c.conj().T@c for c in ops.values())
        expected=np.diag([Evered2023().decay_rate_per_us*(int(i//4==2)+int(i%4==2))
                          for i in range(16)])
        np.testing.assert_allclose(G,expected,atol=1e-13)
    def test_full_channel(self):
        E=construct_channel()
        J=choi_from_channel(E)
        Ks=kraus_from_choi(J)
        self.assertEqual(J.shape,(64,64))
        np.testing.assert_allclose(sum(K.conj().T@K for K in Ks),np.eye(4),atol=1e-7)
        self.assertGreaterEqual(np.linalg.eigvalsh(J)[0],-1e-8)
        U=channel_report()
        self.assertAlmostEqual(U["p_total"],1,places=7)
        self.assertAlmostEqual(sum(U["pauli_weights"].values()),
                               U["p_computational"],places=7)
if __name__=="__main__": unittest.main()
