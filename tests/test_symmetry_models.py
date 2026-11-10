"""Independent three-level/qutrit Hamiltonian symmetry checks.

Run from repository root:
    python -m unittest discover -s tests -p 'test_*.py' -v

Requires numpy. These tests do not use the SAT implementation and
do NOT establish a physical NA-to-LMG reduction.
"""
import itertools
import unittest

import numpy as np

LEVELS = ("0", "1", "r")
INDEX = {a: k for k, a in enumerate(LEVELS)}


def local(a, b, site, n):
    e = np.zeros((3, 3), dtype=complex)
    e[INDEX[a], INDEX[b]] = 1
    ops = [np.eye(3, dtype=complex) for _ in range(n)]
    ops[site] = e
    result = ops[0]
    for op in ops[1:]:
        result = np.kron(result, op)
    return result


def T(a, b, n):
    zero = np.zeros((3**n, 3**n), dtype=complex)
    return sum((local(a, b, i, n) for i in range(n)), zero)


def parity(a, n):
    op = np.eye(3)
    op[INDEX[a], INDEX[a]] = -1
    result = op
    for _ in range(n - 1):
        result = np.kron(result, op)
    return result


def H_lmg(n, eps=(0.1, 0.4, 0.9), lam=0.75):
    if n < 2:
        raise ValueError("N must be >=2 for this LMG normalization")
    h = (-eps[0] * T("0", "0", n) + eps[1] * T("1", "1", n)
         + eps[2] * T("r", "r", n)) / n
    for a, b in itertools.permutations(LEVELS, 2):
        tab = T(a, b, n)
        h -= lam / (n * (n - 1)) * (tab @ tab)
    return h


def H_na(n, omega_hf=1.1, omega_r=1.7, v=0.5):
    h = omega_hf/2 * (T("0", "1", n) + T("1", "0", n))
    h += omega_r/2 * (T("1", "r", n) + T("r", "1", n))
    h += 0.4 * T("1", "1", n) + 0.6 * T("r", "r", n)
    for i in range(n):
        for j in range(i + 1, n):
            h += v * local("r", "r", i, n) @ local("r", "r", j, n)
    return h


def comm_norm(a, b):
    return float(np.linalg.norm(a @ b - b @ a))


class TestPhysicalQutritSymmetries(unittest.TestCase):
    def test_lmg_population_parities_commute(self):
        for n in (2, 3):
            h = H_lmg(n)
            self.assertTrue(np.allclose(h, h.conj().T))
            for level in LEVELS:
                self.assertLess(comm_norm(h, parity(level, n)), 1e-10)

    def test_native_two_drives_break_all_individual_parities(self):
        h = H_na(2)
        for level in LEVELS:
            self.assertGreater(comm_norm(h, parity(level, 2)), 1e-5)

    def test_switching_off_a_drive_restores_a_population_parity(self):
        self.assertLess(comm_norm(H_na(2, omega_r=0), parity("r", 2)), 1e-10)
        self.assertLess(comm_norm(H_na(2, omega_hf=0), parity("0", 2)), 1e-10)

    def test_three_parities_only_have_two_independent_generators(self):
        for n in (2, 3):
            product = parity("0", n) @ parity("1", n) @ parity("r", n)
            self.assertTrue(np.allclose(product, (-1)**n * np.eye(3**n)))

    def test_binary_pauli_representatives_only_on_allowed_states(self):
        # The two-qubit basis |10> is NOT a qutrit state.
        physical_indices = (0, 1, 3)  # |00>, |01>, |11>
        z = np.diag([1, -1])
        eye = np.eye(2)
        reps = {
            "0": -np.kron(eye, z),
            "1": np.kron(z, z),
            "r": np.kron(z, eye),
        }
        signs = {
            "0": [-1, +1, +1],
            "1": [+1, -1, +1],
            "r": [+1, +1, -1],
        }
        for level in LEVELS:
            actual = [int(reps[level][i, i]) for i in physical_indices]
            self.assertEqual(actual, signs[level])


if __name__ == "__main__":
    unittest.main()
