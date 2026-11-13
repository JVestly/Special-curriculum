"""Regression checks for finite SAT/SMT study, run: python -m unittest -v test_rydberg_sat_study.py"""
import unittest, itertools, json
from pathlib import Path
import numpy as np
from rydberg_sat_study import (pauli_to_bits,bits_to_pauli,symp,span,
    enumerate_stabilizers,pair_probabilities,correctable_error_set,ml_success,
    distance,smt,run_z3)
ROOT=Path(__file__).resolve().parent

class TestRydbergSAT(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.groups=enumerate_stabilizers()
        cls.q,cls.w=pair_probabilities(ROOT/'rydberg_pair_pauli_weights.json')
        cls.hard=[pauli_to_bits(x) for x in ('IIII','ZIII','IZII','IIZI','IIIZ')]

    def test_symplectic(self):
        self.assertEqual(symp(pauli_to_bits('XI'),pauli_to_bits('ZI'),n=2),1)
        self.assertEqual(symp(pauli_to_bits('XI'),pauli_to_bits('IZ'),n=2),0)
        for word in itertools.product('IXYZ',repeat=4):
            p=''.join(word)
            self.assertEqual(bits_to_pauli(pauli_to_bits(p)),p)

    def test_enumeration(self):
        self.assertEqual(len(self.groups),11475)
        for g,b in list(self.groups.items())[::137]:
            self.assertEqual(len(g),8)
            self.assertEqual(span(b),g)
            self.assertTrue(all(symp(x,y)==0 for x in g for y in g))

    def test_noise_accounting(self):
        self.assertAlmostEqual(sum(self.q.values()),0.9984648622673,places=7)
        self.assertAlmostEqual(float(self.w.sum()),sum(self.q.values())**2,places=12)
        self.assertTrue(np.all(self.w>=0))

    def test_code_distance_and_correction(self):
        d2=[(g,b) for g,b in self.groups.items() if distance(g,b)>=2]
        self.assertEqual(len(d2),2268)
        strict=[(g,b) for g,b in d2 if correctable_error_set(g,b,self.hard)]
        self.assertEqual(len(strict),768)
        seeded=[(g,b) for g,b in strict if pauli_to_bits('ZZII') in g]
        self.assertEqual(len(seeded),8)
        both=[(g,b) for g,b in strict if all(pauli_to_bits(x) in g for x in ('ZZII','IIZZ'))]
        self.assertFalse(both)

    def test_ml_decoder(self):
        g,b=next(iter(self.groups.items()))
        p,_=ml_success(g,b,self.w)
        self.assertLessEqual(p,self.w.sum()+1e-12)
        self.assertGreaterEqual(p,0)

    def test_hamiltonian_symmetries(self):
        # Exact two-qutrit Rydberg-only H: |1><r| drive plus diagonal terms.
        Q=np.diag([1,-1,-1]);I=np.eye(3)
        drive=np.array([[0,0,0],[0,0,1],[0,1,0]],dtype=complex)
        n_r=np.diag([0,0,1]);n_1=np.diag([0,1,0])
        H=np.kron(drive,I)+np.kron(I,drive)+27*np.kron(n_r,n_r)+3*np.kron(n_1,I)
        for sym in (np.kron(Q,I),np.kron(I,Q),np.kron(Q,Q)):
            self.assertLess(np.linalg.norm(sym@H-H@sym),1e-12)
        hf=np.array([[0,1,0],[1,0,0],[0,0,0]])
        Hhf=H+np.kron(hf,I)
        self.assertGreater(np.linalg.norm(np.kron(Q,I)@Hhf-Hhf@np.kron(Q,I)),1)

    def test_smt_sat(self):
        # Actual Boolean satisfiability, when libz3 is installed.
        status,gens=run_z3(smt(self.hard,(),2))
        if status=='not-installed':self.skipTest('Native Z3 absent; exact enumeration tests still run.')
        self.assertEqual(status,'sat')
        self.assertTrue(correctable_error_set(span(gens),gens,self.hard))
        self.assertGreaterEqual(distance(span(gens),gens),2)
        status,gens=run_z3(smt(self.hard,('ZZII','IIZZ'),2))
        self.assertEqual(status,'unsat')
        self.assertIsNone(gens)

    def test_regression_results(self):
        r=json.loads((ROOT/'SAT_Rydberg/results.json').read_text())
        self.assertEqual(r['unseeded']['hard_kl_feasible'],768)
        self.assertEqual(r['seed_ZZII']['hard_kl_feasible'],8)
        self.assertEqual(r['seed_ZZII_IIZZ']['hard_kl_feasible'],0)
        self.assertAlmostEqual(r['unseeded']['best']['logical_failure_conditional_on_computational_survival'],.0003858895603902651,places=11)
        self.assertAlmostEqual(r['weighted_only_d2']['seed_ZZII_IIZZ']['logical_failure_conditional_on_computational_survival'],.0006444978873057223,places=11)

if __name__=='__main__':unittest.main()