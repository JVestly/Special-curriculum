"""Intermediate-state scattering in a reduced neutral-atom Rydberg channel.

The virtual |e> = 6P_3/2 excited state is adiabatically eliminated at
large single-photon detuning. The effective jump operator is obtained from

   H_e = -Delta |e><e| + (Omega_b exp(+i phi)|1><e| +
                              Omega_red |e><r| + h.c.)/2.

To first order in Omega/Delta, the intermediate-state amplitude is
  <e|psi> ~= [Omega_b exp(-i phi)<1|psi> + Omega_red<r|psi>]/(2 Delta).

For an ASSUMED emission branch e -> |target>, the effective jump is
  C_eff(t)=sqrt(Gamma_e)*|target> [Omega_b exp(-i phi(t))<1| +
                                      Omega_red<r|] / (2 Delta).

No intermediate-state Hamiltonian is evolved directly. The effective 2-photon
Hamiltonian already contains the coherent 1<->r coupling: retaining it and
adding only the scattering correction avoids counting it twice. This model
is a leading-order approximation; experimental 6P multilevel structure,
intermediate-state branching ratios, pulse Stark shifts and laser scattering
from |0> are NOT modeled.

Gamma_e ~ 1/(0.1128 us) from measured Rb 6P3/2 radiative lifetime in
experiments cited by https://doczz.net/doc/8276929/lifetime-measurements-of-the-5d-states-of-rubidium
Omega_420/2pi=237 MHz, Omega_1013/2pi=303 MHz, Delta/2pi=7.8 GHz
from Evered et al. Nature 622 (2023), DOI:10.1038/s41586-023-06481-y.

ASSUMPTION: decay destination default |1>, a Rayleigh-like model. The true
spin-changing/leakage branching ratios are NOT fixed by gamma_e. Alternative
branch destination |0> or |L> is supported as a modeling sensitivity check.
"""
from __future__ import annotations
from dataclasses import dataclass
import json
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
from numpy.polynomial.hermite import hermgauss

from rydberg_qutrit import Evered2023, evered_two_atom_hamiltonian_numpy
from rydberg_open_system import site_operator
from rydberg_full_channel import (
    COMP_INDICES, _basis_input_operators, choi_from_channel, kraus_from_choi,
    reference_unitary, pauli_strings_2q
)
from rydberg_diagnostics import N_R, detuning_quadrature, outcomes


@dataclass(frozen=True)
class ScatteringParameters:
    omega_blue_mhz: float = 237.0    # Omega_420 / 2pi
    omega_red_mhz: float = 303.0     # Omega_1013 / 2pi
    intermediate_detuning_mhz: float = 7800.0  # Delta / 2pi
    intermediate_lifetime_us: float = 0.1128   # 112.8 ns
    target: str = '1'              # ASSUMED branch: 1 (Rayleigh-like), 0, L
    scattering_scale: float = 1.0  # sensitivity parameter, default no correction

    @property
    def decay_rate_us_inv(self):
        return 1 / self.intermediate_lifetime_us

    def validate(self):
        vals = [self.omega_blue_mhz,self.omega_red_mhz,
                self.intermediate_detuning_mhz,self.intermediate_lifetime_us,
                self.scattering_scale]
        if any(not np.isfinite(v) for v in vals):
            raise ValueError('Scattering parameters must be finite')
        if self.intermediate_detuning_mhz == 0 or self.intermediate_lifetime_us <= 0:
            raise ValueError('Detuning must be nonzero; lifetime positive')
        if self.scattering_scale < 0:
            raise ValueError('Scale must be nonnegative')
        if self.target not in {'0','1','L'}:
            raise ValueError('target must be 0, 1, or L')


def effective_local_jump(phi, p=ScatteringParameters()):
    """4x4 site jump sqrt(Gamma_e)|target><e_virtual(phi)|.

    Dimensionless admixtures use ratios Omega_420/(2Delta), Omega_1013/(2Delta),
    because the 2pi factors cancel. Uses phase convention from Hamiltonian.
    """
    p.validate()
    dest = {'0':0,'1':1,'L':3}[p.target]
    j = np.zeros((4,4),complex)
    scale = np.sqrt(p.decay_rate_us_inv*p.scattering_scale)
    j[dest,1] = scale*p.omega_blue_mhz*np.exp(-1j*phi)/(2*p.intermediate_detuning_mhz)
    j[dest,2] = scale*p.omega_red_mhz/(2*p.intermediate_detuning_mhz)
    return j


def dissipators_rydberg_decay(rydberg=Evered2023()):
    """Six Rydberg-branch operators and sum C^dag C."""
    from rydberg_full_channel import all_decay_operators
    cs = list(all_decay_operators(rydberg).values())
    cdc = sum((c.conj().T@c for c in cs),np.zeros((16,16),complex))
    return cs, cdc


def integrate_channel(
    rydberg=Evered2023(), scattering=ScatteringParameters(),
    *, detuning_rad_us=0.0, include_scattering=True,
    rtol=1.5e-8, atol=1e-10,
):
    """Return channel basis E(|a><b|), shape (4,4,16,16).

    Lindblad with the 6 Rydberg decay operators plus 2 additional local
    scattering operators derived from the virtual intermediate state.
    Input is only the 4D computational subspace; output includes |r>, |L>.
    """
    scattering.validate()
    cs_r, cdc_r = dissipators_rydberg_decay(rydberg)
    ys = _basis_input_operators()
    ident = np.eye(4, dtype=complex)
    tau = rydberg.gate_duration_us
    def rhs(t,y):
        batch = y.reshape(16,16,16)
        h=evered_two_atom_hamiltonian_numpy(t,rydberg,include_loss_level=True)
        if detuning_rad_us:
            h=h-detuning_rad_us*N_R
        out = -1j*(h@batch-batch@h)
        for c in cs_r:
            out += c@batch@c.conj().T
        cdc = cdc_r
        if include_scattering:
            jloc=effective_local_jump(rydberg.phase_r_rad(t),scattering)
            c1=np.kron(jloc,ident)
            c2=np.kron(ident,jloc)
            for c in (c1,c2):
                out+=c@batch@c.conj().T
                cdc=cdc + c.conj().T@c
        out -= 0.5*(cdc@batch+batch@cdc)
        return out.reshape(-1)
    sol=solve_ivp(rhs,(0.0,tau),ys.reshape(-1),method='DOP853',rtol=rtol,atol=atol)
    if not sol.success:
        raise RuntimeError(sol.message)
    return sol.y[:,-1].reshape(4,4,16,16)


def averaged_channel(
    rydberg=Evered2023(), scattering=ScatteringParameters(),
    *, include_scattering=True, t2_star_us=3.0, quadrature_order=5,
):
    """Average quasi-static detuning trajectories; None disables detuning noise.

    The detuning is assumed common to the two atoms, as an explicit modeling
    assumption and not an experimentally measured correlation matrix.
    """
    if t2_star_us is None:
        nodes, weights=[0.0],[1.0]
    else:
        nodes, weights=detuning_quadrature(t2_star_us,quadrature_order)
    total=np.zeros((4,4,16,16),complex)
    for delta,w in zip(nodes,weights):
        total+=w*integrate_channel(
            rydberg,scattering,detuning_rad_us=float(delta),
            include_scattering=include_scattering,
        )
    return total


def analyze(E, rydberg=Evered2023(), *, eigen_tol=1e-10):
    J=choi_from_channel(E)
    eig=np.linalg.eigvalsh((J+J.conj().T)/2)
    Ks=kraus_from_choi(J,eigen_tol=eigen_tol)
    eye=np.eye(4)
    completeness=sum((K.conj().T@K for K in Ks),np.zeros((4,4),complex))
    U_ref,_,_=reference_unitary(rydberg)
    A=[U_ref.conj().T @ K[COMP_INDICES,:] for K in Ks]
    weights={name:float(sum(abs(np.trace(P.conj().T@a))**2 for a in A)/16)
             for name,P in pauli_strings_2q().items()}
    rho_avg=sum((k@k.conj().T for k in Ks),np.zeros((16,16),complex))/4
    pops=outcomes(rho_avg)
    reconstructed=sum((np.outer(K.T.ravel(),K.T.ravel().conj()) for K in Ks),
                      np.zeros((64,64),complex))
    return dict(
       J=J, Kraus=np.stack(Ks), pauli=weights, rho_avg=rho_avg, outcomes=pops,
       min_choi_eig=float(eig.min()),
       trace_error=float(np.max(abs(completeness-eye))),
       reconstruction_error=float(np.max(abs(reconstructed-J))),
       normalization=float(sum(weights.values())+pops['rydberg']+pops['any_loss']),
       nonidentity=float(sum(w for name,w in weights.items() if name!='II')),
       reference=U_ref,
    )


def export_results(data, dest, model):
    """Export ALL COMPLEX 16x4 Kraus matrices, Choi and Pauli tables.

    Kraus decompositions are not unique. We export canonical Choi eigenmodes
    with numpy eigh phase conventions; changing numerical libraries may
    rotate nearly degenerate Kraus operators.
    """
    dest=Path(dest); dest.mkdir(parents=True,exist_ok=True)
    np.savez_compressed(dest/'full_kraus_chosen_model.npz',
                        kraus=data['Kraus'],choi=data['J'],
                        reference=data['reference'],
                        pauli_names=np.array(list(data['pauli'])),
                        pauli_weights=np.array(list(data['pauli'].values())))
    with (dest/'kraus_operators.txt').open('w', encoding='utf-8') as f:
        f.write('ALL Kraus matrices (16 rows physical |0,1,r,L>^2 by 4 columns computational |00>,|01>,|10>,|11>).\n')
        f.write('Kraus operators extracted from Choi eigenvectors; these are NOT Pauli errors.\n')
        f.write('Numerical decomposition not unique; no-jump labels are NOT assigned to individual matrices.\n')
        f.write('Each complex number is displayed as real +/- imaginary*1j.\n\n')
        for j,k in enumerate(data['Kraus']):
            f.write(f'K_{j:02d}, Frobenius norm = {np.linalg.norm(k):.12g}:\n')
            f.write(np.array2string(k, precision=12, separator=', ', max_line_width=155, threshold=np.inf))
            f.write('\n\n')
    obj={
        'model':model,
        'shape_choi':list(data['J'].shape),
        'shape_each_kraus':list(data['Kraus'][0].shape),
        'kraus_count':int(len(data['Kraus'])),
        'pauli_probabilities':data['pauli'],
        'outcomes':data['outcomes'],
        'nonidentity_pauli':data['nonidentity'],
        'min_choi_eig':data['min_choi_eig'],
        'trace_error':data['trace_error'],
        'reconstruction_error':data['reconstruction_error'],
        'normalization':data['normalization'],
    }
    (dest/'scattering_channel_results.json').write_text(json.dumps(obj,indent=2))
    return obj
