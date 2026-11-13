"""State-resolved Rydberg loss diagnostics and experimentally informed T2* extension.

Reduced effective model, not the full Evered et al. eight-state calculation.
- 4 states/atom: |0>, |1>, |r>, |L>
- two atoms, phase-shaped Rydberg pulse
- spontaneous Rydberg decay: |r> -> |0>, |1>, |L>
- optional quasi-static Gaussian Rydberg detuning, calibrated by a measured
  T2* assuming Gaussian Ramsey envelope exp[-(t/T2*)**2]. This is NOT an
  additional Lindblad collapse operator: it describes inhomogeneous averaging.

Published parameters: Evered et al., Nature 622 (2023),
https://doi.org/10.1038/s41586-023-06481-y (T2*=3 us, T_r=88 us).
"""
from __future__ import annotations
from dataclasses import dataclass
from typing import Iterable
import numpy as np
from scipy.integrate import solve_ivp
from numpy.polynomial.hermite import hermgauss

from rydberg_qutrit import Evered2023, evered_two_atom_hamiltonian_numpy
from rydberg_full_channel import all_decay_operators, COMP_INDICES, _basis_input_operators, choi_from_channel, kraus_from_choi, reference_unitary, pauli_strings_2q

DIM=16
R=np.diag([0.,0.,1.,0.])
N_R=np.kron(R,np.eye(4)) + np.kron(np.eye(4),R)
P_LL=np.zeros((16,16),complex);P_LL[15,15]=1.0
# disjoint final-outcome sectors, not stochastic trajectory classification
MASK_COMP=np.zeros(16,bool); MASK_COMP[COMP_INDICES]=True
MASK_LOSS=np.array([3 in divmod(i,4) for i in range(16)])
MASK_RYDBERG=(~MASK_COMP)&(~MASK_LOSS)


def detuning_quadrature(t2_star_us:float=3.0, order:int=7):
    """Gaussian angular detuning, assuming exp[-(t/T2*)^2] Ramsey envelope.

    For normal delta~N(0,sigma^2), <e^{-i delta t}>=exp(-sigma^2*t^2/2).
    Thus sigma=sqrt(2)/T2*. Nodes output angular frequencies rad/us.
    Default: common-mode detuning on both atoms, not independent Doppler shifts.
    """
    if t2_star_us<=0 or order<=0: raise ValueError('positive T2* and quadrature order required')
    nodes, weights=hermgauss(order)
    sigma=np.sqrt(2.0)/t2_star_us
    delta=np.sqrt(2.0)*sigma*nodes
    weights=weights/np.sqrt(np.pi)
    return delta, weights


def _initial_kets(labels:Iterable[str]):
    out=[]
    for s in labels:
        if len(s)!=2 or any(x not in ('0','1','r','L') for x in s):
            raise ValueError(f'Invalid basis ket: {s}')
        index=4*('01rL'.index(s[0]))+'01rL'.index(s[1])
        k=np.zeros(16,complex); k[index]=1
        out.append(k)
    return np.array(out)


def integrate_batch(initial_batch:np.ndarray,
                    params:Evered2023|None=None,
                    *, delta_rad_per_us:float=0.0,
                    dephasing_rate_per_us:float=0.0,
                    rtol:float=1.5e-8,atol:float=1.5e-10,
                    grid_points:int|None=None):
    """Integrate batched density matrices or arbitrary basis matrices.

    On-site Markovian dephasing is OPTIONAL and defaults to zero.
    collapse sqrt(2*gamma_phi) n_r has 1-r coherence decay gamma_phi.
    """
    params=params or Evered2023()
    y0=np.asarray(initial_batch,complex)
    if y0.ndim!=3 or y0.shape[1:]!=(16,16): raise ValueError('Expected (batch,16,16)')
    if dephasing_rate_per_us<0:raise ValueError('dephasing rate must be >=0')
    cs=list(all_decay_operators(params).values())
    if dephasing_rate_per_us>0:
        for site in range(2):
            proj=np.kron(R,np.eye(4)) if site==0 else np.kron(np.eye(4),R)
            cs.append(np.sqrt(2*dephasing_rate_per_us)*proj)
    cs=tuple(cs)
    cdc=sum((c.conj().T@c for c in cs),np.zeros((16,16),complex))
    tau=params.gate_duration_us
    t_eval=np.linspace(0,tau,grid_points) if grid_points else None
    def rhs(t,y):
        rho=y.reshape(-1,16,16)
        H=evered_two_atom_hamiltonian_numpy(t,params,include_loss_level=True)
        H=H-delta_rad_per_us*N_R
        out=-1j*(H@rho-rho@H)
        out+=sum((c@rho@c.conj().T for c in cs),np.zeros_like(rho))
        out-=0.5*(cdc@rho+rho@cdc)
        return out.reshape(-1)
    sol=solve_ivp(rhs,(0,tau),y0.reshape(-1),method='DOP853',rtol=rtol,atol=atol,t_eval=t_eval)
    if not sol.success: raise RuntimeError(sol.message)
    return sol.t,sol.y.T.reshape(-1,len(y0),16,16)


def outcomes(rho:np.ndarray):
    """Final state sector probabilities. For 2 atoms, p_LL is both lost."""
    diag=np.diag(rho)
    comp=diag[MASK_COMP].sum().real
    loss=diag[MASK_LOSS].sum().real
    r=diag[MASK_RYDBERG].sum().real
    both=diag[15].real
    return dict(computational=float(comp),rydberg=float(r),
                any_loss=float(loss),both_lost=float(both),
                trace=float(np.trace(rho).real),
                mean_r=float(np.trace(N_R@rho).real))


def state_diagnostics(params=None, labels=('00','01','10','11','rr'),
                      *, t2_star_us:float|None=None, quad_order=7,
                      points=51):
    """Compare pure-state inputs through identical pulse and full decay model.

    Returns table plus per-state time curves. If T2* set, averages Lindblad
    channels over a *common* quasi-static detuning on both atoms.
    """
    p=params or Evered2023()
    ks=_initial_kets(labels)
    initial=np.array([np.outer(k,k.conj()) for k in ks])
    if t2_star_us is None: deltas,weights=[0.0],[1.0]
    else: deltas,weights=detuning_quadrature(t2_star_us,quad_order)
    accum=None;time=None
    for d,w in zip(deltas,weights):
        t,series=integrate_batch(initial,p,delta_rad_per_us=float(d),grid_points=points)
        if accum is None:accum=np.zeros_like(series)
        accum+=w*series;time=t
    result={}
    for j,label in enumerate(labels):
        series=accum[:,j,:,:]
        outcome=[outcomes(rho) for rho in series]
        result[label]={key:np.array([item[key] for item in outcome])
                       for key in outcome[0]}
    return time,result


def full_channel_t2_star(params=None, *, t2_star_us:float|None=3.0,
                         quad_order=7):
    """Channel E(|a><b|), averaging common-mode static detuning trajectories.
    This preserves CPTP (up to solver tolerance) and allows Choi/Kraus analysis.
    """
    p=params or Evered2023()
    if t2_star_us is None: deltas,weights=[0.0],[1.0]
    else: deltas,weights=detuning_quadrature(t2_star_us,quad_order)
    acc=np.zeros((16,16,16),complex)
    for d,w in zip(deltas,weights):
        _,states=integrate_batch(_basis_input_operators(),p,delta_rad_per_us=float(d))
        acc+=w*states[-1]
    return acc.reshape(4,4,16,16)


def report_from_channel(E,params=None):
    p=params or Evered2023()
    J=choi_from_channel(E)
    eig=np.linalg.eigvalsh((J+J.conj().T)/2)
    Ks=kraus_from_choi(J,eigen_tol=1e-9)
    err=np.max(np.abs(sum((k.conj().T@k for k in Ks))-np.eye(4)))
    uref,_,_=reference_unitary(p)
    A=[uref.conj().T@k[COMP_INDICES,:] for k in Ks]
    weights={name:float(np.sum([abs(np.trace(P.conj().T@a))**2 for a in A])/16)
             for name,P in pauli_strings_2q().items()}
    rho_avg=sum((k@k.conj().T for k in Ks),np.zeros((16,16),complex))/4
    pops=outcomes(rho_avg)
    return {'pauli_weights':weights,'outcomes':pops,
            'nonidentity_pauli':sum(w for name,w in weights.items() if name!='II'),
            'kraus_count':len(Ks),'min_choi_eigenvalue':float(eig.min()),
            'trace_error':float(err),'total_normalization':sum(weights.values())+pops['rydberg']+pops['any_loss']}

if __name__=='__main__':
    p=Evered2023()
    for t2 in (None,3.0):
        t,res=state_diagnostics(p,t2_star_us=t2)
        print('T2*=',t2)
        for name,x in res.items():
            print(name, {key:float(x[key][-1]) for key in ('computational','rydberg','any_loss','both_lost','trace')})
