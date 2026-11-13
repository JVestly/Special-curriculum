"""Reduced Evered-2023 two-atom Rydberg decay: Choi, Kraus, Pauli and leakage.

For each atom, |0>,|1>,|r>,|L>. The physical output is 16-dimensional;
the computational input is 4-dimensional. Parameters are model-informed,
not a reproduction of the experimental CZ fidelity.
"""
import numpy as np
from scipy.integrate import solve_ivp
from scipy.linalg import polar
from rydberg_qutrit import Evered2023, evered_two_atom_hamiltonian_numpy
from rydberg_open_system import site_operator

COMP = np.array([0,1,4,5])
def all_decay_operators(params=None):
    p = params or Evered2023()
    branches = [(0,p.branching_to_0),(1,p.branching_to_1),(3,p.branching_to_loss)]
    assert np.isclose(sum(b for _,b in branches),1)
    result = {}
    for dst,ratio in branches:
        local=np.zeros((4,4),complex)
        local[dst,2]=np.sqrt(ratio*p.decay_rate_per_us)
        for i in range(2):
            result[f"atom{i}_r_to_{dst}"]=site_operator(local,i,2)
    return result

def construct_channel(params=None,rtol=2e-8,atol=2e-10):
    """Compute E(|a><b|) for all 16 operator-basis elements at once."""
    p = params or Evered2023()
    operators=list(all_decay_operators(p).values())
    cdc=[c.conj().T@c for c in operators]
    initial=np.zeros((16,16,16),complex)
    for a in range(4):
        for b in range(4):
            initial[4*a+b,COMP[a],COMP[b]]=1
    def rhs(t,y):
        e=y.reshape(16,16,16)
        h=evered_two_atom_hamiltonian_numpy(t,p,include_loss_level=True)
        derivative=-1j*(h@e-e@h)
        for c,d in zip(operators,cdc):
            derivative+=c@e@c.conj().T - (d@e+e@d)/2
        return derivative.reshape(-1)
    sol=solve_ivp(rhs,(0,p.gate_duration_us),initial.reshape(-1),
                  method="DOP853",rtol=rtol,atol=atol)
    if not sol.success: raise RuntimeError(sol.message)
    return sol.y[:,-1].reshape(4,4,16,16)

def choi_from_channel(E):
    """Choi matrix with reference-input index first (64 x 64)."""
    return E.transpose(0,2,1,3).reshape(64,64)

def kraus_from_choi(J,tol=1e-9):
    J=(J+J.conj().T)/2
    values,vectors=np.linalg.eigh(J)
    if values[0]<-tol: raise ValueError("Choi matrix not positive")
    return [np.sqrt(lam)*vectors[:,i].reshape(4,16).T
            for i,lam in enumerate(values) if lam>tol]

def reference_unitary(params=None):
    """Nearest unitary to the noiseless computational block, NOT a CZ target."""
    p=params or Evered2023()
    def rhs(t,y):
        H=evered_two_atom_hamiltonian_numpy(t,p,include_loss_level=True)
        return (-1j*H@y.reshape(16,16)).reshape(-1)
    sol=solve_ivp(rhs,(0,p.gate_duration_us),np.eye(16, dtype=complex).reshape(-1),
                  method="DOP853",rtol=2e-9,atol=2e-11)
    if not sol.success: raise RuntimeError(sol.message)
    U=sol.y[:,-1].reshape(16,16)[np.ix_(COMP,COMP)]
    uref,_=polar(U)
    return uref, 1-min(np.linalg.svd(U,compute_uv=False))**2

def pauli_strings_2q():
    I=np.eye(2, dtype=complex)
    X=np.array([[0,1],[1,0]],complex)
    Y=np.array([[0,-1j],[1j,0]],complex)
    Z=np.diag([1,-1]).astype(complex)
    one=dict(I=I,X=X,Y=Y,Z=Z)
    return {a+b:np.kron(A,B) for a,A in one.items() for b,B in one.items()}

def channel_report(params=None):
    """Unconditional Pauli-twirl weights plus actual population leakage."""
    E=construct_channel(params)
    J=choi_from_channel(E)
    Ks=kraus_from_choi(J)
    uref,residual=reference_unitary(params)
    comp=[uref.conj().T@K[COMP,:] for K in Ks]
    q={name:float(sum(abs(np.trace(P.conj().T@A))**2 for A in comp)/16)
       for name,P in pauli_strings_2q().items()}
    rho=sum(K@K.conj().T for K in Ks)/4
    rydberg=0.
    loss=0.
    for index in range(16):
        i,j=divmod(index,4)
        if 3 in (i,j): loss+=rho[index,index].real
        elif 2 in (i,j): rydberg+=rho[index,index].real
    comp_total=float(rho[np.ix_(COMP,COMP)].trace().real)
    return {
        "pauli_weights":q,
        "p_computational":comp_total,
        "p_rydberg_residual":float(rydberg),
        "p_loss":float(loss),
        "p_total":comp_total+rydberg+loss,
        "kraus_count":len(Ks),
        "choi_shape":list(J.shape),
        "max_trace_error":float(np.max(abs(sum(K.conj().T@K for K in Ks)-np.eye(4)))),
        "min_choi_eigenvalue":float(np.linalg.eigvalsh((J+J.conj().T)/2)[0]),
        "max_noiseless_projection_residual":float(residual),
    }

if __name__=="__main__":
    import json
    print(json.dumps(channel_report(),indent=2))
