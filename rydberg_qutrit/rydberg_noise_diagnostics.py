"""Rydberg state diagnostics and quasi-static T2* noise extension.

Uses rydberg_qutrit.py and rydberg_full_channel.py already in this directory.
Experimental input: Evered et al., Nature 622 (2023), DOI 10.1038/s41586-023-06481-y.
The common-mode Gaussian detuning model for reported T2*=3us is an assumption,
not an experimentally validated Markovian collapse rate.
Run with OPENBLAS_NUM_THREADS=1 for efficient small-system integration.
"""
import numpy as np
from numpy.polynomial.hermite import hermgauss
from scipy.integrate import solve_ivp
from rydberg_qutrit import Evered2023, evered_two_atom_hamiltonian_numpy
from rydberg_full_channel import (
    all_decay_operators, _basis_input_operators, COMP_INDICES,
    choi_from_channel, kraus_from_choi, reference_unitary, pauli_strings_2q
)

R = np.diag([0, 0, 1, 0]).astype(complex)
N_R = np.kron(R, np.eye(4)) + np.kron(np.eye(4), R)
COMP = np.zeros(16, bool)
COMP[COMP_INDICES] = True
LOSS = np.array([3 in divmod(i, 4) for i in range(16)])
RYDBERG = ~(COMP | LOSS)


def detuning_samples(t2star_us=None, order=5):
    """Gaussian shot-to-shot angular detuning, coherent envelope exp[-(t/T2*)^2]."""
    if t2star_us is None:
        return [(0.0, 1.0)]
    if t2star_us <= 0:
        raise ValueError("T2* must be positive")
    x, w = hermgauss(order)
    return list(zip(2 * x / t2star_us, w / np.sqrt(np.pi)))


def evolve(operators, delta_rad_us=0.0, params=None, points=None):
    """Integrate Lindblad evolution of a batch of operator-basis matrices."""
    p = params or Evered2023()
    initial = np.asarray(operators, complex)
    jumps = tuple(all_decay_operators(p).values())
    cdc = sum((c.conj().T @ c for c in jumps), np.zeros((16, 16), complex))

    def rhs(t, y):
        a = y.reshape(-1, 16, 16)
        h = evered_two_atom_hamiltonian_numpy(t, p, include_loss_level=True)
        h = h - delta_rad_us * N_R
        out = -1j * (h @ a - a @ h)
        out -= 0.5 * (cdc @ a + a @ cdc)
        for c in jumps:
            out += c @ a @ c.conj().T
        return out.reshape(-1)

    times = None if points is None else np.linspace(0, p.gate_duration_us, points)
    sol = solve_ivp(
        rhs, (0, p.gate_duration_us), initial.reshape(-1),
        method="DOP853", rtol=1.5e-8, atol=1.5e-10, t_eval=times
    )
    if not sol.success:
        raise RuntimeError(sol.message)
    return sol.t, sol.y.T.reshape(-1, *initial.shape)


def _pure_state_matrices(labels):
    mats = []
    for label in labels:
        if len(label) != 2 or not set(label).issubset(set("01rL")):
            raise ValueError(label)
        idx = 4 * "01rL".index(label[0]) + "01rL".index(label[1])
        v = np.eye(16, dtype=complex)[:, idx]
        mats.append(np.outer(v, v.conj()))
    return np.array(mats)


def final_sectors(rho):
    diagonals = np.real(np.diag(rho))
    return {
        "computational": float(diagonals[COMP].sum()),
        "residual_rydberg": float(diagonals[RYDBERG].sum()),
        "any_L": float(diagonals[LOSS].sum()),
        "both_L": float(diagonals[15]),
        "trace": float(np.real(np.trace(rho))),
    }


def state_diagnostics(labels=("00", "01", "10", "11", "rr"),
                      t2star_us=None, order=5):
    """State-resolved diagnostic, includes noncomputational |rr> input."""
    operators = _pure_state_matrices(labels)
    averaged = np.zeros_like(operators)
    for delta, weight in detuning_samples(t2star_us, order):
        _, sequence = evolve(operators, float(delta))
        averaged += weight * sequence[-1]
    return {name: final_sectors(a) for name, a in zip(labels, averaged)}


def channel_basis(t2star_us=None, order=5):
    """State-independent channel E(|a><b|), 4D input -> 16D output."""
    basis = _basis_input_operators()
    averaged = np.zeros_like(basis)
    for delta, weight in detuning_samples(t2star_us, order):
        _, sequence = evolve(basis, float(delta))
        averaged += weight * sequence[-1]
    return averaged.reshape(4, 4, 16, 16)


def pauli_report(E):
    """Return trace-preserving physical channel and non-TP computational Pauli weights.

    The Pauli weights are relative to the nearest noiseless computational
    unitary, NOT independently verified as ideal CZ errors.
    """
    J = choi_from_channel(E)
    K = kraus_from_choi(J)
    Uref, _, _ = reference_unitary()
    A = [Uref.conj().T @ k[COMP_INDICES, :] for k in K]
    q = {
        s: float(sum(abs(np.trace(P.conj().T @ a))**2 for a in A) / 16)
        for s, P in pauli_strings_2q().items()
    }
    average_rho = sum((k @ k.conj().T for k in K), np.zeros((16,16), complex))/4
    sectors = final_sectors(average_rho)
    return {
        "pauli_weights": q,
        "sectors": sectors,
        "nonidentity_pauli": float(sum(v for name,v in q.items() if name != "II")),
        "normalization": float(sum(q.values()) + sectors["residual_rydberg"] + sectors["any_L"]),
        "choi_min_eigenvalue": float(np.linalg.eigvalsh((J+J.conj().T)/2)[0]),
        "kraus_count": len(K),
        "trace_error": float(np.max(np.abs(sum((k.conj().T@k for k in K))-np.eye(4)))),
    }


if __name__ == "__main__":
    import json
    from pathlib import Path
    outputs = {}
    for name,t2 in (("decay_only",None),("decay_plus_T2star",3.0)):
        print("Running",name,flush=True)
        states = state_diagnostics(t2star_us=t2)
        channel = pauli_report(channel_basis(t2star_us=t2))
        outputs[name] = {"states": states, "channel": channel}
        print(name,channel["nonidentity_pauli"],channel["sectors"],flush=True)
    Path(__file__).with_name("intermediate_t2star_results.json").write_text(
        json.dumps(outputs,indent=2)+"\n"
    )
