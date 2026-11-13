"""First open-system module: irreversible Rydberg loss to |L>.

Local basis: |0>, |1>, |r>, |L> (indices 0,1,2,3).
Units: Hamiltonians H/hbar in rad/us, collapse operators in 1/sqrt(us),
       time in us. Two independent atoms: [C_L tensor I, I tensor C_L].

This module is purposely a *single noise mechanism*, not a full Evered 2023
model. It neglects the 10.6% branching into computational levels, and
intermediate-level decay, laser noise, finite pulse shaping and other errors.

The functions here evolve a chosen rho(0) to VALIDATE loss dynamics.
Constructing a state-independent quantum channel / Choi is a later module.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Callable

import numpy as np
from scipy.integrate import solve_ivp

from rydberg_qutrit import Evered2023, evered_two_atom_hamiltonian_numpy

LOCAL_DIM = 4
R_INDEX = 2
LOSS_INDEX = 3


def site_operator(local_op: np.ndarray, atom: int, n_atoms: int) -> np.ndarray:
    """Tensor a local 4x4 operator into the chosen atom of n atoms."""
    if not (0 <= atom < n_atoms):
        raise ValueError("Atom index out of range")
    if np.shape(local_op) != (4, 4):
        raise ValueError("Expected a 4x4 local operator")
    factors = [np.eye(4, dtype=complex) for _ in range(n_atoms)]
    factors[atom] = np.asarray(local_op, dtype=complex)
    result = factors[0]
    for op in factors[1:]:
        result = np.kron(result, op)
    return result


def loss_collapse_operators(n_atoms: int, gamma_loss_per_us: float) -> list[np.ndarray]:
    """One physical noise mechanism r -> L, independently on each atom.

    C_i = sqrt(gamma_loss) I^tensor(i) x |L><r| x I^tensor(N-i-1).
    Calling this for two atoms returns TWO C_i operators, not one collective
    jump; a single shared jump would introduce unphysical collective decay.
    """
    if not isinstance(n_atoms, int) or not (1 <= n_atoms <= 3):
        raise ValueError("n_atoms must be 1, 2, or 3 for this dense prototype")
    if not np.isfinite(gamma_loss_per_us) or gamma_loss_per_us < 0:
        raise ValueError("gamma_loss_per_us must be finite and nonnegative")
    local_jump = np.zeros((4, 4), dtype=complex)
    local_jump[LOSS_INDEX, R_INDEX] = np.sqrt(gamma_loss_per_us)
    return [site_operator(local_jump, i, n_atoms) for i in range(n_atoms)]


def lindblad_rhs(
    time_us: float,
    rho_vector: np.ndarray,
    hamiltonian_at_time: Callable[[float], np.ndarray],
    collapses: list[np.ndarray],
    cdc_terms: list[np.ndarray] | None = None,
) -> np.ndarray:
    """Lindblad d rho/dt, flattened using NumPy row-major order."""
    dim = collapses[0].shape[0] if collapses else int(round(np.sqrt(rho_vector.size)))
    rho = rho_vector.reshape(dim, dim)
    h = hamiltonian_at_time(time_us)
    derivative = -1j * (h @ rho - rho @ h)
    if cdc_terms is None:
        cdc_terms = [c.conj().T @ c for c in collapses]
    for c, cdc in zip(collapses, cdc_terms):
        derivative += c @ rho @ c.conj().T - 0.5 * (cdc @ rho + rho @ cdc)
    return derivative.reshape(-1)


def basis_density(label: str) -> np.ndarray:
    """Make a pure-state density matrix from e.g. '11', 'r1', or 'rL'.

    This is for *state evolution/testing*, not needed for channel construction.
    """
    table = {"0": 0, "1": 1, "r": 2, "L": 3}
    if not label or any(c not in table for c in label):
        raise ValueError("Use a nonempty string from {0,1,r,L}")
    n = len(label)
    idx = 0
    for letter in label:
        idx = 4 * idx + table[letter]
    ket = np.zeros(4**n, dtype=complex)
    ket[idx] = 1
    return np.outer(ket, ket.conj())


@dataclass
class LossEvolution:
    times_us: np.ndarray
    rho: np.ndarray  # shape (n_times, 4^n, 4^n)
    any_loss_probability: np.ndarray
    expected_number_lost: np.ndarray
    expected_rydberg_population: np.ndarray


def simulate_lindblad(
    hamiltonian_at_time: Callable[[float], np.ndarray],
    collapses: list[np.ndarray],
    rho0: np.ndarray,
    times_us: np.ndarray,
    *,
    rtol: float = 1e-8,
    atol: float = 1e-10,
) -> LossEvolution:
    """Evolve one diagnostic initial state with SciPy (works without QuTiP)."""
    times = np.asarray(times_us, dtype=float)
    if times.ndim != 1 or len(times) < 2 or times[0] < 0 or np.any(np.diff(times) <= 0):
        raise ValueError("times_us must be increasing and nonnegative")
    dim = rho0.shape[0]
    if rho0.shape != (dim, dim):
        raise ValueError("rho0 must be square")
    n_atoms = int(round(np.log(dim) / np.log(4)))
    if 4**n_atoms != dim:
        raise ValueError("dimension must be 4**n_atoms")
    if any(c.shape != (dim, dim) for c in collapses):
        raise ValueError("collapse operator dimension mismatch")
    if not np.allclose(rho0, rho0.conj().T) or not np.isclose(np.trace(rho0), 1):
        raise ValueError("rho0 must be Hermitian and trace one")

    cdc = [c.conj().T @ c for c in collapses]
    sol = solve_ivp(
        lambda t, y: lindblad_rhs(t, y, hamiltonian_at_time, collapses, cdc),
        (float(times[0]), float(times[-1])),
        rho0.reshape(-1),
        t_eval=times,
        method="DOP853", rtol=rtol, atol=atol,
    )
    if not sol.success:
        raise RuntimeError(f"Lindblad integration failed: {sol.message}")
    rho = sol.y.T.reshape(-1, dim, dim)
    loss_projector = np.zeros((4, 4), dtype=complex)
    loss_projector[LOSS_INDEX, LOSS_INDEX] = 1
    rydberg_projector = np.zeros((4, 4), dtype=complex)
    rydberg_projector[R_INDEX, R_INDEX] = 1
    lost_expect = np.zeros(len(times))
    rydberg_expect = np.zeros(len(times))
    for i in range(n_atoms):
        lost_expect += np.einsum("ab,tba->t", site_operator(loss_projector, i, n_atoms), rho).real
        rydberg_expect += np.einsum("ab,tba->t", site_operator(rydberg_projector, i, n_atoms), rho).real
    no_loss_local = np.diag([1, 1, 1, 0])
    no_loss_global = no_loss_local.astype(complex)
    for _ in range(n_atoms - 1):
        no_loss_global = np.kron(no_loss_global, no_loss_local)
    p_no_loss = np.einsum("ab,tba->t", no_loss_global, rho).real
    return LossEvolution(
        times_us=times,
        rho=rho,
        any_loss_probability=1 - p_no_loss,
        expected_number_lost=lost_expect,
        expected_rydberg_population=rydberg_expect,
    )


def simulate_evered_loss(
    params: Evered2023 | None = None,
    initial_state: str = "11",
    num_times: int = 101,
) -> LossEvolution:
    """Two-atom driven pulse, only Rydberg-to-loss dissipative mechanism."""
    p = params if params is not None else Evered2023()
    gamma_loss = p.decay_rate_per_us * p.branching_to_loss
    collapses = loss_collapse_operators(2, gamma_loss)
    rho0 = basis_density(initial_state)
    if rho0.shape != (16, 16):
        raise ValueError("initial_state must describe two atoms")
    times = np.linspace(0.0, p.gate_duration_us, num_times)
    return simulate_lindblad(
        lambda t: evered_two_atom_hamiltonian_numpy(
            t, p, include_loss_level=True
        ),
        collapses,
        rho0,
        times,
    )


def loss_collapse_operators_qutip(n_atoms: int, gamma_loss_per_us: float):
    """Optional QuTiP wrappers matching the NumPy collapse matrices."""
    try:
        import qutip as qt
    except ImportError as exc:
        raise ImportError("Install QuTiP 5: pip install 'qutip>=5'") from exc
    return [qt.Qobj(c, dims=[[4]*n_atoms, [4]*n_atoms])
            for c in loss_collapse_operators(n_atoms, gamma_loss_per_us)]


if __name__ == "__main__":
    params = Evered2023()
    result = simulate_evered_loss(params)
    print("Rydberg loss branching ratio:", params.branching_to_loss)
    print(f"Gamma_loss: {params.decay_rate_per_us * params.branching_to_loss:.8f} /us")
    print(f"Pulse duration: {params.gate_duration_us:.6f} us")
    print(f"P(at least one loss) after pulse: {result.any_loss_probability[-1]:.8f}")
    print(f"Expected lost atoms after pulse: {result.expected_number_lost[-1]:.8f}")
    print(f"Max trace drift: {np.max(np.abs(np.trace(result.rho, axis1=1, axis2=2)-1)):.2e}")