"""Bloqade-convention three-level neutral-atom Rydberg Hamiltonian.

Basis per atom: |0> = [1,0,0], |1> = [0,1,0], |r> = [0,0,1].
Tensor basis has atom 0 on the left: |00>, |01>, |0r>, |10>, ... |rr>.

Units:
    positions_um: micrometers
    omega_*_mhz, delta_*_mhz, C6_mhz_um6: cycles per microsecond (MHz)
    phase_*_rad: radians
    resulting H / hbar: angular frequency in rad / microsecond
    time_us: microseconds

The experimental parameters in `Evered2023` are a reduced effective three-level
model, NOT the full eight-level atom model in the publication. In particular,
this file currently implements only the HAMILTONIAN, not dissipation or leakage.

References:
- Bloqade.jl three-level conventions: https://queracomputing.github.io/Bloqade.jl/dev/3-level/
- Evered et al., Nature 622, 268-272 (2023), DOI:10.1038/s41586-023-06481-y
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

import numpy as np

TWO_PI = 2.0 * np.pi


@dataclass(frozen=True)
class Evered2023:
    """Published / approximate parameters for 87Rb, n=53, as labelled.

    omega_r_mhz: reported effective two-photon Rabi frequency Omega/2pi.
    blockade_mhz: approximate measured V/2pi at separation_um.
    rydberg_lifetime_us: assumed lifetime in the paper's model.
    C6 is INFERRED from V(R) at 2 um assuming isotropic V=C6/R^6.
    The paper's two-atom gate uses a more detailed 8-level error model.
    """
    omega_r_mhz: float = 4.6
    separation_um: float = 2.0
    blockade_mhz: float = 450.0
    rydberg_lifetime_us: float = 88.0
    branching_to_0: float = 0.053
    branching_to_1: float = 0.053
    branching_to_loss: float = 0.894
    pulse_amplitude_rad: float = TWO_PI * 0.1122
    pulse_omega_ratio: float = 1.0431
    pulse_phase_offset_rad: float = -0.7318
    dimensionless_pulse_time: float = 1.215  # Omega*T/(2pi)

    @property
    def inferred_c6_mhz_um6(self) -> float:
        return self.blockade_mhz * self.separation_um**6

    @property
    def gate_duration_us(self) -> float:
        return self.dimensionless_pulse_time / self.omega_r_mhz

    @property
    def decay_rate_per_us(self) -> float:
        return 1.0 / self.rydberg_lifetime_us

    def phase_r_rad(self, time_us: float) -> float:
        """Fixed-amplitude gate phase profile from Evered et al. (2023).

        phi(t) = A*cos(omega*t - phi0) + delta0*t, delta0=0.
        Applied with the Bloqade sign convention e^(+i phi)|1><r|.
        """
        omega_r_rad_per_us = TWO_PI * self.omega_r_mhz
        return float(self.pulse_amplitude_rad * np.cos(
            self.pulse_omega_ratio * omega_r_rad_per_us * time_us
            - self.pulse_phase_offset_rad
        ))


def _site_values(value: float | Sequence[float], n: int, name: str) -> np.ndarray:
    arr = np.asarray(value, dtype=float)
    if arr.ndim == 0:
        return np.full(n, float(arr), dtype=float)
    if arr.shape != (n,):
        raise ValueError(f"{name} must be scalar or have one value per atom ({n})")
    if not np.all(np.isfinite(arr)):
        raise ValueError(f"{name} contains non-finite values")
    return arr


def _site_operator(op: np.ndarray, i: int, n: int) -> np.ndarray:
    terms = [np.eye(3, dtype=complex) for _ in range(n)]
    terms[i] = op
    result = terms[0]
    for item in terms[1:]:
        result = np.kron(result, item)
    return result


def rydberg_hamiltonian_numpy(
    positions_um: Sequence[float] | Sequence[Sequence[float]],
    *,
    omega_r_mhz: float | Sequence[float] = 0.0,
    delta_r_mhz: float | Sequence[float] = 0.0,
    phase_r_rad: float | Sequence[float] = 0.0,
    omega_hf_mhz: float | Sequence[float] = 0.0,
    delta_hf_mhz: float | Sequence[float] = 0.0,
    phase_hf_rad: float | Sequence[float] = 0.0,
    c6_mhz_um6: float = 0.0,
    include_loss_level: bool = False,
) -> np.ndarray:
    """Construct H/hbar for N qutrits in angular rad/us (dense NumPy matrix).

    Matches Bloqade 3-level convention including the -Delta_hf on |r>:

      sum_i [ Omega_hf/2 * (e^iphi_hf|0><1| + h.c.)
              - Delta_hf * |1><1|
              + Omega_r/2 * (e^iphi_r|1><r| + h.c.)
              - (Delta_hf + Delta_r) * |r><r| ]
       + sum_{i<j} C6 / |r_i-r_j|^6 * n^r_i n^r_j.

    NOTE: all MHz values are converted to angular rad/us with factor 2*pi.
    This dense implementation is for small groups of atoms (e.g. 1-3), NOT
    a full BB code or a large many-qutrit Lindblad simulation.
    """
    positions = np.asarray(positions_um, dtype=float)
    if positions.ndim == 1:
        positions = positions.reshape(-1, 1)
    if positions.ndim != 2 or positions.shape[1] not in (1, 2, 3):
        raise ValueError("positions_um must be (N,), (N,1), (N,2), or (N,3)")
    n = len(positions)
    if not (1 <= n <= 5):
        raise ValueError("Use 1-5 atoms with this dense Hamiltonian implementation")
    if not np.all(np.isfinite(positions)):
        raise ValueError("positions must be finite")
    if not np.isfinite(c6_mhz_um6):
        raise ValueError("c6_mhz_um6 must be finite")

    wr = _site_values(omega_r_mhz, n, "omega_r_mhz")
    dr = _site_values(delta_r_mhz, n, "delta_r_mhz")
    pr = _site_values(phase_r_rad, n, "phase_r_rad")
    wh = _site_values(omega_hf_mhz, n, "omega_hf_mhz")
    dh = _site_values(delta_hf_mhz, n, "delta_hf_mhz")
    ph = _site_values(phase_hf_rad, n, "phase_hf_rad")

    local_dim = 4 if include_loss_level else 3
    e0 = np.eye(local_dim, dtype=complex)[:, 0]
    e1 = np.eye(local_dim, dtype=complex)[:, 1]
    er = np.eye(local_dim, dtype=complex)[:, 2]
    op01 = np.outer(e0, e1.conj())
    op1r = np.outer(e1, er.conj())
    n1 = np.outer(e1, e1.conj())
    nr = np.outer(er, er.conj())

    dim = local_dim**n
    H = np.zeros((dim, dim), dtype=complex)
    nr_sites = []
    for i in range(n):
        h_local = (
            (wh[i] / 2) * (np.exp(1j * ph[i]) * op01 + np.exp(-1j * ph[i]) * op01.conj().T)
            - dh[i] * n1
            + (wr[i] / 2) * (np.exp(1j * pr[i]) * op1r + np.exp(-1j * pr[i]) * op1r.conj().T)
            - (dh[i] + dr[i]) * nr
        )
        H += TWO_PI * _site_operator(h_local, i, n)
        nr_sites.append(_site_operator(nr, i, n))

    for i in range(n):
        for j in range(i + 1, n):
            R = float(np.linalg.norm(positions[i] - positions[j]))
            if R <= 0:
                raise ValueError("Atoms must occupy distinct positions")
            V_mhz = c6_mhz_um6 / R**6
            H += TWO_PI * V_mhz * (nr_sites[i] @ nr_sites[j])

    return H


def rydberg_hamiltonian(*args, **kwargs):
    """Return the same Hamiltonian as a QuTiP Qobj, with qutrit tensor dims."""
    try:
        import qutip as qt
    except ImportError as exc:
        raise ImportError("Install QuTiP to get a Qobj: pip install 'qutip>=5'") from exc
    matrix = rydberg_hamiltonian_numpy(*args, **kwargs)
    local_dim = 4 if kwargs.get('include_loss_level', False) else 3
    n = int(round(np.log(matrix.shape[0]) / np.log(local_dim)))
    return qt.Qobj(matrix, dims=[[local_dim] * n, [local_dim] * n], isherm=True)


def evered_two_atom_hamiltonian(time_us: float, params: Evered2023 | None = None, *, include_loss_level: bool = False):
    """QuTiP Hamiltonian for the *reduced* two-qutrit Evered-inspired gate pulse.

    No full gate benchmarking is implied. Missing intermediate levels, noise,
    loss, finite pulse turn-on/off, and additional electronic levels.
    """
    p = params if params is not None else Evered2023()
    return rydberg_hamiltonian(
        [0.0, p.separation_um],
        omega_r_mhz=p.omega_r_mhz,
        phase_r_rad=p.phase_r_rad(time_us),
        c6_mhz_um6=p.inferred_c6_mhz_um6,
        include_loss_level=include_loss_level,
    )


def evered_two_atom_hamiltonian_numpy(time_us: float, params: Evered2023 | None = None, *, include_loss_level: bool = False):
    """NumPy version for independent verification when QuTiP isn't installed."""
    p = params if params is not None else Evered2023()
    return rydberg_hamiltonian_numpy(
        [0.0, p.separation_um],
        omega_r_mhz=p.omega_r_mhz,
        phase_r_rad=p.phase_r_rad(time_us),
        c6_mhz_um6=p.inferred_c6_mhz_um6,
        include_loss_level=include_loss_level,
    )


if __name__ == "__main__":
    p = Evered2023()
    H = evered_two_atom_hamiltonian_numpy(0.0, p)
    print("Basis: |00>, |01>, |0r>, |10>, |11>, |1r>, |r0>, |r1>, |rr>")
    print("Shape of 2-atom H:", H.shape)
    print("Gate duration [us]:", p.gate_duration_us)
    print("Rabi Omega/2pi [MHz]:", p.omega_r_mhz)
    print("V_12/2pi [MHz]:", H[8, 8].real / TWO_PI)
    print("Hermiticity error:", np.max(np.abs(H - H.conj().T)))