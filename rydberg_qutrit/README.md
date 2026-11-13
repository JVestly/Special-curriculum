# Three-level Rydberg Hamiltonian — thesis prototype

The first physical-model module of the hardware-informed QEC compiler.

## Install

```bash
pip install numpy scipy 'qutip>=5'
```

## Example

```python
from rydberg_qutrit import Evered2023, evered_two_atom_hamiltonian, rydberg_hamiltonian
p = Evered2023()
H = evered_two_atom_hamiltonian(0.0)
print(H.shape)                # (9, 9)
print(p.gate_duration_us)     # about 0.264 microseconds

# Static two-atom Hamiltonian, zero drive phases:
H_static = rydberg_hamiltonian(
    [0, 2], omega_r_mhz=4.6, c6_mhz_um6=p.inferred_c6_mhz_um6
)

# A time-dependent physical drive, not a logical program:
import qutip as qt
H_t = qt.QobjEvo(lambda t: evered_two_atom_hamiltonian(t, p))
```

The Hamiltonian implements Bloqade's 3-level convention:

\[
H/\hbar = \sum_i \left[ \frac{\Omega_i^{hf}}2(e^{i\phi_i^{hf}}|0_i\rangle\langle1_i|+h.c.)
-\Delta_i^{hf}|1_i\rangle\langle1_i| +
\frac{\Omega_i^r}2(e^{i\phi_i^r}|1_i\rangle\langle r_i|+h.c.)
-(\Delta_i^{hf}+\Delta_i^r)|r_i\rangle\langle r_i|\right]
+\sum_{i<j}(C_6/R_{ij}^6)|r_ir_j\rangle\langle r_ir_j|.
\]

Values entered in MHz are converted to angular rad/us by 2*pi; t is in us.

## Experimental grounding and limitations

Parameter values sourced from Evered et al., Nature (2023):
- Effective 2-photon Rabi drive \(\Omega_r/2\pi=4.6\) MHz.
- Approximate blockade \(V/2\pi=450\) MHz at \(R=2\) um.
- Fixed-amplitude two-qubit pulse \(\phi(t)=A\cos(\omega t-\phi_0)\), \(A=2\pi(0.1122)\), \(\omega=1.0431\Omega\), \(\phi_0=-0.7318\), \(\Omega T/2\pi=1.215\).
- Rydberg lifetime of 88 us and decay branching fractions 0.053 to |0>, 0.053 to |1>, 0.894 to other levels (for future Lindblad work).

The inferred value \(C_6/2\pi \approx 28800\) MHz·um^6 is only an *effective* van der Waals approximation inferred from the blockade at 2 um. It is not the Bloqade default value (which corresponds to different atomic choices).

**This does not reproduce the published gate fidelity.** It omits intermediate excited states, other Rydberg levels, scattering, dephasing, nonideal turn-on/off and leakage. Most measured decay branches leave the 3-level model; a trace-preserving Lindblad model must include an explicit sink level or a documented effective treatment. Do not incorrectly renormalize the 0.053+0.053 branches to 1.

## Verification

```bash
python -m unittest -v test_rydberg_qutrit.py
```

The tests run without QuTiP (the NumPy Hamiltonian is the source of truth). The QuTiP wrapper requires QuTiP 5.

## Next

Introduce Rydberg decay with correct branching (including the leakage sink), then construct the Lindblad channel and derive Pauli weights after projecting into the 2-level computational subspace.

References: https://queracomputing.github.io/Bloqade.jl/dev/3-level/ ; https://doi.org/10.1038/s41586-023-06481-y

## Full three-branch Rydberg-decay channel (intermediate results: 2026-11-13)

Added `rydberg_full_channel.py`, `test_rydberg_full_channel.py`, and
`demo_full_channel.py`. Reuses the Hamiltonian above with **all** Rydberg
relaxation branching fractions from Evered et al. (2023):

- `r -> 0`: 0.053 of the total 1/(88 us) decay rate;
- `r -> 1`: 0.053;
- `r -> L`: 0.894.

Six site-local Lindblad collapse operators (three outcomes on each of two
atoms) implement a **trace-preserving** 16D physical evolution. This is
complete **only with respect to the three chosen Rydberg decay branches**;
intermediate-state scattering, laser noise, atom temperature, SPAM, crosstalk
and other hardware imperfections are **not** included.

The map under study is a rectangular CPTP channel from the four-dimensional
computational INPUT to the 16-dimensional physical OUTPUT. To construct it
without selecting a particular initial state, we evolve all 16 input matrix
units `|a><b|`; the resulting Choi matrix is 64 x 64. Diagonalizing the Choi
matrix produces rectangular Kraus operators of shape (16,4).

The reduced computational Kraus operators are `V^dag K_j` (shape 4x4). For a
physical-operation error model we left-compose with the inverse of the
**closest unitary to the noiseless computational block**. This reference
unitary is a computational convention, not a validation against the published
experimental CZ gate. It can hide coherent control errors relative to a
specified CZ target. The Pauli-twirled weights are

    q_P = (1/16) sum_j |Tr[P^dag U_ref^dag V^dag K_j]|^2.

These weights are *unconditional* with respect to leakage, i.e. they sum to
computational-subspace survival rather than one. A self-consistent partition
for the input rho=I/4 is

    sum_P q_P + p_residual_Rydberg + p_loss_sink = 1.

At the pulse duration 0.2641304348 us, the first numerical result is:

| Outcome | Probability |
|---|---:|
| Computational Pauli weights (sum) | 99.8893886931% |
| Rydberg residual (no loss) | 0.0073246994% |
| Sink-level population (>=1 L) | 0.1032866072% |
| Total | 99.9999999997% |

Within the computational contribution, `q_II = 99.8834478261%` and
`sum_(P != II) q_P = 0.0059408671%`. The computational nonidentity Pauli
weights are strongly nonuniform; examples include `q_IX = 0.0006089318%`,
`q_IZ = 0.0004873800%`, and `q_ZZ = 0.0001050747%`.

Checks: 15 unit tests pass, including a constructed isometric identity channel,
Choi-to-Kraus reconstruction and numerical physical-channel CPTP checks.
Maximum trace-preservation error is approximately 1.1e-11; minimum Choi
eigenvalue approximately -8.1e-16 (within roundoff). The computed noiseless
projection residual is 0.007207%. This is a numerical diagnostic only.

**Interpretation:** leakage is the largest modeled departure from the
computational space. Pauli-only error correction would miss most of the
modeled departures. An erasure/leakage-aware decoder may be relevant, but
whether loss is detectable must be specified rather than assumed.

To reproduce:

```bash
pip install numpy scipy matplotlib
cd rydberg_qutrit
OPENBLAS_NUM_THREADS=1 python -m unittest discover -s . -p 'test_*.py' -v
OPENBLAS_NUM_THREADS=1 python demo_full_channel.py
```

Results are written to `channel_results.json` and
`rydberg_full_channel_probabilities.png` (also SVG). The older
`rydberg_loss_probability.png` and `rydberg_occupation.png` remain useful
intermediate single-input diagnostics, not state-independent error models.