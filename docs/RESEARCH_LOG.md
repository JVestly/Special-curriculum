# Special curriculum research log

## 2026-11-10 — Initial repository and symmetry checks

Repo initially contained `PhysicsInformedSAT_Solver.py` and `rydberg_interaction.ipynb`. The solver implements binary symplectic functions, Pauli-string commutation, a GF(2) nullspace and Pauli Hamiltonian symmetry search; `solve()` still raises `NotImplementedError`.

An equivalent five-test direct-matrix check (numpy, independent of SAT) passed locally for LMG population parity, native three-level two-drive parity breaking, drive-off restoration, parity dependence, and restricted binary-register parity representatives. The tests were added under `tests/test_symmetry_models.py` for verification on the branch.

**These tests do not demonstrate a physical LMG realization on neutral-atom hardware or a valid H-seeded code.**

### Next actions
Define precisely what operators may be used as stabilizer seeds on the claimed physical Hilbert space; then implement and independently verify SAT discovery for a known error-correcting code.

Thesis schedule and compiler-specific tasks are managed separately from this paper.
