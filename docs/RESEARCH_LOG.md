# Research log

## 2026-11-10 — Actual first session
Repository audit: the root contains PhysicsInformedSAT_Solver.py and rydberg_interaction.ipynb. The solver has binary symplectic utilities, Pauli commutation, a GF(2) nullspace and symmetry search; its solve() method currently raises NotImplementedError. No standalone tests were present in the root listing before this PR.

**Independent small-matrix check:** five unittest test methods were executed locally with numpy for 2/3-atom three-level LMG and native driven-NA toy Hamiltonians. All five passed. Tests assess population parity commuting for LMG, broken parity with two native drives, drive-off restorations, parity dependence and restricted binary parity representations. This does NOT establish a NA→LMG reduction or a valid neutral-atom H-code.

Scientific scope: thesis primary all weekdays (15–20h/week), paper 3 fixed sessions/week (5–8h) ending Dec 18; same tested core contributes to both. The project charter and detailed schedule are proposed in the planning PR, not merged into main.

Next experiment: separately verify the physical-subspace encoding and exact correctability formulation before treating any candidate G as a physical QECC stabilizer.