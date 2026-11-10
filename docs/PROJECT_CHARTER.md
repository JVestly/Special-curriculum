# Special curriculum paper — research specification

**Submission deadline:** Friday, 18 December 2026.
**Scope:** this repository contains the paper's experiments, SAT code-search implementation, mathematical checks, and manuscript support. The separate master's thesis has its own repository and task list. The paper receives **three dedicated sessions per week** (Tuesday, Thursday, Friday), with an optional fourth only if necessary.

## Research question

For a fixed number of physical/encoded registers, logical dimension, and declared error set, how do **Hamiltonian-derived Pauli symmetry constraints** affect SAT stabilizer-code feasibility, runtime and code properties compared with the same search without imposed symmetry seeds?

This is a study of **symmetry-informed code discovery**, not proof of a universal quantum compiler and not a claim that Pauli symmetries automatically give physically realizable neutral-atom QEC codes.

## Code-discovery interface

```
CodeRequest(n, k, pauli_errors, initial_stabilizers, constraints)
  -> CodeResult(stabilizers, status, correctness_report, solver_metadata)
```

- SAT constructs mutually commuting, independent signed stabilizer generators; does not synthesize encoding circuits.
- Every pair E_a†E_b of target errors must satisfy the correctability condition, including degenerate cases.
- Symmetry candidates are drawn from the Pauli-string commutant of a specified Hamiltonian; a GF(2) nullspace basis need not be mutually commuting.
- Commutation with the Hamiltonian for a few seed stabilizers does **not** mean the complete selected codespace is invariant under its dynamics.
- A three-level/qutrit model can be represented in a larger binary Hilbert space for computation. That does **not** add physical neutral-atom qubits. Distinguish physical qutrit constraints from virtual-register code parameters.

## Experimental comparison

1. Verify the Hamiltonian symmetry conditions using analytical identities and small physical-qutrit matrices.
2. Discover and independently validate a known small stabilizer code as the SAT baseline.
3. Under fixed (n,k,errors,solver settings), compare SAT search *without seeds* and *with justified compatible seeds*. Record SAT/UNSAT/timeout, candidate properties and runtime.
4. Define a small, documented neutral-atom-inspired Pauli error set/channel using published information, explicitly excluding or separately modeling atom loss and leakage.
5. Produce reproducible plots/tables. A result showing that seeding makes the problem infeasible or gives no advantage is valid if reported accurately.

## Models and physical interpretation

- Three-level neutral-atom Hamiltonian: population dynamics driven by 0↔1 and 1↔r transitions plus Rydberg interactions.
- Generalized three-level LMG Hamiltonian: pair-transition benchmark with conserved level-population parities.
- **No physical neutral-atom → three-level LMG reduction has been demonstrated.** Their mathematical results must remain separate until such a derivation is established.

## Success by December 18

A reproducible SAT-discovered known code and independent algebraic verification; well-defined seeded/unseeded experiment(s) and limitations; documented model assumptions and at least one paper-quality comparison; manuscript reviewed with coauthor. No physical compiler, decoding circuit or thesis-specific deliverables are required for this paper.

The paper's mathematical utilities can be reused later, but implementation and documentation for the master's compiler belong in its own repository.
