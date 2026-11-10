# Project charter — review draft (2026-11-10)

## Priority and purpose
This is a **thesis-led research programme** with a separate special-curriculum paper due Friday 18 December 2026. Thesis work continues every weekday; the paper receives **three sessions each week (Tuesday, Thursday, Friday)**, 90–150 minutes each, and an optional fourth only for blockers. Planned load: thesis 15–20 focused hours/week; paper 5–8. Shared verified code is reused, but the two projects have different scientific claims.

### Thesis objective
Build a *hardware-agnostic mathematical core* for Hamiltonian-informed stabilizer-code discovery, then validate a **restricted** compilation/QEC pipeline for arbitrary *supported* small logical Clifford inputs. Hardware-agnostic means a common interface works with at least two different Hamiltonian models; it does **not** mean arbitrary Hamiltonians yield useful QEC codes or that every device supports the same physical gates.

### Paper objective
For an identical (n,k,error set), test whether **fixed Hamiltonian-derived symmetry seeds** change SAT-code feasibility, search properties and code characteristics compared with unseeded SAT. Include an experimentally informed neutral-atom Pauli error example with labeled approximations. A three-level LMG model is an algebraic benchmark until a physical NA reduction is demonstrated.

## Precise interfaces
- HamiltonianModel: physical Hilbert space, Pauli decomposition or specified embedding, coefficients and units, assumptions, physical-subspace projector.
- PhysicalEncoding: maps physical degrees of freedom to binary registers, including unused states. **Three qutrits encoded in six registers are not six physical neutral-atom qubits.**
- NoiseModel: error operators, frequencies if known, physical provenance, application level (data/gate/measurement); do not disguise loss or leakage as Pauli errors.
- CodeRequest: n, k, correctable error set, optional *signed, commuting independent* Hamiltonian seeds, explicit constraints.
- Code: stabilizers, logical operators, independent correctness report; the SAT output is a code, **not** a quantum encoder.
- LogicalProgram: general-input supported Clifford instructions, clear unsupported-gate errors.
- HardwareBackend: computational-qubit definition, native gates, allowed interactions/connectivity, measurement, reset, gate noise and costs. A Hamiltonian alone is not a hardware backend.
- CompilationResult: encoding circuit, supported logical circuit, syndrome circuit, decoder/recovery, numerical verification and costs.

## Mathematical standards
1. Individual Pauli symmetries of an expanded H are found from the GF(2) symplectic nullspace, then *commuting independent* seeds are selected. A nullspace basis need not commute mutually.
2. A Hamiltonian-*seeded* code does not guarantee all SAT-generated stabilizers commute with H. Do not assert dynamical codespace invariance without independently checking [H,Pi_code]=0.
3. Correctability tests use all **pair products E_a† E_b**, accounting for degeneracy. Signs and -I exclusion are explicit.
4. Symmetry operators found in a qutrit embedding may be valid only on the physical subspace; don't report virtual code parameters as physical-qubit parameters.
5. Every simulation has a fixed configuration, reference, solver version and independent mathematical verifier. Timeouts/UNSAT and negative results are recorded without embellishment.

## Exit criteria
**Paper, Dec 18:** one working and independently verified SAT code discovery, fixed-condition seeded/unseeded comparison, documented NA-inspired Pauli noise assumptions, reproducible results, coauthored manuscript. H-seeding is a benchmark if not physically realizable.

**Thesis, spring 2027 (submission date TBC):** two distinct Hamiltonian examples through the same mathematical core, supported general-input logical circuit, correct small Clifford encoder for arbitrary input, syndrome and one recovery cycle, a decoder, at least one declared physical gate-level mapping, measured resources, full tests. NOT a universal fault-tolerant compiler.

## Gate decisions
- Nov 13: contracts and small-matrix symmetry checks verified.
- Nov 20: first SAT-discovered known code or scope reduced to an independently validated tiny example.
- Nov 27: paper comparison defined and initial seeded/unseeded results recorded.
- Dec 4: freeze paper experiments.
- Dec 11: full coauthor draft.
- Dec 18: submission and reproduction snapshot.

## Excluded until baseline works
Pulse-level Lindblad, universal physical-device compilation, arbitrary logical gates, full qutrit QECC search, native hardware benchmarking, unverified NA-to-LMG equivalence, and exhaustive SAT-solver competitions.