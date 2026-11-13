# Rydberg-informed SAT: H-seeded and unseeded stabilizer codes

**Date:** 2026-11-13. **Status:** complete finite-size SAT feasibility and exact probability-weighted decoding **within stated assumptions**. It is **not** yet a scalable hardware-aware QEC benchmark or a circuit-level logical error simulation.

## Research question

Can the same hardware-derived Pauli error probabilities lead to different codes when a Rydberg Hamiltonian symmetry is required as a stabilizer? Compare **unseeded ε-code searches** and **Hamiltonian-seeded ε-code searches** under identical code size, correction requirements, and an exact maximum-likelihood Pauli decoder.

## Physical input and assumption

The input is the previously computed 16-string Pauli projection of a **two-atom Rydberg pulse** using the corrected Evered-inspired decay and intermediate-state-scattering model plus assumed Gaussian common detuning fluctuations. The sum of the two-qubit Pauli weights is **0.9984648622673**, with the rest in residual Rydberg and effective loss states. These 16 values are stored in `rydberg_pair_pauli_weights.json`.

For the four-qubit example, **assume two identical independent disjoint pair operations**, on (1,2) and (3,4):

\[
q_{P_1P_2P_3P_4}=q_{P_1P_2}\,q_{P_3P_4}.
\]

This provides all **256** Pauli-string weights, with total computational survival probability **0.99693208118245**. It is not a calibrated four-atom channel; possible cross-pair correlations and propagation of errors through subsequent gates are excluded. **Loss and residual Rydberg** remain outside the Pauli decoder; their combined modeled probability is about **0.306792%** for the two pairs, not a measured logical failure rate.

## The H-code seed is a genuine symmetry of the Rydberg-only Hamiltonian

For each **three-level** atom, define

\[
Q=|0\rangle\langle0|-|1\rangle\langle1|-|r\rangle\langle r|=\operatorname{diag}(1,-1,-1).
\]

During a Rydberg-only drive with **no hyperfine** \(|0\rangle\leftrightarrow|1\rangle\) pulse, the drive connects \(|1\rangle\) and \(|r\rangle\), which both have eigenvalue \(-1\) under \(Q\). Detunings and the Rydberg blockade are diagonal. Thus \([Q_i,H(t)]=0\) for each atom, at all times in this restricted pulse. We checked the commutator numerically against the Bloqade-convention qutrit Hamiltonian. The pair symmetry \(Q_1Q_2\) restricts to **\(ZZII\)** on the computational states, and \(Q_3Q_4\) restricts to **\(IIZZ\)**.

**Limitations:** Switching on a hyperfine drive breaks the symmetry. Scattering and Rydberg decay jumps can also break it. Furthermore, only the seeded generators are required to commute with \(H\); the other discovered stabilizers are **not** necessarily Hamiltonian symmetries. Therefore these are **H-seeded code-design** experiments, not a proof that the complete QEC code is dynamically invariant during the pulse.

## SAT constraints: full Pauli stabilizers, not X-only

With \(n=4,k=1\) we search for \(r=n-k=3\) independent, commuting Pauli generators \(g_1,g_2,g_3\). Each uses binary symplectic variables \((x_i|z_i)\), so arbitrary \(I,X,Y,Z\) entries are allowed.

**Hard constraints:**

1. The three generators commute pairwise: \(\langle g_i,g_j\rangle_s=0\).
2. Their GF(2) rank is exactly 3 (no dependent generator).
3. The distance is at least 2: every weight-1 Pauli either is a stabilizer or anticommutes with a generator; equivalently there is no weight-1 nontrivial logical Pauli.
4. For the strict **ε-code** mode, all errors in \(\mathcal A=\{IIII,ZIII,IZII,IIZI,IIIZ\}\) must be **correctable**: for every pair \(E_a,E_b\), their product (up to phase) must either belong to the stabilizer group or anticommute with at least one generator. This is the **Pauli Knill–Laflamme** condition including degeneracy.
5. In seeded runs, require \(ZZII\), or both \(ZZII\) and \(IIZZ\), as generators. Since independent generators can be changed within their span, the exact-search filter checks **group membership** rather than a particular basis.

We constructed SMT-LIB2 **Boolean** satisfiability formulas and solved them with Z3 through its native library. The same finite-size search was **independently verified by enumerating all 11,475 commuting independent rank-3 stabilizer groups** on four qubits. This makes the optimum certified for the chosen \([[4,1]]\) search, rather than a solver heuristic.

## What the probabilities optimize

The SAT solver establishes which codes satisfy the hard constraints. Among *all* feasible codes, rank by **exact maximum-likelihood decoder success** under all 256 Pauli weights. For each syndrome, errors are grouped into stabilizer-equivalent classes; the decoder chooses the class with largest total probability. This includes degeneracy and avoids treating mere detection as correction:

\[
p_{\rm succ}=\sum_{s}\max_{L\text{ class}}\sum_{E\in(s,L)}q_E,\qquad
p_{\rm fail\mid survive}=\frac{p_{\rm survive}-p_{\rm succ}}{p_{\rm survive}}.
\]

No stochastic sampling error: every Pauli string is summed exactly. No assumption that the Pauli-twirled model is an exact representation of the underlying coherent/leakage channel.

## Results A: strict ε-code (correct every single-qubit Z; d≥2)

| Seed constraint | All stabilizer groups containing seed | KL + distance feasible | SAT result | Min. conditional logical error |
|:--|--:|--:|:--|--:|
| None | 11,475 | **768** | SAT | **0.038588956%** |
| \(ZZII\) | 315 | **8** | SAT | **0.038590422%** |
| \(ZZII\) and \(IIZZ\) | 15 | **0** | UNSAT | N/A |

- **Unseeded optimum:** \(S=\langle XXII,YZXX,IIZY\rangle\), **\([[4,1,2]]\)**.
- **Single-H-seed optimum:** \(S=\langle IIXX,ZZII,XXYZ\rangle\), **\([[4,1,2]]\)**.
- **Both seeds:** no code can satisfy the above constraints for this \([[4,1]]\) problem.

The unseeded code wins by **0.000001466 percentage points**, a very small difference in this simplified model. This experiment shows **no improvement from one H seed**. Requiring both H seeds is **incompatible with correction of all single-Z errors** under the target size/distance.

## Results B: probability-prioritized code (d≥2, no mandatory correction of all single-Z errors)

| Seed constraint | d≥2 candidates | SAT result | Min. conditional logical error |
|:--|--:|:--|--:|
| None | 2,268 | SAT | **0.037866637%** |
| \(ZZII\) | 36 | SAT | **0.037868033%** |
| Both seeds | 4 | SAT | **0.064449789%** |

The weighted-only mode can outperform the strict ε-code mode slightly by **not guaranteeing correction of every single-Z error**. Both-seed codes exist in this relaxed problem, but their optimum under our noise is considerably worse. Do **not** compare these modes as if they impose the same correction guarantees. Within each mode, seeding versus not seeding is matched.

## Important correction to the earlier preliminary trial

An earlier \([[4,1]]\) code with conditional error **0.0341%** had **distance 1**. It was a phase-error-specific code, not a code with guaranteed detection of arbitrary single-qubit errors. It is **superseded** by the distance-2 results above; it must not be used as the main result.

## What the results do and do not show

**Established here:** A genuine general-Pauli Boolean SAT/SMT feasibility formulation; exact finite-size certification; optional Hamiltonian seed constraints; Pauli Knill–Laflamme correction requirements; distance-2 enforcement; and probability-based selection using all 256 errors and optimal syndrome decoding.

**Not established:** a physical implementation of syndrome extraction, gates or decoders; leakage correction; code performance under the original non-Pauli Kraus channel; cross-pair correlated noise; an advantage over a hardware-unaware baseline; scalability beyond this finite search; or full invariance of the encoded code under the driven Hamiltonian. Hardware-unaware matched comparisons and logical-error-versus-noise-strength studies remain next experiments.

## Reproduce

Install Python with NumPy. From the `Results` directory, run `python rydberg_sat_study.py` to enumerate and optimize; if the native Z3 library is present, the SAT cases are solved and cross-checked. Run `python -m unittest -v test_rydberg_sat_study.py` for regression checks. SMT-LIB2 cases are included and may also be solved with a standalone Z3 executable. All input Pauli weights and numerical outputs are committed.

**Related sources:** DalFavero et al., *Satisfying Quantum Codes: Physics-Informed and Hardware-Aware Code Design with SAT Solvers*, arXiv:2608.23460 (2026); Evered et al., *Nature* **622**, 268–272 (2023); Bloqade.jl, *Hamiltonian in a 3-level system*, https://queracomputing.github.io/Bloqade.jl/dev/3-level/.