# Intermediate Rydberg noise results — 13 November 2026

These are **intermediate numerical results**, not experimental gate fidelities or a complete hardware noise budget. Code: [rydberg_noise_diagnostics.py](rydberg_noise_diagnostics.py). The previous Hamiltonian, three-branch Lindblad noise model, and Choi/Kraus/Pauli extraction remain intact.

## Physics and assumptions

Reduced two-atom model with local |0>, |1>, |r>, |L> levels. Both atoms have three independently modeled Rydberg decay destinations |r>→|0>, |r>→|1>, |r>→|L>. We take total Rydberg lifetime T_r=88 μs and branching fractions (0.053, 0.053, 0.894) from Evered et al., *Nature* **622**, 268–272 (2023), https://doi.org/10.1038/s41586-023-06481-y. Gate duration ≈0.264 μs, with a reduced Rydberg Hamiltonian and phase-shaped drive.

The publication reports ground–Rydberg T₂* ≈3 μs. Our optional extension **assumes**, rather than experimentally establishes, *common-mode quasi-static Gaussian detuning* shared by both atoms. We set the angular-frequency spread sigma=sqrt(2)/T₂*, giving a Gaussian free-induction envelope exp[-(t/T₂*)²], and average Lindblad channels over a five-point Gauss–Hermite quadrature. This is **not** implemented by adding a Lindblad T₂* collapse operator. Other spatial correlation assumptions would produce different correlated Pauli weights.

## Diagnostic outcomes with three decay branches

For one pulse, starting from pure basis states, the probability of at least one atom reaching the effective out-of-manifold |L> sink is:

| Initial state | P(≥1 in L) | P(residual r) | P(both L) |
|---|---:|---:|---:|
| |00> | 0.000000% | 0.000000% | 0.000000% |
| |01> | 0.137063% | 0.005876% | 0.000000% |
| |10> | 0.137063% | 0.005876% | 0.000000% |
| |11> | 0.139021% | 0.017547% | 0.000086% |
| |rr> | 0.535078% | 99.440643% | 0.000382% |

|rr> is a useful diagnostic but *not* a computational input. 89.4% is the *conditional decay branching fraction*, not the total loss probability during a pulse.

## State-independent channel results

For the 4D computational input channel to the 16D physical output, derive Choi (64×64), rectangular Kraus operators (16×4), then Pauli-twirl the computational restriction after removing the nearest noiseless computational unitary. This reference is not independently validated as an ideal CZ, and the effective Pauli weights do not fully preserve coherent noise or leakage. The output-sector probabilities below are evaluated on the maximally mixed computational input.

| Probability | Decay only | Decay + optional T₂* detuning model |
|---|---:|---:|
| II | 99.883448% | 99.790875% |
| All nonidentity Paulis | 0.005941% | 0.082488% |
| Residual r | 0.007325% | 0.023364% |
| Effective L sink | 0.103287% | 0.103273% |
| Total | ~100% | ~100% |

Selected Pauli-twirled weights:

| String | Decay only | Decay + T₂* detuning |
|---|---:|---:|
| IZ | 0.000487% | 0.026487% |
| ZI | 0.000487% | 0.026487% |
| ZZ | 0.000105% | 0.024653% |
| IX | 0.000609% | 0.000609% |

The nonidentity Pauli weight increases ≈13.9 times under this model, whereas Rydberg sink loss barely changes. This demonstrates sensitivity to *modeled coherence noise* and its spatial correlations, **not** an experimentally validated hardware advantage or a logical-error-rate result. Normalization is satisfied within ~10⁻¹¹ and the Choi matrix is positive up to floating-point tolerance.

## Interpretation and next steps

- The same Rydberg decay rates produce different state-dependent loss probabilities because Rydberg occupation depends on the initial state and pulse.
- Common-mode detuning increases IZ, ZI, and ZZ weights substantially; the ZZ growth reflects an **assumption of correlated detuning**, which needs to be tested against independent-site fluctuations.
- To approach the actual Evered et al. hardware noise model, add intermediate-level scattering using explicit intermediate levels and measured parameters. The publication's more complete calculation contains eight atomic levels per atom, so this reduced model must not be presented as a replication.
- For SAT-based code discovery, evaluate uniform-noise vs physics-informed error sets under the *same* independent evaluation model and treat the |L> sink separately from ordinary Pauli correction.

This file documents interim code and observations for the special-curriculum repository, without integrating the Lindblad derivation into the paper's core SAT claims.
