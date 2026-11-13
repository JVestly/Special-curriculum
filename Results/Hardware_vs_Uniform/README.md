# Hardware-informed vs uniform Pauli SAT comparison

Date: 13 November 2026. This is an intermediate numerical **[[4,1,2]] code-capacity study**, not experimental QEC.

See [results_Hardware_vs_Uniform.ipynb](../results_Hardware_vs_Uniform.ipynb) for assumptions, formulas, reproducible code cells, plots, and the results table.

| Design | Generators | Conditional error under the Rydberg model |
|---|---|---:|
| Rydberg-aware | `XXXX, ZYII, XXYZ` | 0.03786664% |
| Rydberg-aware with `ZZII` H-seed | `XXXX, ZZII, XXYZ` | 0.03786803% |
| Uniform 15-error two-qubit Pauli baseline, no H-seed | `XXII, IIXX, ZZZZ` | 0.11681243% |

The 67.6% reduction applies **only against this specified uninformed representative**. The uniform model has **324 exactly tied optimal codes**: their true Rydberg logical-error rates range from **0.03786664% to 0.11681243%**, median **0.03859042%**. **64 uniform-optimal codes are also Rydberg-optimal.** Thus there is no unconditional guarantee of improvement; physical probabilities provide a principled way to distinguish codes that uniform noise cannot rank.

Both design models match the same total Pauli error rate and computational survival, and both are evaluated on the same Rydberg Pauli model on two independent disjoint atom pairs. Every candidate has distance `>=2`; no hardware-specific mandatory single-Z correction constraints are imposed in the primary comparison. SAT/SMT checks validity; exhaustive enumeration of 2,268 distance-2 codes plus exact syndrome maximum-likelihood calculation selects optima. This is **not weighted MaxSAT**.

The total modeled probability of leaving the computational subspace across both pairs is **0.306792%**. Counting all such outcomes as failures (without attempting erasure correction) yields **0.344542%** versus **0.423246%**, a **pessimistic proxy**, not measured leakage-aware logical failure. A separate scaled-noise sweep does **not** establish a pseudothreshold.

**Reproduction:** run all cells in the notebook with NumPy, Matplotlib, and the sibling `Results/rydberg_sat_study.py` and `Results/rydberg_pair_pauli_weights.json`. The original SAT study and physical model are documented in `Results/SAT_Rydberg/README.md`.
