# Special curriculum paper — session plan to 18 December 2026

**Status:** proposed schedule. Three paper sessions a week, normally **Tuesday, Thursday and Friday**, 90–150 minutes each (approximately 5–8 focused hours/week). Thesis work is managed in the separate master's repository.

| Date | Main paper task | Definition of done |
|---|---|---|
| Tue 10 Nov | Derive/test three-level LMG and native Rydberg parity results | Small-matrix commutators reproduced; limits recorded |
| Thu 12 Nov | Evaluate symmetry generators as legitimate SAT seeds | Physical qutrit vs virtual binary representation explicitly distinguished |
| Fri 13 Nov | Finalize Hamiltonian symmetry derivation and tests | Hermiticity, parity relations, documented counterexamples |
| Tue 17 Nov | SAT variables and generator commutation | Boolean formulation and tests |
| Thu 19 Nov | Stabilizer independence and error-pair correctability | Independent checker for Knill–Laflamme conditions |
| Fri 20 Nov | First reproducibly discovered small code | Known correctable code verified; UNSAT test; blockers recorded |
| Tue 24 Nov | Document neutral-atom error mechanisms | Traceable source parameters and explicit approximations |
| Thu 26 Nov | Validate symmetry seed compatibility | Signed commuting independent seeds; allowed physical space |
| Fri 27 Nov | Fixed seeded/unseeded SAT comparison | Initial table at identical n,k,error criteria, solver settings |
| Tue 1 Dec | Run controlled experiments | Reproducible parameters and raw results |
| Thu 3 Dec | Ablations, negative cases and runtime evidence | SAT/UNSAT/timeout and quantitative comparison |
| Fri 4 Dec | Freeze paper experimental settings | Scripts/configuration locked, figures drafted |
| Tue 8 Dec | Write method and assumptions | Mathematical formulation and evidence documented |
| Thu 10 Dec | Verify data, captions and limitations | Final figures and tables reproducible |
| Fri 11 Dec | Complete full manuscript draft | Draft to coauthor for review |
| Tue 15 Dec | Mathematical/citation audit | Independent verification and source checks |
| Thu 17 Dec | Incorporate feedback | Submission-ready paper, no new unverified claims |
| Fri 18 Dec | Submit and archive experiments | Final manuscript and reproducibility snapshot |

## Weekly checkpoints
- **Nov 10–13:** correct physical interpretation of LMG and Rydberg symmetry calculations.
- **Nov 16–20:** functioning SAT code discovery and independent correctness checks.
- **Nov 23–27:** documented physical error assumptions and comparable H-seeded/unseeded runs.
- **Nov 30–Dec 4:** completed controlled experiments and frozen configuration.
- **Dec 7–11:** coauthor draft.
- **Dec 14–18:** review, validation and submission.

### Decision points
- If a candidate symmetry cannot seed a valid QEC code on the claimed physical degrees of freedom, report it as a mathematical benchmark rather than claiming physical realizability.
- If full single-qubit-Pauli correction is too hard to formulate in time, first demonstrate correct small error-specific code design; distinguish its distance/coverage honestly.
- If seeded SAT is slower, UNSAT, or not better, retain that negative result rather than changing experimental conditions to create a positive one.

### Friday review
List observed results, mathematical assumptions, tests passed, open failures, elapsed paper hours and next week's first reproducible task.
