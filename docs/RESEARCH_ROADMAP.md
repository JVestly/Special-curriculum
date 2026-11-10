# Daily dual-track plan: 10 November–18 December 2026

**Priority:** master's thesis every weekday, 15–20 focused hours/week. Paper only Tue/Thu/Fri 90–150 minutes each; fourth session is optional, not expected. Separate paper and thesis status should be recorded at end of each session.

**Completion rule:** each task results in a committed note, code/test artifact, explicit model decision or documented blocker with reproducer. No plan entry is automatically "completed" because its date passed.

| Date | Thesis main block | Paper additional block |
|---|---|---|
| Tue 10 Nov | Thesis: define MVP, interfaces and hard limits | Paper: start 3-level LMG vs native NA parity derivation; run independent tests |
| Wed 11 Nov | Thesis: HamiltonianModel/NoiseModel/PhysicalEncoding contracts and examples | — |
| Thu 12 Nov | Thesis: GF(2) and direct-matrix test environment | Paper: determine precisely which parity operators could seed a physical-qubit code |
| Fri 13 Nov | Thesis: arbitrary-input small Clifford LogicalProgram contract | Paper: finish qutrit/native comparison tests; document limitations |
| Mon 16 Nov | Thesis: typed Pauli/stabilizer objects + tests | — |
| Tue 17 Nov | Thesis: define CodeRequest/Code interface | Paper: SAT variables and commutation clauses |
| Wed 18 Nov | Thesis: derive logical X/Z and independent code verifier | — |
| Thu 19 Nov | Thesis: specify arbitrary logical-state encoder checks | Paper: rank/independence and all pair-product error conditions |
| Fri 20 Nov | Thesis: connect known code to logical operators | Paper: discover first known code; [[5,1,3]] stretch |
| Mon 23 Nov | Thesis: encoder prototype using a known code | — |
| Tue 24 Nov | Thesis: decoder interface and syndrome conventions | Paper: NA-inspired error-model evidence/provenance |
| Wed 25 Nov | Thesis: verify encoding on |0>,|1>,|+> and generic state | — |
| Thu 26 Nov | Thesis: second Hamiltonian adapter with same interface | Paper: verify H-seed eligibility and parameter assumptions |
| Fri 27 Nov | Thesis: integrate code object and decoder contract | Paper: fixed n,k,error seeded vs unseeded comparison |
| Mon 30 Nov | Thesis: first encode→error→syndrome→recovery demonstration | — |
| Tue 1 Dec | Thesis: stable shared mathematical API | Paper: benchmark runs/configurations/raw results |
| Wed 2 Dec | Thesis: second model accepted without core rewrites | — |
| Thu 3 Dec | Thesis: decoder logical-failure verification harness | Paper: ablation, SAT/UNSAT and timeout table |
| Fri 4 Dec | Thesis: write draft thesis architecture Methods | Paper: freeze experiments; figures and provenance |
| Mon 7 Dec | Thesis: two distinct supported Clifford programs | — |
| Tue 8 Dec | Thesis: syndrome-measurement circuit prototype | Paper: Methods and model-limitations section |
| Wed 9 Dec | Thesis: HardwareBackend gates/connectivity/cost spec | — |
| Thu 10 Dec | Thesis: decoder ↔ syndrome-bit order/Pauli frame | Paper: final figures, captions and reproducibility review |
| Fri 11 Dec | Thesis: circuits/decoder integration test | Paper: complete full coauthor manuscript draft |
| Mon 14 Dec | Thesis: document one working end-to-end demonstration | — |
| Tue 15 Dec | Thesis: native gate compilation/routing interface | Paper: mathematical proof and citation audit |
| Wed 16 Dec | Thesis: error paths/unsupported input tests | — |
| Thu 17 Dec | Thesis: reproduce representative thesis benchmark | Paper: build, coauthor review and final corrections |
| Fri 18 Dec | Thesis: end-of-term progress and January backlog | Paper: submit and archive reproduction snapshot |

## Weekly gates (paper)
- Week Nov 10–13: physical qutrit symmetries independently tested, seed limitations written.
- Nov 16–20: independently verified tiny SAT-discovered QECC, or explicit blocker and reduced task.
- Nov 23–27: comparable H-seeded/unseeded SAT experiments under identical error constraints.
- Nov 30–Dec 4: freeze experimental settings and record raw data with reproducible scripts.
- Dec 7–11: full coauthor draft, figs and source checks.
- Dec 14–18: corrections only; submitted manuscript and archived run configuration.

## Weekly/monthly thesis plan
- **December remainder:** demonstrated basic encoding/recovery and general logical-input contract; second Hamiltonian adapter.
- **January 2027:** W1 freeze interfaces; W2 synthesize verified Clifford encoder; W3 two general-input logical circuits; W4 integration tests/first performance study.
- **February:** W1 stabilizer measurement circuits; W2 coset-aware decoder; W3 complete error correction cycle; W4 code-capacity and first noisy circuit tests.
- **March:** W1 choose supported gate backend; W2 native gate mapping and resource accounting; W3 circuit-level noisy comparison; W4 second Hamiltonian-family validation and limitations.
- **April:** W1 frozen benchmarks; W2 ablations; W3 novelty/limitations write-up; W4 full draft and independent verification.
- **May:** thesis figures, review, reproducibility and defense demo. Adjust once the actual university submission date is confirmed.

## Friday review questions
What was actually verified this week? What changed in the thesis and paper deliverables? What is still conjecture? How many hours did each track take? What is the first reproducible step for Monday?