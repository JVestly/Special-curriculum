# Agent-assisted research rules

Agents accelerate implementation and literature audits. **Scientific assumptions and conclusions require human review.** Each task specifies mathematics, input/output contract, files in scope, expected independent tests and a falsifiable definition of done.

1. No broad "build a quantum compiler" tickets. One reviewable module or falsification task per PR.
2. Independent oracle wherever possible: GF(2) vs direct matrices, SAT witness vs tiny exhaustive search, decoder output vs logical-coset enumeration.
3. Never equate a virtual two-qubit encoding with two physical qubits. Represent unsupported features explicitly.
4. Track proven / tested / assumed / hypothetical separately. Cite papers and actual noise-model provenance.
5. Reproducibility: fixed configurations, dependency versions, seeds, commands, raw data, result plots and negative findings.
6. Plan updates are done after weekly Friday review. Do not silently enlarge scope when an experiment fails.
7. Suggested parallel technical tracks: Hamiltonian adapters; SAT code search; independent correctness verifier; decoder and compiler circuit stage. Freeze their interface before parallel development.