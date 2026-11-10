# Research and development standards

This project uses short, verifiable implementation and research tasks. **Scientific assumptions and conclusions require review.** Each task specifies the mathematical basis, input/output contract, files in scope, independent tests and a falsifiable definition of done.

1. Avoid broad, undefined development tasks. Make each pull request a reviewable module or a narrowly framed test of a claim.
2. Use independent checks wherever possible: GF(2) calculations against direct matrices, SAT solutions against small exhaustive searches, decoder outputs against logical-coset enumeration.
3. Never equate virtual registers in a qutrit embedding with independently controllable physical qubits. Represent unsupported features explicitly.
4. Distinguish proven, experimentally tested, assumed and hypothetical statements. Cite sources and document noise-model provenance.
5. Preserve reproducibility: fixed configurations, dependency versions, random seeds, commands, raw results and negative findings.
6. Review the plan every Friday. When a calculation fails, document the issue and revise the corresponding task rather than silently expanding the scope.
7. Keep Hamiltonian representations, SAT code search, correctness verification, decoding and circuit construction as separate modules with stable, documented interfaces.
