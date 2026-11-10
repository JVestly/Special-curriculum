# Three-level Hamiltonian symmetry investigation — 2026-11-10

**Status: a controlled physical-qutrit comparison, not a claimed NA→LMG reduction.**

For each atom use the basis (|0>,|1>,|r>). Define local transitions E_ab^(i)=|a_i><b_i| and collective T_ab=Σ_i E_ab^(i), on the **3^N-dimensional physical Hilbert space**.

## Two models, not one
- Driven three-level neutral-atom model: H_NA = (Omega_hf/2)(T_01+T_10) + (Omega_r/2)(T_1r+T_r1) + Σ_a e_a T_aa + Σ_(i<j) V_ij n_(r,i)n_(r,j). The two single-atom drives generally break the population parities of their connected levels.
- Generalized three-level LMG benchmark: H_LMG = (-eps_0 T_00 + eps_1 T_11 + eps_r T_rr)/N - lambda/[N(N-1)] Σ_(a≠b) T_ab², for N≥2. The individual level energies are a generalization of a specific published symmetric choice.

## Derive the symmetry; do not guess the generators
Population number N_c=T_cc obeys [N_c,T_ab]=(δ_ca-δ_cb)T_ab. Hence [N_c,T_ab²]=2(δ_ca-δ_cb)T_ab².

Define Pi_c = exp(iπ N_c). Then Pi_c T_ab² Pi_c† = exp(2iπ(δ_ca-δ_cb)) T_ab² = T_ab². Diagonal energy terms also commute with each Pi_c. Thus *all three LMG population parities commute with H_LMG*.

Since N_0+N_1+N_r=N I, their product equals (-1)^N I; **only two are independent** modulo overall phase. This does NOT show that these are all Hamiltonian symmetries: collective-model permutations and non-Pauli operators remain possible.

For the driven NA Hamiltonian, a 0↔1 drive flips the parity of the 0 and 1 occupations and a 1↔r drive flips 1 and r. With both drives active, the individual population parities generally do not commute with H_NA; with a drive switched off, a corresponding parity may return.

## Optional qubit representation; not new physical atoms
Encoding |0>→|00>, |1>→|01>, |r>→|11> reproduces local parity on these three allowed states using -IZ, ZZ, ZI respectively. For N=3, parity group members include weight-three IZIZIZ and ZIZIZI. This is a **six-virtual-register** representation of three qutrits, NOT an independent six-qubit physical code. The action on unused |10> is an extension choice; test commutation in the chosen embedding separately.

## What is now verified?
The accompanying independent direct-matrix tests (not using the SAT implementation) check N=2 and N=3 LMG parities, native two-drive parity breaking, drive-off cases, the product identity, and the restricted binary parity representatives. Passing the tests validates **these facts**, not the physical plausibility of LMG pair transitions on a NA processor.

## H-code gate before SAT
A candidate G must: act on the claimed physical information carriers; have the correct symmetry-sector sign; be mutually commuting and independent with other seeds; define an allowed nonempty physical code subspace; be compatible with the target *error-pair* correctability conditions. Extra SAT-selected stabilizers need not preserve the whole Hamiltonian dynamics.

**Next question:** which *actual* physical computational-qubit Hamiltonian and error set permit an implementable H-seeded code? Keep LMG benchmark results separate from NA-hardware claims.