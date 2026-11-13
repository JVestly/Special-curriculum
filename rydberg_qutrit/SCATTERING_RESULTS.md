# Intermediate-state scattering: validated numerical extension (13 Nov 2026)

**Status: intermediate model, NOT an experimentally calibrated two-qubit gate prediction.**

## What was added

The existing model has a two-atom effective Rydberg Hamiltonian in the local
basis `|0>, |1>, |r>, |L>` and six Rydberg decay channels (three destinations
per atom). In addition, the noise model averages over a **common-mode**
quasi-static Gaussian detuning distribution, using a measured ground-Rydberg
`T2*=3 microseconds` and an *assumed* Gaussian Ramsey envelope.

We added intermediate-state scattering via *adiabatic elimination* of a
short-lived `|e>` state. In the reference paper, `|e>` is a manifold of
6P3/2 states, excited by 420 nm light and coupled to Rydberg by 1013 nm light.
Our reduced single-intermediate-state approximation takes

```
(H_e / hbar) = -Delta |e><e|
 + 1/2 [Omega_b exp(+i phi)|1><e| + Omega_red|e><r| + h.c.]
```

At large detuning `|Delta| >> |Omega_b|, |Omega_red|`, the virtual intermediate
amplitude is approximately

```
<e|psi> = [Omega_b exp(-i phi)<1|psi> + Omega_red<r|psi>]/(2 Delta).
```

Given a decay from `|e>` to an **assumed** destination `|a>`, a leading-order
effective jump operator is

```
C_eff,a(t) = sqrt(Gamma_e) |a> [Omega_b exp(-i phi(t))<1| + Omega_red<r|]/(2 Delta).
```

For two atoms, the same local process acts independently at each site.
The coherent Raman coupling `Omega_eff = Omega_b Omega_red/(2 Delta)` is
already included in the original Rydberg Hamiltonian, so we do not add it
again. We also do not add a separate AC Stark shift; the original calibrated
phase/detuning convention remains the reference. This is an approximation.

### Parameters and provenance

| Parameter | Value | Source |
|---|---:|---|
| `Omega_420 / 2pi` | 237 MHz | Evered et al. (2023) |
| `Omega_1013 / 2pi` | 303 MHz | Evered et al. (2023) |
| `Delta / 2pi` | 7,800 MHz | Evered et al. (2023) |
| Effective `Omega_r / 2pi` | 4.6 MHz | Evered et al. (2023) |
| `T_r` | 88 microseconds | Evered et al. Rydberg noise model |
| `T2*` | 3 microseconds | Evered et al., modeled as assumed Gaussian common-mode detuning |
| `T_e` | 112.8 ns | Rb 6P3/2 lifetime used in Sheng et al., PRA 78, 062506 (2008), referring to earlier measurements |
| e-scattering branch destination | `|1>` | **ASSUMED, not calibrated** |

The complete intermediate-state branching probabilities are **not** known
from these parameters. Changing the branch target to `|0>` (Raman-like) or
`|L>` (loss-like) is supported as a *sensitivity check*, but neither is a
measured branching distribution. The actual experimental system has several
intermediate hyperfine/Zeeman levels, and has laser/motion effects omitted here.

## Channel results

All reported probabilities are computed from a CPTP channel for **arbitrary
input on the computational subspace** (4 input dimensions). The physical output
has 16 dimensions, including Rydberg and loss states. The final outcome
partition is averaged over a maximally mixed computational input `I_4/4`.
The Pauli weights are extracted from the computational block of the Kraus
operators after removing the closest unitary of the *noiseless reduced pulse*
(not necessarily the exact experimental CZ gate). Thus the Pauli weights are
those of the Pauli-twirled **surviving computational channel**, not literal
physical Pauli jump probabilities. The identity term is not the same as
'no photon scattered'.

| Channel outcome | Decay + T2* | Plus intermediate scattering into `|1>` |
|---|---:|---:|
| Identity Pauli weight | 99.7908754% | 99.7638516% |
| All nonidentity Pauli weights | 0.0824876% | 0.0934635% |
| Residual Rydberg | 0.0233637% | 0.0394094% |
| Output in loss state | 0.1032733% | 0.1032754% |
| Total | ~100% | ~100% |

The model produces a change of **+0.0109759 percentage points** in the
nonidentity Pauli-twirled weight (about a 13.3% relative increase). Loss into
`|L>` is nearly unchanged because the **assumed** e-scattering destination is
`|1>`, not `|L>`.

This does not demonstrate that the physical scattering branch is actually
dominated by return to `|1>`! In the decay-only sensitivity analysis (no T2*
fluctuations), putting all e-scattering into `|L>` instead changes the
input-averaged loss from ~0.10329% to ~0.13826%, and changes the Pauli weights.
This is proof that *a measured lifetime alone is not sufficient to predict
all channel error probabilities*.

## Kraus operators: full matrices

For the chosen (assumed-return-to-|1>) model:

- Choi matrix: `64 x 64` (reference dimension 4 × output dimension 16).
- Kraus matrices: **25 operators**, each `16 x 4`.
- `scattering_results/kraus_operators.txt`: **all** 25 complex matrices.
- `scattering_results/full_kraus_chosen_model.npz`: complete machine-readable
  Kraus array, Choi matrix, reference unitary and all 16 Pauli weights.
- `scattering_results/scattering_channel_results.json`: structured summary.
- Choi eigenvalues and Kraus count depend on numerical precision/tolerance.
  Individual Kraus matrices are **not unique**; their ensemble is the channel.

Numerical verification:

- `max |sum_j K_j^dagger K_j - I_4| = 1.11e-15`
- `max |J - sum_j vec(K_j) vec(K_j)^dagger| = 7.78e-16`
- `sum_P q_P + p_rydberg + p_loss ≈ 1` to 1e-15.
- All six new unit tests pass.

### Reproduce

```
python -m pip install numpy scipy matplotlib
cd rydberg_qutrit
python -m unittest test_rydberg_scattering -v
python run_scattering_results.py
python plot_scattering.py
```

The full input-averaged channel extraction evaluates every 16 computational
operator basis elements, not a single initial `rho(0)`. The saved `.npz`
contains concrete matrices ready to use in the next analysis step.

## Comparison and limitations before SAT

The physically grounded quantities are the experimentally characterized
laser frequencies, detuning, and atomic lifetimes. The **effective scattering
branch and common-mode fluctuation model are modeling assumptions**. Thus these
numbers should be described as *model-derived predictions under explicitly
stated assumptions*, not experimental channel tomography.

Before treating this as a QEC noise model, we should (1) obtain 6P scattering
branching ratios for the specific drive polarization and manifold, (2) include
intermediate-state substructure or validate elimination against an explicit
small-system e-level model, (3) distinguish detected loss from undetected
leakage, and (4) define a physically realizable error channel for syndrome
extraction or idle storage instead of silently extrapolating a gate-pulse
channel to an entire QEC cycle.

## References

- S. J. Evered et al., *High-fidelity parallel entangling gates on a neutral-atom quantum computer*, Nature 622, 268–272 (2023), https://doi.org/10.1038/s41586-023-06481-y.
- D. Sheng et al., *Lifetime measurements of the 5d states of rubidium*, Physical Review A 78, 062506 (2008), https://doi.org/10.1103/PhysRevA.78.062506 (uses a published average 6P3/2 lifetime 112.8±1.7 ns from prior measurements).
