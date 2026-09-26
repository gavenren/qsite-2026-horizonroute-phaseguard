# PhaseGuard: test the circuit before you trust the phase map

## Result

Circuit design controls how much quantum order remains visible under noise. Our 28-CNOT preparation retains 67 of 70 ordered grid labels at p = 0.01. A 254-CNOT exact-state loader retains none on its separate grid. At p = 0.05, both lose every strong-order label. This is a detection result. It does not show that the Hamiltonian changed its equilibrium phases.

We supply a 41 × 41 clean map at N = 12, a 15 × 15 noisy VQE map at N = 8, and a 21 × 21 exact-loader comparison. Further checks use N = 16, fidelity susceptibility, error mitigation, and finite measurement samples.

## Model and phase signals

The periodic ANNNI Hamiltonian is:

H = −Σ Zi Zi+1 + κ Σ Zi Zi+2 − h Σ Xi.

We set J1 = 1 and scan 0 ≤ κ ≤ 1, 0 ≤ h ≤ 2. Neighbor coupling favors uniform order. Frustration favors the period-four antiphase. The transverse field favors paramagnetic order [1, 2].

Sparse diagonalization gives the clean reference. We use the exact even global-X symmetry sector. This prevents false fidelity spikes from numerical mixing of nearly degenerate even and odd states. At h = 0, we select the equal positive superposition of basis ground states. We do not assign a unique phase to the degenerate point κ = 0.5, h = 0.

Raw magnetization can vanish in an ordered finite system. We instead measure C(d), the site average of ⟨Zi Zi+d⟩, and S(q) = Σd C(d) cos(qd)/N. Our normalized signals are:

F = [N S(0) − 1]/(N − 1); A = [N S(π/2) − 1]/(N/2 − 1).

Each signal is one in its perfect ordered state and zero for uncorrelated Z outcomes. A threshold fixed at κ = 0, h = 1 gives 0.562884 for N = 8 and 0.520118 for N = 12. The larger signal sets the ordered label. Below both thresholds, we report weak order or a para-like signal. The noisy data do not set the threshold.

---PAGE---

## Main noise experiment: 28-CNOT VQE

Four layers use 32 RY gates and 28 chain CNOTs. Clean optimization minimizes energy plus 0.5(1−⟨X⊗N⟩). The known symmetry constraint does not use target-state fidelity. Parameters remain fixed for p = 0, 0.01, and 0.05.

Each CNOT has target depolarization immediately after it. The target is unchanged with probability 1−p. Each X, Y, or Z error has probability p/3. One-qubit gates and direct Z measurement are ideal. Every noisy circuit starts in |0…0⟩. We simulate the full preparation, not a final noise channel on an ideal input.

| Detected region | p = 0 | p = 0.01 | p = 0.05 |
|---|---:|---:|---:|
| Ferro-like points | 34 | 33 | 0 |
| Antiphase-like points | 36 | 34 | 0 |
| Weak-order points | 155 | 158 | 225 |

At p = 0.01, ferro-label retention is 97.1%; antiphase-label retention is 94.4%. Mean positive-field fidelity falls from 0.9852 to 0.7981. At p = 0.05 it is 0.3453. Thus a phase label and state fidelity measure different failures.

At κ = 0, the interpolated order boundary moves from h ≈ 1.00 to 0.92 at p = 0.01. At κ ≈ 0.286, it moves from 0.50 to 0.47. These are linear threshold estimates on a grid with Δh ≈ 0.143, not precise thermodynamic boundaries. At p = 0.05 no crossing remains; we report an unresolved boundary, not a transition at zero field.

Clean labels agree with exact N = 8 labels at 221 of 225 points. Preparation errors remain: 12 positive-field points have fidelity below 0.95; the minimum is 0.665. Open circles mark these points. Mean and maximum energy errors are 0.0040 and 0.0190 per site. The selected optimizer converged at 196 points. All stop records and parameters are saved. At h = 0, we also record ground-space fidelity because the ground state is degenerate.

## Clean size check

At κ = 0.3, the exact order boundary moves from 0.478 to 0.460 to 0.453 for N = 8, 12, 16. The reference is 0.442. The N = 16 fidelity peak is approximately 0.43. Line-cut spacing is 0.02. Antiphase size effects are larger. At κ = 0.8, the order estimates are 0.635, 0.522, and 0.446. We do not treat these as converged infinite-chain results.

---PAGE---

## Why the circuit matters

PennyLane Möttönen preparation uses 254 CNOTs per exact-loaded state [4]. On its 21 × 21 grid, 127 points pass the clean order threshold. None pass at p = 0.01 or 0.05. Three same-point controls isolate circuit cost: at p = 0.01, short-circuit fidelities are 0.809, 0.794, and 0.820. Exact loading gives 0.189, 0.149, and 0.178.

The antiphase is slightly less robust by grid-label retention in our short circuit. This is not a universal ranking. Circuit structure and the selected state matter. ZZ signals also miss phase errors that reduce fidelity. A weak-order label can hide a highly mixed state.

Linear zero-noise extrapolation uses Cest(0) = 2C(p) − C(2p) [5]. For exact loading at p = 0.01, it recovers 99 of 127 ordered labels with no false ordered labels. Overall label agreement rises from 71.2% to 93.7%. Correlation error falls from 0.1423 to 0.0932. At p = 0.05, no ordered label returns. Numerical noise scaling is used; this is not hardware mitigation. Equal independent measurement variances would increase fivefold.

## Checks and limits

PennyLane checks of saved VQE points agree within 1.13 × 10⁻¹² for correlations. Independent random-circuit matrix checks also pass. The largest N = 12 eigenstate residual is 2.42 × 10⁻¹⁰. In 2,000 batches of 1,024 Z samples, both selected short-circuit ordered controls keep their correct labels at p = 0.01. Shared bitstring samples preserve covariance between pair estimates.

We do not claim a confirmed floating phase. Intermediate momenta also occur in a modulated paramagnet. Larger-system gap and correlation-decay scaling are needed [3]. The supplied lower reference is PT, not BKT; the upper line is KT/BKT. We also fix the Ising limit hI(0) = 1. Finite size, threshold choice, preparation error, and the idealized noise model limit all conclusions.

**AI disclosure.** This project has substantial AI assistance for code, analysis, and text. All results come from executed local simulations. No QPU result or quantum speedup is claimed.

## Sources

[1] [QSITE challenge](https://github.com/benmcdonough20/QSITE-2026-QuantumCoalition)

[2] [PennyLane ANNNI demo](https://pennylane.ai/demos/tutorial_annni)

[3] [Cea et al., ANNNI phase diagram](https://arxiv.org/abs/2402.11022)

[4] [PennyLane Möttönen preparation](https://docs.pennylane.ai/en/stable/code/api/pennylane.MottonenStatePreparation.html)

[5] [Temme et al., error mitigation](https://arxiv.org/abs/1612.02058)
