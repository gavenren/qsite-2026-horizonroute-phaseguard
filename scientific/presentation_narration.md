# PhaseGuard presentation narration

## Slide 1 — The question

A quantum phase map can look wrong for two different reasons. The model can change its ground state. Or the circuit can lose the signal that identifies that state. PhaseGuard separates these effects. We study the ANNNI chain from the challenge. We map clean states, prepare them with quantum circuits, and add the required gate noise. Our main result is that circuit design controls detectable order. A circuit with twenty-eight CNOT gates keeps most ordered labels at one-percent noise. A much longer exact-state loader loses them all.

## Slide 2 — Make the physics checkable

The ANNNI model has three competing terms. Neighboring spins prefer to align. Frustration favors the up, up, down, down pattern. The transverse field favors alignment along a different direction. We scan frustration from zero to one and field from zero to two. Our clean reference has twelve spins and forty-one points on each axis. We also compare eight, twelve, and sixteen spins. The solver uses the exact even symmetry sector. This prevents numerical mixing of nearly equal states from creating false transition signals. The multicritical point has many ground states, so we do not assign it a unique phase.

## Slide 3 — Detect order without breaking symmetry

A finite ordered state can have zero average magnetization. Therefore, we use spin correlations. One signal detects uniform order. The other detects the period-four antiphase. Both signals have a known ordered limit and an uncorrelated limit. We fix the decision threshold at the known Ising critical point. We keep it when noise is added. Fidelity susceptibility gives a second phase signal. At frustration zero point three, the exact order estimate moves toward the reference as the system grows. The antiphase side has larger size effects. These checks prevent us from treating one small-system map as an exact infinite-chain result.

## Slide 4 — The main noisy preparation

Our main noisy map has eight spins and fifteen points on each axis. A variational circuit has four layers, thirty-two rotations, and twenty-eight CNOT gates. Clean training minimizes energy with a known symmetry constraint. It does not optimize target-state fidelity. We then keep those parameters fixed. After every CNOT, we apply the required depolarizing channel to its target. Every circuit starts in the zero state. We simulate the complete noisy preparation. A fast density-matrix simulator agrees with PennyLane checks to about twelve decimal places. All parameters and optimizer stop records are saved.

## Slide 5 — What noise changes

The clean map has thirty-four ferromagnetic labels and thirty-six antiphase labels. At one-percent gate noise, these become thirty-three and thirty-four. Thus most ordered points remain visible. At the Ising end, the interpolated detection boundary moves from field one point zero zero to about zero point nine two. The grid is coarse, so these are approximate threshold estimates. At five-percent noise, no point passes the ordered threshold. That means the detector cannot resolve order. It does not mean a measured thermodynamic transition moved to zero field. Mean state fidelity falls much more than label agreement, which shows why both checks matter.

## Slide 6 — Compare the circuit cost

For comparison, PennyLane builds an exact state loader with two hundred and fifty-four CNOT gates. We test it on a twenty-one by twenty-one grid. At one-percent noise, every strong-order label disappears. Three same-point controls show the difference clearly. The short circuit gives state fidelities near eighty percent. Exact loading gives only about fifteen to nineteen percent. Both use the same error probability after each CNOT. This is evidence for a lower gate budget. It is not a universal rule that one phase is more fragile. Some phase errors reduce fidelity while leaving the measured spin correlations almost unchanged.

## Slide 7 — Mitigation has a useful range

We repeat the exact-loader experiment at twice the noise and use linear zero-noise extrapolation. At one-percent noise, it restores ninety-nine of the one hundred and twenty-seven clean ordered labels. It adds no false ordered labels. Overall label agreement improves from seventy-one point two percent to ninety-three point seven percent. The correlation error also decreases. At five-percent noise, no ordered label returns. The accumulated error is too large for this simple correction. This is a simulated noise-scaling test, not a hardware result. With equal independent measurement variance, the extrapolation would multiply that variance by five.

## Slide 8 — Keep the claim within the evidence

The clean variational labels match exact labels at two hundred and twenty-one of two hundred and twenty-five points. However, twelve positive-field points have clean fidelity below ninety-five percent. We mark them as preparation-uncertain. We also test finite measurement samples and compare system sizes. We do not claim a confirmed floating phase. Intermediate correlation momenta are insufficient evidence. The package includes raw arrays, code, saved parameters, an executed notebook, figures, and checks. This project has substantial AI assistance for code, analysis, and text. All reported results are local simulations. The practical lesson is to check the physics, report the circuit, and test the noise before trusting the map.
