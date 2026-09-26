"""Create the submission notebook and execute its saved-data analysis cells."""
import json
import nbformat as nbf
from nbclient import NotebookClient
from phaseguard import ROOT

cells=[]
def md(s):cells.append(nbf.v4.new_markdown_cell(s))
def code(s):cells.append(nbf.v4.new_code_cell(s))
md("""# PhaseGuard
## Circuit cost can hide quantum order

This notebook maps the 1D ANNNI model. It compares clean ground states with noisy state preparation.

**Main claim:** a phase detector can lose order because of its circuit. That loss does not prove a thermodynamic phase change.

The noise maps use N=8 and target depolarization after every preparation CNOT. We compare a 254-CNOT exact loader with a 28-CNOT variational circuit. A separate clean map uses N=12. Line cuts compare N=8, 12, and 16.

This work has substantial AI assistance. All results below come from executed local simulations. No QPU result is claimed.
""")
code("""from pathlib import Path
import json, numpy as np
from IPython.display import display, Image, Markdown, Code
from phaseguard import ANNNI, reference, order, classify, attenuation, crossing
ROOT = Path.cwd()
summary = json.loads((ROOT/'data/summary.json').read_text())
display(Markdown('**Saved result:** ' + str(summary['N8_grid_shape']) + ' noisy grid; ' + str(summary['N12_grid_shape']) + ' clean grid.'))
""")
md("""## 1. Hamiltonian and exact checks

H = −Σ ZiZi+1 + κ Σ ZiZi+2 − h Σ Xi, with periodic boundaries and J1=1.

The solver works in the exact even global-spin-flip sector. This prevents false fidelity peaks from numerical mixing of nearly degenerate even and odd states. At h=0 it uses the equal positive superposition of basis ground states. The multicritical point has a ground-state manifold, not a unique selected phase.

All implementation files are in this package. The code is imported below; the final section displays the full source. Use the commands in README.md to repeat all simulations.
""")
code("""checks = json.loads((ROOT/'data/validation.json').read_text())
density_checks = json.loads((ROOT/'data/density_validation.json').read_text())
assert checks['hamiltonian_max_error'] < 1e-10
assert max(r['matrix_max_error'] for r in density_checks) < 1e-9
print('Hamiltonian error:', checks['hamiltonian_max_error'])
print('PennyLane density-matrix error:', max(r['matrix_max_error'] for r in density_checks))
print('N=12 largest ground-state residual:', summary['N12_max_residual'])
""")
md("""## 2. Symmetry-safe phase signals

We measure C(d), the site average of ⟨Zi Zi+d⟩. We use S(q)=Σd C(d) cos(qd)/N.

The normalized signals are F=[N S(0)−1]/(N−1) and A=[N S(π/2)−1]/(N/2−1). Both equal one in their perfect ordered limit and zero for uncorrelated Z measurements. Raw ⟨Z⟩ is not used.

Each size uses a threshold fixed at the exact Ising point κ=0, h=1. The larger signal sets the ordered label. If neither signal passes the threshold, the label is weak order / para-like. The threshold is not fitted to the reference curves or noisy data.

Solid curves are the challenge Ising and upper KT references. The lower dashed curve is PT. These are qualitative finite-size references. We do not label an intermediate correlation momentum as a confirmed floating phase.
""")
code("""print('Fixed N=8 threshold:', summary['N8_threshold'])
print('Fixed N=12 threshold:', summary['N12_threshold'])
display(Image(filename=str(ROOT/'figures/clean_N12_diagnostics.png')))
display(Image(filename=str(ROOT/'figures/finite_size.png')))
""")
md("""## 3. Actual gate noise in state preparation

For each N=8 ground state, PennyLane constructs a Möttönen circuit. The circuit starts in |0…0⟩. The main grid runs the full gate list. A depolarizing channel follows every CNOT on its target: unchanged with probability 1−p, and each X/Y/Z error with probability p/3.

One-qubit gates and final direct Z measurements are ideal. A fast exact density-matrix implementation applies the same PennyLane gate list. It is checked against `default.mixed`. The clean p=0 state is obtained from exact diagonalization after independent checks of the compiled clean circuit.

These are detected-order maps of noisy outputs. Weak order is not evidence that the Hamiltonian entered a different equilibrium phase.
""")
code("""display(Image(filename=str(ROOT/'figures/preparation_triptych.png')))
print('Preparation CNOT range:', summary['preparation_CNOT_range'])
print(json.dumps(summary['noise'], indent=2))
""")
md("""## 4. Boundary estimates and resolution

Boundary values use linear interpolation of the first downward threshold crossing. Grid spacing is Δh=0.1 for the noisy map. A missing boundary is reported as `None`. It means the selected order signal never passes the fixed threshold; it is not a measured transition at h=0.

The clean map has Δh=0.05. Line cuts have Δh=0.02. Finite-size and threshold error exceed interpolation precision near the narrow floating region.
""")
code("""print(json.dumps(summary['preparation_boundaries'], indent=2))
print('Threshold sensitivity:', summary['threshold_sensitivity'])
print('N=8 reference agreement outside the floating strip:',
      summary['reference_agreement_N8_clean_outside_floating_band'])
""")
md("""## Short-circuit and finite-shot controls

Four layers of RY rotations and nearest-neighbor CNOTs use 28 CNOTs. Clean optimization uses energy plus 0.5(1−⟨X⊗N⟩). The symmetry penalty does not use target-state fidelity. The same clean parameters are used at each noise strength. All parameters and optimizer outcomes are saved in `variational_control/`.

The finite-shot check draws complete Z bitstrings. Every pair estimate uses the same batch, so the simulation retains covariance between observables. We use 2,000 independent batches of 1,024 shots at each control. Saved percentile intervals describe this repeated-sampling test.
""")
md("""## Full 15 × 15 variational map

The short-circuit map optimizes clean parameters at every grid point. It holds those parameters fixed for p=0, 0.01, and 0.05. Energy and fidelity checks compare every clean point with exact diagonalization. Open circles mark positive-field points with clean fidelity below 0.95. These points have preparation uncertainty before gate noise is added.

At h=0, the ground state is degenerate. The data also store ground-space fidelity. We exclude h=0 from single-state fidelity summaries. A low energy error alone does not prove high fidelity when the energy gap is small.
""")
code("""vsummary=json.loads((ROOT/'data/variational_summary.json').read_text())
display(Image(filename=str(ROOT/'figures/variational_triptych.png')))
display(Image(filename=str(ROOT/'figures/variational_quality.png')))
display(Image(filename=str(ROOT/'figures/variational_boundaries.png')))
print(json.dumps(vsummary,indent=2))
""")
code("""display(Image(filename=str(ROOT/'figures/circuit_comparison.png')))
shot_check=json.loads((ROOT/'data/finite_shots.json').read_text())
print('Shots per batch:',shot_check['shots_per_repetition'])
print('Number of batches:',shot_check['repetitions'])
for row in shot_check['rows']:
    if row['p']==.01:
        print(row['phase'],row['preparation'],'correct fraction:',row['correct_label_fraction'])
""")
md("""## 5. Mitigation test

We also simulate 2p. Linear zero-noise extrapolation uses Cest(0)=2C(p)−C(2p). It removes the first-order term when a smooth low-noise expansion is useful. It can retain a large bias when many gates accumulate error.

This test changes the numerical channel strength. It is not a hardware folding experiment. For independent estimates with the same variance, the coefficients 2 and −1 multiply variance by five. Saved maps use exact expectations and have no shot variance.
""")
code("""display(Image(filename=str(ROOT/'figures/mitigation.png')))
for row in summary['noise']:
    print('p:',row['p'],'correlation MAE:',row['correlation_MAE'],
          'ZNE MAE:',row['zne_correlation_MAE'])
""")
md("""## 6. Readout control and physical meaning

Direct Z measurement obtains every ZZ correlation without a CNOT. Thus this readout has no gate noise under the task model.

For a separate circuit test, a routed parity measurement over distance d uses 3(d−1)+1 CNOTs. Its exact signal is C(d)(1−4p/3)^(2d−1). The adjoint propagation and PennyLane checks agree. Analytic inversion recovers the mean only under the exact known model. It increases variance.

This optional test shows why a phase-robustness claim must specify the circuit. It is not used in the main preparation-noise maps.
""")
code("display(Image(filename=str(ROOT/'figures/readout_control.png')))")
md("""## 7. Limits and sources

N≤16 cannot settle the thermodynamic floating phase. Intermediate wavevectors can also occur in a gapped modulated paramagnet. The preparation circuit is an expensive exact-loading baseline. Its noise sensitivity does not establish an unavoidable limit for short VQE circuits. The main classifier is deliberately fixed before noise is applied.

Sources: [challenge](https://github.com/benmcdonough20/QSITE-2026-QuantumCoalition), [PennyLane ANNNI](https://pennylane.ai/demos/tutorial_annni), [noise convention](https://pennylane.ai/challenges/heisenberg_model), [Cea et al.](https://arxiv.org/abs/2402.11022), [Möttönen](https://docs.pennylane.ai/en/stable/code/api/pennylane.MottonenStatePreparation.html), [zero-noise extrapolation](https://arxiv.org/abs/1612.02058).
""")
md("## 8. Full source\n\nThe following cells display the exact implementation used for the data and figures.")
for name in ["phaseguard.py","noisy_preparation.py","analyze.py","analyze_variational.py","finite_shots.py","variational_control/run_control.py","variational_map/run_map.py"]:
    code(f"display(Code(filename=str(ROOT/{name!r}), language='python'))")
nb=nbf.v4.new_notebook(cells=cells,metadata={"kernelspec":{"name":"python3","display_name":"Python 3","language":"python"}})
path=ROOT/"PhaseGuard.ipynb"
nbf.write(nb,path)
NotebookClient(nb,timeout=120,kernel_name="python3",resources={"metadata":{"path":str(ROOT)}}).execute()
nbf.write(nb,path)
print(path)
