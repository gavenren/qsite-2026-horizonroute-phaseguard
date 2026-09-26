"""Shallow PennyLane variational preparation for three ANNNI control points.

This uses 4 layers of 8 RY rotations and 7 nearest-neighbor CNOT gates.
Clean parameters are optimized once, then held fixed for all noise values.
No reference phase curves or labels enter the objective.
"""
from __future__ import annotations

import argparse
import json
import platform
import sys
import time
from pathlib import Path

import numpy as np
import pennylane as qml
from pennylane import numpy as pnp
from scipy.optimize import minimize

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
from phaseguard import ANNNI, order, structure

N = 8
LAYERS = 4
SEED = 250925


def ansatz(params, p=0.0, noisy=False):
    params = params.reshape(LAYERS,N)
    for row in params:
        for q in range(N):
            qml.RY(row[q],wires=q)
        for q in range(N-1):
            qml.CNOT(wires=[q,q+1])
            if noisy:
                qml.DepolarizingChannel(p,wires=q+1)


def seed_parameters(name):
    x = np.zeros((LAYERS,N))
    if name == 'ferromagnetic':
        x[-1,0] = np.pi/2
    elif name == 'paramagnetic':
        x[-1,:] = np.pi/2
    else:
        # A product antiphase in the final output. The inverse of a chain
        # CNOT layer gives the required input bits for that last layer.
        bits = np.array([0,0,1,1,0,0,1,1])
        before = bits ^ np.r_[0,bits[:-1]]
        x[-1,:] = np.pi*before
    return x.ravel()


def run(symmetry_weight=0.0, output='results.json', iterations=400):
    t0 = time.monotonic()
    rng = np.random.default_rng(SEED)
    model = ANNNI(N)
    records = []
    pure_dev = qml.device('lightning.qubit',wires=N)

    @qml.qnode(qml.device('default.qubit',wires=N))
    def pure_state(params):
        ansatz(params)
        return qml.state()

    for name,k,h in [('ferromagnetic',.1,.2),('antiphase',.9,.1),('paramagnetic',.2,1.8)]:
        exact_energy, exact_state, residual = model.ground(k,h)
        hamiltonian = sum([-qml.Z(q)@qml.Z((q+1)%N)
                           +k*qml.Z(q)@qml.Z((q+2)%N)-h*qml.X(q)
                           for q in range(N)])
        parity = qml.prod(*(qml.X(q) for q in range(N)))
        objective_hamiltonian = hamiltonian + symmetry_weight*(qml.Identity(0)-parity)

        @qml.qnode(pure_dev,interface='autograd',diff_method='adjoint')
        def energy(params):
            ansatz(params)
            return qml.expval(objective_hamiltonian)

        gradient = qml.grad(energy)

        def objective(x):
            x = pnp.array(x,requires_grad=True)
            return float(energy(x)), np.asarray(gradient(x),dtype=float)

        physical_seed = seed_parameters(name)
        starts = [physical_seed,
                  physical_seed+rng.normal(0,.15,LAYERS*N),
                  physical_seed+rng.normal(0,.4,LAYERS*N),
                  rng.uniform(-np.pi,np.pi,LAYERS*N),
                  rng.normal(0,.5,LAYERS*N),
                  np.full(LAYERS*N,np.pi/4)]
        if name == 'antiphase' and symmetry_weight > 0:
            symmetric_anti = np.zeros((LAYERS,N))
            symmetric_anti[-2,1] = -np.pi/2
            symmetric_anti[-1,0] = np.pi/2
            symmetric_anti[-1,[2,4,6]] = np.pi
            starts[0] = symmetric_anti.ravel()
        trials = []
        best = None
        for trial, x0 in enumerate(starts):
            fit = minimize(objective,x0,jac=True,method='L-BFGS-B',
                           options={'maxiter':iterations,'ftol':1e-13,'gtol':1e-7,'maxls':30})
            candidate = dict(trial=trial,objective=float(fit.fun),
                             iterations=int(fit.nit),evaluations=int(fit.nfev),
                             success=bool(fit.success),message=str(fit.message))
            trials.append(candidate)
            if best is None or fit.fun < best.fun:
                best = fit
            print(name,candidate,'elapsed',round(time.monotonic()-t0,1),flush=True)

        params = np.asarray(best.x)
        state = np.asarray(pure_state(params))
        clean_fidelity = float(abs(np.vdot(exact_state,state))**2)
        classical_energy = float(np.vdot(state,model.hamiltonian(k,h)@state).real)
        parity_value = float(np.vdot(state,state[::-1]).real)
        expected_objective = classical_energy + symmetry_weight*(1-parity_value)
        assert abs(expected_objective-best.fun) < 1e-9
        assert classical_energy >= exact_energy-1e-9
        noisy_records = []
        for p in [0,.01,.05]:
            @qml.qnode(qml.device('default.mixed',wires=N))
            def noisy_state():
                ansatz(params,p,noisy=True)
                return qml.density_matrix(wires=range(N))

            rho = np.asarray(noisy_state())
            corr = model.corr_diag @ np.diag(rho).real
            f,a,sq = order(corr)
            e = float(np.trace(model.hamiltonian(k,h)@rho).real)
            fidelity = float(np.vdot(exact_state,rho@exact_state).real)
            rec = dict(p=p,energy=e,energy_error_per_site=(e-exact_energy)/N,
                       fidelity_to_exact=fidelity,ferro_order=float(f),anti_order=float(a),
                       structure_factor_q0=float(sq[0]),
                       structure_factor_qpi2=float(sq[N//4]),
                       correlations=corr.tolist(),trace=float(np.trace(rho).real),
                       purity=float(np.trace(rho@rho).real))
            assert abs(rec['trace']-1) < 1e-9
            if p == 0:
                assert abs(e-classical_energy) < 1e-9
                assert abs(fidelity-clean_fidelity) < 1e-9
            noisy_records.append(rec)
            print(name,'noise',p,'energy',e,'fidelity',fidelity,flush=True)
        records.append(dict(phase_point=name,kappa=k,h=h,n=N,layers=LAYERS,
                            cnot_count=LAYERS*(N-1),rotation_count=LAYERS*N,
                            exact_energy=exact_energy,exact_residual=residual,
                            best_clean_energy=classical_energy,
                            clean_energy_error_per_site=(classical_energy-exact_energy)/N,
                            clean_fidelity=clean_fidelity,global_x_parity=parity_value,
                            symmetry_penalty_weight=symmetry_weight,optimizer_trials=trials,
                            parameters=params.tolist(),noise=noisy_records))
        result = dict(description='Four-layer real variational preparation; fixed clean parameters under target depolarizing noise.',
                      seed=SEED,pennylane_version=qml.__version__,python_version=platform.python_version(),
                      elapsed_seconds=time.monotonic()-t0,records=records)
        (HERE/output).write_text(json.dumps(result,indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--symmetry-weight',type=float,default=0.0)
    parser.add_argument('--output',default='results.json')
    parser.add_argument('--iterations',type=int,default=400)
    args = parser.parse_args()
    run(args.symmetry_weight,args.output,args.iterations)
