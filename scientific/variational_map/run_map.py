"""Bounded 15x15 shallow VQE map with explicit target-CNOT noise.

Optimize the clean energy plus a known global-spin-flip symmetry penalty.
Freeze those parameters for p=0, 0.01, 0.05. Compare every point with ED.
"""
from __future__ import annotations
import json
import sys
import time
from pathlib import Path

import numpy as np
import pennylane as qml
from pennylane import numpy as pnp
from scipy.optimize import minimize

HERE=Path(__file__).resolve().parent
BASE=HERE.parent
sys.path.insert(0,str(BASE))
from phaseguard import ANNNI, classify
from noisy_preparation import RealDensity

N=8
LAYERS=4
POINTS=15
SEED=719834
THRESHOLD=.5628839496768857
SYMMETRY_WEIGHT=.5
MAXITER=200
RESCUE_ERROR_PER_SITE=.0025


def ansatz(x):
    for row in x.reshape(LAYERS,N):
        for q in range(N):
            qml.RY(row[q],wires=q)
        for q in range(N-1):
            qml.CNOT(wires=[q,q+1])


def main():
    t0=time.monotonic()
    rng=np.random.default_rng(SEED)
    model=ANNNI(N)
    engine=RealDensity(N)
    ks=np.linspace(0,1,POINTS)
    hs=np.linspace(0,2,POINTS)
    ps=np.array([0,.01,.05])
    shape=(POINTS,POINTS)
    complete=np.zeros(shape,bool)
    energy=np.full((3,)+shape,np.nan)
    corr=np.full((3,)+shape+(N,),np.nan)
    fidelity=np.full((3,)+shape,np.nan)
    ground_space_fidelity=np.full_like(fidelity,np.nan)
    purity=np.full_like(fidelity,np.nan)
    trace=np.full_like(fidelity,np.nan)
    exact_energy=np.full(shape,np.nan)
    exact_corr=np.full(shape+(N,),np.nan)
    parameters=np.full(shape+(LAYERS*N,),np.nan)
    objective_value=np.full(shape,np.nan)
    parity_value=np.full(shape,np.nan)
    error_per_site=np.full(shape,np.nan)
    optimizer_success=np.zeros(shape,bool)
    iteration_count=np.zeros(shape,int)
    evaluation_count=np.zeros(shape,int)
    trial_count=np.zeros(shape,int)
    residual=np.full(shape,np.nan)
    trials=[]

    control=json.loads((BASE/'variational_control/symmetry_results.json').read_text())
    fixed_seeds=[np.array(r['parameters']) for r in control['records']]
    exact_ferro=np.zeros((LAYERS,N));exact_ferro[-1,0]=np.pi/2
    exact_anti=np.zeros((LAYERS,N));exact_anti[-2,1]=-np.pi/2
    exact_anti[-1,0]=np.pi/2;exact_anti[-1,[2,4,6]]=np.pi
    exact_para=np.zeros((LAYERS,N));exact_para[-1,:]=np.pi/2
    fixed_seeds += [exact_ferro.ravel(),exact_anti.ravel(),exact_para.ravel()]
    parity=qml.prod(*(qml.X(q) for q in range(N)))
    dev=qml.device('lightning.qubit',wires=N)

    @qml.qnode(qml.device('default.qubit',wires=N))
    def state_circuit(x):
        ansatz(x)
        return qml.state()

    def save(final=False):
        name='variational_grid_N8.npz' if final else 'variational_partial.npz'
        np.savez_compressed(HERE/name,kappa=ks,h=hs,p=ps,corr=corr,
                            energy=energy,fidelity=fidelity,
                            ground_space_fidelity=ground_space_fidelity,
                            purity=purity,trace=trace,exact_energy=exact_energy,
                            exact_corr=exact_corr,parameters=parameters,
                            objective_value=objective_value,global_x_parity=parity_value,
                            clean_energy_error_per_site=error_per_site,
                            optimizer_success=optimizer_success,
                            iteration_count=iteration_count,evaluation_count=evaluation_count,
                            trial_count=trial_count,ed_residual=residual,complete=complete,
                            cnot_count=np.full(shape,28),threshold=THRESHOLD)
        meta=dict(seed=SEED,n=N,layers=LAYERS,points=POINTS,noise_values=ps.tolist(),
                  cnot_count=28,rotation_count=32,max_iterations_per_trial=MAXITER,
                  max_rescue_trials=2,rescue_energy_error_per_site=RESCUE_ERROR_PER_SITE,
                  symmetry_penalty_weight=SYMMETRY_WEIGHT,
                  objective='Energy + 0.5*(1 - expectation(X on all 8 sites))',
                  clean_parameters_frozen_under_noise=True,
                  selected_point_count=int(complete.sum()),
                  elapsed_seconds=time.monotonic()-t0,pennylane_version=qml.__version__,
                  classification_threshold=THRESHOLD,all_trials=trials)
        (HERE/'metadata.json').write_text(json.dumps(meta,indent=2))

    previous_params=None
    for ik,k in enumerate(ks):
        h_indices=range(POINTS-1,-1,-1) if ik%2==0 else range(POINTS)
        for ih in h_indices:
            h=hs[ih]
            e,ground,r=model.ground(k,h)
            exact_energy[ih,ik]=e
            exact_corr[ih,ik]=model.observables(ground)[0]
            residual[ih,ik]=r
            ham=sum([-qml.Z(q)@qml.Z((q+1)%N)+k*qml.Z(q)@qml.Z((q+2)%N)-h*qml.X(q)
                     for q in range(N)])
            target=ham+SYMMETRY_WEIGHT*(qml.Identity(0)-parity)

            @qml.qnode(dev,interface='autograd',diff_method='adjoint')
            def objective_circuit(x):
                ansatz(x)
                return qml.expval(target)

            gradient=qml.grad(objective_circuit)
            def objective(x):
                x=pnp.array(x,requires_grad=True)
                return float(objective_circuit(x)),np.asarray(gradient(x),dtype=float)

            seeds=list(fixed_seeds)
            if previous_params is not None:
                seeds.append(previous_params)
            if ik and complete[ih,ik-1]:
                seeds.append(parameters[ih,ik-1])
            seed_values=[float(objective_circuit(x)) for x in seeds]
            ranked=np.argsort(seed_values)
            best=None
            selected_trial=None
            point_trials=[]
            for trial in range(3):
                if trial==0:
                    x0=seeds[int(ranked[0])]
                elif trial==1:
                    x0=seeds[int(ranked[1])]
                else:
                    x0=best.x+rng.normal(0,.15,LAYERS*N)
                fit=minimize(objective,x0,jac=True,method='L-BFGS-B',
                             options={'maxiter':MAXITER,'ftol':1e-11,'gtol':2e-6,'maxls':25})
                record=dict(trial=trial,objective=float(fit.fun),iterations=int(fit.nit),
                            evaluations=int(fit.nfev),success=bool(fit.success),message=str(fit.message))
                point_trials.append(record)
                if best is None or fit.fun<best.fun:
                    best=fit
                    selected_trial=trial
                # The penalty is nonnegative. The objective excess therefore
                # gives a safe stopping bound on the physical energy excess.
                if (best.fun-e)/N <= RESCUE_ERROR_PER_SITE:
                    break

            x=np.asarray(best.x)
            psi=np.asarray(state_circuit(x))
            assert np.max(abs(psi.imag))<1e-12
            psi=psi.real
            sparse_h=model.hamiltonian(k,h)
            pure_energy=float(psi@(sparse_h@psi))
            pure_parity=float(psi@psi[::-1])
            assert abs(pure_energy+SYMMETRY_WEIGHT*(1-pure_parity)-best.fun)<1e-8
            assert pure_energy>=e-1e-8
            parameters[ih,ik]=x
            objective_value[ih,ik]=best.fun
            parity_value[ih,ik]=pure_parity
            error_per_site[ih,ik]=(pure_energy-e)/N
            optimizer_success[ih,ik]=best.success
            iteration_count[ih,ik]=sum(v['iterations'] for v in point_trials)
            evaluation_count[ih,ik]=sum(v['evaluations'] for v in point_trials)
            trial_count[ih,ik]=len(point_trials)
            corr[0,ih,ik]=model.observables(psi)[0]
            energy[0,ih,ik]=pure_energy
            fidelity[0,ih,ik]=float(abs(ground@psi)**2)
            purity[0,ih,ik]=trace[0,ih,ik]=1
            support=np.isclose(model.nn+k*model.nnn,e,atol=1e-10) if h==0 else None
            ground_space_fidelity[0,ih,ik]=float(np.sum(psi[support]**2)) if h==0 else fidelity[0,ih,ik]
            with qml.tape.QuantumTape() as tape:
                ansatz(x)
            gates=tape.operations
            assert sum(g.name=='CNOT' for g in gates)==28
            for ip,p in enumerate(ps[1:],1):
                rho=engine.run(gates,p)
                corr[ip,ih,ik]=model.corr_diag@np.diag(rho)
                energy[ip,ih,ik]=float(np.trace(sparse_h@rho))
                fidelity[ip,ih,ik]=float(ground@rho@ground)
                ground_space_fidelity[ip,ih,ik]=float(np.diag(rho)[support].sum()) if h==0 else fidelity[ip,ih,ik]
                purity[ip,ih,ik]=float(np.sum(rho*rho))
                trace[ip,ih,ik]=float(np.trace(rho))
                assert abs(trace[ip,ih,ik]-1)<1e-9
            complete[ih,ik]=True
            previous_params=x
            trials.append(dict(ik=ik,ih=ih,kappa=float(k),h=float(h),selected_trial=selected_trial,trials=point_trials))
        save()
        done=complete.sum()
        print(f'column {ik+1}/{POINTS}; points={done}; elapsed={time.monotonic()-t0:.1f}s; '
              f'mean clean E/site error={np.nanmean(error_per_site):.6g}; '
              f'max error={np.nanmax(error_per_site):.6g}',flush=True)
    assert complete.all()
    save(final=True)
    (HERE/'variational_partial.npz').unlink()
    result=dict(points=int(complete.sum()),elapsed_seconds=time.monotonic()-t0,
                mean_energy_error_per_site=float(error_per_site.mean()),
                max_energy_error_per_site=float(error_per_site.max()),
                mean_clean_fidelity_h_positive=float(fidelity[0,1:,:].mean()),
                min_clean_fidelity_h_positive=float(fidelity[0,1:,:].min()),
                optimizer_converged_points=int(optimizer_success.sum()),
                total_trials=int(trial_count.sum()),
                class_counts=[np.bincount(classify(corr[ip],THRESHOLD).ravel(),minlength=3).tolist() for ip in range(3)])
    (HERE/'summary.json').write_text(json.dumps(result,indent=2))
    print(json.dumps(result,indent=2))


if __name__=='__main__':
    main()
