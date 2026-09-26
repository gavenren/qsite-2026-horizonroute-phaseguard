"""One bounded neighbor pass at points with large clean preparation error."""
import json
import sys
import time
from pathlib import Path
import numpy as np
import pennylane as qml
from pennylane import numpy as pnp
from scipy.optimize import minimize

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent))
from phaseguard import ANNNI,classify
from noisy_preparation import RealDensity
from run_map import ansatz,THRESHOLD,SYMMETRY_WEIGHT,MAXITER


def main():
    t0=time.monotonic()
    path=HERE/'variational_grid_N8.npz'
    with np.load(path) as source:
        d={k:source[k].copy() for k in source.files}
    assert d['complete'].all()
    np.savez_compressed(HERE/'initial_grid_N8.npz',**d)
    model=ANNNI(8);engine=RealDensity(8)
    points=d['complete'].shape[0]
    target_mask=d['clean_energy_error_per_site']>.003
    target_mask[1:] |= d['fidelity'][0,1:]<.95
    targets=[tuple(int(x) for x in v) for v in np.argwhere(target_mask)]
    targets.sort(key=lambda pair:d['clean_energy_error_per_site'][pair],reverse=True)
    records=[]
    dev=qml.device('lightning.qubit',wires=8)
    parity=qml.prod(*(qml.X(q) for q in range(8)))
    @qml.qnode(qml.device('default.qubit',wires=8))
    def pure_state(x):
        ansatz(x)
        return qml.state()
    for index,(ih,ik) in enumerate(targets):
        k,h=d['kappa'][ik],d['h'][ih]
        ham=sum([-qml.Z(q)@qml.Z((q+1)%8)+k*qml.Z(q)@qml.Z((q+2)%8)-h*qml.X(q) for q in range(8)])
        target=ham+SYMMETRY_WEIGHT*(qml.Identity(0)-parity)
        @qml.qnode(dev,interface='autograd',diff_method='adjoint')
        def cost(x):
            ansatz(x)
            return qml.expval(target)
        grad=qml.grad(cost)
        def objective(x):
            x=pnp.array(x,requires_grad=True)
            return float(cost(x)),np.asarray(grad(x),float)
        neighbor_indices=[(ih,ik)]+[(a,b) for a,b in [(ih-1,ik),(ih+1,ik),(ih,ik-1),(ih,ik+1)] if 0<=a<points and 0<=b<points]
        seeds=[d['parameters'][a,b].copy() for a,b in neighbor_indices]
        values=[float(cost(x)) for x in seeds]
        seed=seeds[int(np.argmin(values))]
        fit=minimize(objective,seed,jac=True,method='L-BFGS-B',options={'maxiter':MAXITER,'ftol':1e-11,'gtol':2e-6,'maxls':25})
        before=float(d['objective_value'][ih,ik])
        changed=fit.fun<before-1e-11
        rec=dict(ih=ih,ik=ik,kappa=float(k),h=float(h),before=before,after=float(fit.fun),
                 changed=bool(changed),iterations=int(fit.nit),evaluations=int(fit.nfev),
                 success=bool(fit.success),message=str(fit.message),seed_neighbor=neighbor_indices[int(np.argmin(values))])
        records.append(rec)
        d['iteration_count'][ih,ik]+=fit.nit;d['evaluation_count'][ih,ik]+=fit.nfev;d['trial_count'][ih,ik]+=1
        if changed:
            x=np.asarray(fit.x)
            e,ground,res=model.ground(k,h)
            psi=np.asarray(pure_state(x)).real
            sparse=model.hamiltonian(k,h)
            energy=float(psi@(sparse@psi));par=float(psi@psi[::-1])
            assert abs(energy+SYMMETRY_WEIGHT*(1-par)-fit.fun)<1e-8
            d['parameters'][ih,ik]=x;d['objective_value'][ih,ik]=fit.fun
            d['global_x_parity'][ih,ik]=par;d['clean_energy_error_per_site'][ih,ik]=(energy-e)/8
            d['optimizer_success'][ih,ik]=fit.success
            d['corr'][0,ih,ik]=model.observables(psi)[0];d['energy'][0,ih,ik]=energy
            d['fidelity'][0,ih,ik]=float(abs(ground@psi)**2)
            support=np.isclose(model.nn+k*model.nnn,e,atol=1e-10) if h==0 else None
            d['ground_space_fidelity'][0,ih,ik]=float(np.sum(psi[support]**2)) if h==0 else d['fidelity'][0,ih,ik]
            with qml.tape.QuantumTape() as tape:ansatz(x)
            for ip,p in enumerate(d['p'][1:],1):
                rho=engine.run(tape.operations,p)
                d['corr'][ip,ih,ik]=model.corr_diag@np.diag(rho)
                d['energy'][ip,ih,ik]=float(np.trace(sparse@rho))
                d['fidelity'][ip,ih,ik]=float(ground@rho@ground)
                d['ground_space_fidelity'][ip,ih,ik]=float(np.diag(rho)[support].sum()) if h==0 else d['fidelity'][ip,ih,ik]
                d['purity'][ip,ih,ik]=float(np.sum(rho*rho));d['trace'][ip,ih,ik]=float(np.trace(rho))
        print(f'refine {index+1}/{len(targets)}: ({k:.3f},{h:.3f}) changed={changed}; {time.monotonic()-t0:.1f}s',flush=True)
    np.savez_compressed(path,**d)
    (HERE/'refinement.json').write_text(json.dumps(dict(max_iterations_per_point=MAXITER,elapsed_seconds=time.monotonic()-t0,records=records),indent=2))
    result=dict(points=int(d['complete'].sum()),refined_points=len(records),
                mean_energy_error_per_site=float(d['clean_energy_error_per_site'].mean()),
                max_energy_error_per_site=float(d['clean_energy_error_per_site'].max()),
                mean_clean_fidelity_h_positive=float(d['fidelity'][0,1:,:].mean()),
                min_clean_fidelity_h_positive=float(d['fidelity'][0,1:,:].min()),
                optimizer_converged_points=int(d['optimizer_success'].sum()),
                total_trials=int(d['trial_count'].sum()),
                class_counts=[np.bincount(classify(d['corr'][ip],THRESHOLD).ravel(),minlength=3).tolist() for ip in range(3)])
    (HERE/'summary.json').write_text(json.dumps(result,indent=2))
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
