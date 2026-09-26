"""Finite-shot stress test at three saved control points.

Each repetition samples the full Z bitstring distribution. All pair estimates
therefore have the covariance of a shared measurement batch.
"""
import json
import numpy as np
import pennylane as qml
from phaseguard import ROOT, ANNNI, order
from noisy_preparation import RealDensity,compile_state


def run(shots=1024,repeats=2000):
    rng=np.random.default_rng(20260925)
    model=ANNNI(8);engine=RealDensity(8)
    f_values,a_values,_=order(model.corr_diag.T)
    threshold=.5628839496768858
    rows=[]
    controls=json.loads((ROOT/'variational_control/symmetry_results.json').read_text())['records']
    for record in controls:
        _,state,_=model.ground(record['kappa'],record['h'])
        full,_,_=compile_state(state,8)
        with qml.tape.QuantumTape() as tape:
            for layer in np.array(record['parameters']).reshape(4,8):
                for q,theta in enumerate(layer):qml.RY(theta,wires=q)
                for q in range(7):qml.CNOT(wires=[q,q+1])
        for kind,gates in [('shallow',tape.operations),('exact_loader',full)]:
            for p in [0,.01,.05]:
                rho=engine.run(gates,p)
                prob=np.maximum(0,np.diag(rho));prob/=prob.sum()
                counts=rng.multinomial(shots,prob,size=repeats)
                ff=counts@f_values/shots;aa=counts@a_values/shots
                label=np.full(repeats,2);label[(ff>=aa)&(ff>=threshold)]=0
                label[(aa>ff)&(aa>=threshold)]=1
                target={'ferromagnetic':0,'antiphase':1,'paramagnetic':2}[record['phase_point']]
                rows.append(dict(phase=record['phase_point'],preparation=kind,p=p,
                                 correct_label_fraction=float(np.mean(label==target)),
                                 ferro_interval_95=np.quantile(ff,[.025,.975]).tolist(),
                                 anti_interval_95=np.quantile(aa,[.025,.975]).tolist(),
                                 ferro_exact=float(prob@f_values),anti_exact=float(prob@a_values)))
    result=dict(seed=20260925,shots_per_repetition=shots,repetitions=repeats,
                interpretation='Empirical 2.5 and 97.5 percentiles of repeated samples; ideal direct Z measurement.',rows=rows)
    (ROOT/'data/finite_shots.json').write_text(json.dumps(result,indent=2))
    for row in rows: print(row['phase'],row['preparation'],row['p'],row['correct_label_fraction'])


if __name__=='__main__':run()
