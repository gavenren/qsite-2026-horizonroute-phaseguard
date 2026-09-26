"""Check the saved map against its parameters and PennyLane density matrices."""
import json
import sys
from pathlib import Path
import numpy as np
import pennylane as qml

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent))
from phaseguard import ANNNI,classify
from run_map import ansatz,THRESHOLD


def main():
    data=np.load(HERE/'variational_grid_N8.npz')
    assert data['complete'].all()
    assert data['corr'].shape==(3,15,15,8)
    for key in ['corr','energy','fidelity','parameters','exact_energy','exact_corr']:
        assert np.isfinite(data[key]).all(),key
    assert np.max(abs(data['trace']-1))<1e-9
    assert np.max(abs(data['corr'][...,0]-1))<1e-9
    assert np.min(data['energy']-data['exact_energy'])>-1e-8
    assert np.min(data['fidelity'])>-1e-10 and np.max(data['fidelity'])<1+1e-9
    assert np.min(data['purity'])>=1/256-1e-10
    model=ANNNI(8)
    records=[]
    for ih,ik in [(2,6),(9,2),(2,12),(14,14)]:
        params=data['parameters'][ih,ik]
        k,h=data['kappa'][ik],data['h'][ih]
        _,ground,_=model.ground(k,h)
        for ip,p in enumerate(data['p']):
            @qml.qnode(qml.device('default.mixed',wires=8))
            def circuit():
                for row in params.reshape(4,8):
                    for q in range(8):
                        qml.RY(row[q],wires=q)
                    for q in range(7):
                        qml.CNOT(wires=[q,q+1])
                        qml.DepolarizingChannel(p,wires=q+1)
                return qml.density_matrix(wires=range(8))
            rho=np.asarray(circuit())
            corr=model.corr_diag@np.diag(rho).real
            energy=float(np.trace(model.hamiltonian(k,h)@rho).real)
            fidelity=float(np.vdot(ground,rho@ground).real)
            record=dict(ih=ih,ik=ik,p=float(p),
                        correlation_error=float(np.max(abs(corr-data['corr'][ip,ih,ik]))),
                        energy_error=abs(energy-data['energy'][ip,ih,ik]),
                        fidelity_error=abs(fidelity-data['fidelity'][ip,ih,ik]))
            assert max(record[key] for key in ['correlation_error','energy_error','fidelity_error'])<1e-8
            records.append(record)
    exact_labels=classify(data['exact_corr'],THRESHOLD)
    shallow_labels=classify(data['corr'][0],THRESHOLD)
    confusion=np.zeros((3,3),int)
    for a,b in zip(exact_labels.ravel(),shallow_labels.ravel()):
        confusion[a,b]+=1
    low_fidelity=np.argwhere(data['fidelity'][0,1:,:]<.95)
    weak=[]
    for ih0,ik in low_fidelity:
        ih=ih0+1
        weak.append(dict(kappa=float(data['kappa'][ik]),h=float(data['h'][ih]),
                         fidelity=float(data['fidelity'][0,ih,ik]),
                         energy_error_per_site=float(data['clean_energy_error_per_site'][ih,ik])))
    result=dict(pennylane_map_spot_checks=records,
                max_correlation_error=max(r['correlation_error'] for r in records),
                max_energy_error=max(r['energy_error'] for r in records),
                max_fidelity_error=max(r['fidelity_error'] for r in records),
                clean_label_agreement=float(np.mean(exact_labels==shallow_labels)),
                clean_confusion_rows_exact_columns_vqe=confusion.tolist(),
                positive_field_points_below_095_fidelity=weak,
                positive_field_point_count=210,
                max_trace_error=float(np.max(abs(data['trace']-1))))
    (HERE/'audit.json').write_text(json.dumps(result,indent=2))
    print(json.dumps(result,indent=2))


if __name__=='__main__':
    main()
