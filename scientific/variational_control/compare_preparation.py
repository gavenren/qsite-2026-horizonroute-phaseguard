"""Compare shallow and full-amplitude preparation at the same three points."""
import json
import sys
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent))
from phaseguard import ANNNI, order
from noisy_preparation import compile_state, RealDensity


def main():
    source = json.loads((HERE/'symmetry_results.json').read_text())
    assert len(source['records']) == 3
    model=ANNNI(8)
    engine=RealDensity(8)
    rows=[]
    for record in source['records']:
        k,h=record['kappa'],record['h']
        exact_energy,state,residual=model.ground(k,h)
        gates,state,truncation=compile_state(state,8)
        count=sum(g.name=='CNOT' for g in gates)
        for shallow in record['noise']:
            p=shallow['p']
            rho=engine.run(gates,p)
            corr=model.corr_diag@np.diag(rho)
            f,a,sq=order(corr)
            full=dict(p=p,cnot_count=count,fidelity_to_exact=float(state@rho@state),
                      energy=float(np.trace(model.hamiltonian(k,h)@rho)),
                      ferro_order=float(f),anti_order=float(a),
                      structure_factor_q0=float(sq[0]),structure_factor_qpi2=float(sq[2]),
                      correlations=corr.tolist())
            if p == 0:
                assert full['fidelity_to_exact'] > 1-1e-10
            rows.append(dict(phase_point=record['phase_point'],kappa=k,h=h,p=p,
                             exact_energy=exact_energy,
                             shallow_cnot_count=record['cnot_count'],
                             shallow=shallow,full_amplitude=full))
    result=dict(description='Same-point comparison; all target CNOT channels use the same p.',rows=rows)
    (HERE/'comparison.json').write_text(json.dumps(result,indent=2))
    for row in rows:
        print(row['phase_point'],row['p'],
              'shallow fidelity',round(row['shallow']['fidelity_to_exact'],6),
              'full fidelity',round(row['full_amplitude']['fidelity_to_exact'],6),
              'CNOTs',row['shallow_cnot_count'],row['full_amplitude']['cnot_count'])


if __name__ == '__main__':
    main()
