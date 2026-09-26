"""Small-N control: noise in decomposed exact amplitude preparation.

This is separate from the main readout-noise experiment. No preparation gates
are silently treated as noise-free. Native RY/RZ are ideal by the task model.
"""
import json
import time
import numpy as np
import pennylane as qml
from phaseguard import ANNNI, ROOT, order


def run():
    n = 8
    model = ANNNI(n)
    records = []
    t0 = time.monotonic()
    for name,k,h in [("ferromagnetic",.15,.25),("antiphase",.85,.2),("paramagnetic",.2,1.5)]:
        _,state,_ = model.ground(k,h)
        gates = qml.MottonenStatePreparation(state,wires=range(n)).decomposition()
        # Decomposition may contain a global phase. Applying it is harmless.
        cnot_count = sum(g.name == "CNOT" for g in gates)
        assert all(g.name in {"RY","RZ","CNOT","GlobalPhase"} for g in gates)
        for p in [0,.01,.05]:
            @qml.qnode(qml.device("default.mixed",wires=n))
            def circuit():
                for gate in gates:
                    qml.apply(gate)
                    if gate.name == "CNOT":
                        qml.DepolarizingChannel(p,wires=gate.wires[1])
                return qml.density_matrix(wires=range(n))
            rho = np.asarray(circuit())
            corr = model.corr_diag @ np.diag(rho).real
            f,a,_ = order(corr)
            rec = dict(phase=name,kappa=k,h=h,n=n,p=p,cnot_count=cnot_count,
                       total_operations=len(gates),ferro=float(f),anti=float(a),
                       fidelity=float((state @ rho @ state).real),trace=float(np.trace(rho).real))
            records.append(rec)
            print(rec, f"elapsed={time.monotonic()-t0:.1f}s",flush=True)
            if p == 0:
                assert rec["fidelity"] > 1-1e-8
    (ROOT/"data/preparation_control.json").write_text(json.dumps(records,indent=2))


if __name__ == "__main__": run()
