"""Exact density-matrix simulation of PennyLane's compiled Möttönen circuits.

The same qml operation list is used by the reference default.mixed simulator.
The specialized implementation is real valued; this ground-state problem has
real nonnegative amplitudes. It rejects nonzero complex rotations.
"""
from __future__ import annotations
import argparse
import json
import time
import numpy as np
import pennylane as qml
from phaseguard import ANNNI, ROOT, SEED, order, classify


class RealDensity:
    def __init__(self,n):
        self.n=n
        self.dim=2**n
        ids=np.arange(self.dim)
        self.flip={t: ids ^ (1<<(n-1-t)) for t in range(n)}
        self.zero={t: ids[((ids>>(n-1-t))&1)==0] for t in range(n)}
        self.same={t: (((ids[:,None]^ids[None,:])>>(n-1-t))&1)==0 for t in range(n)}
        self.perm={(c,t):ids ^ (((ids>>(n-1-c))&1) << (n-1-t))
                   for c in range(n) for t in range(n) if c!=t}

    def run(self,gates,p):
        rho=np.zeros((self.dim,self.dim))
        rho[0,0]=1
        for gate in gates:
            if gate.name == "RY":
                angle=float(gate.data[0]); wire=int(gate.wires[0])
                c,s=np.cos(angle/2),np.sin(angle/2)
                a=self.zero[wire]; b=self.flip[wire][a]
                ra,rb=rho[a,:].copy(),rho[b,:].copy()
                rho[a,:]=c*ra-s*rb; rho[b,:]=s*ra+c*rb
                ra,rb=rho[:,a].copy(),rho[:,b].copy()
                rho[:,a]=c*ra-s*rb; rho[:,b]=s*ra+c*rb
            elif gate.name == "CNOT":
                c,t=map(int,gate.wires); perm=self.perm[c,t]
                rho=rho[np.ix_(perm,perm)]
                if p:
                    flipped=rho[np.ix_(self.flip[t],self.flip[t])]
                    same=self.same[t]
                    # Blocks diagonal in target Z mix with their flipped block.
                    # Off-diagonal target blocks shrink by 1-4p/3.
                    rho=np.where(same,(1-2*p/3)*rho+(2*p/3)*flipped,(1-4*p/3)*rho)
            elif gate.name == "GlobalPhase":
                pass
            elif gate.name == "RZ" and abs(float(gate.data[0]))<1e-12:
                pass
            else:
                raise ValueError(f"Unsupported gate: {gate}")
        return rho


def compile_state(state,n):
    # For h>0 the stoquastic ground state is positive. Remove tiny eigensolver
    # roundoff in tails only; the norm of this change is checked and recorded.
    if np.min(state)<-1e-9:
        raise ValueError("Ground state is not nonnegative.")
    positive=np.maximum(state,0)
    positive/=np.linalg.norm(positive)
    error=float(np.linalg.norm(positive-state))
    assert error<1e-8
    return qml.MottonenStatePreparation(positive,wires=range(n)).decomposition(),positive,error


def validate():
    records=[]
    for n in [4,8]:
        model=ANNNI(n); engine=RealDensity(n)
        for k,h in [(.15,.25),(.85,.2),(.2,1.5)]:
            _,state,_=model.ground(k,h)
            gates,state,error=compile_state(state,n)
            for p in [0,.01,.05]:
                rho=engine.run(gates,p)
                @qml.qnode(qml.device("default.mixed",wires=n))
                def reference():
                    for gate in gates:
                        qml.apply(gate)
                        if gate.name=="CNOT":
                            qml.DepolarizingChannel(p,wires=gate.wires[1])
                    return qml.density_matrix(wires=range(n))
                ref=np.asarray(reference())
                err=float(np.max(abs(rho-ref)))
                trace=abs(float(np.trace(rho))-1)
                minimum=float(np.linalg.eigvalsh(rho).min())
                assert err<1e-9 and trace<1e-10 and minimum>-1e-10
                if p==0: assert abs(float(state@rho@state)-1)<1e-10
                records.append(dict(n=n,kappa=k,h=h,p=p,matrix_max_error=err,
                                    trace_error=trace,minimum_eigenvalue=minimum))
                print(records[-1],flush=True)
    (ROOT/"data/density_validation.json").write_text(json.dumps(records,indent=2))


def scan(points=21,mitigate=False):
    n=8; model=ANNNI(n); engine=RealDensity(n)
    ks=np.linspace(0,1,points); hs=np.linspace(0,2,points)
    ps=np.array([0,.02,.10] if mitigate else [0,.01,.05])
    prefix="mitigation" if mitigate else "preparation"
    shape=(len(ps),points,points)
    corr=np.empty(shape+(n,)); fidelity=np.empty(shape); count=np.empty((points,points),int)
    residual=np.empty((points,points)); truncation=np.empty_like(residual)
    purity=np.empty(shape); trace=np.empty(shape)
    t0=time.monotonic()
    for ik,k in enumerate(ks):
        previous=None
        for ih,h in enumerate(hs):
            _,state,r=model.ground(k,h,previous); previous=state
            gates,state,err=compile_state(state,n)
            residual[ih,ik]=r; truncation[ih,ik]=err
            count[ih,ik]=sum(g.name=="CNOT" for g in gates)
            # Exact p=0 state is used after independent gate-level validation.
            corr[0,ih,ik]=model.observables(state)[0]
            fidelity[0,ih,ik]=1; purity[0,ih,ik]=1; trace[0,ih,ik]=1
            for ip,p in enumerate(ps[1:],1):
                rho=engine.run(gates,p)
                corr[ip,ih,ik]=model.corr_diag@np.diag(rho)
                fidelity[ip,ih,ik]=state@rho@state
                purity[ip,ih,ik]=np.sum(rho*rho)
                trace[ip,ih,ik]=np.trace(rho)
        print(f"prep map {ik+1}/{points}: {time.monotonic()-t0:.1f}s",flush=True)
        # Keep an explicit partial checkpoint, never expose it as a full map.
        np.savez_compressed(ROOT/f"data/{prefix}_partial.npz",completed_columns=ik+1,
                            kappa=ks,h=hs,p=ps,corr=corr[:,:,:ik+1],fidelity=fidelity[:,:,:ik+1])
    assert np.max(abs(trace-1))<1e-9
    np.savez_compressed(ROOT/f"data/{prefix}_grid_N8.npz",kappa=ks,h=hs,p=ps,
                        corr=corr,fidelity=fidelity,purity=purity,trace=trace,
                        cnot_count=count,residual=residual,truncation=truncation)
    (ROOT/f"data/{prefix}_partial.npz").unlink()


if __name__=="__main__":
    parser=argparse.ArgumentParser(); parser.add_argument("action",choices=["validate","scan"])
    parser.add_argument("--points",type=int,default=21)
    parser.add_argument("--mitigate",action="store_true"); args=parser.parse_args()
    if args.action=="validate": validate()
    else: scan(args.points,args.mitigate)
