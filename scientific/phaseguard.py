"""PhaseGuard: finite ANNNI ground states and gate-noisy correlation readout.

All noise maps use ideal state injection. The circuit is a destructive parity
readout on a ring. This is not a noisy ground-state preparation algorithm.
"""
from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np
from scipy.sparse import csr_matrix, diags
from scipy.sparse.linalg import eigsh

ROOT = Path(__file__).resolve().parent
SEED = 20260925


class ANNNI:
    def __init__(self, n):
        if n % 4:
            raise ValueError("Use N divisible by four to fit the period-four state.")
        self.n = n
        self.ids = np.arange(2**n, dtype=np.int64)
        self.z = 1.0 - 2.0 * ((self.ids[:, None] >> np.arange(n-1, -1, -1)) & 1)
        self.nn = -(self.z * np.roll(self.z, -1, axis=1)).sum(axis=1)
        self.nnn = (self.z * np.roll(self.z, -2, axis=1)).sum(axis=1)
        rows = np.repeat(self.ids, n)
        cols = (self.ids[:, None] ^ (1 << np.arange(n))).ravel()
        self.flip = csr_matrix((np.ones(len(rows)), (rows, cols)), shape=(2**n, 2**n))
        # Exact even global-spin-flip sector. At low h the even/odd splitting
        # falls below machine precision; a full-space solver can mix them and
        # create false fidelity spikes even while its residual is small.
        half = 2**(n-1)
        reps = np.arange(half)
        flipped = reps[:,None] ^ (1 << np.arange(n))
        even_cols = np.minimum(flipped, 2**n-1-flipped).ravel()
        self.even_flip = csr_matrix((np.ones(half*n),(np.repeat(reps,n),even_cols)),shape=(half,half))
        self.corr_diag = np.array([(self.z * np.roll(self.z, -d, axis=1)).mean(axis=1)
                                   for d in range(n)])

    def hamiltonian(self, kappa, h):
        return diags(self.nn + kappa*self.nnn, format="csr") - h*self.flip

    def ground(self, kappa, h, initial=None):
        diag = self.nn + kappa*self.nnn
        if h == 0:
            support = np.isclose(diag, diag.min(), atol=1e-12)
            state = support.astype(float) / np.sqrt(support.sum())
            return float(diag.min()), state, 0.0
        half = len(diag)//2
        ham = diags(diag[:half],format="csr") - h*self.even_flip
        if initial is None:
            initial = np.ones(half) / np.sqrt(half)
        else:
            initial = (initial[:half]+initial[half:][::-1])/np.sqrt(2)
            initial /= np.linalg.norm(initial)
        # For h>0, Perron-Frobenius fixes the ground state in the even sector.
        vals, states = eigsh(ham, k=1, which="SA", v0=initial, tol=1e-11, maxiter=5000)
        sector_state = states[:,0]
        state = np.concatenate([sector_state,sector_state[::-1]])/np.sqrt(2)
        if state.sum() < 0:
            state = -state
        residual = float(np.linalg.norm(self.hamiltonian(kappa,h) @ state - vals[0]*state))
        return float(vals[0]), state, residual

    def observables(self, state):
        prob = np.abs(state)**2
        corr = self.corr_diag @ prob
        mx = float(state @ (self.flip @ state) / self.n)
        return corr, mx


def reference(k):
    k = np.asarray(k, float)
    ising = np.full_like(k, np.nan)
    mask = k < .5
    # Algebraically rationalized expression: exact limit h_I(0)=1.
    ising[mask] = 2*(1-2*k[mask]) / (1+np.sqrt((1-3*k[mask]+4*k[mask]**2)/(1-k[mask])))
    lower = np.where(k > .5, 1.05*(k-.5), np.nan)
    upper = np.where(k > .5, 1.05*np.sqrt(np.maximum(0, (k-.5)*(k-.1))), np.nan)
    return ising, lower, upper


def readout_gates(n, i, j):
    """Route i along the shortest ring path, then read ZZ as Z on j."""
    clockwise = (j-i) % n
    step = 1 if clockwise <= n/2 else -1
    distance = min(clockwise, n-clockwise)
    if distance == 0:
        return []
    route = [(i+step*s) % n for s in range(distance+1)]
    gates = []
    for a, b in zip(route[:-2], route[1:-1]):
        gates.extend([(a,b), (b,a), (a,b)])
    gates.append((route[-2], j))
    return gates


def heisenberg_readout(n, i, j):
    """Exact adjoint Pauli propagation. Counts only channels in the light cone."""
    support = {j}
    exponent = 0
    for control, target in reversed(readout_gates(n, i, j)):
        if target in support:
            exponent += 1  # D_p^*(Z) = (1-4p/3) Z
            support.symmetric_difference_update({control})
    return support, exponent


def attenuation(n, p, routed=True):
    d = np.minimum(np.arange(n), n-np.arange(n))
    exponent = np.where(d == 0, 0, 2*d-1 if routed else 1)
    return (1-4*p/3)**exponent


def structure(corr):
    """S(q)=N^-1 sum_d C(d) cos(qd); the self term is retained."""
    n = corr.shape[-1]
    q = 2*np.pi*np.arange(n//2+1)/n
    return np.einsum("...d,qd->...q", corr, np.cos(q[:,None]*np.arange(n))) / n


def order(corr):
    n = corr.shape[-1]
    sq = structure(corr)
    # Subtract the uncorrelated 1/N floor, normalize perfect order to 1.
    ferro = (n*sq[..., 0]-1)/(n-1)
    anti = (n*sq[..., n//4]-1)/(n/2-1)
    return ferro, anti, sq


def classify(corr, threshold):
    f, a, sq = order(corr)
    labels = np.full(f.shape, 2, dtype=int)  # paramagnetic-like
    labels[(f >= a) & (f >= threshold)] = 0
    labels[(a > f) & (a >= threshold)] = 1
    # Do not turn an intermediate wavevector into a floating-phase label.
    # A gapped, modulated paramagnet can have the same finite-ring peak.
    return labels


def crossing(hs, values, threshold):
    idx = np.flatnonzero((values[:-1] >= threshold) & (values[1:] < threshold))
    if len(idx) == 0:
        return None
    j = idx[0]
    return float(hs[j]+(hs[j+1]-hs[j])*(values[j]-threshold)/(values[j]-values[j+1]))


def scan(n=12, points=41):
    model = ANNNI(n)
    ks = np.linspace(0,1,points)
    hs = np.linspace(0,2,points)
    corr = np.empty((points,points,n))
    energy = np.empty((points,points))
    mx = np.empty_like(energy)
    chi = np.full_like(energy, np.nan)
    residual = np.empty_like(energy)
    t0 = time.monotonic()
    for ik,k in enumerate(ks):
        previous = None
        for ih,h in enumerate(hs):
            e, psi, r = model.ground(k,h,previous)
            energy[ih,ik], residual[ih,ik] = e,r
            corr[ih,ik],mx[ih,ik] = model.observables(psi)
            if previous is not None and ih > 1:
                # Convention: chi_F = -2 ln |<psi(h-dh)|psi(h)>| / dh^2.
                chi[ih,ik] = -2*np.log(max(abs(previous @ psi), 1e-15))/(hs[1]-hs[0])**2
            previous = psi
        print(f"N={n}: column {ik+1}/{points}, {time.monotonic()-t0:.1f}s", flush=True)
    np.savez_compressed(ROOT/f"data/grid_N{n}.npz", kappa=ks,h=hs,corr=corr,
                        energy=energy,mx=mx,chi=chi,residual=residual)


def validate():
    import pennylane as qml
    rng = np.random.default_rng(SEED)
    checks = {}
    model = ANNNI(8)
    terms = [-qml.Z(i)@qml.Z((i+1)%8) for i in range(8)]
    terms += [.37*qml.Z(i)@qml.Z((i+2)%8) for i in range(8)]
    terms += [-.63*qml.X(i) for i in range(8)]
    matrix = qml.matrix(sum(terms),wire_order=range(8))
    checks["hamiltonian_max_error"] = float(np.max(abs(matrix-model.hamiltonian(.37,.63).toarray())))
    for k,h in [(0,0),(.8,0),(0,.7),(.4,.5),(.8,.6)]:
        e,psi,r = model.ground(k,h)
        exact = np.linalg.eigvalsh(model.hamiltonian(k,h).toarray())[0]
        checks[f"energy_error_{k}_{h}"] = abs(e-exact)
    checks["ising_limit_error"] = abs(float(reference(np.array([0]))[0][0])-1)
    e,psi,_ = model.ground(0,0)
    checks["ferro_limit_error"] = abs(order(model.observables(psi)[0])[0]-1)
    e,psi,_ = model.ground(1,0)
    checks["antiphase_limit_error"] = abs(order(model.observables(psi)[0])[1]-1)
    for n in [4,8,12,16]:
        for i in range(n):
            for j in range(n):
                if i != j:
                    support, exp = heisenberg_readout(n,i,j)
                    assert support == {i,j}
                    assert exp == 2*min((j-i)%n,(i-j)%n)-1
    # Random complex states test a stronger case than a single ground state.
    noise_errors = []
    for n in [4,6]:
        state = rng.normal(size=2**n)+1j*rng.normal(size=2**n)
        state /= np.linalg.norm(state)
        dev = qml.device("default.mixed",wires=n)
        for i,j in [(0,1),(0,n//2),(n-1,1)]:
            for p in [0,.01,.05]:
                @qml.qnode(dev)
                def circuit():
                    qml.StatePrep(state,wires=range(n))
                    for c,t in readout_gates(n,i,j):
                        qml.CNOT(wires=[c,t])
                        qml.DepolarizingChannel(p,wires=t)
                    return qml.expval(qml.Z(j))
                exact = np.vdot(state,qml.matrix(qml.Z(i)@qml.Z(j),wire_order=range(n))@state).real
                _, exp = heisenberg_readout(n,i,j)
                noise_errors.append(abs(float(circuit())-exact*(1-4*p/3)**exp))
    checks["pennylane_noise_max_error"] = max(noise_errors)
    checks["pennylane_noise_circuit_count"] = len(noise_errors)
    assert checks["hamiltonian_max_error"] < 1e-12
    assert max(v for k,v in checks.items() if "error" in k) < 1e-8
    (ROOT/"data/validation.json").write_text(json.dumps(checks,indent=2))
    print(json.dumps(checks,indent=2))


def cuts():
    records = []
    for n in [8,12,16]:
        model = ANNNI(n)
        hs = np.linspace(.02,1.7,85)
        for k in [0,.3,.6,.8,1.]:
            previous = None
            corr = []
            chi = []
            energies = []
            for h in hs:
                e,psi,r = model.ground(k,h,previous)
                corr.append(model.observables(psi)[0])
                chi.append(None if previous is None else -2*np.log(max(abs(previous@psi),1e-15))/.02**2)
                energies.append(e)
                previous = psi
            records.append(dict(n=n,kappa=k,h=hs.tolist(),corr=np.array(corr).tolist(),chi=chi,energy=energies))
            print(f"line cut N={n}, kappa={k}",flush=True)
    (ROOT/"data/line_cuts.json").write_text(json.dumps(records,allow_nan=False))


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("action", choices=["scan","validate","cuts"])
    p.add_argument("--n",type=int,default=12)
    p.add_argument("--points",type=int,default=41)
    a = p.parse_args()
    (ROOT/"data").mkdir(exist_ok=True)
    if a.action == "scan": scan(a.n,a.points)
    elif a.action == "validate": validate()
    elif a.action == "cuts": cuts()
