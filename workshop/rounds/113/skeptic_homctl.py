"""Round 113 (skeptic), T10: planted controls for a Hom_K(T,T[m]) tilting test (engine: rounds/050/skeptic_tilt.py, imported unchanged).
Usage (repo root): skeptic_homctl.py N [MAXALG]   (LNAs of length N and their relation duals, first MAXALG of each)
Controls (truth known by construction, not by any tilting engine):
 P1 zero differential: T = A_{i!=k} + P_k + (+P_h)[-1].  P_h and P_h[-1] are both summands, so id_{P_h} gives Hom(T,T[+-1]) != 0: truth NO.
 P2 duplicated summand: T = sum_{i!=k} P_i + P_j (j != k twice, P_k missing): Hom(T,T[m]) = 0, not a generator (K0 det 0): truth NO, Hom-only test says YES.
 P3 right-T, wrong child: End(T_k) Cartan vs Cartan of the child at another vertex k' != k (labelled): truth 'does not match' whenever Cartans differ.
 R  real steps T_k (S=-1): count Hom-tilting, and K0 det (must be +-1)."""
import sys, importlib.util, itertools
sys.path.insert(0, '.')
import numpy as np
spec = importlib.util.spec_from_file_location('skt', 'workshop/rounds/050/skeptic_tilt.py'); skt = importlib.util.module_from_spec(spec); spec.loader.exec_module(skt)
from quivermutation import nakayama as nk, pathAlgebra, arrowPaths as ap, mutation, reduction, invariants
N = int(sys.argv[1]); MAXA = int(sys.argv[2]) if len(sys.argv) > 2 else 10**9

def tilt_ok(A, T, V):
    z = True
    for m in (-1, 1):
        if any(skt.homdim(A, T[i], T[j], m) for i in V for j in V): z = False
    return z
def k0det(T, V, alg):
    idx = {v: n for n, v in enumerate(V)}; M = np.zeros((len(V), len(V)))
    for r, i in enumerate(V):
        for d, vs in T[i][0].items():
            for v in vs: M[r, idx[v]] += (-1) ** d
    return round(abs(np.linalg.det(M)))
res = dict(alg=0, P1=0, P1no=0, P2=0, P2homyes=0, P2k0no=0, P3=0, P3no=0, R=0, Rtilt=0, Rk0=0, Rcart=0)
algs = []
for lna in list(nk.LinearNakayamaAlgebra.allOfLength(N))[:MAXA]:
    algs += [lna, pathAlgebra.dualPathAlgebra(lna)]
for alg in algs:
    A = skt.Alg(alg); V = sorted(alg.quiver.nodes); res['alg'] += 1
    for k in V:
        outs = ap.arrowsOutOf(alg.quiver, k)
        if not outs: continue
        T = skt.complexT(A, k, V, -1)
        res['R'] += 1; res['Rtilt'] += tilt_ok(A, T, V); res['Rk0'] += (k0det(T, V, alg) == 1)
        # P1
        Tz = dict(T); Tz[k] = ({-1: [k], 0: [b[1] for b in outs]}, {-1: {}})
        Tz2 = {i: ({0: [i]}, {}) for i in V if i != k}
        Tz2['k'] = Tz[k]; Tz2['h'] = ({1: [b[1] for b in outs]}, {})   # explicit split of the cone: P_k at 0? keep P_k and (+P_h)[-1] as stalk summands
        Tz2['k'] = ({0: [k]}, {}); Tz2['h'] = ({1: [b[1] for b in outs]}, {})
        Vz = list(Tz2)
        res['P1'] += 1; res['P1no'] += (not tilt_ok(A, Tz2, Vz))
        # P2
        js = [j for j in V if j != k]
        if js:
            j = js[0]; T2 = {i: ({0: [i]}, {}) for i in V if i != k}; T2['dup'] = ({0: [j]}, {})
            V2 = list(T2); res['P2'] += 1
            res['P2homyes'] += tilt_ok(A, T2, V2)
            M = np.zeros((len(V), len(V))); idx = {v: n for n, v in enumerate(V)}
            for r, i in enumerate(V2): M[r, idx[T2[i][0][0][0]]] += 1
            res['P2k0no'] += (round(abs(np.linalg.det(M))) == 0)
        # P3: Cartan match of End(T_k) vs child at k'
        if tilt_ok(A, T, V):
            H = np.array([[skt.homdim(A, T[i], T[j], 0) for j in V] for i in V])
            for k2 in V:
                if k2 == k or not ap.arrowsOutOf(alg.quiver, k2): continue
                try:
                    raw = mutation.quiverMutationAtVertex(alg, k2); c = reduction.reducePathAlgebra(raw)
                    if sorted(c.vertices()) != V: continue
                    C = np.array(invariants.cartanMatrix(c, exact=True).tolist(), dtype=int)
                except Exception: continue
                res['P3'] += 1; res['P3no'] += (not (H == C.T).all())
print(N, res)
