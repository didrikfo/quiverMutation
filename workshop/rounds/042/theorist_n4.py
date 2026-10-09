"""Round 042 (theorist): the n = 4 non-LNA example of E-141 through the same lens: C', Z = C_A|V\\v, u, F = Z Z^-T, orbit relation F^s e_w = e_i, moments c_k, Q.
Usage: theorist_n4.py (repo root)"""
import sys
import numpy as np, sympy as sp
sys.path.insert(0, '.'); ARGV = sys.argv; sys.argv = ['x']
src = open('workshop/rounds/033/experimentalist_bothdie.py').read().split("mode = _a[1]")[0]
exec(compile(src, 'bd', 'exec'))
A = build([(1,2),(1,3),(2,3),(2,4),(3,4)], [[[1,2,3,4],[1,3,4]]], [1,2,3,4]); v = 3
raw = mutation.quiverMutationAtVertex(A, v); ch = reduction.reducePathAlgebra(raw)
V = sorted(A.quiver.nodes); n = len(V); k = V.index(v)
CA = np.array(invariants.cartanMatrix(A), dtype=int); CB = np.array(invariants.cartanMatrix(ch), dtype=int)
out = [w for (a, w) in A.quiver.edges() if a == v]; J = perI(A, v)
rm = np.eye(n, dtype=int); rm[k, k] = -1
for w in out: rm[k, V.index(w)] += 1
Cp = rm @ CA @ rm.T; keep = [a for a in range(n) if a != k]
print('out', out, 'J', J); print("C'\n", Cp); print('C_B - C\'\n', CB - Cp)
Z = Cp[np.ix_(keep, keep)]; u = Cp[k, keep]; print('Z\n', Z, '\nu', u, ' column v of C\' is e_v:', all(Cp[a, k] == 0 for a in keep))
Y = np.rint(np.linalg.inv(Z)).astype(int); F = Z @ Y.T; print('F\n', F)
x = sp.symbols('x'); P = lambda C: sp.expand((x * sp.Matrix(C.tolist()) + sp.Matrix(C.T.tolist())).det())
print('Q =', sp.factor(P(CB) - P(Cp)))
