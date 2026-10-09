"""Round 042 (theorist): the J != 0 steps that do NOT have the shape (C' e_v = e_v, |out v| = 1, u = e_w - e_i): their u, J, C_B row v, Q(x) coefficients.
Usage: theorist_exceptions.py dump.pkl"""
import sys, pickle
import numpy as np, sympy as sp
from collections import Counter
D = pickle.load(open(sys.argv[1], 'rb')); n = D['n']; x = sp.symbols('x')
P = lambda C: sp.Poly((x * sp.Matrix(C.tolist()) + sp.Matrix(C.T.tolist())).det(), x)
c = Counter(); N = 0
for r in D['recs']:
    Jn = {i: j for i, j in r['J'].items() if j}
    if not Jn: continue
    V, v = r['V'], r['v']; k = V.index(v)
    rm = np.eye(n, dtype=int); rm[k, k] = -1
    for w in r['out']: rm[k, V.index(w)] += 1
    Cp = rm @ r['CA'] @ rm.T; keep = [a for a in range(n) if a != k]
    colv = all(Cp[a, k] == 0 for a in keep)
    if len(r['out']) == 1 and len(Jn) == 1 and colv:
        w = keep.index(V.index(r['out'][0])); i = keep.index(V.index(list(Jn)[0])); e = np.eye(n - 1, dtype=int)
        if (Cp[k, keep] == e[w] - e[i]).all(): continue
    N += 1
    q = (P(r['CB']) - P(Cp)).all_coeffs()[::-1]; low = next((j for j, a in enumerate(q) if a), None)
    c[('colv', colv, 'u', tuple(int(a) for a in Cp[k, keep]), 'dimJ', tuple(sorted(Jn.values())), 'CB row v', tuple(int(a) for a in r['CB'][k, keep]), 'Q low', low, int(q[low]) if low is not None else None)] += 1
print('exceptions', N)
for kx, vv in sorted(c.items(), key=lambda t: -t[1])[:15]: print(vv, kx)
