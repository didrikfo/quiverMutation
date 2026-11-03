"""Round 111 toolsmith, item 2 (continued): congruence invariants of the four E-147 class-2 key-kept steps (rebuilt as in toolsmith_c2_exact.py):
Smith forms of S = C+C^T and K = C-C^T, rank of radical of S, det(xC+C^T) equality (D = 0), Smith form of C+tC^T... a difference in any is a proved Z-congruence NO.
usage (repo root): timeout 10m .venv/bin/python workshop/rounds/111/toolsmith_c2_inv.py"""
import sys, pickle, ast
sys.path.insert(0, '.'); sys.path.insert(0, 'workshop/rounds/111')
import numpy as np, sympy as sp
from sympy.matrices.normalforms import smith_normal_form
import toolsmith_c1_13_exact as X
g = X.__dict__
def mk(edges, reprr, nodes):
    return g['build']([(a, b) for a, b, k in edges], [[[p[0][0]] + [a[1] for a in p] for p in d] for d in ast.literal_eval(reprr)], nodes)
snf = lambda M: [int(abs(x)) for x in smith_normal_form(sp.Matrix(M.tolist()), domain=sp.ZZ).diagonal()]
x = sp.symbols('x')
for t, z in enumerate(pickle.load(open('workshop/rounds/043/skeptic_zero_n7c2.pkl', 'rb'))):
    N = sorted({q for e in z['pe'] for q in e[:2]}); v = z['v']
    A = mk(z['pe'], z['pr'], N); ch = g['reduction'].reducePathAlgebra(g['mutation'].quiverMutationAtVertex(A, v))
    CA = np.array(g['cartan'](A), dtype=np.int64); CB = np.array(g['cartan'](ch), dtype=np.int64)
    out = {}
    for nm, f in (('S', lambda C: C + C.T), ('K', lambda C: C - C.T), ('C+2Ct', lambda C: C + 2 * C.T), ('C-2Ct', lambda C: C - 2 * C.T), ('C+3Ct', lambda C: C + 3 * C.T), ('2C+Ct', lambda C: 2*C + C.T)):
        a, b = snf(f(CA)), snf(f(CB)); out[nm] = (a == b, a, b)
    pa = sp.expand((x * sp.Matrix(CA.tolist()) + sp.Matrix(CA.T.tolist())).det()); pb = sp.expand((x * sp.Matrix(CB.tolist()) + sp.Matrix(CB.T.tolist())).det())
    print('step', t, 'det(xC+Ct) equal', sp.expand(pa - pb) == 0, 'poly', sp.factor(pa))
    for k, (e, a, b) in out.items(): print('   SNF', k, 'equal' if e else 'DIFFER', a, b)
