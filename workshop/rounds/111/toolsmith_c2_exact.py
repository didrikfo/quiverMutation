"""Round 111 toolsmith, answer to referee item 2: E-147's four pickled class-2 (n=7, key class index 2) key-kept J != 0 steps
(workshop/rounds/043/skeptic_zero_n7c2.pkl, written by `skeptic_c2.py 7 2 3 500 30 guard`).  Rebuild parent/child from arrows+relations as
043/skeptic_zero.py, take the library Cartan matrices, print signature/det/Phi charpoly, equality, and the exact Z-congruence search of toolsmith_c1_13_exact.py.
usage (repo root): timeout 10m .venv/bin/python workshop/rounds/111/toolsmith_c2_exact.py"""
import sys, pickle, ast
sys.path.insert(0, '.'); ARGV = sys.argv
import numpy as np, sympy as sp
sys.path.insert(0, 'workshop/rounds/111')
import toolsmith_c1_13_exact as X   # also execs bothdie (nk, search, ... in X's globals)
g = X.__dict__
def mk(edges, reprr, nodes):
    arrows = [(a, b) for a, b, k in edges]; rels = []
    for d in ast.literal_eval(reprr):
        rels.append([[p[0][0]] + [a[1] for a in p] for p in d])
    return g['build'](arrows, rels, nodes)
for t, z in enumerate(pickle.load(open('workshop/rounds/043/skeptic_zero_n7c2.pkl', 'rb'))):
    N = sorted({x for e in z['pe'] for x in e[:2]}); v = z['v']
    A = mk(z['pe'], z['pr'], N); ch = g['reduction'].reducePathAlgebra(g['mutation'].quiverMutationAtVertex(A, v))
    CA = np.array(g['cartan'](A), dtype=np.int64); CB = np.array(g['cartan'](ch), dtype=np.int64)
    ka, kb = g['search']._coxeterKeyOrNone(A), g['search']._coxeterKeyOrNone(ch)
    SA, SB = CA + CA.T, CB + CB.T
    sig = lambda S: (int((np.linalg.eigvalsh(S.astype(float)) > 1e-9).sum()), int((np.linalg.eigvalsh(S.astype(float)) < -1e-9).sum()))
    cp = lambda C: sp.factor(sp.Matrix(X.phi(C).tolist()).charpoly().as_expr())
    print('step', t, 'v', v, 'key equal', ka == kb, 'CA==CB', bool((CA == CB).all()), 'det S', round(np.linalg.det(SA.astype(float))), round(np.linalg.det(SB.astype(float))),
          'sig', sig(SA), sig(SB), 'Phi charpoly', cp(CA), '|', cp(CB))
    print('  CA rows', CA.tolist()); print('  CB rows', CB.tolist())
    for sl in (0.5, 5):
        try: print('  exact search slack', sl, X.congruent(CA, CB, sl), flush=True)
        except Exception as e: print('  exact search slack', sl, 'FAILED', type(e).__name__, e, flush=True)
