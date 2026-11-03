"""Round 111 toolsmith, item 2 (continued): YES-only box search for P in GL_7(Z), P C_A P^T = C_B, on the four E-147 class-2 steps (S positive semidefinite, rank 6;
the Salem-based exhaustive search does not apply).  Rows p_i of P in [-b,b]^7 with p C_A p^T = 1 (diagonal of C_B), depth-first with the pairwise conditions
p_i C_A p_j^T = C_B[i,j] both ways, det = +-1 at the end.  A hit proves Z-congruence; a miss proves nothing.
usage (repo root): timeout 10m .venv/bin/python workshop/rounds/111/toolsmith_c2_box.py [b]"""
import sys, pickle, ast, itertools, time
sys.path.insert(0, '.'); ARGV = sys.argv; sys.path.insert(0, 'workshop/rounds/111')
import numpy as np
import toolsmith_c1_13_exact as X
g = X.__dict__
b = int(ARGV[1]) if len(ARGV) > 1 else 2
def mk(edges, reprr, nodes):
    return g['build']([(a, c) for a, c, k in edges], [[[p[0][0]] + [a[1] for a in p] for p in d] for d in ast.literal_eval(reprr)], nodes)
V = np.array(list(itertools.product(range(-b, b + 1), repeat=7)), dtype=np.int64)
for t, z in enumerate(pickle.load(open('workshop/rounds/043/skeptic_zero_n7c2.pkl', 'rb'))):
    N = sorted({q for e in z['pe'] for q in e[:2]}); v = z['v']
    A = mk(z['pe'], z['pr'], N); ch = g['reduction'].reducePathAlgebra(g['mutation'].quiverMutationAtVertex(A, v))
    CA = np.array(g['cartan'](A), dtype=np.int64); CB = np.array(g['cartan'](ch), dtype=np.int64)
    q = np.einsum('ni,ij,nj->n', V, CA, V); cand = V[q == 1]; t0 = time.time()
    W = cand @ CA; sol = []
    def rec(rows):
        if sol: return
        i = len(rows)
        if i == 7:
            P = np.array([cand[r] for r in rows])
            if abs(round(np.linalg.det(P.astype(float)))) == 1 and (P @ CA @ P.T == CB).all(): sol.append(P)
            return
        ok = np.ones(len(cand), bool)
        for j, r in enumerate(rows):
            ok &= (cand @ CA @ cand[r] == CB[i, j]) & (cand[r] @ CA @ cand.T == CB[j, i])
        for k in np.nonzero(ok)[0]:
            rec(rows + [k])
            if sol: return
    rec([])
    print('step', t, 'box', b, 'candidates', len(cand), 'FOUND P (congruent)' if sol else 'no P in box', '%.0fs' % (time.time() - t0), flush=True)
    if sol: print(sol[0].tolist())
