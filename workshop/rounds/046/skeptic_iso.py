"""Round 046 (skeptic): search for P in GL_n(Z) with P C_A P^T = C_B (columns f_i of P^T have entries in [-B, B]) by backtracking:
f_i must satisfy f_i^T C_B f_j = C_A[i][j] for all i, j (chi_B(f_i, f_j) = C_A[i,j]).  Finding P proves the Cartan forms congruent (necessary for derived
equivalence, NOT sufficient); not finding one in the box proves nothing.  Also tests the transposed form C_B^T.
Usage: skeptic_iso.py in.pkl BOX [kind]"""
import sys, pickle, itertools, time
import numpy as np
recs = pickle.load(open(sys.argv[1], 'rb')); B = int(sys.argv[2]); kind = sys.argv[3] if len(sys.argv) > 3 else 'fail'
def search(CA, CB, B, tlimit=60):
    n = len(CA); CA = np.array(CA); CB = np.array(CB)
    g = np.indices((2 * B + 1,) * n).reshape(n, -1).T.astype(np.int64) - B
    q = np.einsum('ij,jk,ik->i', g, CB, g); cand = g[q == CA[0, 0]]     # all diag are 1
    t0 = time.time(); sol = []
    def rec(i, chosen):
        if time.time() - t0 > tlimit: return None
        if i == n:
            P = np.array(chosen)
            return P if abs(round(np.linalg.det(P.astype(float)))) == 1 else False
        c = cand
        for j, f in enumerate(chosen):
            c = c[(c @ CB @ f == CA[i, j]) & (f @ CB @ c.T == CA[j, i])]
            if len(c) == 0: return False
        for f in c:
            r = rec(i + 1, chosen + [f])
            if r is None or (r is not False): return r
        return False
    return rec(0, []), len(cand)
for k, r in enumerate(x for x in recs if x['kind'] == kind):
    out = []
    for name, CB in (('C_B', r['CB']), ('C_B^T', np.array(r['CB']).T.tolist())):
        res, nc = search(r['CA'], CB, B)
        out.append((name, 'timeout' if res is None else ('FOUND' if res is not False else 'none'), nc))
    print(k, 'depth', r['depth'], out, flush=True)
